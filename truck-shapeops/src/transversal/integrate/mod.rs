use crate::alternative::Alternative;

use super::*;
use truck_geometry::prelude::*;
use truck_meshalgo::prelude::*;
use truck_topology::*;

mod hinted;
use hinted::HintedIntersectionCurve;

/// Only solids consisting of faces whose surface is implemented this trait can be used for set operations.
pub trait ShapeOpsSurface:
    ParametricSurface3D
    + ParameterDivision2D
    + SearchParameter<D2, Point = Point3>
    + SearchNearestParameter<D2, Point = Point3>
    + Invertible
    + Send
    + Sync {
}
impl<S> ShapeOpsSurface for S where S: ParametricSurface3D
        + ParameterDivision2D
        + SearchParameter<D2, Point = Point3>
        + SearchNearestParameter<D2, Point = Point3>
        + Invertible
        + Send
        + Sync
{
}

/// Only solids consisting of edges whose curve is implemented this trait can be used for set operations.
pub trait ShapeOpsCurve<S: ShapeOpsSurface>:
    ParametricCurve3D
    + ParameterDivision1D<Point = Point3>
    + Cut
    + Invertible
    + From<IntersectionCurve<BSplineCurve<Point3>, S, S>>
    + SearchParameter<D1, Point = Point3>
    + SearchNearestParameter<D1, Point = Point3>
    + Send
    + Sync {
}
impl<C, S: ShapeOpsSurface> ShapeOpsCurve<S> for C where C: ParametricCurve3D
        + ParameterDivision1D<Point = Point3>
        + Cut
        + Invertible
        + From<IntersectionCurve<BSplineCurve<Point3>, S, S>>
        + SearchParameter<D1, Point = Point3>
        + SearchNearestParameter<D1, Point = Point3>
        + Send
        + Sync
{
}

type AltCurveShell<C, S> =
    Shell<Point3, Alternative<C, IntersectionCurve<PolylineCurve<Point3>, S, S>>, S>;

fn altshell_to_shell<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    altshell: &AltCurveShell<C, S>,
    tol: f64,
) -> Option<Shell<Point3, C, S>> {
    altshell.try_mapped(
        |p| Some(*p),
        |c| match c {
            Alternative::FirstType(c) => Some(c.clone()),
            Alternative::SecondType(ic) => {
                let hinted = HintedIntersectionCurve::new(ic, tol)?;
                // A curve of no length can't be approximated (its knot
                // vector would be scaled by 0).
                let (t0, t1) = hinted.range_tuple();
                if !(t1 - t0 > 0.0) {
                    return None;
                }
                let bsp =
                    BSplineCurve::quadratic_approximation(&hinted, hinted.range_tuple(), tol, 100)?;
                Some(
                    IntersectionCurve::new(ic.surface0().clone(), ic.surface1().clone(), bsp)
                        .into(),
                )
            }
        },
        |s| Some(s.clone()),
    )
}

fn process_one_pair_of_shells<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell0: &Shell<Point3, C, S>,
    shell1: &Shell<Point3, C, S>,
    tol: f64,
) -> Option<[Shell<Point3, C, S>; 2]> {
    nonpositive_tolerance!(tol);
    let poly_shell0 = shell0.triangulation(tol);
    let poly_shell1 = shell1.triangulation(tol);
    let altshell0: AltCurveShell<C, S> =
        shell0.mapped(|x| *x, |c| Alternative::FirstType(c.clone()), Clone::clone);
    let altshell1: AltCurveShell<C, S> =
        shell1.mapped(|x| *x, |c| Alternative::FirstType(c.clone()), Clone::clone);
    let loops_store::LoopsStoreQuadruple {
        geom_loops_store0: loops_store0,
        geom_loops_store1: loops_store1,
        ..
    } = loops_store::create_loops_stores(&altshell0, &poly_shell0, &altshell1, &poly_shell1)?;
    let mut cls0 = divide_face::divide_faces(&altshell0, &loops_store0, tol)?;
    cls0.integrate_by_component();
    let mut cls1 = divide_face::divide_faces(&altshell1, &loops_store1, tol)?;
    cls1.integrate_by_component();
    let [mut and0, mut or0, unknown0] = cls0.and_or_unknown();
    let [mut and1, mut or1, unknown1] = cls1.and_or_unknown();
    classify_unknown(unknown0, &poly_shell1, &mut and0, &mut or0)?;
    classify_unknown(unknown1, &poly_shell0, &mut and1, &mut or1)?;
    and0.append(&mut and1);
    or0.append(&mut or1);
    Some([
        altshell_to_shell(&and0, tol)?,
        altshell_to_shell(&or0, tol)?,
    ])
}

/// Sorts the faces no intersection loop touched (`unknown`) into those
/// inside the other shell (`and`) and those outside (`or`), by casting a
/// ray from a vertex of each face through the other shell's mesh. A ray
/// grazing the mesh (through an edge or a vertex of its triangles) can
/// miscount, so the faces vote: a connected piece of such faces lies all
/// on one side, the side three quarters of its faces' rays find (`None` if
/// neither side gets that many).
fn classify_unknown<C, S>(
    unknown: AltCurveShell<C, S>,
    other: &Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>>,
    and: &mut AltCurveShell<C, S>,
    or: &mut AltCurveShell<C, S>,
) -> Option<()>
where
    C: ShapeOpsCurve<S>,
    S: ShapeOpsSurface,
{
    for component in unknown.connected_components() {
        let (mut inside, mut outside) = (0_usize, 0_usize);
        for face in component.face_iter() {
            let pt = face.boundaries()[0].vertex_iter().next()?.point();
            let dir = hash::take_one_unit(pt);
            let count = other.iter().try_fold(0, |count, face| {
                let poly = face.surface()?;
                Some(count + poly.signed_crossing_faces(pt, dir))
            })?;
            match count >= 1 {
                true => inside += 1,
                false => outside += 1,
            }
        }
        // Three quarters of the rays must agree: rays that disagree more
        // than a stray graze does mean the meshes can't tell.
        let total = inside + outside;
        if 4 * inside >= 3 * total {
            and.extend(component);
        } else if 4 * outside >= 3 * total {
            or.extend(component);
        } else {
            return None;
        }
    }
    Some(())
}

/// AND operation between two solids.
pub fn and<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> Option<Solid<Point3, C, S>> {
    let mut iter0 = solid0.boundaries().iter();
    let mut iter1 = solid1.boundaries().iter();
    let shell0 = iter0.next()?;
    let shell1 = iter1.next()?;
    let [mut and_shell, _] = process_one_pair_of_shells(shell0, shell1, tol)?;
    for shell in iter0 {
        let [res, _] = process_one_pair_of_shells(&and_shell, shell, tol)?;
        and_shell = res;
    }
    for shell in iter1 {
        let [res, _] = process_one_pair_of_shells(&and_shell, shell, tol)?;
        and_shell = res;
    }
    let boundaries = and_shell.connected_components();
    Solid::try_new(boundaries).ok()
}

/// OR operation between two solids.
pub fn or<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> Option<Solid<Point3, C, S>> {
    let mut iter0 = solid0.boundaries().iter();
    let mut iter1 = solid1.boundaries().iter();
    let shell0 = iter0.next()?;
    let shell1 = iter1.next()?;
    let [_, mut or_shell] = process_one_pair_of_shells(shell0, shell1, tol)?;
    for shell in iter0 {
        let [_, res] = process_one_pair_of_shells(&or_shell, shell, tol)?;
        or_shell = res;
    }
    for shell in iter1 {
        let [_, res] = process_one_pair_of_shells(&or_shell, shell, tol)?;
        or_shell = res;
    }
    let boundaries = or_shell.connected_components();
    Solid::try_new(boundaries).ok()
}

#[cfg(test)]
mod tests;
