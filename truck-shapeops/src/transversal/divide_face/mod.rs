#![allow(clippy::many_single_char_names)]

use super::faces_classification::FacesClassification;
use super::loops_store::*;
use rustc_hash::FxHashMap as HashMap;
use std::ops::Deref;
use truck_meshalgo::prelude::*;
use truck_topology::*;

fn create_parameter_boundary<P, C, S>(
    face: &Face<P, C, S>,
    wire: &Wire<P, C>,
    polys: &mut HashMap<EdgeID<C>, PolylineCurve<P>>,
    tol: f64,
) -> Option<PolylineCurve<Point2>>
where
    P: Copy,
    C: BoundedCurve<Point = P> + ParameterDivision1D<Point = P>,
    S: Clone + SearchParameter<D2, Point = P>,
{
    let surface = face.surface();
    let pt = wire.front_vertex().unwrap().point();
    let p: Point2 = surface.search_parameter(pt, None, 100)?.into();
    let vec = wire.edge_iter().try_fold(vec![p], |mut vec, edge| {
        let poly = polys.entry(edge.id()).or_insert_with(|| {
            let curve = edge.curve();
            let div = curve.parameter_division(curve.range_tuple(), tol).1;
            PolylineCurve(div)
        });
        let mut p = *vec.last().unwrap();
        let closure = |q: &P| -> Option<Point2> {
            p = surface.search_parameter(*q, Some(p.into()), 100)?.into();
            Some(p)
        };
        let add: Option<Vec<Point2>> = match edge.orientation() {
            true => poly.iter().skip(1).map(closure).collect(),
            false => poly.iter().rev().skip(1).map(closure).collect(),
        };
        vec.append(&mut add?);
        Some(vec)
    })?;
    Some(PolylineCurve(vec))
}

#[derive(Clone, Debug)]
struct WireChunk<C> {
    poly: PolylineCurve<Point2>,
    wire: BoundaryWire<Point3, C>,
}

/// `wire` (closed) as simple loops: where it comes back to a vertex it
/// passed (an intersection loop touching the face's boundary at a point),
/// the part in between is a loop of its own. A wire through a vertex twice
/// can't bound a face.
fn simple_loops<P, C>(wire: &Wire<P, C>) -> Vec<Wire<P, C>> {
    let mut loops = Vec::new();
    let mut path: Vec<Edge<P, C>> = Vec::new();
    for edge in wire.edge_iter() {
        path.push(edge.clone());
        let back = edge.back().id();
        if let Some(k) = path.iter().position(|e| e.front().id() == back) {
            if k > 0 || path.len() < wire.len() {
                let looped: Wire<P, C> = path.drain(k..).collect();
                loops.push(looped);
            }
        }
    }
    if !path.is_empty() {
        loops.push(path.into_iter().collect());
    }
    loops
}

type FaceWithShapesOpStatus<C, S> = (Face<Point3, C, S>, ShapesOpStatus);
fn divide_one_face<C, S>(
    face: &Face<Point3, C, S>,
    loops: &Loops<Point3, C>,
    tol: f64,
) -> Option<Vec<FaceWithShapesOpStatus<C, S>>>
where
    C: BoundedCurve<Point = Point3> + ParameterDivision1D<Point = Point3>,
    S: Clone + SearchParameter<D2, Point = Point3>,
{
    let (mut pre_faces, mut negative_wires) = (Vec::new(), Vec::new());
    let mut map = HashMap::default();
    let wires = loops.iter().flat_map(|wire| {
        simple_loops(wire)
            .into_iter()
            .map(|w| BoundaryWire::new(w, wire.status()))
    });
    wires.into_iter().try_for_each(|wire| {
        let poly = create_parameter_boundary(face, &wire, &mut map, tol)?;
        match poly.area() > 0.0 {
            true => pre_faces.push(vec![WireChunk { poly, wire }]),
            false => negative_wires.push(WireChunk { poly, wire }),
        }
        Some(())
    })?;
    // A hole belongs to the innermost outer loop around it: with nested
    // loops (a ring-shaped tool crossing the face), the first one found
    // can be an outer loop around another one, which then gets no hole.
    // Around it: around most of its points (one can lie on the outer
    // loop, where the inclusion test can go either way). A hole of no area
    // (a sliver between coincident edges) holds nothing.
    negative_wires.into_iter().try_for_each(|chunk| {
        let hole_area = chunk.poly.area();
        let bdb: BoundingBox<Point2> = chunk.poly.iter().collect();
        if hole_area.abs() <= 1.0e-9 * bdb.diagonal().magnitude2() {
            return Some(());
        }
        let points: Vec<Point2> = chunk.poly.iter().copied().collect();
        let op = pre_faces
            .iter_mut()
            .filter(|face| {
                // Not its own twin (the same loop the other way round,
                // which "includes" the hole's points on its boundary).
                let area = face[0].poly.area();
                let inside = points.iter().filter(|p| face[0].poly.include(**p)).count();
                area > -hole_area * (1.0 + 1e-6) && 2 * inside > points.len()
            })
            .min_by(|a, b| a[0].poly.area().total_cmp(&b[0].poly.area()))?;
        op.push(chunk);
        Some(())
    })?;
    let vec = pre_faces
        .into_iter()
        .map(|pre_face| {
            let surface = face.surface();
            let op = pre_face
                .iter()
                .find(|chunk| chunk.wire.status() != ShapesOpStatus::Unknown);
            let status = match op {
                Some(chunk) => chunk.wire.status(),
                None => ShapesOpStatus::Unknown,
            };
            let wires: Vec<Wire<Point3, C>> = pre_face
                .into_iter()
                .map(|chunk| chunk.wire.deref().clone())
                .collect();
            // Degenerate splits (e.g. when the other solid removes the whole
            // face) give invalid wires: fail the operation, don't panic.
            let mut new_face = Face::try_new(wires, surface).ok()?;
            if !face.orientation() {
                new_face.invert();
            }
            Some((new_face, status))
        })
        .collect::<Option<Vec<_>>>()?;
    Some(vec)
}

pub fn divide_faces<C, S>(
    shell: &Shell<Point3, C, S>,
    loops_store: &LoopsStore<Point3, C>,
    tol: f64,
) -> Option<FacesClassification<Point3, C, S>>
where
    C: BoundedCurve<Point = Point3> + ParameterDivision1D<Point = Point3>,
    S: Clone + SearchParameter<D2, Point = Point3>,
{
    let mut res = FacesClassification::<Point3, C, S>::default();
    shell
        .iter()
        .zip(loops_store)
        .try_for_each(|(face, loops)| {
            if loops
                .iter()
                .all(|wire| wire.status() == ShapesOpStatus::Unknown)
            {
                res.push(face.clone(), ShapesOpStatus::Unknown);
            } else {
                let vec = divide_one_face(face, loops, tol)?;
                vec.into_iter()
                    .for_each(|(face, status)| res.push(face, status));
            }
            Some(())
        })?;
    Some(res)
}

#[cfg(test)]
mod tests;
