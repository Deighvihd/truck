use truck_geometry::prelude::*;

// A surface curved one way and straight the other (a long, shallow
// channel: a parabola swept 100 along y) is split only across. Splitting
// cells the way they were longest cut it into near-square cells along its
// whole straight length.
#[test]
fn a_surface_straight_one_way_is_split_only_across() {
    let surface = BSplineSurface::new(
        (KnotVec::bezier_knot(2), KnotVec::bezier_knot(1)),
        vec![
            vec![Point3::new(-1.0, 0.0, 0.0), Point3::new(-1.0, 100.0, 0.0)],
            vec![Point3::new(0.0, 0.0, 1.0), Point3::new(0.0, 100.0, 1.0)],
            vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 100.0, 0.0)],
        ],
    );
    let (udiv, vdiv) = surface.parameter_division(((0.0, 1.0), (0.0, 1.0)), 0.001);
    assert!(udiv.len() > 8, "{}", udiv.len());
    assert_eq!(vdiv.len(), 2, "{vdiv:?}");
}
