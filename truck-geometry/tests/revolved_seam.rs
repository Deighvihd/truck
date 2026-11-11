use truck_geometry::prelude::*;

/// Points on the profile itself (angle 0) are found, as at any other angle.
#[test]
fn points_on_the_profile_are_found() {
    // 98.04249999999999, not 98.0425: one rounding off, as a scaled model's
    // coordinates are.
    let r = 39.217 * (100.0 / 40.0);
    let line = Line(Point3::new(r, 1.25, 0.0), Point3::new(r, 47.3225, 0.0));
    let rc = RevolutedCurve::by_revolution(line, Point3::origin(), Vector3::unit_y());
    for p in [
        Point3::new(r, 47.3225, 0.0),
        Point3::new(r, 1.25, 0.0),
        Point3::new(r, 20.0, 0.0),
        Point3::new(0.0, 1.25, -r),
    ] {
        let found = rc.search_parameter(p, None, 100);
        let (u, v) = found.unwrap_or_else(|| panic!("{p:?} not found"));
        assert!(rc.subs(u, v).near(&p), "{p:?}: {u} {v}");
    }
}
