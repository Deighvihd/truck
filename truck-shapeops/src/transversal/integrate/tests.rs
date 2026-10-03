use truck_meshalgo::prelude::*;
use truck_modeling::*;

#[test]
fn punched_cube() {
    let v = builder::vertex(Point3::origin());
    let e = builder::tsweep(&v, Vector3::unit_x());
    let f = builder::tsweep(&e, Vector3::unit_y());
    let cube: Solid = builder::tsweep(&f, Vector3::unit_z());

    let v = builder::vertex(Point3::new(0.5, 0.25, -0.5));
    let w = builder::rsweep(
        &v,
        Point3::new(0.5, 0.5, 0.0),
        Vector3::unit_z(),
        Rad(7.0),
        3,
    );
    let f = builder::try_attach_plane(&[w]).unwrap();
    let mut cylinder = builder::tsweep(&f, Vector3::unit_z() * 2.0);
    cylinder.not();
    let and = crate::and(&cube, &cylinder, 0.05).unwrap();

    let poly = and.triangulation(0.01).to_polygon();
    let file = std::fs::File::create("punched-cube.obj").unwrap();
    obj::write(&poly, file).unwrap();
}

// A tool larger than the whole solid: the operation must fail with `None`,
// not panic building an unclosed solid.
#[test]
fn tool_engulfs_solid() {
    let v = builder::vertex(Point3::origin());
    let e = builder::tsweep(&v, Vector3::unit_x());
    let f = builder::tsweep(&e, Vector3::unit_y());
    let cube: Solid = builder::tsweep(&f, Vector3::unit_z());

    let v = builder::vertex(Point3::new(0.5, -1.5, -0.5));
    let w = builder::rsweep(
        &v,
        Point3::new(0.5, 0.5, 0.0),
        Vector3::unit_z(),
        Rad(7.0),
        3,
    );
    let f = builder::try_attach_plane(&[w]).unwrap();
    let mut cylinder = builder::tsweep(&f, Vector3::unit_z() * 2.0);
    cylinder.not();
    assert!(crate::and(&cube, &cylinder, 0.05).is_none());
}

// A tool engulfing the solid with a face on the solid's face: splitting the
// faces gave invalid wires and `Face::debug_new` panicked ("This wire is not
// simple.") in debug builds, or built an invalid face in release.
#[test]
fn engulfing_tool_with_coincident_cap() {
    let v = builder::vertex(Point3::origin());
    let e = builder::tsweep(&v, Vector3::unit_x() * 100.0);
    let f = builder::tsweep(&e, Vector3::unit_y() * 50.0);
    let plate: Solid = builder::tsweep(&f, Vector3::unit_z() * 25.0);

    let v0 = builder::vertex(Point3::new(200.0, 25.0, -1.0));
    let v1 = builder::vertex(Point3::new(-100.0, 25.0, -1.0));
    let arc0 = builder::circle_arc(&v0, &v1, Point3::new(50.0, 175.0, -1.0));
    let arc1 = builder::circle_arc(&v1, &v0, Point3::new(50.0, -125.0, -1.0));
    let disk: Face = builder::try_attach_plane(&[wire![arc0, arc1]]).unwrap();
    let mut cylinder: Solid = builder::tsweep(&disk, Vector3::unit_z() * 26.0);
    cylinder.not();
    assert!(crate::and(&plate, &cylinder, 0.125).is_none());
}

// A pocket cut from a plate's bottom face, its cap just past the face. The
// leader of each intersection curve was a quadratic B-spline fitted at
// uniform parameters of an unevenly spaced polyline, so it overshot and ran
// backwards in places (within tolerance of the curve). The tessellated
// boundary then crossed itself, and the pocket wall's mesh folded over.
#[test]
fn intersection_curves_do_not_backtrack() {
    let v = builder::vertex(Point3::origin());
    let e = builder::tsweep(&v, Vector3::unit_x() * 100.0);
    let f = builder::tsweep(&e, Vector3::unit_y() * 50.0);
    let plate: Solid = builder::tsweep(&f, Vector3::unit_z() * 25.0);

    let v0 = builder::vertex(Point3::new(60.0, 25.0, -1.0));
    let v1 = builder::vertex(Point3::new(40.0, 25.0, -1.0));
    let arc0 = builder::circle_arc(&v0, &v1, Point3::new(50.0, 35.0, -1.0));
    let arc1 = builder::circle_arc(&v1, &v0, Point3::new(50.0, 15.0, -1.0));
    let disk: Face = builder::try_attach_plane(&[wire![arc0, arc1]]).unwrap();
    let mut cylinder: Solid = builder::tsweep(&disk, Vector3::unit_z() * 13.5);
    cylinder.not();
    let pocketed = crate::and(&plate, &cylinder, 0.125).unwrap();

    let mut checked = 0;
    for edge in pocketed.edge_iter() {
        let curve = edge.curve();
        if !matches!(curve, Curve::IntersectionCurve(_)) {
            continue;
        }
        checked += 1;
        let (t0, t1) = curve.range_tuple();
        const N: usize = 1000;
        let points: Vec<Point3> = (0..=N)
            .map(|i| curve.subs(t0 + (t1 - t0) * i as f64 / N as f64))
            .collect();
        let backtracks = points
            .windows(3)
            .filter(|w| (w[1] - w[0]).dot(w[2] - w[1]) < 0.0)
            .count();
        assert_eq!(backtracks, 0, "{edge:?}");
    }
    assert!(checked > 0);

    // The pocket is 10 units in radius and 12.5 deep.
    let poly = pocketed.triangulation(0.025).to_polygon();
    let volume = poly.volume();
    let expected = 100.0 * 50.0 * 25.0 - std::f64::consts::PI * 100.0 * 12.5;
    assert!(
        (volume - expected).abs() < expected * 1e-4,
        "{volume} vs {expected}"
    );
}
