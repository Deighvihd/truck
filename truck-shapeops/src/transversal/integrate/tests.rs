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
