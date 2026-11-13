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

fn square_wire(x0: f64, y0: f64, x1: f64, y1: f64, z: f64) -> Wire {
    let v = [
        builder::vertex(Point3::new(x0, y0, z)),
        builder::vertex(Point3::new(x1, y0, z)),
        builder::vertex(Point3::new(x1, y1, z)),
        builder::vertex(Point3::new(x0, y1, z)),
    ];
    (0..4)
        .map(|i| builder::line(&v[i], &v[(i + 1) % 4]))
        .collect()
}

// A tool with a hole through it (a frame) crossing the plate's top face:
// the face is split along two loops, one inside the other. The hole loop
// was put into the first outer loop containing it (the face's own
// boundary) rather than the innermost (the frame's outline), and a loop's
// twin "contained" it too, so the split faces had non-simple wires.
#[test]
fn frame_shaped_tool() {
    let v = builder::vertex(Point3::origin());
    let e = builder::tsweep(&v, Vector3::unit_x() * 100.0);
    let f = builder::tsweep(&e, Vector3::unit_y() * 50.0);
    let plate: Solid = builder::tsweep(&f, Vector3::unit_z() * 25.0);

    let outer = square_wire(10.0, 10.0, 90.0, 40.0, 10.0);
    let inner = square_wire(40.0, 20.0, 60.0, 30.0, 10.0).inverse();
    let frame: Face = builder::try_attach_plane(&[outer, inner]).unwrap();
    let tool: Solid = builder::tsweep(&frame, Vector3::unit_z() * 20.0);

    let joined = crate::or(&plate, &tool, 0.05).unwrap();
    assert_eq!(joined.boundaries().len(), 1);
    let mut hole = tool.clone();
    hole.not();
    let cut = crate::and(&plate, &hole, 0.05).unwrap();
    assert_eq!(cut.boundaries().len(), 1);
    // The plate's top: the rim round the pocket, and the island's top.
    let tops = cut.boundaries()[0]
        .face_iter()
        .filter(|f| {
            let p = f.boundaries()[0].front_vertex().unwrap().point();
            f.boundaries().iter().all(|w| {
                w.vertex_iter()
                    .all(|v| (v.point().z - 25.0).abs() < 1e-6 && p.z == v.point().z)
            })
        })
        .count();
    assert_eq!(tops, 2);
}

// A box joined to a cylinder, reaching out past its wall near the top:
// the wall keeps a notch where the box's sides meet it. Undivided along
// its height (it's straight that way), the wall's mesh had only the top
// and bottom rims and the notch's corners, and triangles spanned a third
// of the way round, cutting deep inside: the join measured less than the
// cylinder alone.
#[test]
fn notched_wall_meshes_on_its_surface() {
    let v = builder::vertex(Point3::new(9.5, 0.0, 0.0));
    let w = builder::rsweep(&v, Point3::origin(), Vector3::unit_z(), Rad(7.0), 3);
    let f = builder::try_attach_plane(&[w]).unwrap();
    let cylinder: Solid = builder::tsweep(&f, Vector3::unit_z() * 27.7);

    let top = square_wire(-5.5, 0.2, 5.2, 9.74, 24.0);
    let face: Face = builder::try_attach_plane(&[top]).unwrap();
    let tool: Solid = builder::tsweep(&face, Vector3::unit_z() * 3.0);

    let joined = crate::or(&cylinder, &tool, 0.05).unwrap();
    let volume = |s: &Solid| s.triangulation(0.01).to_polygon().volume();
    let (before, after) = (volume(&cylinder), volume(&joined));
    assert!(after > before, "{before} → {after}");
}
