use super::*;

type CompressedSolid = truck_topology::compress::CompressedSolid<Point3, Curve, Surface>;

macro_rules! dir ( () => { concat!(env!("CARGO_MANIFEST_DIR"), "/../resources/shape/") });

const SHAPE_JSONS: &[&str] = &[
    concat!(dir!(), "bottle.json"),
    concat!(dir!(), "punched-cube.json"),
    concat!(dir!(), "torus-punched-cube.json"),
    concat!(dir!(), "sphere.json"),
    concat!(dir!(), "torus.json"),
];

fn read_jsons() -> Vec<Vec<u8>> {
    let closure = |json| std::fs::read(json).unwrap();
    SHAPE_JSONS.iter().map(closure).collect()
}

#[test]
fn solid_is_closed() {
    for (i, json) in read_jsons().into_iter().enumerate() {
        let solid: Solid = serde_json::from_reader(json.as_slice()).unwrap();
        let mut poly = solid.triangulation(0.01).to_polygon();
        poly.put_together_same_attrs(TOLERANCE * 2.0)
            .remove_degenerate_faces()
            .remove_unused_attrs();
        assert_eq!(
            poly.shell_condition(),
            ShellCondition::Closed,
            "not closed: file no. {i}"
        );
    }
}

#[test]
fn csolid_is_closed() {
    for (i, json) in read_jsons().into_iter().enumerate() {
        let solid: CompressedSolid = serde_json::from_reader(json.as_slice()).unwrap();
        let mut poly = solid.triangulation(0.01).to_polygon();
        poly.put_together_same_attrs(TOLERANCE * 2.0)
            .remove_degenerate_faces()
            .remove_unused_attrs();
        assert_eq!(
            poly.shell_condition(),
            ShellCondition::Closed,
            "not closed: file no. {i}"
        );
    }
}

#[test]
fn compare_occt_mesh() {
    let jsons = read_jsons();
    let solid: Solid = serde_json::from_slice(jsons[2].as_slice()).unwrap();
    let res = solid.triangulation(0.005).to_polygon();
    let path = concat!(dir!(), "../obj/by_occt.obj");
    let ans = obj::read(std::fs::read(path).unwrap().as_slice()).unwrap();
    assert!(res.is_clung_to_by(ans.positions(), 0.05));
    assert!(ans.is_clung_to_by(res.positions(), 0.05));
}

#[test]
fn compare_occt_mesh_csolid() {
    let jsons = read_jsons();
    let solid: CompressedSolid = serde_json::from_slice(jsons[2].as_slice()).unwrap();
    let res = solid.triangulation(0.005).to_polygon();
    let path = concat!(dir!(), "../obj/by_occt.obj");
    let ans = obj::read(std::fs::read(path).unwrap().as_slice()).unwrap();
    assert!(res.is_clung_to_by(ans.positions(), 0.05));
    assert!(ans.is_clung_to_by(res.positions(), 0.05));
}

#[test]
fn large_number_meshing() {
    let json = std::fs::read(concat!(dir!(), "large-torus.json")).unwrap();
    let torus: Solid = serde_json::from_slice(json.as_slice()).unwrap();
    let _ = torus.triangulation(1.0).to_polygon();
}

fn special_cylinder_model() -> Shell {
    let v0 = builder::vertex(Point3::new(0.0, 1.0, 0.0));
    let v1 = builder::vertex(Point3::new(0.0, 1.0, 1.0));
    let v2 = builder::vertex(Point3::new(0.0, -1.0, 0.0));
    let v3 = builder::vertex(Point3::new(0.0, -1.0, 1.0));

    let edge0 = builder::line(&v0, &v1);
    let edge1 = builder::line(&v2, &v3);
    let edge2 = builder::circle_arc(&v0, &v2, Point3::new(-1.0, 0.0, 0.0));
    let edge3 = builder::circle_arc(&v2, &v0, Point3::new(1.0, 0.0, 0.0));
    let edge4 = builder::circle_arc(&v1, &v3, Point3::new(-1.0, 0.0, 1.0));
    let edge5 = builder::circle_arc(&v3, &v1, Point3::new(1.0, 0.0, 1.0));

    let face0 =
        builder::try_attach_plane(&[vec![edge2.inverse(), edge3.inverse()].into()]).unwrap();
    let face1 = builder::try_attach_plane(&[vec![edge4.clone(), edge5.clone()].into()]).unwrap();

    let surface_row = RevolutedCurve::<Curve>::by_revolution(
        Line(Point3::new(1.0, 0.0, 1.0), Point3::new(1.0, 0.0, 0.0)).into(),
        Point3::origin(),
        Vector3::unit_z(),
    );
    let surface: Surface = Processor::new(surface_row).into();

    let face2 = Face::new(
        vec![vec![edge2, edge1.clone(), edge4.inverse(), edge0.inverse()].into()],
        surface.clone(),
    );
    let face3 = Face::new(
        vec![vec![edge3, edge0, edge5.inverse(), edge1.inverse()].into()],
        surface,
    );

    vec![face0, face1, face2, face3].into()
}

#[test]
fn special_cylinder() {
    let shell = special_cylinder_model();
    let mut mesh = shell.triangulation(0.01).to_polygon();
    mesh.put_together_same_attrs(TOLERANCE)
        .remove_degenerate_faces()
        .remove_unused_attrs();
    assert_eq!(mesh.shell_condition(), ShellCondition::Closed);
}

#[test]
fn special_cylinder_csolid() {
    let shell = special_cylinder_model().compress();
    let mut mesh = shell.triangulation(0.01).to_polygon();
    mesh.put_together_same_attrs(TOLERANCE)
        .remove_degenerate_faces()
        .remove_unused_attrs();
    assert_eq!(mesh.shell_condition(), ShellCondition::Closed);
}

#[test]
fn robust_closed() {
    let cube: Solid = {
        let v = builder::vertex(Point3::origin());
        let e = builder::tsweep(&v, Vector3::unit_x());
        let f = builder::tsweep(&e, Vector3::unit_y());
        builder::tsweep(&f, Vector3::unit_z())
    };

    let o = Point3::new(0.5, 0.5, 0.5);
    cube.edge_iter().for_each(|edge| {
        let curve = edge.curve();

        if let Curve::Line(line) = curve {
            let m = line.subs(0.5);
            let p = m + 0.2 * (o - m);
            let bsp = BSplineCurve::new(KnotVec::bezier_knot(2), vec![line.0, p, line.1]);
            edge.set_curve(bsp.into());
        }
    });

    assert!(cube
        .triangulation(0.01)
        .face_iter()
        .all(|face| face.surface().is_none()));

    let mut mesh = cube.robust_triangulation(0.01).to_polygon();
    mesh.put_together_same_attrs(TOLERANCE2)
        .remove_unused_attrs();
    assert_eq!(mesh.shell_condition(), ShellCondition::Closed);
}

/// The segments of a face's boundary that its mesh leaves out. The face
/// is JSON: its surface, orientation, and wires of polylines (one per
/// edge, oriented as the wire runs).
fn missing_boundary(json: &str) -> usize {
    type TEdge = truck_topology::Edge<Point3, PolylineCurve<Point3>>;
    type TWire = truck_topology::Wire<Point3, PolylineCurve<Point3>>;
    type TFace = truck_topology::Face<Point3, PolylineCurve<Point3>, Surface>;
    type TShell = truck_topology::Shell<Point3, PolylineCurve<Point3>, Surface>;
    let data: serde_json::Value = serde_json::from_str(json).unwrap();
    let surface: Surface = serde_json::from_value(data["surface"].clone()).unwrap();
    let mut wires: Vec<Vec<Vec<[f64; 3]>>> = serde_json::from_value(data["wires"].clone()).unwrap();
    // As the face runs (an inverted face's wires backwards): as the
    // surface runs.
    let inverted = data["orientation"] == false;
    if inverted {
        wires.iter_mut().for_each(|wire| {
            wire.reverse();
            wire.iter_mut().for_each(|e| e.reverse());
        });
    }
    let wires: Vec<TWire> = wires
        .iter()
        .map(|wire| {
            let points =
                |e: &Vec<[f64; 3]>| -> Vec<Point3> { e.iter().map(|&p| p.into()).collect() };
            let vertices: Vec<_> = wire
                .iter()
                .map(|e| truck_topology::Vertex::new(points(e)[0]))
                .collect();
            wire.iter()
                .enumerate()
                .map(|(i, e)| {
                    let next = &vertices[(i + 1) % vertices.len()];
                    TEdge::new(&vertices[i], next, PolylineCurve(points(e)))
                })
                .collect()
        })
        .collect();
    let mut face = TFace::new(wires, surface);
    if inverted {
        face.invert();
    }
    let shell: TShell = vec![face].into_iter().collect();
    let meshed = shell.triangulation(0.025);
    let face = &meshed[0];
    let polygon = face.surface().unwrap();
    let positions = polygon.positions();
    // Positions to a billionth (consecutive edges' polylines can end and
    // start a rounding error apart).
    let key = |p: Point3| [p.x, p.y, p.z].map(|x| (x * 1.0e9).round() as i64);
    let mut edges = std::collections::HashSet::new();
    for tri in polygon.faces().triangle_iter() {
        for k in 0..3 {
            let (a, b) = (
                key(positions[tri[k].pos]),
                key(positions[tri[(k + 1) % 3].pos]),
            );
            edges.insert((a, b));
            edges.insert((b, a));
        }
    }
    face.boundaries()
        .iter()
        .flat_map(|w| w.edge_iter())
        .map(|e| {
            let c = e.oriented_curve();
            c.windows(2)
                .filter(|w| key(w[0]) != key(w[1]) && !edges.contains(&(key(w[0]), key(w[1]))))
                .count()
        })
        .sum()
}

/// Faces from booleans whose meshes left boundary segments out: slivers
/// along the boundary whose center lies within rounding of it (taken as
/// outside), and a boundary vertex whose segment can't be added, so a
/// chord past it is (a sliver's side is not the boundary's).
#[test]
fn trimmed_faces_mesh_their_whole_boundary() {
    for name in ["sliver-plane", "sliver-nurbs", "chord-plane"] {
        let path = format!(
            "{}/tests/tessellation/trimming/{name}.json",
            env!("CARGO_MANIFEST_DIR")
        );
        let json = std::fs::read_to_string(path).unwrap();
        assert_eq!(missing_boundary(&json), 0, "{name}");
    }
}

/// A face bounded by two edges along the same line (a sliver a boolean can
/// leave between coincident edges) has no area: it meshes to nothing, not
/// to the whole parameter range of its surface.
#[test]
fn flat_face_meshes_to_nothing() {
    let plane = Plane::new(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    );
    let v = builder::vertices([Point3::new(-2.0, 0.5, 0.0), Point3::new(2.0, 0.5, 0.0)]);
    let wire: Wire = vec![builder::line(&v[0], &v[1]), builder::line(&v[1], &v[0])].into();
    let face = Face::new(vec![wire], Surface::Plane(plane));
    let shell: Shell = vec![face].into();
    let meshed = shell.triangulation(0.01);
    let polygon = meshed[0].surface().unwrap();
    assert_eq!(polygon.faces().len(), 0);
}

/// A flat loop along a whole side of its surface's parameter range is a
/// seam of a surface closed that way (a full revolution's band): the face
/// is the whole range, as before.
#[test]
fn flat_loop_along_a_side_is_a_seam() {
    let plane = Plane::new(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    );
    let v = builder::vertices([Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)]);
    let wire: Wire = vec![builder::line(&v[0], &v[1]), builder::line(&v[1], &v[0])].into();
    let face = Face::new(vec![wire], Surface::Plane(plane));
    let shell: Shell = vec![face].into();
    let meshed = shell.triangulation(0.01);
    let polygon = meshed[0].surface().unwrap();
    assert!(polygon.faces().len() > 0);
}

/// A circle smaller than the tolerance, in two halves: each half within
/// the tolerance of the one chord between its ends would mesh as that
/// chord, there and back, a disk of no area. Bent edges get two segments.
#[test]
fn a_circle_smaller_than_the_tolerance_keeps_its_area() {
    let r = 0.01;
    let plane = Plane::new(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    );
    let v = builder::vertices([Point3::new(r, 0.0, 0.0), Point3::new(-r, 0.0, 0.0)]);
    let wire: Wire = vec![
        builder::circle_arc(&v[0], &v[1], Point3::new(0.0, r, 0.0)),
        builder::circle_arc(&v[1], &v[0], Point3::new(0.0, -r, 0.0)),
    ]
    .into();
    let face = Face::new(vec![wire], Surface::Plane(plane));
    let shell: Shell = vec![face].into();
    let meshed = shell.triangulation(0.05);
    let polygon = meshed[0].surface().unwrap();
    let positions = polygon.positions();
    let area: f64 = polygon
        .faces()
        .triangle_iter()
        .map(|t| {
            let [a, b, c] = t.map(|v| positions[v.pos]);
            (b - a).cross(c - a).magnitude() / 2.0
        })
        .sum();
    // A square inscribed in the circle: 2r².
    assert!(area > r * r, "area {area}");
}
