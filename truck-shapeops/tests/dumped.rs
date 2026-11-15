//! Booleans dumped from rcad (`RCAD_DUMP_BOOLEANS`), as JSON: the
//! operation and the two solids. `DUMPS=dir cargo test -p truck-shapeops
//! --test dumped -- --ignored --nocapture` runs a directory of them.
use truck_modeling::*;
use truck_topology::compress::CompressedSolid;

/// Runs a dumped boolean.
fn run(json: &str) -> Option<Solid> {
    let json: serde_json::Value = serde_json::from_str(json).unwrap();
    let load = |k: &str| -> Solid {
        let c: CompressedSolid<Point3, Curve, Surface> =
            serde_json::from_value(json[k].clone()).unwrap();
        Solid::extract(c).unwrap()
    };
    let (target, mut tool) = (load("target"), load("tool"));
    let tol = 0.125;
    match json["op"].as_str().unwrap() {
        "Join" => truck_shapeops::or(&target, &tool, tol),
        "Intersect" => truck_shapeops::and(&target, &tool, tol),
        _ => {
            tool.not();
            truck_shapeops::and(&target, &tool, tol)
        }
    }
}

/// A ring joined to a rod it doesn't touch (their boxes overlap): no
/// intersection, so every face is sorted by a ray from one of its
/// vertices, and one ray grazed the rod's mesh, counting the ring's face
/// inside. The ring's faces vote now.
#[test]
fn faces_apart_vote_on_their_side() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/booleans/apart-grazing-ray.json"
    );
    let result = run(&std::fs::read_to_string(path).unwrap()).unwrap();
    assert_eq!(result.boundaries().len(), 2);
}

#[test]
#[ignore]
fn dumped_booleans() {
    let dir = std::env::var("DUMPS").unwrap();
    let mut paths: Vec<_> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .collect();
    paths.sort();
    for path in paths {
        let result = run(&std::fs::read_to_string(&path).unwrap());
        println!(
            "{}: {}",
            path.file_name().unwrap().to_string_lossy(),
            match &result {
                Some(s) => format!("{} shells", s.boundaries().len()),
                None => "NONE".into(),
            }
        );
    }
}

/// A join whose splitting of a face gave a loop through one vertex twice
/// (not a face's boundary) and a hole whose first point lay on the loop
/// around it (so no loop seemed to hold it): no result. Loops are split
/// where they come back to a vertex, and a hole goes in the loop around
/// most of its points.
#[test]
fn face_splits_through_a_vertex_twice() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/booleans/split-face-loops.json"
    );
    let result = run(&std::fs::read_to_string(path).unwrap()).unwrap();
    assert_eq!(result.boundaries().len(), 1);
}

/// A hole about four tolerances wide (radius 0.48) cut up through a tall
/// block from below its bottom face. The tool's lower half crosses the
/// bottom face along a half circle, and finding that curve's points
/// trusted any point within 10 tolerances of its leader: the far side of
/// the hole too. The edge's second half ran round the upper side, so the
/// bottom face had the upper half circle twice (its mesh doubled). Each
/// edge along the hole now stays on its own side of the hole's axis.
#[test]
fn small_hole_edges_keep_their_side() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/booleans/thin-hole-from-bottom.json"
    );
    let result = run(&std::fs::read_to_string(path).unwrap()).unwrap();
    assert_eq!(result.boundaries().len(), 1);
    // The hole's axis and radius (working units, as rcad dumped them).
    let (cx, cy, r) = (4.0101875, -1.9310656, 0.4796875);
    let tol = 0.125;
    let mut arcs = 0;
    for edge in result.edge_iter() {
        let curve = edge.curve();
        if matches!(curve, Curve::Line(_)) {
            continue;
        }
        arcs += 1;
        let (t0, t1) = curve.range_tuple();
        let side = (curve.subs((t0 + t1) / 2.0).y - cy).signum();
        for i in 0..=32 {
            let p = curve.subs(t0 + (t1 - t0) * i as f64 / 32.0);
            let off = ((p.x - cx).powi(2) + (p.y - cy).powi(2)).sqrt() - r;
            assert!(off.abs() < tol, "{p:?} is {off} off the hole");
            assert!((p.y - cy) * side > -tol, "{p:?} is across the hole");
        }
    }
    // Two half circles at each end of the hole, each seen from both of
    // its faces.
    assert_eq!(arcs, 8);
}

/// A prism joined through a block with a round hole, one of its long
/// edges grazing the hole's wall: the edge crosses the wall twice, close
/// together. The second intersection curve ending on that edge started
/// from the mesh's point between the crossings, and Newton's method ran
/// to the first crossing (already a vertex, the end of the edge as cut),
/// so the edge couldn't be cut there: no result. The curve's end is now
/// the crossing inside the edge nearest the mesh's point.
#[test]
fn an_edge_grazing_a_hole_crosses_it_twice() {
    use truck_meshalgo::prelude::*;
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/booleans/edge-grazing-a-hole.json"
    );
    let json = std::fs::read_to_string(path).unwrap();
    let result = run(&json).unwrap();
    assert_eq!(result.boundaries().len(), 1);
    // The tool sticks out of the block (into the hole and past the far
    // face), so the join holds more than the block, and less than both.
    let value: serde_json::Value = serde_json::from_str(&json).unwrap();
    let load = |k: &str| -> Solid {
        let c: CompressedSolid<Point3, Curve, Surface> =
            serde_json::from_value(value[k].clone()).unwrap();
        Solid::extract(c).unwrap()
    };
    let volume = |s: &Solid| s.triangulation(0.01).to_polygon().volume();
    let (block, tool) = (volume(&load("target")), volume(&load("tool")));
    let joined = volume(&result);
    assert!(
        block < joined && joined < block + tool,
        "{block} {tool} {joined}"
    );
}
