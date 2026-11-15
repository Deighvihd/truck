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
