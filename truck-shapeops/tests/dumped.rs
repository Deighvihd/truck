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
