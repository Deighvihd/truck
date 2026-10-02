use truck_topology::{shell::ShellCondition, *};

// Cases where truck version <= 0.4.0 fails
#[test]
fn singular_vertices_in_certain_boundaries() {
    let v = Vertex::news([(); 12]);
    let edge = [
        Edge::new(&v[0], &v[1], ()),
        Edge::new(&v[1], &v[2], ()),
        Edge::new(&v[0], &v[7], ()),
        Edge::new(&v[1], &v[8], ()),
        Edge::new(&v[2], &v[9], ()),
        Edge::new(&v[3], &v[4], ()),
        Edge::new(&v[4], &v[3], ()),
        Edge::new(&v[5], &v[6], ()),
        Edge::new(&v[6], &v[5], ()),
        Edge::new(&v[7], &v[8], ()),
        Edge::new(&v[8], &v[9], ()),
        Edge::new(&v[7], &v[10], ()),
        Edge::new(&v[9], &v[11], ()),
        Edge::new(&v[10], &v[11], ()),
    ];
    let shell: Shell<(), (), ()> = vec![
        Face::new(
            vec![
                wire![
                    edge[9].inverse(),
                    edge[2].inverse(),
                    edge[0].clone(),
                    edge[3].clone(),
                ],
                wire![edge[5].clone(), edge[6].clone()],
            ],
            (),
        ),
        Face::new(
            vec![
                wire![
                    edge[3].inverse(),
                    edge[1].clone(),
                    edge[4].clone(),
                    edge[10].inverse(),
                ],
                wire![edge[7].clone(), edge[8].clone()],
            ],
            (),
        ),
        Face::new(
            vec![wire![
                edge[10].clone(),
                edge[12].clone(),
                edge[13].inverse(),
                edge[11].inverse(),
                edge[9].clone(),
            ]],
            (),
        ),
    ]
    .into();
    assert_eq!(shell.shell_condition(), ShellCondition::Oriented);
    assert!(shell.singular_vertices().is_empty());
}

// `connected_components` used to return faces (and components) in hash-map
// order, which depends on memory addresses, so results varied run to run.
#[test]
fn connected_components_keep_shell_order() {
    // Two strips of triangles, their faces interleaved in the shell.
    let strip = || {
        let v = Vertex::news([(); 5]);
        let e = |i: usize, j: usize| Edge::new(&v[i], &v[j], ());
        let (e01, e12, e23, e34) = (e(0, 1), e(1, 2), e(2, 3), e(3, 4));
        let (e02, e13, e24) = (e(0, 2), e(1, 3), e(2, 4));
        [
            Face::new(vec![wire![e01.clone(), e12.clone(), e02.inverse()]], ()),
            Face::new(vec![wire![e12.inverse(), e13.clone(), e23.inverse()]], ()),
            Face::new(vec![wire![e23.clone(), e34.clone(), e24.inverse()]], ()),
        ]
    };
    for _ in 0..20 {
        let (a, b) = (strip(), strip());
        let shell: Shell<(), (), ()> = vec![
            a[0].clone(),
            b[0].clone(),
            a[1].clone(),
            b[1].clone(),
            a[2].clone(),
            b[2].clone(),
        ]
        .into();
        let ids = |faces: &[Face<(), (), ()>]| faces.iter().map(|f| f.id()).collect::<Vec<_>>();
        let components = shell.connected_components();
        assert_eq!(components.len(), 2);
        assert_eq!(
            ids(&components[0].face_iter().cloned().collect::<Vec<_>>()),
            ids(&a)
        );
        assert_eq!(
            ids(&components[1].face_iter().cloned().collect::<Vec<_>>()),
            ids(&b)
        );
    }
}
