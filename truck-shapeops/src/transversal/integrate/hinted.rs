use std::ops::Bound;
use truck_geometry::prelude::*;
use truck_meshalgo::prelude::PolylineCurve;

/// An intersection curve whose points are found starting from the surface
/// parameters of its polyline leader's vertices, instead of searching both
/// surfaces from scratch for every point.
///
/// Approximating an intersection curve by a B-spline samples it many times,
/// and each unhinted sample searches a grid over both surfaces, which made
/// booleans with curved faces take seconds.
///
/// The parameter is arc length along the leader. The leader's vertices are
/// unevenly spaced, and a quadratic B-spline interpolating samples taken at
/// uniform leader parameters overshoots between them, running backwards in
/// places.
#[derive(Clone, Debug)]
pub(super) struct HintedIntersectionCurve<'a, S> {
    curve: &'a IntersectionCurve<PolylineCurve<Point3>, S, S>,
    /// Parameters on (surface0, surface1) of each leader vertex.
    hints: Vec<(Point2, Point2)>,
    /// Arc length along the leader at each of its vertices.
    lengths: Vec<f64>,
    /// How far a hinted result may be from the leader before it is distrusted.
    bound: f64,
}

const TRIALS: usize = 100;

impl<'a, S> HintedIntersectionCurve<'a, S>
where S: ParametricSurface3D + SearchNearestParameter<D2, Point = Point3>
{
    /// Walks the leader's vertices, each starting from the previous one's
    /// parameters. Falls back to a full search when a hinted step fails or
    /// lands far from the leader (e.g. converged to another branch).
    pub(super) fn new(
        curve: &'a IntersectionCurve<PolylineCurve<Point3>, S, S>,
        tol: f64,
    ) -> Option<Self> {
        let bound = 10.0 * tol;
        let leader = curve.leader();
        let mut hints = Vec::with_capacity(leader.len());
        let mut previous: Option<(Point2, Point2)> = None;
        for (i, target) in leader.iter().enumerate() {
            let t = i as f64;
            let walked = previous.and_then(|(uv0, uv1)| {
                curve
                    .search_triple_with_hints(t, uv0.into(), uv1.into(), TRIALS)
                    .filter(|(p, _, _)| p.distance(*target) < bound)
            });
            let (_, uv0, uv1) = match walked {
                Some(found) => found,
                None => curve.search_triple(t, TRIALS)?,
            };
            hints.push((uv0, uv1));
            previous = Some((uv0, uv1));
        }
        let mut lengths = Vec::with_capacity(leader.len());
        let mut length = 0.0;
        for (i, point) in leader.iter().enumerate() {
            if i > 0 {
                length += point.distance(leader[i - 1]);
            }
            lengths.push(length);
        }
        Some(Self {
            curve,
            hints,
            lengths,
            bound,
        })
    }

    /// The leader parameter at arc length `s`, and its rate of change.
    fn leader_parameter(&self, s: f64) -> (f64, f64) {
        let last = self.lengths.len().saturating_sub(1);
        // The segment containing `s`: the last vertex at or before it.
        let i = self
            .lengths
            .partition_point(|l| *l <= s)
            .saturating_sub(1)
            .min(last.saturating_sub(1));
        let (l0, l1) = (self.lengths[i], self.lengths[(i + 1).min(last)]);
        match l1 - l0 > 0.0 {
            true => (
                i as f64 + ((s - l0) / (l1 - l0)).clamp(0.0, 1.0),
                1.0 / (l1 - l0),
            ),
            false => (i as f64, 0.0),
        }
    }
}

impl<S> ParametricCurve for HintedIntersectionCurve<'_, S>
where S: ParametricSurface3D + SearchNearestParameter<D2, Point = Point3>
{
    type Point = Point3;
    type Vector = Vector3;
    fn subs(&self, s: f64) -> Point3 {
        let (t, _) = self.leader_parameter(s);
        let last = self.hints.len().saturating_sub(1);
        let i = (t.max(0.0) as usize).min(last.saturating_sub(1));
        let s = (t - i as f64).clamp(0.0, 1.0);
        let (a0, a1) = self.hints[i];
        let (b0, b1) = self.hints[(i + 1).min(last)];
        let hint0 = a0 + (b0 - a0) * s;
        let hint1 = a1 + (b1 - a1) * s;
        let target = self.curve.leader().subs(t);
        self.curve
            .search_triple_with_hints(t, hint0.into(), hint1.into(), TRIALS)
            .filter(|(p, _, _)| p.distance(target) < self.bound)
            .map(|(p, _, _)| p)
            .unwrap_or_else(|| self.curve.subs(t))
    }
    fn der(&self, s: f64) -> Vector3 {
        let (t, dt) = self.leader_parameter(s);
        self.curve.der(t) * dt
    }
    fn der2(&self, s: f64) -> Vector3 {
        let (t, dt) = self.leader_parameter(s);
        self.curve.der2(t) * dt * dt
    }
    fn der_n(&self, n: usize, s: f64) -> Vector3 {
        let (t, dt) = self.leader_parameter(s);
        self.curve.der_n(n, t) * dt.powi(n as i32)
    }
    fn parameter_range(&self) -> ParameterRange {
        let length = self.lengths.last().copied().unwrap_or(0.0);
        (Bound::Included(0.0), Bound::Included(length))
    }
}

impl<S> BoundedCurve for HintedIntersectionCurve<'_, S> where S: ParametricSurface3D + SearchNearestParameter<D2, Point = Point3> {}
