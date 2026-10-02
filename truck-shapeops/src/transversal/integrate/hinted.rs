use truck_geometry::prelude::*;
use truck_meshalgo::prelude::PolylineCurve;

/// An intersection curve whose points are found starting from the surface
/// parameters of its polyline leader's vertices, instead of searching both
/// surfaces from scratch for every point.
///
/// Approximating an intersection curve by a B-spline samples it many times,
/// and each unhinted sample searches a grid over both surfaces, which made
/// booleans with curved faces take seconds.
#[derive(Clone, Debug)]
pub(super) struct HintedIntersectionCurve<'a, S> {
    curve: &'a IntersectionCurve<PolylineCurve<Point3>, S, S>,
    /// Parameters on (surface0, surface1) of each leader vertex.
    hints: Vec<(Point2, Point2)>,
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
        Some(Self {
            curve,
            hints,
            bound,
        })
    }
}

impl<S> ParametricCurve for HintedIntersectionCurve<'_, S>
where S: ParametricSurface3D + SearchNearestParameter<D2, Point = Point3>
{
    type Point = Point3;
    type Vector = Vector3;
    fn subs(&self, t: f64) -> Point3 {
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
    fn der(&self, t: f64) -> Vector3 { self.curve.der(t) }
    fn der2(&self, t: f64) -> Vector3 { self.curve.der2(t) }
    fn der_n(&self, n: usize, t: f64) -> Vector3 { self.curve.der_n(n, t) }
    fn parameter_range(&self) -> ParameterRange { self.curve.parameter_range() }
}
