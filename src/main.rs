//! The kentron adapter: answers the B-rep protocol's requests on stdin and stdout, so the
//! conformance suite, or any other client, can drive kentron.
//!
//! Every function is answered by kentron's one call for it, or refused as unsupported with the
//! reason in `supports`. What the adapter itself does is translation: the suite's profiles,
//! frames and paths as kentron's curves, faces, frames and spines, a point as the edge or face it
//! lies on through kentron's distance query at `10·ε`.

use std::f64::consts::TAU;
use std::process::ExitCode;

use brep_protocol::{Kernel, Refusal, model, model::Primitive, profile};
use kentron::geometry::curve::{BSplineCurve, Circle, Curve, TrimmedCurve};
use kentron::geometry::surface::{Plane, Surface};
use kentron::math::{
    AffineTransform, Axis, Direction, Frame, Matrix3, Point, Transform, Transformable, Vector,
};
use kentron::topology::{Model, NodeKind, Shape};
use kentron::{
    blend, boolean, defeature, distance, draft, io, loft, make, meshing, modify, offset,
    hidden_lines, primitives, properties, raycast, sew, sweep, tolerance, validity,
};

/// The kentron release this adapter builds against, as Cargo.toml pins it (kentron's crate
/// exports no version).
const VERSION: &str = "0.21.0";

/// The contract's `10·ε` within which a point names an edge or a face.
const NAMING: f64 = 10.0 * 1e-7;

/// The functions kentron answers in one call. The rest of the protocol's, with why not:
///
/// - `betti`, `section_area`: kentron reaches each only through other calls (counts of shapes;
///   a section whose curves bound an area).
const SUPPORTED: &[&str] = &[
    "box",
    "slab",
    "cylinder",
    "sphere",
    "torus",
    "cone",
    "extrude",
    "revolve",
    "screw",
    "loft",
    "loft_through",
    "sweep",
    "sweep_path",
    "union",
    "intersection",
    "difference",
    "union_all",
    "difference_all",
    "fillet",
    "fillet_radii",
    "fillet_variable",
    "fillet_faces",
    "full_round",
    "fillet_setback",
    "fillet_setback_g2",
    "fillet_conic",
    "fillet_g2",
    "project",
    "chamfer",
    "chamfer_two",
    "chamfer_distances",
    "offset",
    "offset_faces",
    "delete_faces",
    "move_faces",
    "thicken",
    "shell",
    "draft",
    "write_step",
    "sew",
    "read_step",
    "volume",
    "area",
    "classify",
    "normal",
    "curvature",
    "closest",
    "check",
    "centroid",
    "inertia",
    "distance",
    "interference",
    "section",
    "bounding_box",
    "split",
    "mesh",
    "display",
    "raycast",
    "transform",
    "affine",
];

/// kentron's session: one model holding every solid the adapter has handed out.
struct Kentron {
    model: Model,
}

/// kentron's error as the protocol's refusal, its kind kentron's own.
fn refusal(e: kentron::Error) -> Refusal {
    match e {
        kentron::Error::InvalidInput(m) => Refusal::invalid(m),
        kentron::Error::Failed(m) => Refusal::failed(m),
    }
}

fn point(a: model::Vec3) -> Point {
    Point::new(a[0], a[1], a[2])
}

fn vector(a: model::Vec3) -> Vector {
    Vector::new(a[0], a[1], a[2])
}

fn direction(a: model::Vec3) -> Result<Direction, Refusal> {
    Direction::try_new(a[0], a[1], a[2])
        .ok_or_else(|| Refusal::invalid(format!("no direction {a:?}")))
}

fn xyz(p: Point) -> model::Vec3 {
    [p.x, p.y, p.z]
}

/// The suite's frame as kentron's: its z, and its x as the hint, so that kentron's y is the
/// frame's.
fn frame(f: &model::Frame) -> Result<Frame, Refusal> {
    Ok(Frame::new(
        point(f.origin),
        direction(f.z)?,
        direction(f.x)?,
    ))
}

/// A plane a profile is drawn in: `(u, v)` is `o + u·x + v·y`.
struct Sketch {
    o: model::Vec3,
    x: model::Vec3,
    y: model::Vec3,
}

impl Sketch {
    fn at(&self, uv: [f64; 2]) -> Point {
        let [u, v] = uv;
        Point::new(
            self.o[0] + u * self.x[0] + v * self.y[0],
            self.o[1] + u * self.x[1] + v * self.y[1],
            self.o[2] + u * self.x[2] + v * self.y[2],
        )
    }

    /// The plane's normal, `x × y`.
    fn normal(&self) -> Result<Direction, Refusal> {
        let (x, y) = (self.x, self.y);
        direction([
            x[1] * y[2] - x[2] * y[1],
            x[2] * y[0] - x[0] * y[2],
            x[0] * y[1] - x[1] * y[0],
        ])
    }

    /// kentron's frame of the plane: its normal as z and x as x, so that its y is the sketch's.
    fn frame_at(&self, o: Point) -> Result<Frame, Refusal> {
        Ok(Frame::new(o, self.normal()?, direction(self.x)?))
    }
}

/// A circle of the sketch about `center`: counter-clockwise in the sketch for `ccw`, its
/// parameter the sketch's angle (clockwise, the angle's negative).
fn circle(sketch: &Sketch, center: [f64; 2], radius: f64, ccw: bool) -> Result<Circle, Refusal> {
    let normal = sketch.normal()?;
    let z = if ccw { normal } else { -normal };
    Ok(Circle::new(
        Frame::new(sketch.at(center), z, direction(sketch.x)?),
        radius,
    ))
}

/// The periodic B-spline of one period, as kentron takes it: the protocol's first `n` of its
/// `n + degree` points, and the knots of that period with their multiplicities.
fn periodic(sp: &model::Spline, sketch: &Sketch) -> Result<Curve, Refusal> {
    let degree = sp.degree;
    let n = sp
        .ctrl
        .len()
        .checked_sub(degree)
        .filter(|&n| n > 0)
        .ok_or_else(|| Refusal::invalid("a periodic spline needs more points than its degree"))?;
    let poles: Vec<Point> = sp.ctrl[..n].iter().map(|&c| sketch.at(c)).collect();
    let period = sp
        .knots
        .get(degree..degree + n + 1)
        .ok_or_else(|| Refusal::invalid("a periodic spline with too few knots"))?;
    let (knots, mults) = distinct(period);
    let curve = BSplineCurve::try_new_periodic(poles, knots, mults, degree).map_err(refusal)?;
    Ok(Curve::BSpline(curve))
}

/// Knot values with their multiplicities, from a list that repeats a knot once per
/// multiplicity.
fn distinct(flat: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let (mut knots, mut mults): (Vec<f64>, Vec<usize>) = (Vec::new(), Vec::new());
    for &k in flat {
        if knots.last() == Some(&k) {
            *mults.last_mut().expect("a knot before") += 1;
        } else {
            knots.push(k);
            mults.push(1);
        }
    }
    (knots, mults)
}

impl Kentron {
    /// A loop of the profile as kentron's wire, in the sketch's plane.
    fn wire(&mut self, segments: &[profile::Seg], sketch: &Sketch) -> Result<Shape, Refusal> {
        let m = &mut self.model;
        let mut edges = Vec::with_capacity(segments.len());
        for s in segments {
            let edge = match s {
                profile::Seg::Line { a, b } => make::edge_between(m, sketch.at(*a), sketch.at(*b)),
                profile::Seg::Arc {
                    center,
                    radius,
                    start,
                    sweep,
                } => {
                    let ccw = *sweep > 0.0;
                    let c = Curve::Circle(circle(sketch, *center, *radius, ccw)?);
                    // On the clockwise circle, the sketch's angle `t` is the parameter `-t`.
                    let first = if ccw { *start } else { -*start };
                    if sweep.abs() >= TAU {
                        make::edge_from_curve_range(m, c, first, first)
                    } else {
                        make::edge_from_curve_range(m, c, first, first + sweep.abs())
                    }
                }
                profile::Seg::Conic { a, apex, b, w } => {
                    let poles = vec![sketch.at(*a), sketch.at(*apex), sketch.at(*b)];
                    let bezier = BSplineCurve::try_new(poles, vec![0.0, 1.0], vec![3, 3], 2)
                        .map_err(refusal)?;
                    if w.is_nan() || *w <= 0.0 {
                        return Err(Refusal::invalid(format!("a conic of weight {w}")));
                    }
                    make::edge_from_curve(
                        m,
                        Curve::BSpline(bezier.with_weights(vec![1.0, *w, 1.0])),
                    )
                }
                profile::Seg::Spline(sp) => make::edge_from_curve(m, periodic(sp, sketch)?),
                profile::Seg::Gn { .. } => {
                    return Err(Refusal::unsupported(
                        "a blend-section (gn) segment: kentron has no such curve",
                    ));
                }
                _ => return Err(Refusal::unsupported("a segment of an unknown kind")),
            };
            edges.push(edge.map_err(refusal)?);
        }
        make::wire(m, &edges).map_err(refusal)
    }

    /// The profile's face in the sketch's plane: the outer loop, then the holes.
    fn face(&mut self, p: &profile::Profile, sketch: &Sketch) -> Result<Shape, Refusal> {
        if p.loops.is_empty() {
            return Err(Refusal::invalid("a profile with no loop"));
        }
        let wires = p
            .loops
            .iter()
            .map(|l| self.wire(l, sketch))
            .collect::<Result<Vec<_>, _>>()?;
        let plane = Surface::Plane(Plane::from_frame(sketch.frame_at(point(sketch.o))?));
        make::face_bounded_on(&mut self.model, plane, &wires).map_err(refusal)
    }

    /// The profile's one loop as a wire, for a call that takes sections: holes are refused, not
    /// dropped.
    fn section_wire(&mut self, p: &profile::Profile, sketch: &Sketch) -> Result<Shape, Refusal> {
        match &p.loops[..] {
            [one] => self.wire(one, sketch),
            [] => Err(Refusal::invalid("a profile with no loop")),
            _ => Err(Refusal::unsupported(
                "a section with holes: kentron's loft takes wires",
            )),
        }
    }

    /// The one edge, or face, of `s` within `10·ε` of each point, by kentron's distance from the
    /// point to each. A point near several is not settled here: unsupported. A seam, where a face
    /// meets itself, is no edge of the protocol's (it lists seams apart from edges), so it names
    /// none.
    fn entities(
        &mut self,
        s: Shape,
        points: &[model::Vec3],
        kind: NodeKind,
    ) -> Result<Vec<Shape>, Refusal> {
        let what = if kind == NodeKind::Edge {
            "edge"
        } else {
            "face"
        };
        let mut seams = std::collections::HashSet::new();
        if kind == NodeKind::Edge {
            for face in self.model.unique_subshapes(s, NodeKind::Face) {
                let Ok(face_id) = kentron::topology::FaceId::try_from(face.node) else { continue };
                for edge in self.model.unique_subshapes(face, NodeKind::Edge) {
                    let Ok(edge_id) = kentron::topology::EdgeId::try_from(edge.node) else { continue };
                    if self.model.is_seam_on_face(edge_id, face_id) {
                        seams.insert(edge_id);
                    }
                }
            }
        }
        let candidates: Vec<Shape> = self
            .model
            .unique_subshapes(s, kind)
            .into_iter()
            .filter(|&e| kind != NodeKind::Edge || !self.model.is_degenerate(e))
            .filter(|&e| kentron::topology::EdgeId::try_from(e.node).is_ok_and(|id| !seams.contains(&id)) || kind != NodeKind::Edge)
            .collect();
        let mut found = Vec::with_capacity(points.len());
        for &p in points {
            let vertex = make::vertex(&mut self.model, point(p));
            let mut near = Vec::new();
            for &c in &candidates {
                let d = distance::distance(&mut self.model, vertex, c).map_err(refusal)?;
                if d.value <= NAMING {
                    near.push(c);
                }
            }
            match near[..] {
                [one] => found.push(one),
                [] => return Err(Refusal::invalid(format!("no {what} passes through {p:?}"))),
                _ => {
                    return Err(Refusal::unsupported(format!(
                        "{p:?} lies on {} {what}s: kentron takes one",
                        near.len()
                    )));
                }
            }
        }
        Ok(found)
    }

    /// The empty solid, as kentron's booleans test for it (`is_empty` in `boolean/mod.rs`, which
    /// is private): a shape with no face, edge or vertex.
    fn is_empty(&self, s: Shape) -> bool {
        [NodeKind::Face, NodeKind::Edge, NodeKind::Vertex].iter().all(|&k| self.model.explore(s, k).next().is_none())
    }

    /// An unbounded planar face: kentron's plane through `origin` square to `normal`.
    fn plane(&mut self, origin: model::Vec3, normal: model::Vec3) -> Result<Shape, Refusal> {
        let face = self.model.new_empty(NodeKind::Face);
        let plane = Plane::new(point(origin), direction(normal)?);
        self.model
            .set_face_surface(face, Surface::Plane(plane), tolerance::LINEAR)
            .map_err(refusal)?;
        Ok(face)
    }
}

impl Kernel for Kentron {
    type Solid = Shape;

    fn name(&self) -> String {
        "kentron".into()
    }

    fn version(&self) -> String {
        VERSION.into()
    }

    fn supports(&self) -> Vec<String> {
        SUPPORTED.iter().map(|s| s.to_string()).collect()
    }

    fn tolerance(&self) -> f64 {
        tolerance::LINEAR
    }

    fn make(&mut self, p: &Primitive) -> Result<Shape, Refusal> {
        let f = p.frame();
        let e = refusal;
        let m = &mut self.model;
        match p {
            Primitive::Box { size, .. } => primitives::Cuboid::new(size[0], size[1], size[2])
                .with_frame(frame(&f)?)
                .build(m)
                .map_err(e),
            Primitive::Cylinder { radius, height, .. } => {
                primitives::Cylinder::new(*radius, *height)
                    .with_frame(frame(&f)?)
                    .build(m)
                    .map_err(e)
            }
            Primitive::Sphere { radius, .. } => primitives::Sphere::new(*radius)
                .with_frame(frame(&f)?)
                .build(m)
                .map_err(e),
            Primitive::Torus { major, minor, .. } => primitives::Torus::new(*major, *minor)
                .with_frame(frame(&f)?)
                .build(m)
                .map_err(e),
            Primitive::Cone {
                bottom,
                top,
                height,
                ..
            } => primitives::Cone::new(*bottom, *top, *height)
                .with_frame(frame(&f)?)
                .build(m)
                .map_err(e),
            Primitive::Extrude {
                profile, height, ..
            } => {
                let face = self.face(
                    profile,
                    &Sketch {
                        o: f.origin,
                        x: f.x,
                        y: f.y,
                    },
                )?;
                let v = vector(f.z) * *height;
                Ok(sweep::extrude(&mut self.model, face, v).map_err(e)?.shape)
            }
            Primitive::Revolve { profile, angle, .. } => {
                let face = self.face(
                    profile,
                    &Sketch {
                        o: f.origin,
                        x: f.x,
                        y: f.z,
                    },
                )?;
                let axis = Axis::new(point(f.origin), direction(f.z)?);
                Ok(sweep::revolve(&mut self.model, face, axis, *angle)
                    .map_err(e)?
                    .shape)
            }
            Primitive::Screw {
                profile,
                rise,
                angle,
                ..
            } => {
                // A turn rises 2π·rise, and the whole screw rise·angle.
                let face = self.face(
                    profile,
                    &Sketch {
                        o: f.origin,
                        x: f.x,
                        y: f.z,
                    },
                )?;
                let axis = Axis::new(point(f.origin), direction(f.z)?);
                Ok(
                    sweep::screw_extrude(&mut self.model, face, axis, TAU * *rise, *rise * *angle)
                        .map_err(e)?
                        .shape,
                )
            }
            Primitive::Sweep { profile, path, .. } => {
                let face = self.face(
                    profile,
                    &Sketch {
                        o: f.origin,
                        x: f.x,
                        y: f.y,
                    },
                )?;
                let spine = self.spine(&f, path)?;
                Ok(sweep::pipe(&mut self.model, face, spine).map_err(e)?.shape)
            }
            Primitive::SweepPath {
                profile,
                path,
                twist,
                ..
            } => {
                let face = self.face(
                    profile,
                    &Sketch {
                        o: f.origin,
                        x: f.x,
                        y: f.y,
                    },
                )?;
                let poles: Vec<Point> = path
                    .ctrl
                    .iter()
                    .map(|c| {
                        Sketch {
                            o: f.origin,
                            x: f.x,
                            y: f.y,
                        }
                        .at([c[0], c[1]])
                            + vector(f.z) * c[2]
                    })
                    .collect();
                let (knots, mults) = distinct(&path.knots);
                let curve = BSplineCurve::try_new(poles, knots, mults, path.degree).map_err(e)?;
                let edge =
                    make::edge_from_curve(&mut self.model, Curve::BSpline(curve)).map_err(e)?;
                let spine = make::wire(&mut self.model, &[edge]).map_err(e)?;
                // The pipe's frame carried by the least rotations from tangent to tangent: the
                // rotation-minimizing frame the function asks for.
                let options = sweep::PipeOptions {
                    guide: Some(sweep::Guide::Discrete),
                    twist: *twist,
                    ..Default::default()
                };
                Ok(sweep::pipe_with(&mut self.model, face, spine, &options)
                    .map_err(e)?
                    .shape)
            }
            Primitive::Loft {
                profile,
                height,
                top,
                ..
            } => {
                let base = self.section_wire(
                    profile,
                    &Sketch {
                        o: f.origin,
                        x: f.x,
                        y: f.y,
                    },
                )?;
                let up = xyz(point(f.origin) + vector(f.z) * *height);
                let top = self.section_wire(
                    top,
                    &Sketch {
                        o: up,
                        x: f.x,
                        y: f.y,
                    },
                )?;
                // Each point of the profile joined to its image in the scaled copy: the two
                // sections' edges matched as they are given.
                let options = loft::LoftOptions {
                    solid: true,
                    ruled: true,
                    matching: loft::Matching::AsGiven,
                    ..Default::default()
                };
                Ok(loft::loft(&mut self.model, &[base, top], &options)
                    .map_err(e)?
                    .shape)
            }
            Primitive::Slab {
                size,
                knots_u,
                knots_v,
                heights,
                ..
            } => primitives::Slab::new(size[0], size[1], knots_u.clone(), knots_v.clone(), heights.clone())
                .with_frame(frame(&f)?)
                .build(m)
                .map_err(e),
        }
    }

    fn union(&mut self, a: &Shape, b: &Shape) -> Result<Shape, Refusal> {
        boolean::union(&mut self.model, *a, *b).map_err(refusal)
    }

    fn intersection(&mut self, a: &Shape, b: &Shape) -> Result<Shape, Refusal> {
        boolean::intersection(&mut self.model, *a, *b).map_err(refusal)
    }

    fn difference(&mut self, a: &Shape, b: &Shape) -> Result<Shape, Refusal> {
        boolean::difference(&mut self.model, *a, *b).map_err(refusal)
    }

    fn union_all(&mut self, solids: &[&Shape]) -> Result<Shape, Refusal> {
        let Some((first, rest)) = solids.split_first() else {
            return Err(Refusal::invalid("union_all of no solids"));
        };
        let tools: Vec<Shape> = rest.iter().map(|&&s| s).collect();
        boolean::combine(
            &mut self.model,
            &[**first],
            &tools,
            boolean::Combination::Union,
        )
        .map_err(refusal)
    }

    fn difference_all(&mut self, target: &Shape, tools: &[&Shape]) -> Result<Shape, Refusal> {
        let tools: Vec<Shape> = tools.iter().map(|&&s| s).collect();
        boolean::combine(
            &mut self.model,
            &[*target],
            &tools,
            boolean::Combination::Difference,
        )
        .map_err(refusal)
    }

    fn volume(&mut self, s: &Shape) -> Result<f64, Refusal> {
        Ok(properties::volume(&self.model, *s).mass)
    }

    fn area(&mut self, s: &Shape) -> Result<f64, Refusal> {
        Ok(properties::area(&self.model, *s).mass)
    }

    fn check(&mut self, s: &Shape) -> Result<Vec<String>, Refusal> {
        let report = validity::check(&self.model, *s);
        Ok(report
            .errors()
            .map(|p| format!("{:?}: {}", p.kind, p.description))
            .collect())
    }

    fn classify(&mut self, s: &Shape, points: &[model::Vec3]) -> Result<Vec<model::Side>, Refusal> {
        // An empty solid holds no point.
        if self.is_empty(*s) {
            return Ok(vec![model::Side::Out; points.len()]);
        }
        let classifier = boolean::PointClassifier::new(&self.model, *s).map_err(refusal)?;
        Ok(points
            .iter()
            .map(|&p| match classifier.locate(point(p)) {
                boolean::PointLocation::Inside => model::Side::In,
                boolean::PointLocation::Boundary => model::Side::On,
                boolean::PointLocation::Outside => model::Side::Out,
            })
            .collect())
    }

    fn closest(&mut self, s: &Shape, points: &[model::Vec3]) -> Result<Vec<model::Vec3>, Refusal> {
        points
            .iter()
            .map(|&p| {
                let vertex = make::vertex(&mut self.model, point(p));
                let d = distance::distance(&mut self.model, vertex, *s).map_err(refusal)?;
                Ok(xyz(d.on_second))
            })
            .collect()
    }

    fn centroid(&mut self, s: &Shape) -> Result<Option<model::Vec3>, Refusal> {
        if self.is_empty(*s) {
            return Ok(None);
        }
        Ok(Some(
            xyz(properties::volume(&self.model, *s).center_of_mass),
        ))
    }

    fn inertia(&mut self, s: &Shape) -> Result<[[f64; 3]; 3], Refusal> {
        if self.is_empty(*s) {
            return Ok([[0.0; 3]; 3]);
        }
        let m = properties::volume(&self.model, *s).inertia;
        Ok(std::array::from_fn(|i| std::array::from_fn(|j| m[(i, j)])))
    }

    fn distance(&mut self, a: &Shape, b: &Shape) -> Result<f64, Refusal> {
        Ok(distance::distance(&mut self.model, *a, *b)
            .map_err(refusal)?
            .value)
    }

    fn bounding_box(&mut self, s: &Shape) -> Result<[model::Vec3; 2], Refusal> {
        let b = properties::tight_bounding_box(&self.model, *s);
        Ok([xyz(b.min), xyz(b.max)])
    }

    fn fillet(&mut self, s: &Shape, edges: &[model::Vec3], radius: f64) -> Result<Shape, Refusal> {
        let edges = self.entities(*s, edges, NodeKind::Edge)?;
        blend::fillet(&mut self.model, *s, &edges, radius).map_err(refusal)
    }

    fn fillet_radii(
        &mut self,
        s: &Shape,
        edges: &[model::Vec3],
        radii: &[f64],
    ) -> Result<Shape, Refusal> {
        if edges.len() != radii.len() {
            return Err(Refusal::invalid(format!(
                "{} edges but {} radii",
                edges.len(),
                radii.len()
            )));
        }
        let edges = self.entities(*s, edges, NodeKind::Edge)?;
        let pairs: Vec<(Shape, f64)> = edges.into_iter().zip(radii.iter().copied()).collect();
        blend::fillet_edges(&mut self.model, *s, &pairs).map_err(refusal)
    }

    fn chamfer(
        &mut self,
        s: &Shape,
        edges: &[model::Vec3],
        distance: f64,
    ) -> Result<Shape, Refusal> {
        let edges = self.entities(*s, edges, NodeKind::Edge)?;
        blend::chamfer(&mut self.model, *s, &edges, distance).map_err(refusal)
    }

    fn chamfer_distances(
        &mut self,
        s: &Shape,
        edges: &[model::Vec3],
        distances: &[f64],
    ) -> Result<Shape, Refusal> {
        if edges.len() != distances.len() {
            return Err(Refusal::invalid(format!(
                "{} edges but {} distances",
                edges.len(),
                distances.len()
            )));
        }
        let edges = self.entities(*s, edges, NodeKind::Edge)?;
        let chamfered: Vec<blend::ChamferedEdge> = edges
            .into_iter()
            .zip(distances)
            .map(|(edge, &d)| blend::ChamferedEdge {
                edge,
                size: blend::ChamferSize::Distances(d, d),
                face: None,
            })
            .collect();
        blend::chamfer_edges(&mut self.model, *s, &chamfered).map_err(refusal)
    }

    fn chamfer_two(
        &mut self,
        s: &Shape,
        edge: model::Vec3,
        face: model::Vec3,
        d1: f64,
        d2: f64,
    ) -> Result<Shape, Refusal> {
        let edge = self.entities(*s, &[edge], NodeKind::Edge)?[0];
        let face = self.entities(*s, &[face], NodeKind::Face)?[0];
        let chamfered = blend::ChamferedEdge {
            edge,
            size: blend::ChamferSize::Distances(d1, d2),
            face: Some(face),
        };
        blend::chamfer_edges(&mut self.model, *s, &[chamfered]).map_err(refusal)
    }

    fn offset(&mut self, s: &Shape, distance: f64) -> Result<Shape, Refusal> {
        offset::shape(&mut self.model, *s, distance, offset::Join::Arc).map_err(refusal)
    }

    fn offset_faces(
        &mut self,
        s: &Shape,
        faces: &[model::Vec3],
        distance: f64,
    ) -> Result<Shape, Refusal> {
        let faces = self.entities(*s, faces, NodeKind::Face)?;
        let moved: Vec<(Shape, f64)> = faces.into_iter().map(|f| (f, distance)).collect();
        offset::shape_with_faces(&mut self.model, *s, 0.0, &moved, offset::Join::Intersection)
            .map_err(refusal)
    }

    fn delete_faces(&mut self, s: &Shape, faces: &[model::Vec3]) -> Result<Shape, Refusal> {
        let faces = self.entities(*s, faces, NodeKind::Face)?;
        Ok(defeature::remove_faces(&mut self.model, *s, &faces)
            .map_err(refusal)?
            .shape)
    }

    fn thicken(
        &mut self,
        s: &Shape,
        faces: &[model::Vec3],
        thickness: f64,
    ) -> Result<Shape, Refusal> {
        let faces = self.entities(*s, faces, NodeKind::Face)?;
        let [face] = faces[..] else {
            return Err(Refusal::unsupported(format!(
                "thickening {} faces in one call: kentron thickens one sheet",
                faces.len()
            )));
        };
        offset::thicken(&mut self.model, face, thickness, offset::Join::Arc).map_err(refusal)
    }

    fn shell(&mut self, s: &Shape, open: &[model::Vec3], thickness: f64) -> Result<Shape, Refusal> {
        let open = self.entities(*s, open, NodeKind::Face)?;
        // Walls inward, the moved faces extended or trimmed to meet: joined by intersection.
        offset::thick_solid(
            &mut self.model,
            *s,
            &open,
            -thickness,
            offset::Join::Intersection,
        )
        .map_err(refusal)
    }

    fn draft(
        &mut self,
        s: &Shape,
        faces: &[model::Vec3],
        origin: model::Vec3,
        pull: model::Vec3,
        angle: f64,
    ) -> Result<Shape, Refusal> {
        let faces = self.entities(*s, faces, NodeKind::Face)?;
        let pull = direction(pull)?;
        let neutral = Plane::new(point(origin), pull);
        let drafted: Vec<draft::DraftedFace> = faces
            .into_iter()
            .map(|face| draft::DraftedFace {
                face,
                angle,
                neutral,
            })
            .collect();
        draft::draft(&mut self.model, *s, pull, &drafted).map_err(refusal)
    }

    fn write_step(&mut self, s: &Shape) -> Result<String, Refusal> {
        let mut out = Vec::new();
        io::step::write(&self.model, *s, "kentron solid", None, &mut out).map_err(refusal)?;
        String::from_utf8(out)
            .map_err(|e| Refusal::failed(format!("the STEP text is not UTF-8: {e}")))
    }

    fn read_step(&mut self, text: &str) -> Result<Shape, Refusal> {
        io::step::read(&mut self.model, &mut text.as_bytes()).map_err(refusal)
    }

    fn sew(&mut self, text: &str, tolerance: f64) -> Result<Shape, Refusal> {
        let faces = io::step::read(&mut self.model, &mut text.as_bytes()).map_err(refusal)?;
        let options = sew::SewOptions {
            tolerance,
            ..Default::default()
        };
        let sewn = sew::sew(&mut self.model, &[faces], &options).map_err(refusal)?;
        // kentron's sewing gives shells; the closed one is the solid's boundary.
        let shells: Vec<Shape> = if sewn.kind() == NodeKind::Shell {
            vec![sewn]
        } else {
            self.model.explore(sewn, NodeKind::Shell).collect()
        };
        let [shell] = shells[..] else {
            return Err(Refusal::failed(format!(
                "sewing gave {} shells",
                shells.len()
            )));
        };
        let solid = self.model.new_empty(NodeKind::Solid);
        self.model.add(solid, shell).map_err(refusal)?;
        Ok(solid)
    }

    fn transform(&mut self, s: &Shape, m: [[f64; 4]; 3]) -> Result<Shape, Refusal> {
        let t = Transform::from_rows(m).ok_or_else(|| Refusal::invalid("a singular matrix"))?;
        modify::transform(&mut self.model, *s, &t, Default::default()).map_err(refusal)
    }

    fn affine(&mut self, s: &Shape, m: [[f64; 4]; 3]) -> Result<Shape, Refusal> {
        let matrix = Matrix3::from_rows(m.map(|r| [r[0], r[1], r[2]]));
        let t = AffineTransform::new(matrix, Vector::new(m[0][3], m[1][3], m[2][3]));
        modify::transform_affine(&mut self.model, *s, &t).map_err(refusal)
    }

    fn section(
        &mut self,
        s: &Shape,
        origin: model::Vec3,
        normal: model::Vec3,
        tolerance: f64,
    ) -> Result<Vec<Vec<model::Vec3>>, Refusal> {
        let plane = self.plane(origin, normal)?;
        let cut = boolean::section(&mut self.model, *s, plane).map_err(refusal)?;
        let mut pieces = Vec::new();
        for edge in self.model.unique_subshapes(cut, NodeKind::Edge) {
            let Some((curve, first, last)) = self.model.edge_curve(edge) else {
                continue;
            };
            let trimmed = Curve::Trimmed(TrimmedCurve::new(curve, first, last));
            let sampled = kentron::geometry::sample::quasi_uniform_deflection(&trimmed, tolerance)
                .ok_or_else(|| Refusal::failed("a section edge could not be sampled"))?;
            let mut points: Vec<model::Vec3> = sampled.points.into_iter().map(xyz).collect();
            if edge.orientation == kentron::topology::Orientation::Reversed {
                points.reverse();
            }
            pieces.push(points);
        }
        Ok(pieces)
    }

    fn mesh(&mut self, s: &Shape, tolerance: f64) -> Result<brep_protocol::Mesh, Refusal> {
        let options = meshing::MeshOptions {
            deflection: tolerance,
            ..Default::default()
        };
        let (points, triangles) =
            meshing::triangulate(&mut self.model, *s, options).map_err(refusal)?;
        Ok(brep_protocol::Mesh::new(
            points.into_iter().map(xyz).collect(),
            triangles,
        ))
    }

    fn loft_through(
        &mut self,
        f: &model::Frame,
        profiles: &[profile::Profile],
        heights: &[f64],
    ) -> Result<Shape, Refusal> {
        if profiles.len() != heights.len() {
            return Err(Refusal::invalid(format!(
                "{} profiles at {} heights",
                profiles.len(),
                heights.len()
            )));
        }
        let sections = profiles
            .iter()
            .zip(heights)
            .map(|(p, &h)| {
                let at = xyz(point(f.origin) + vector(f.z) * h);
                self.section_wire(
                    p,
                    &Sketch {
                        o: at,
                        x: f.x,
                        y: f.y,
                    },
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let options = loft::LoftOptions {
            solid: true,
            ruled: false,
            ..Default::default()
        };
        Ok(loft::loft(&mut self.model, &sections, &options)
            .map_err(refusal)?
            .shape)
    }

    fn normal(&mut self, s: &Shape, points: &[model::Vec3]) -> Result<Vec<Option<model::Vec3>>, Refusal> {
        let points: Vec<Point> = points.iter().map(|&p| point(p)).collect();
        let normals = properties::normals(&mut self.model, *s, &points).map_err(refusal)?;
        Ok(normals.into_iter().map(|n| n.map(|d| [d.x(), d.y(), d.z()])).collect())
    }

    fn curvature(&mut self, s: &Shape, points: &[model::Vec3]) -> Result<Vec<Option<[f64; 2]>>, Refusal> {
        let points: Vec<Point> = points.iter().map(|&p| point(p)).collect();
        let found = properties::curvatures(&mut self.model, *s, &points).map_err(refusal)?;
        Ok(found.into_iter().map(|k| k.map(|k| [k.max, k.min])).collect())
    }

    fn interference(&mut self, a: &Shape, b: &Shape) -> Result<brep_protocol::Meeting, Refusal> {
        let found = boolean::interference(&mut self.model, *a, *b).map_err(refusal)?;
        Ok(brep_protocol::Meeting { volume: found.volume, gap: found.gap })
    }

    fn fillet_variable(&mut self, s: &Shape, edge: model::Vec3, from: model::Vec3, r_from: f64, r_to: f64) -> Result<Shape, Refusal> {
        let edge = self.entities(*s, &[edge], NodeKind::Edge)?[0];
        // kentron's radii run from the edge's first vertex to its last: `r_from` at whichever
        // end is nearer `from`.
        let (first, last) = self.model.edge_vertices(edge);
        let (first, last) = (self.model.point(first), self.model.point(last));
        let (to_first, to_last) = (first.distance(point(from)), last.distance(point(from)));
        let radii = if to_first < to_last {
            (r_from, r_to)
        } else if to_last < to_first {
            (r_to, r_from)
        } else {
            return Err(Refusal::invalid("`from` is as near one end of the edge as the other"));
        };
        blend::fillet_linear(&mut self.model, *s, &[(edge, radii.0, radii.1)]).map_err(refusal)
    }

    fn fillet_faces(&mut self, s: &Shape, face_a: model::Vec3, face_b: model::Vec3, radius: f64) -> Result<Shape, Refusal> {
        let faces = self.entities(*s, &[face_a, face_b], NodeKind::Face)?;
        blend::fillet_faces(&mut self.model, *s, faces[0], faces[1], radius).map_err(refusal)
    }

    fn full_round(&mut self, s: &Shape, side_a: model::Vec3, middle: model::Vec3, side_b: model::Vec3) -> Result<Shape, Refusal> {
        let faces = self.entities(*s, &[side_a, middle, side_b], NodeKind::Face)?;
        blend::full_round(&mut self.model, *s, faces[0], faces[1], faces[2]).map_err(refusal)
    }

    fn fillet_setback(&mut self, s: &Shape, edges: &[model::Vec3], radius: f64, setback: f64) -> Result<Shape, Refusal> {
        let edges = self.entities(*s, edges, NodeKind::Edge)?;
        blend::fillet_setback(&mut self.model, *s, &edges, radius, setback).map_err(refusal)
    }

    fn fillet_setback_g2(&mut self, s: &Shape, edges: &[model::Vec3], radius: f64, setback: f64) -> Result<Shape, Refusal> {
        let edges = self.entities(*s, edges, NodeKind::Edge)?;
        blend::fillet_setback_g2(&mut self.model, *s, &edges, radius, setback).map_err(refusal)
    }

    fn fillet_conic(&mut self, s: &Shape, edges: &[model::Vec3], setback: f64, rho: f64) -> Result<Shape, Refusal> {
        let edges = self.entities(*s, edges, NodeKind::Edge)?;
        blend::fillet_conic(&mut self.model, *s, &edges, setback, rho).map_err(refusal)
    }

    fn fillet_g2(&mut self, s: &Shape, edges: &[model::Vec3], offset: f64, center: f64) -> Result<Shape, Refusal> {
        let edges = self.entities(*s, edges, NodeKind::Edge)?;
        blend::fillet_g2(&mut self.model, *s, &edges, offset, center).map_err(refusal)
    }

    fn project(&mut self, s: &Shape, d: model::Vec3, tolerance: f64) -> Result<brep_protocol::Projection, Refusal> {
        // The protocol's view axes as kentron's view frame: its x across, its z toward the
        // viewer, so that its y is the protocol's axis up.
        let toward = direction([-d[0], -d[1], -d[2]])?;
        let (u, _) = brep_protocol::view_axes(d);
        let view = Frame::new(Point::ORIGIN, toward, direction(u)?);
        let drawn = hidden_lines::hidden_lines(&self.model, *s, &view, tolerance).map_err(refusal)?;
        let flat = |lines: Vec<hidden_lines::DrawnLine>| lines.into_iter().map(|l| l.points.into_iter().map(|p| [p.x, p.y]).collect()).collect();
        Ok(brep_protocol::Projection::new(flat(drawn.visible), flat(drawn.hidden)))
    }

    fn split(&mut self, s: &Shape, origin: model::Vec3, normal: model::Vec3) -> Result<(Shape, Shape), Refusal> {
        let n = direction(normal)?;
        let plane = make::whole_face(&mut self.model, Surface::Plane(Plane::new(point(origin), n))).map_err(refusal)?;
        let cells = boolean::cells(&mut self.model, &[*s, plane]).map_err(refusal)?;
        // The solid's pieces, sorted by side: each lies wholly on one side of the plane, its
        // centroid with it.
        let (ahead, behind) = (self.model.new_empty(NodeKind::Compound), self.model.new_empty(NodeKind::Compound));
        for cell in cells.iter().filter(|c| c.inside[0] && c.shape.kind() == NodeKind::Solid) {
            let centre = properties::volume(&self.model, cell.shape).center_of_mass;
            let side = if (centre - point(origin)).dot(n.into()) > 0.0 { ahead } else { behind };
            self.model.add(side, cell.shape).map_err(refusal)?;
        }
        Ok((ahead, behind))
    }

    fn move_faces(&mut self, s: &Shape, faces: &[model::Vec3], matrix: [[f64; 4]; 3]) -> Result<Shape, Refusal> {
        let faces = self.entities(*s, faces, NodeKind::Face)?;
        let motion = Transform::from_rows(matrix).ok_or_else(|| Refusal::invalid("a singular matrix"))?;
        draft::move_faces(&mut self.model, *s, &faces, &motion).map_err(refusal)
    }

    fn display(&mut self, s: &Shape, tolerance: f64) -> Result<brep_protocol::Display, Refusal> {
        let d = meshing::display(&mut self.model, *s, tolerance).map_err(refusal)?;
        let dir = |n: Direction| [n.x(), n.y(), n.z()];
        let faces = d
            .faces
            .into_iter()
            .map(|f| brep_protocol::DisplayFace::new(f.positions.into_iter().map(xyz).collect(), f.normals.into_iter().map(dir).collect(), f.triangles))
            .collect();
        let edges = d
            .edges
            .into_iter()
            .map(|e| brep_protocol::DisplayEdge::new(e.points.into_iter().map(xyz).collect(), e.faces, e.turn, e.split))
            .collect();
        let seams = d.seams.into_iter().map(|p| p.into_iter().map(xyz).collect()).collect();
        Ok(brep_protocol::Display::new(faces, edges, seams))
    }

    fn raycast(&mut self, s: &Shape, origin: model::Vec3, direction: model::Vec3) -> Result<Vec<f64>, Refusal> {
        raycast::raycast(&self.model, *s, point(origin), vector(direction)).map_err(refusal)
    }

    fn raycasts(&mut self, s: &Shape, rays: &[(model::Vec3, model::Vec3)]) -> Result<Vec<Vec<f64>>, Refusal> {
        let rays: Vec<(Point, Vector)> = rays.iter().map(|&(o, d)| (point(o), vector(d))).collect();
        raycast::raycast_all(&self.model, *s, &rays).map_err(refusal)
    }

    fn reset(&mut self) {
        self.model = Model::new();
    }
}

impl Kentron {
    /// The sweep's path as kentron's spine: from the frame's origin along its z, each line straight
    /// on and each arc bending toward the current x about an axis along the frame's y, which the
    /// bends leave as it is.
    fn spine(&mut self, f: &model::Frame, path: &[model::PathSeg]) -> Result<Shape, Refusal> {
        let e = refusal;
        let y = direction(f.y)?;
        let (mut at, mut x, mut z) = (point(f.origin), vector(f.x), vector(f.z));
        let mut edges = Vec::with_capacity(path.len());
        for piece in path {
            match *piece {
                model::PathSeg::Line { length } => {
                    let end = at + z * length;
                    edges.push(make::edge_between(&mut self.model, at, end).map_err(e)?);
                    at = end;
                }
                model::PathSeg::Arc { radius, angle } => {
                    // About the axis through `at + radius·x` along y: from `at`, the circle's
                    // x is `−x`, and its parameter runs toward `z`.
                    let centre = at + x * radius;
                    let back = Direction::try_from_vector(-x)
                        .ok_or_else(|| Refusal::invalid("a bend with no direction"))?;
                    let c = Circle::new(Frame::new(centre, y, back), radius);
                    edges.push(
                        make::edge_from_curve_range(&mut self.model, Curve::Circle(c), 0.0, angle)
                            .map_err(e)?,
                    );
                    let turn = Transform::rotation(Axis::new(centre, y), angle);
                    at = at.transformed(&turn);
                    x = x.transformed(&turn);
                    z = z.transformed(&turn);
                }
            }
        }
        make::wire(&mut self.model, &edges).map_err(e)
    }
}

fn main() -> ExitCode {
    let mut kernel = Kentron {
        model: Model::new(),
    };
    match brep_protocol::serve(&mut kernel) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{e}");
            ExitCode::FAILURE
        }
    }
}
