use std::{path::PathBuf, str::FromStr};

use serde::Deserialize;

pub use super::compartment_combinations::Expression as CombinationExpression;
use crate::core::config::{Axes, CompartmentID, Dimensions, defaults};

impl<'de> Deserialize<'de> for CombinationExpression {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        CombinationExpression::from_str(&s).map_err(serde::de::Error::custom)
    }
}

/// Avogadro's number (per mol).
const N_A: f64 = 6.0221415e23;

// TODO: I think it's cursed that we store these defaults here. I'd like to create a file
// collecting all of these consts, one day.
fn bead_radius_default() -> f32 {
    defaults::BEAD_RADIUS as f32 // nm
}

fn max_tries_mult_default() -> u64 {
    defaults::MAX_TRIES_MULT
}

fn max_tries_rot_div_default() -> u64 {
    defaults::MAX_TRIES_ROT_DIV
}

#[derive(Deserialize)]
pub struct General {
    pub seed: Option<u64>,
    #[serde(default = "bead_radius_default")]
    pub bead_radius: f32,
    #[serde(default = "max_tries_mult_default")]
    pub max_tries_mult: u64,
    #[serde(default = "max_tries_rot_div_default")]
    pub max_tries_rot_div: u64,
}

impl Default for General {
    fn default() -> Self {
        Self {
            seed: Default::default(),
            bead_radius: bead_radius_default(),
            max_tries_mult: max_tries_mult_default(),
            max_tries_rot_div: max_tries_rot_div_default(),
        }
    }
}

#[derive(Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Shape {
    Spherical,
    Cuboid,
    None,
}

impl std::fmt::Display for Shape {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Shape::Spherical => "spherical",
            Shape::Cuboid => "cuboid",
            Shape::None => "empty ('none')",
        }
        .fmt(f)
    }
}

#[derive(Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Mask {
    Shape(Shape),
    Analytical { shape: Shape, center: Option<[f32; 3]>, radius: Option<f32> },
    Voxels { path: PathBuf },
    Combination(CombinationExpression),
}

#[derive(Deserialize)]
pub struct Compartment {
    pub id: CompartmentID,
    #[serde(flatten)]
    pub mask: Mask,
}

impl Compartment {
    pub fn is_predefined(&self) -> bool {
        match &self.mask {
            Mask::Shape(_) | Mask::Analytical { .. } | Mask::Voxels { .. } => true,
            Mask::Combination(_) => false,
        }
    }
}

pub(crate) fn true_by_default() -> bool {
    true
}

#[derive(Deserialize)]
pub struct Space {
    pub size: Dimensions,
    pub resolution: f32,
    pub compartments: Vec<Compartment>,
    #[serde(default = "true_by_default")]
    pub periodic: bool,
    // TODO: constraint system (satisfied _somewhat_ by the notion of a rule).
}

#[derive(Deserialize)]
#[serde(untagged)]
pub enum RuleExpression {
    Rule(String),
    Or(Vec<RuleExpression>),
}

fn parse_axes<'de, D>(deserializer: D) -> Result<Axes, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    s.parse().map_err(serde::de::Error::custom)
}

#[derive(Deserialize)]
pub struct Segment {
    pub name: String,
    pub tag: Option<String>,
    #[serde(flatten)]
    pub quantity: Quantity,
    pub path: PathBuf,
    pub compartments: Vec<CompartmentID>,
    #[serde(default)]
    pub rules: Vec<RuleExpression>,
    #[serde(default, deserialize_with = "parse_axes")]
    pub rotation_axes: Axes,
    #[serde(default)]
    pub initial_rotation: [f32; 3],
    // TODO: center?
}

#[derive(Deserialize, Clone, Copy)]
#[serde(rename_all = "lowercase")]
pub enum Quantity {
    Number(usize),
    /// Concentration in mol/L.
    Concentration(f64),
}

impl Quantity {
    /// Determine the number of segments that is implied by this [`Quantity`].
    ///
    /// In case this `Quantity` is a [`Quantity::Concentration`], the number of segments is
    /// lazily determined from the provided `volume`, and rounded.
    ///
    /// The value returned by `volume` must be in cubic nanometers (nm³).
    pub fn bake<F: Fn() -> f64>(&self, volume: F) -> usize {
        match *self {
            Quantity::Number(n) => n,
            Quantity::Concentration(c) => {
                // n = N_A * c * V
                let v = volume() * 1e-24; // From nm³ to L.
                let n = N_A * c * v;
                f64::round(n) as usize
            }
        }
    }

    /// Returns whether the contained value can be interpreted as resulting in zero placements.
    ///
    /// When the quantity is a `Number(0)` or `Concentration(0.0)`, the baked number is certainly
    /// zero. When `Number(n)` for `n > 0`, the baked number is certainly not zero.
    ///
    /// But, in case of a positive concentration, whether the final number is zero or not depends
    /// on the associated volume.
    ///
    /// If the concentration is smaller than zero, it is treated as a zero.
    pub fn is_zero(&self) -> bool {
        match *self {
            Quantity::Number(n) => n == 0,
            Quantity::Concentration(c) => c <= 0.0,
        }
    }
}

impl std::fmt::Display for Quantity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Quantity::Number(n) => write!(f, "{n} instances"),
            Quantity::Concentration(c) => write!(f, "{c} mol/L"),
        }
    }
}

pub type TopolIncludes = Vec<String>;

#[derive(Deserialize)]
pub struct Output {
    pub title: String,
    pub topol_includes: Option<TopolIncludes>,
}

#[derive(Deserialize)]
pub struct Config {
    #[serde(default)]
    pub general: General,
    pub space: Space,
    pub segments: Vec<Segment>,
    pub output: Output,
}

mod convert {
    use std::collections::{HashMap, HashSet};

    use super::*;
    use crate::core::config;
    use crate::core::config::legacy::convert::rule::parse_rule;

    impl From<CombinationExpression> for config::Expr<CompartmentID> {
        fn from(ce: CombinationExpression) -> Self {
            type Expr = config::Expr<CompartmentID>;
            fn unflatten(
                expressions: Vec<CombinationExpression>,
                binary: impl Fn(Box<Expr>, Box<Expr>) -> Expr,
            ) -> Expr {
                // We recursively convert the children first.
                let exprs: Vec<config::Expr<CompartmentID>> =
                    expressions.into_iter().map(|expression| expression.into()).collect();
                // Then, we stitch them together into a tree.
                let flat = exprs;
                flat.into_iter()
                    .reduce(|acc, item| binary(Box::new(acc), Box::new(item)))
                    .expect("TODO")
            }

            type CE = CombinationExpression;
            match ce {
                CE::Id(id) => config::Expr::Term(id),
                CE::Not(expression) => config::Expr::Not(Box::new((*expression).into())),
                CE::Union(expressions) => unflatten(expressions, config::Expr::Or),
                CE::Intersect(expressions) => unflatten(expressions, config::Expr::And),
            }
        }
    }

    impl From<Mask> for config::Mask {
        fn from(mask: Mask) -> Self {
            match mask {
                Mask::Shape(shape) => match shape {
                    // TODO: This sucks but is true.
                    Shape::Spherical => panic!("a sphere without a radius is an undefined shape"),
                    Shape::Cuboid | Shape::None => config::Mask::All,
                },
                Mask::Analytical { shape: Shape::Spherical, center, radius } => {
                    config::Mask::Shape(config::Shape::Sphere {
                        center: match center {
                            None => config::Anchor::Center,
                            Some(center) => config::Anchor::Point(center),
                        },
                        // TODO: This sucks but is true.
                        radius: radius.expect("a sphere without a radius is an undefined shape"),
                    })
                }
                Mask::Analytical {
                    shape: Shape::Cuboid | Shape::None,
                    // Cursed that these could be set, but we'll just keep quit about it.
                    center: _,
                    radius: _,
                } => config::Mask::All,
                Mask::Voxels { path } => config::Mask::Voxels(path),
                Mask::Combination(expression) => {
                    // TODO: This conversion sucks and will be changed.
                    config::Mask::Combination(expression.into())
                }
            }
        }
    }

    impl From<Compartment> for config::Compartment {
        fn from(compartment: Compartment) -> Self {
            let Compartment { id, mask } = compartment;
            config::Compartment { id, mask: mask.into() }
        }
    }

    impl From<Segment> for config::Segment {
        fn from(segment: Segment) -> Self {
            let Segment {
                name,
                tag,
                quantity,
                path,
                compartments,
                rules,
                rotation_axes,
                initial_rotation,
            } = segment;

            // As you can judge from the comments in the upcoming section, some non-trivial things
            // are happening here. For context, we need to do some juggling between the different
            // data types. This is a complex affair that is coordinated over the other conversion
            // functions as well. So, to understand this, make sure you also grasp how canonical
            // ids are used in the root From<Config> section as well.

            // We don't support this, anymore. It is still possible to rotate the contents of the
            // structure file, of course.
            if initial_rotation != <[f32; 3]>::default() {
                unimplemented!("segment initial rotation is deprecated")
            }

            // This is a bit complicated. We need to
            // (1) come up with a label for each unique RuleExpression in the file,
            // (2) those rule expressions need to be formulated as a constraint and
            //     pushed to the constraints list. This is done in the conversion of
            //     legacy::Config to config::Config.
            // Note that this is inefficient. That does not matter, since this just serves to
            // convert a legacy format to the new format. It lets us write these conversions in
            // such a way that the functions don't have to coordinate with each other which
            // makes them nice and seperate. That is preferred here.
            let compartment_ids_from_rules = rules
                .iter()
                .map(|re| {
                    let legacy_rule = parse_rule(re).expect("TODO");
                    rule::canonical_id_compartment(&legacy_rule)
                })
                // Only retain unique rules.
                .collect::<HashSet<_>>();

            // Same goes for the rotation_axes, if non-default. This notion is now conceived of as
            // a rule, instead of as a per-segment property. We also come up with an injective id,
            // here, and we assign it as the sole rule in the rules field.
            let axes_rule = if rotation_axes != Default::default() {
                let id = rule::canonical_id_rotation_axes(rotation_axes);
                vec![id].into_boxed_slice()
            } else {
                Default::default()
            };

            // The list of compartments represents a union of masks that can be merged together to
            // provide the accessible space for a segment. In the legacy format, where rules were
            // applied per segment, the presence of rules meant that the union of compartments was
            // intersected with the masks distilled from the rules. In other words, the rules were
            // _applied_ to the compartments.
            // Therefore, it would be incorrect to simply merge the compartment_ids_from_rules with
            // the provided compartment ids, since that would create a union of the explicit
            // compartment masks and the virtual masks distilled from the rules.
            // Instead, if any legacy rules are present, we need to create a special compartment
            // not only for the rule, but also for the union of compartment ids intersected with
            // the new legacy rule compartment.
            // From legacy::Segment { compartments: [a, b, c], rules: within 10 of d }, we go to
            //      config::Segment { compartment_ids: [{apply/{and/{or/a'b'c}'{win/d'10}}}] }.
            // In From<Config>, this {apply/...} rule is actually added to the compartments
            // section.
            let compartment_ids = if compartment_ids_from_rules.is_empty() {
                // No rules, so we can just take the union of compartments directly.
                compartments.into_boxed_slice()
            } else {
                let rule_ids = compartment_ids_from_rules.into_iter().collect::<Vec<_>>();
                let id = rule::canonical_id_apply(&compartments, &rule_ids);
                vec![id].into_boxed_slice()
            };

            // Pfhew... that sucked. But we're here now.
            config::Segment {
                name,
                tag,
                quantity: quantity.into(),
                path,
                compartment_ids,
                rules: axes_rule,
            }
        }
    }

    impl From<Quantity> for config::Quantity {
        fn from(quantity: Quantity) -> Self {
            match quantity {
                Quantity::Number(n) => config::Quantity::Number(n as u64),
                Quantity::Concentration(c) => config::Quantity::Concentration(c),
            }
        }
    }

    impl From<Config> for config::Config {
        fn from(config: Config) -> Self {
            let Config {
                general: General { seed, bead_radius, max_tries_mult, max_tries_rot_div },
                space: Space { size, resolution, compartments, periodic },
                segments,
                output: Output { title, topol_includes },
            } = config;

            // More cursed canonical rule id logic.
            // Together with the juggling in From<Segment>, the following constitutes a certified
            // 'tricky bit'. We make the assumption here that we created id-constraint pairs in an
            // injective manner. That is, one id has a single, unique rule associated with it, and
            // vice versa. The management of the rule keys throughout this conversion code aims to
            // uphold that.
            // TODO: This can all be made less.. unclear if we just move it all together into one
            //       happy monolithic family. Consider.

            // The legacy rules are now considered compartment declarations. We go through the
            // segments and collect all of their rules into a set and for each of those items we
            // come up with the canonical id the From<Segment> implementation also generated, along
            // with the compartment mask definition that fits the bill.
            // TODO: Ordering concerns from using HashMaps here?
            let compartments = {
                let mut legacy_rules = HashMap::new();
                let mut rule_applications = HashMap::new();
                for segment in &segments {
                    let mut segment_legacy_rules = HashMap::new();
                    for re in &segment.rules {
                        // First, we do a legacy parse of the rule.
                        let legacy_rule = parse_rule(re).expect("TODO");
                        // Come up with a unique id for this rule that will match how a segment
                        // converts the rule into the same id.
                        let id = rule::canonical_id_compartment(&legacy_rule);
                        segment_legacy_rules.insert(id, legacy_rule);
                    }

                    // Here comes another tricky bit.
                    if segment_legacy_rules.is_empty() {
                        // If there are no legacy rules, we can simply take the segment's list of
                        // compartments, like described in From<Segment>. That means that we don't
                        // have to do anything special here.
                        continue;
                    }

                    let segment_legacy_rule_ids =
                        segment_legacy_rules.keys().cloned().collect::<Vec<_>>();
                    let rule_application_id =
                        rule::canonical_id_apply(&segment.compartments, &segment_legacy_rule_ids);
                    let rule_application = (segment.compartments.clone(), segment_legacy_rule_ids);

                    legacy_rules.extend(segment_legacy_rules);
                    rule_applications.insert(rule_application_id, rule_application);
                }

                // Now we make actual compartments from the hashmaps.
                let legacy_rules = legacy_rules.into_iter().map(|(id, legacy_rule)| {
                    // And convert the legacy Rule into a Mask.
                    let mask = legacy_rule.into_mask();
                    config::Compartment { id, mask }
                });
                let rule_applications = rule_applications.into_iter().map(|(id, application)| {
                    // And convert the legacy rule applications into a Mask.
                    use config::Expr;
                    let (compartment_ids, legacy_rule_ids) = application;
                    let compartments = compartment_ids
                        .into_iter()
                        .map(Expr::Term)
                        .reduce(|acc, r| Expr::Or(Box::new(acc), Box::new(r)))
                        .expect("a rules combination expression cannot be empty");
                    let rules = legacy_rule_ids
                        .into_iter()
                        .map(Expr::Term)
                        .reduce(|acc, r| Expr::And(Box::new(acc), Box::new(r)))
                        .expect("a rules combination expression cannot be empty");
                    let mask = config::Mask::Combination(Expr::And(
                        Box::new(compartments),
                        Box::new(rules),
                    ));
                    config::Compartment { id, mask }
                });

                compartments
                    .into_iter()
                    .map(Into::into)
                    .chain(legacy_rules)
                    .chain(rule_applications)
                    .collect()
            };

            // Currently, constraints can only be rotation axes declarations. Look through all the
            // segments to see if there are non-default rotation axes set. For each of those,
            // create a constraint with its canonical id that will also be generated by
            // From<Segment> and the appropriate Rule. We collect them in a set to automatically
            // retain only unique constraints.
            let constraints = {
                let mut constraints = HashSet::new();
                for segment in &segments {
                    let rotation_axes = segment.rotation_axes;
                    if rotation_axes != Default::default() {
                        let id = rule::canonical_id_rotation_axes(rotation_axes);
                        let constraint = config::Constraint {
                            id,
                            rule: config::Rule::RotationAxes(rotation_axes),
                        };
                        constraints.insert(constraint);
                    }
                }
                constraints.into_iter().collect()
            };

            // Run through some defaults for the general and space sections.
            let bead_radius = if bead_radius != defaults::BEAD_RADIUS as f32 {
                Some(bead_radius as f64)
            } else {
                None
            };
            let max_tries_mult = if max_tries_mult != defaults::MAX_TRIES_MULT {
                Some(max_tries_mult)
            } else {
                None
            };
            let max_tries_rot_div = if max_tries_rot_div != defaults::MAX_TRIES_ROT_DIV {
                Some(max_tries_rot_div)
            } else {
                None
            };
            // Bit silly, but let's match this already ugly structure for clarity.
            let periodic = if periodic != defaults::PERIODIC { Some(periodic) } else { None };

            config::Config {
                general: config::General {
                    title: if title.is_empty() { None } else { Some(title) },
                    seed,
                    bead_radius,
                    max_tries_mult,
                    max_tries_rot_div,
                    rearrange_method: None,
                },
                space: config::Space {
                    dimensions: Some(size),
                    resolution: Some(resolution as f64),
                    periodic,
                },
                includes: topol_includes.unwrap_or_default().into_iter().map(Into::into).collect(),
                constraints,
                compartments,
                segments: segments.into_iter().map(Into::into).collect(),
            }
        }
    }

    // We're vendoring the Rule stuff until that can all be refactored out.
    mod rule {
        use std::num::ParseFloatError;
        use std::str::FromStr;

        use crate::core::config::{self, Axis, Limit, Op};

        use super::{CompartmentID, RuleExpression};

        // TODO: This should be part of the config parsing.
        pub fn parse_rule(expr: &RuleExpression) -> Result<Rule, ParseRuleError> {
            match expr {
                RuleExpression::Rule(s) => Rule::from_str(s),
                RuleExpression::Or(exprs) => {
                    Ok(Rule::Or(exprs.iter().map(parse_rule).collect::<Result<_, _>>()?))
                }
            }
        }

        pub fn canonical_id_compartment(rule: &Rule) -> String {
            match rule {
                Rule::Position(Limit { axis, op, value }) => {
                    let op = match op {
                        Op::LessThan => "lt",
                        Op::GreaterThan => "gt",
                    };
                    format!("{{lim/{axis}'{op}'{value}}}")
                }
                Rule::IsCloser(id, distance) => {
                    format!("{{win/{id}'{distance}}}")
                }
                Rule::Or(rules) => {
                    let ids =
                        rules.iter().map(canonical_id_compartment).collect::<Vec<_>>().join("'");
                    format!("{{or/{ids}}}")
                }
            }
        }

        pub fn canonical_id_rotation_axes(axes: config::Axes) -> String {
            let axes = axes.list().iter().map(ToString::to_string).collect::<String>();
            format!("{{axes/{axes}}}")
        }

        pub fn canonical_id_apply(
            compartment_ids: &[CompartmentID],
            rule_ids: &[CompartmentID],
        ) -> String {
            assert!(!compartment_ids.is_empty(), "expected at least one compartment");
            assert!(!rule_ids.is_empty(), "expected at least one rule");

            let rids = rule_ids.join("'");
            let cids = compartment_ids.join("'");
            format!("{{apply/{rids}/to/{cids}}}")
        }

        #[derive(Debug, Clone, PartialEq)]
        pub enum Rule {
            Position(Limit),
            IsCloser(CompartmentID, f32),

            /// A set of rules where any of them can be true for this [`Rule`] to apply.
            Or(Vec<Rule>),
        }

        impl Rule {
            pub(crate) fn into_mask(self) -> config::Mask {
                use config::Expr;
                match self {
                    Rule::Position(limit) => config::Mask::Limits(config::Expr::Term(limit)),
                    Rule::IsCloser(id, distance) => config::Mask::Within { distance, id },
                    Rule::Or(rules) => config::Mask::Combination(
                        rules
                            .into_iter()
                            .map(|rule| Expr::Term(canonical_id_compartment(&rule)))
                            .reduce(|acc, r| Expr::Or(Box::new(acc), Box::new(r)))
                            .expect("a rules combination expression cannot be empty"),
                    ),
                }
            }
        }

        impl FromStr for Rule {
            type Err = ParseRuleError;

            fn from_str(s: &str) -> Result<Self, Self::Err> {
                let trimmed = s.trim();
                let mut words = trimmed.split_whitespace();
                let keyword = words.next().ok_or(ParseRuleError::Empty)?;
                match keyword {
                    kind @ ("less_than" | "greater_than") => {
                        let axis = words
                            .next()
                            .ok_or(ParseRuleError::SyntaxError("expected axis".to_string()))?
                            .parse()
                            .map_err(ParseRuleError::ParseAxisError)?;
                        let value = words
                            .next()
                            .ok_or(ParseRuleError::SyntaxError(
                                "expected scalar value".to_string(),
                            ))?
                            .parse()
                            .map_err(ParseRuleError::ParseScalarError)?;

                        let poscon = match kind {
                            "greater_than" => Limit { axis, op: Op::GreaterThan, value },
                            "less_than" => Limit { axis, op: Op::LessThan, value },
                            _ => unreachable!(), // By virtue of this branch's pattern.
                        };
                        Ok(Rule::Position(poscon))
                    }
                    "is_closer_to" => {
                        let compartment_id = words.next().ok_or(ParseRuleError::SyntaxError(
                            "expected compartment id".to_string(),
                        ))?;
                        let distance = words
                            .next()
                            .ok_or(ParseRuleError::SyntaxError(
                                "expected scalar value".to_string(),
                            ))?
                            .parse()
                            .map_err(ParseRuleError::ParseScalarError)?;

                        Ok(Rule::IsCloser(compartment_id.to_string(), distance))
                    }
                    unknown => Err(ParseRuleError::UnknownKeyword(unknown.to_string())),
                }
            }
        }

        #[derive(Debug, Clone)]
        pub enum ParseRuleError {
            Empty,
            UnknownKeyword(String),
            SyntaxError(String),
            ParseScalarError(ParseFloatError),
            ParseAxisError(String),
        }

        impl std::fmt::Display for ParseRuleError {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                match self {
                    ParseRuleError::Empty => write!(f, "no rule keyword was provided"),
                    ParseRuleError::UnknownKeyword(unknown) => {
                        write!(f, "encountered an unknown keyword: {unknown:?}")
                    }
                    ParseRuleError::SyntaxError(err) => write!(f, "syntax error: {err}"),
                    ParseRuleError::ParseScalarError(err) => {
                        write!(f, "could not parse float: {err}")
                    }
                    ParseRuleError::ParseAxisError(err) => write!(f, "could not parse axis: {err}"),
                }
            }
        }

        impl std::error::Error for ParseRuleError {}

        impl FromStr for Axis {
            type Err = String;

            fn from_str(s: &str) -> Result<Self, Self::Err> {
                match s {
                    "x" => Ok(Self::X),
                    "y" => Ok(Self::Y),
                    "z" => Ok(Self::Z),
                    weird => Err(format!("expected one of 'x', 'y', or 'z', but found {weird:?}")),
                }
            }
        }
    }
}
