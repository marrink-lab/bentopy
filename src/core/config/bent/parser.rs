use std::path::PathBuf;

use chumsky::input::BorrowInput;
use chumsky::prelude::*;

use crate::core::config::bent::lexer::{Number, Token};
use crate::core::config::{
    Anchor, Axes, Axis, Center, Compartment, Config, Constraint, Dimensions, Expr, General, Limit,
    Mask, Op, Quantity, RearrangeMethod, Rule, Segment, Shape, Space,
};

/// A shorter alias for this otherwise rather unwieldy type.
type E<'src, 'tokens> = extra::Err<Rich<'tokens, Token<'src>>>;

mod components {
    use super::*;

    pub(crate) fn ident<'ts, 's: 'ts, I>(
        word: &'static str,
    ) -> impl Parser<'ts, I, Token<'s>, E<'s, 'ts>>
    where
        I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
    {
        just(Token::Ident(word)).boxed()
    }

    pub(crate) fn axis<'ts, 's: 'ts, I>() -> impl Parser<'ts, I, Axis, E<'s, 'ts>>
    where
        I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
    {
        select! {
            Token::Ident("x") => Axis::X,
            Token::Ident("y") => Axis::Y,
            Token::Ident("z") => Axis::Z,
        }
        .labelled("axis")
    }

    pub(crate) fn limit<'ts, 's: 'ts, I>() -> impl Parser<'ts, I, Limit, E<'s, 'ts>>
    where
        I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
    {
        let op = select! {
            Token::GreaterThan => Op::GreaterThan,
            Token::LessThan => Op::LessThan,
        }
        .labelled("inequality operator");
        let value = select! { Token::Number(n) => n.as_float() }.labelled("number");

        // For example, x > 10.
        let axis_op_value = group((axis(), op, value))
            .map(|(axis, op, value)| Limit { axis, op, value })
            .boxed();
        // For example, 10 < x.
        let value_op_axis = group((value, op, axis()))
            .map(|(value, rev_op, axis)| Limit {
                axis,
                op: rev_op.reverse(),
                value,
            })
            .boxed();

        choice((axis_op_value, value_op_axis)).labelled("limit")
    }

    pub fn expr_parser<'ts, 's: 'ts, I, T: 'ts>(
        atom: impl Parser<'ts, I, T, E<'s, 'ts>> + 'ts,
    ) -> impl Parser<'ts, I, Expr<T>, E<'s, 'ts>>
    where
        I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
    {
        use chumsky::pratt::*;

        recursive(|expr| {
            let atom = atom.map(Expr::Term).boxed();
            let group = expr
                .delimited_by(just(Token::ParenOpen), just(Token::ParenClose))
                .boxed();
            choice((atom, group))
                .pratt((
                    prefix(3, just(Token::Not), |_op, rhs, _| Expr::Not(Box::new(rhs))).boxed(),
                    infix(left(2), just(Token::And), |lhs, _op, rhs, _| {
                        Expr::And(Box::new(lhs), Box::new(rhs))
                    })
                    .boxed(),
                    infix(left(1), just(Token::Or), |lhs, _op, rhs, _| {
                        Expr::Or(Box::new(lhs), Box::new(rhs))
                    })
                    .boxed(),
                ))
                .labelled("expression")
                .as_context()
                .boxed()
        })
    }
}

pub fn terms_parser<'ts, 's: 'ts, I>() -> impl Parser<'ts, I, Expr<String>, E<'s, 'ts>>
where
    I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
{
    components::expr_parser(select! { Token::Ident(s) => s.to_owned() })
        .labelled("combination expression")
        .as_context()
}

pub fn limits_parser<'ts, 's: 'ts, I>() -> impl Parser<'ts, I, Expr<Limit>, E<'s, 'ts>>
where
    I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
{
    components::expr_parser(components::limit())
        .labelled("limits expression")
        .as_context()
}

pub fn include_parser<'ts, 's: 'ts, I>() -> impl Parser<'ts, I, PathBuf, E<'s, 'ts>>
where
    I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
{
    select! { Token::String(s) => s.into() }
}

pub(crate) fn compartment_parser<'ts, 's: 'ts, I>() -> impl Parser<'ts, I, Compartment, E<'s, 'ts>>
where
    I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
{
    let number = select! { Token::Number(n) => n.as_float() as f32 }.boxed();
    let point = select! { Token::Point(point) => point }.boxed();
    let id = select! { Token::Ident(id) => id.to_string() };

    let sphere = {
        let sphere_center = point
            .clone()
            .map(Center::Point)
            .or(components::ident("center").to(Center::Center));
        let diameter = components::ident("diameter")
            .ignore_then(number.clone())
            .map(|d| d * 0.5);
        let radius = components::ident("radius").ignore_then(number.clone());
        components::ident("sphere")
            .ignore_then(components::ident("at"))
            .ignore_then(sphere_center)
            .then_ignore(components::ident("with"))
            .then(choice((diameter, radius)))
            .map(|(center, radius)| Shape::Sphere { center, radius })
            .boxed()
    };
    let cuboid = {
        let anchor = choice((
            point.map(Anchor::Point),
            components::ident("start").to(Anchor::Start),
            components::ident("center").to(Anchor::Center),
            components::ident("end").to(Anchor::End),
        ))
        .boxed();
        components::ident("cuboid")
            .ignore_then(components::ident("from"))
            .ignore_then(anchor.clone())
            .then_ignore(components::ident("to"))
            .then(anchor.clone())
            .map(|(start, end)| Shape::Cuboid { start, end })
            .boxed()
    };
    let voxels = components::ident("from")
        .ignore_then(select! { Token::String(s) => s.into() }.labelled("quoted path"))
        .map(Mask::Voxels)
        .boxed();
    let shape = components::ident("as")
        .ignore_then(choice((sphere, cuboid)))
        .map(Mask::Shape)
        .boxed();
    let limits = components::ident("where")
        .ignore_then(limits_parser())
        .map(Mask::Limits)
        .boxed();
    let within = components::ident("within")
        .ignore_then(select! { Token::Number(n) => n.as_float() }.labelled("distance"))
        .then_ignore(components::ident("of"))
        .then(id.labelled("compartment id"))
        .map(|(distance, id)| Mask::Within {
            distance: distance as f32,
            id,
        })
        .boxed();
    let combination = components::ident("combines")
        .ignore_then(terms_parser())
        .map(Mask::Combination)
        .labelled("compartment combination expression")
        .boxed();
    let all = group(["is", "all"].map(components::ident))
        .to(Mask::All)
        .boxed();

    id.then(choice((voxels, shape, limits, within, combination, all)))
        .map(|(id, mask)| Compartment { id, mask })
        .labelled("compartment declaration")
        .as_context()
}

pub(crate) fn constraint_parser<'ts, 's: 'ts, I>() -> impl Parser<'ts, I, Constraint, E<'s, 'ts>>
where
    I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
{
    let id = select! { Token::Ident(id) => id.to_string() };
    let rotation_axes = components::ident("rotates")
        .ignore_then(
            components::axis()
                .separated_by(just(Token::Comma))
                .collect::<Vec<_>>()
                .map(|axes| Axes {
                    x: axes.contains(&Axis::X),
                    y: axes.contains(&Axis::Y),
                    z: axes.contains(&Axis::Z),
                }),
        )
        .map(Rule::RotationAxes)
        .boxed();

    id.then(rotation_axes)
        .map(|(id, rule)| Constraint { id, rule })
        .labelled("constraint declaration")
        .as_context()
}

pub(crate) fn segment_parser<'ts, 's: 'ts, I>() -> impl Parser<'ts, I, Segment, E<'s, 'ts>>
where
    I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
{
    let identifier = select! { Token::Ident(ident) => ident.to_owned() };
    let name = identifier;
    let tag = identifier;
    let id = name.then(just(Token::Colon).ignore_then(tag).or_not());
    let number = select! { Token::Number(Number::Integer(n)) => n }.map(Quantity::Number);
    let molarity = select! { Token::Concentration(c) => c }.map(Quantity::Concentration);
    let quantity = molarity.or(number);
    let from = components::ident("from");
    let path = from.ignore_then(select! { Token::String(s) => s.into() });
    let ids = select! { Token::Ident(id) => id.to_owned() }
        .separated_by(just(Token::Comma))
        .collect::<Vec<_>>()
        .map(Into::into);
    let compartments = components::ident("in").ignore_then(ids.clone());
    let satisfies = components::ident("satisfies");
    let rules = satisfies
        .ignore_then(ids.labelled("rule names").as_context())
        .or_not()
        .map(Option::unwrap_or_default);

    group((id, quantity, path, compartments, rules))
        .map(
            |((name, tag), quantity, path, compartment_ids, rules)| Segment {
                name,
                tag,
                quantity,
                path,
                compartment_ids,
                rules,
            },
        )
        .labelled("segment declaration")
        .as_context()
}

pub fn section_parser<'ts, 's: 'ts, I, T: 'ts>(
    header: &'static str,
    item: impl Parser<'ts, I, T, E<'s, 'ts>>,
) -> impl Parser<'ts, I, Vec<T>, E<'s, 'ts>>
where
    I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
{
    just(Token::Section(header))
        .ignore_then(item.repeated().collect())
        .labelled(format!("{header} section"))
        .as_context()
}

pub fn general_parser<'ts, 's: 'ts, I>() -> impl Parser<'ts, I, General, E<'s, 'ts>>
where
    I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
{
    enum Field {
        Title(String),
        Seed(u64),
        BeadRadius(f64),
        RearrangeMethod(RearrangeMethod),
        MaxTriesMult(u64),
        MaxTriesRotDiv(u64),
    }

    let title = components::ident("title")
        .ignore_then(select! { Token::String(s) => s.to_string() }.labelled("quoted string"))
        .map(Field::Title);
    let seed = components::ident("seed")
        .ignore_then(select! { Token::Number(Number::Integer(n)) => n })
        .map(Field::Seed);
    let bead_radius = components::ident("bead-radius")
        .ignore_then(select! { Token::Number(n) => n.as_float() })
        .map(Field::BeadRadius);
    let rearrange_method = components::ident("rearrange")
        .ignore_then(choice((
            components::ident("moment").to(RearrangeMethod::Moment),
            components::ident("volume").to(RearrangeMethod::Volume),
            components::ident("bounding-sphere").to(RearrangeMethod::BoundingSphere),
            components::ident("none").to(RearrangeMethod::None),
        )))
        .map(Field::RearrangeMethod);
    let max_tries_mult = components::ident("max-tries-mult")
        .ignore_then(select! { Token::Number(Number::Integer(n)) => n })
        .map(Field::MaxTriesMult);
    let max_tries_rot_div = components::ident("max-tries-rot-div")
        .ignore_then(select! { Token::Number(Number::Integer(n)) => n })
        .map(Field::MaxTriesRotDiv);

    let fields = choice((
        title,
        seed,
        bead_radius,
        rearrange_method,
        max_tries_mult,
        max_tries_rot_div,
    ))
    .labelled("general declaration")
    .as_context()
    .repeated()
    .collect::<Vec<_>>()
    .map(|fields| {
        let mut general = General::default();
        for field in fields {
            match field {
                Field::Title(t) => general.title = Some(t),
                Field::Seed(s) => general.seed = Some(s),
                Field::BeadRadius(r) => general.bead_radius = Some(r),
                Field::RearrangeMethod(m) => general.rearrange_method = Some(m),
                Field::MaxTriesMult(n) => general.max_tries_mult = Some(n),
                Field::MaxTriesRotDiv(n) => general.max_tries_rot_div = Some(n),
            }
        }
        general
    });

    just(Token::Section("general"))
        .ignore_then(fields)
        .labelled("general section")
        .as_context()
}

pub fn space_parser<'ts, 's: 'ts, I>() -> impl Parser<'ts, I, Space, E<'s, 'ts>>
where
    I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
{
    enum Field {
        Dimensions(Dimensions),
        Resolution(f64),
        Periodic(bool),
    }

    let dimensions = components::ident("dimensions")
        .ignore_then(select! { Token::Point(point) => point })
        .map(Field::Dimensions);
    let resolution = components::ident("resolution")
        .ignore_then(select! { Token::Number(n) => n.as_float() })
        .map(Field::Resolution);
    let periodic = components::ident("periodic")
        .ignore_then(choice((
            components::ident("true").to(true),
            components::ident("false").to(false),
        )))
        .map(Field::Periodic);

    let fields = choice((dimensions, resolution, periodic))
        .labelled("space declaration")
        .as_context()
        .repeated()
        .collect::<Vec<_>>()
        .map(|fields| {
            let mut space = Space::default();
            for field in fields {
                match field {
                    Field::Dimensions(d) => space.dimensions = Some(d),
                    Field::Resolution(r) => space.resolution = Some(r),
                    Field::Periodic(b) => space.periodic = Some(b),
                }
            }

            space
        });

    just(Token::Section("space"))
        .ignore_then(fields)
        .labelled("space section")
        .as_context()
}

pub fn includes_parser<'ts, 's: 'ts, I>() -> impl Parser<'ts, I, Vec<PathBuf>, E<'s, 'ts>>
where
    I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
{
    section_parser("includes", include_parser())
}

pub fn compartments_parser<'ts, 's: 'ts, I>() -> impl Parser<'ts, I, Vec<Compartment>, E<'s, 'ts>>
where
    I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
{
    section_parser("compartments", compartment_parser())
}

pub fn constraints_parser<'ts, 's: 'ts, I>() -> impl Parser<'ts, I, Vec<Constraint>, E<'s, 'ts>>
where
    I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
{
    section_parser("constraints", constraint_parser())
}

pub fn segments_parser<'ts, 's: 'ts, I>() -> impl Parser<'ts, I, Vec<Segment>, E<'s, 'ts>>
where
    I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
{
    section_parser("segments", segment_parser())
}

pub fn parser<'ts, 's: 'ts, I>() -> impl Parser<'ts, I, Config, E<'s, 'ts>>
where
    I: BorrowInput<'ts, Token = Token<'s>, Span = SimpleSpan>,
{
    enum Section {
        General(General),
        Space(Space),
        Includes(Vec<PathBuf>),
        Constraints(Vec<Constraint>),
        Compartments(Vec<Compartment>),
        Segments(Vec<Segment>),
    }

    let general = general_parser().map(Section::General);
    let space = space_parser().map(Section::Space);
    let includes = includes_parser().map(Section::Includes);
    let constraints = constraints_parser().map(Section::Constraints);
    let compartments = compartments_parser().map(Section::Compartments);
    let segments = segments_parser().map(Section::Segments);

    choice((
        general,
        space,
        includes,
        constraints,
        compartments,
        segments,
    ))
    .repeated()
    .collect::<Vec<_>>()
    .map(|sections| {
        let mut general = General::default();
        let mut space = Space::default();
        let mut includes = Vec::new();
        let mut constraints = Vec::new();
        let mut compartments = Vec::new();
        let mut segments = Vec::new();

        for section in sections {
            match section {
                Section::General(g) => {
                    // TODO: Find something for this update pattern.
                    general = General {
                        title: g.title.or(general.title),
                        seed: g.seed.or(general.seed),
                        bead_radius: g.bead_radius.or(general.bead_radius),
                        rearrange_method: g.rearrange_method.or(general.rearrange_method),
                        max_tries_mult: g.max_tries_mult.or(general.max_tries_mult),
                        max_tries_rot_div: g.max_tries_rot_div.or(general.max_tries_rot_div),
                    }
                }
                Section::Space(s) => {
                    // TODO: Find something for this update pattern.
                    space = Space {
                        dimensions: s.dimensions.or(space.dimensions),
                        resolution: s.resolution.or(space.resolution),
                        periodic: s.periodic.or(space.periodic),
                    }
                }
                Section::Includes(items) => includes.extend(items),
                Section::Constraints(items) => constraints.extend(items),
                Section::Compartments(items) => compartments.extend(items),
                Section::Segments(items) => segments.extend(items),
            }
        }

        Config {
            general,
            space,
            includes,
            constraints,
            compartments,
            segments,
        }
    })
}

#[cfg(test)]
mod tests {
    use chumsky::Parser;

    use crate::core::config::bent::{lexer::lexer, make_input};

    use super::*;

    macro_rules! lex_and_parse {
        ($parser:path => $src:ident) => {{
            let src = $src;
            // HACK: This sucks, but it allows us to move on and just write the nice tests.
            let tokens = Box::new(lexer().parse(src).into_result().expect("correct tokens")).leak();
            let tokens = make_input((0..src.len()).into(), tokens);
            $parser().parse(tokens).into_result()
        }};
    }

    mod components {
        use super::*;

        // Some helpers to ease the pain of boxing everything.
        type E = Expr<Limit>;

        fn limit(axis: Axis, op: Op, value: f64) -> E {
            term(Limit { axis, op, value })
        }

        fn term(t: Limit) -> E {
            E::Term(t)
        }

        fn not(e: E) -> E {
            E::Not(Box::new(e))
        }

        fn or(lhs: E, rhs: E) -> E {
            E::Or(Box::new(lhs), Box::new(rhs))
        }

        fn and(lhs: E, rhs: E) -> E {
            E::And(Box::new(lhs), Box::new(rhs))
        }

        #[test]
        fn terms() {
            let src = "not ((a or b)) and (not c)";
            let res = lex_and_parse!(terms_parser => src);
            assert_eq!(
                res,
                Ok(Expr::And(
                    Box::new(Expr::Not(Box::new(Expr::Or(
                        Box::new(Expr::Term("a".to_string())),
                        Box::new(Expr::Term("b".to_string()))
                    )))),
                    Box::new(Expr::Not(Box::new(Expr::Term("c".to_string()))))
                ))
            );
        }

        #[test]
        fn limits_parens() {
            let src = "not ((x < 10 or (y > 3.1))) and (not 40 > z)";
            let res = lex_and_parse!(limits_parser => src);
            assert_eq!(
                res,
                Ok(and(
                    not(or(
                        limit(Axis::X, Op::LessThan, 10.0),
                        limit(Axis::Y, Op::GreaterThan, 3.1)
                    )),
                    not(limit(Axis::Z, Op::LessThan, 40.0))
                ))
            );
        }

        #[test]
        fn limits_many_parens() {
            let src = "(not ((x < 10 or (y > 3.1))) and (((not 40 > z))))";
            let res = lex_and_parse!(limits_parser => src);
            assert_eq!(
                res,
                Ok(and(
                    not(or(
                        limit(Axis::X, Op::LessThan, 10.0),
                        limit(Axis::Y, Op::GreaterThan, 3.1)
                    )),
                    not(limit(Axis::Z, Op::LessThan, 40.0))
                ))
            );
        }

        /// This expression can be handled nicely by the Pratt parser.
        #[test]
        fn limits_pratt() {
            let src = "not x < 10 or y > 3.1 and not 40 > z";
            let res = lex_and_parse!(limits_parser => src);
            assert_eq!(
                res,
                Ok(or(
                    not(limit(Axis::X, Op::LessThan, 10.0)),
                    and(
                        limit(Axis::Y, Op::GreaterThan, 3.1),
                        not(limit(Axis::Z, Op::LessThan, 40.0))
                    )
                ))
            );
        }
    }

    mod general {
        use super::*;

        #[test]
        fn empty() {
            let src = "[ general ]
";
            let res = lex_and_parse!(general_parser => src);
            assert_eq!(res, Ok(General::default()))
        }

        #[test]
        fn some() {
            let src = r#"[ general ]
title "abc"
seed 1234"#;
            let res = lex_and_parse!(general_parser => src);
            assert_eq!(
                res,
                Ok(General {
                    title: Some("abc".to_string()),
                    seed: Some(1234),
                    ..Default::default()
                })
            )
        }

        #[test]
        fn all() {
            let src = r#"[ general ]
title       "abc"
seed        1234
bead-radius 0.3

   rearrange moment
max-tries-mult 
  1000
max-tries-rot-div
  100
"#;
            let res = lex_and_parse!(general_parser => src);
            assert_eq!(
                res,
                Ok(General {
                    title: Some("abc".to_string()),
                    seed: Some(1234),
                    bead_radius: Some(0.3),
                    rearrange_method: Some(RearrangeMethod::Moment),
                    max_tries_mult: Some(1000),
                    max_tries_rot_div: Some(100),
                })
            )
        }
    }

    mod space {
        use super::*;

        #[test]
        fn empty() {
            let src = "[ space ]
";
            let res = lex_and_parse!(space_parser => src);
            assert_eq!(res, Ok(Space::default()))
        }

        #[test]
        fn dimensions() {
            let src = "[ space ]
dimensions 10,10,10.0";
            let res = lex_and_parse!(space_parser => src);
            assert_eq!(
                res,
                Ok(Space {
                    dimensions: Some([10.0; 3]),
                    resolution: None,
                    periodic: None
                })
            )
        }

        #[test]
        fn all() {
            let src = "[ space ]
dimensions 10,10,10.0
resolution 0.5
periodic false
";
            let res = lex_and_parse!(space_parser => src);
            assert_eq!(
                res,
                Ok(Space {
                    dimensions: Some([10.0; 3]),
                    resolution: Some(0.5),
                    periodic: Some(false)
                })
            )
        }
    }

    mod includes {
        use super::*;

        #[test]
        fn include() {
            let src = r#""forcefield/martini.itp""#;
            let res = lex_and_parse!(include_parser => src);
            assert_eq!(res, Ok("forcefield/martini.itp".into()))
        }

        #[test]
        fn single() {
            let src = r#"[ includes ]
"forcefield/martini.itp""#;
            let res = lex_and_parse!(includes_parser => src);
            assert_eq!(res, Ok(vec!["forcefield/martini.itp".into()]))
        }

        #[test]
        fn multiple() {
            let src = r#"[ includes ]
"forcefield/martini.itp"
"structures/*.itp""#;
            let res = lex_and_parse!(includes_parser => src);
            assert_eq!(
                res,
                Ok(vec![
                    "forcefield/martini.itp".into(),
                    "structures/*.itp".into()
                ])
            )
        }

        #[test]
        fn single_with_trailing_space() {
            let src = "[ includes ]\t\t
\"forcefield/martini.itp\"       ";
            let res = lex_and_parse!(includes_parser => src);
            assert_eq!(res, Ok(vec!["forcefield/martini.itp".into()]))
        }

        #[test]
        fn quoted_path_with_spaces() {
            let src = r#"[ includes ]   
 "quoted-without-spaces"
 "quoted with spaces"
"#;
            let res = lex_and_parse!(includes_parser => src);
            assert_eq!(
                res,
                Ok(vec![
                    "quoted-without-spaces".into(),
                    "quoted with spaces".into(),
                ])
            )
        }
    }

    mod compartment {
        use super::*;

        #[test]
        fn all() {
            let src = "space is all";
            let res = lex_and_parse!(compartment_parser => src);
            assert_eq!(
                res,
                Ok(Compartment {
                    id: "space".to_string(),
                    mask: Mask::All,
                })
            )
        }

        #[test]
        fn sphere_center() {
            let src = "space as sphere at center with radius 10.0";
            let res = lex_and_parse!(compartment_parser => src);
            assert_eq!(
                res,
                Ok(Compartment {
                    id: "space".to_string(),
                    mask: Mask::Shape(Shape::Sphere {
                        center: Center::Center,
                        radius: 10.0
                    }),
                })
            )
        }

        #[test]
        fn sphere_point() {
            let src = "space as sphere at 1,2.0,3 with radius 10.0";
            let res = lex_and_parse!(compartment_parser => src);
            assert_eq!(
                res,
                Ok(Compartment {
                    id: "space".to_string(),
                    mask: Mask::Shape(Shape::Sphere {
                        center: Center::Point([1.0, 2.0, 3.0]),
                        radius: 10.0
                    }),
                })
            )
        }

        #[test]
        fn sphere_diameter() {
            let src = "space as sphere at center with diameter 20.0";
            let res = lex_and_parse!(compartment_parser => src);
            assert_eq!(
                res,
                Ok(Compartment {
                    id: "space".to_string(),
                    mask: Mask::Shape(Shape::Sphere {
                        center: Center::Center,
                        radius: 10.0
                    }),
                })
            )
        }

        #[test]
        fn cuboid_anchored() {
            let src = "space as cuboid from center to end";
            let res = lex_and_parse!(compartment_parser => src);
            assert_eq!(
                res,
                Ok(Compartment {
                    id: "space".to_string(),
                    mask: Mask::Shape(Shape::Cuboid {
                        start: Anchor::Center,
                        end: Anchor::End
                    })
                })
            )
        }

        #[test]
        fn cuboid_point() {
            let src = "space as cuboid from 3,14,1.5 to end";
            let res = lex_and_parse!(compartment_parser => src);
            assert_eq!(
                res,
                Ok(Compartment {
                    id: "space".to_string(),
                    mask: Mask::Shape(Shape::Cuboid {
                        start: Anchor::Point([3.0, 14.0, 1.5]),
                        end: Anchor::End
                    })
                })
            )
        }

        #[test]
        fn close() {
            let src = "close within 5 of mask";
            let res = lex_and_parse!(compartment_parser => src);
            assert_eq!(
                res,
                Ok(Compartment {
                    id: "close".to_string(),
                    mask: Mask::Within {
                        distance: 5.0,
                        id: "mask".to_string()
                    }
                })
            )
        }

        #[test]
        fn limits() {
            let src = "limited where not x > 10.0 and 3 < y";
            let res = lex_and_parse!(compartment_parser => src);
            assert_eq!(
                res,
                Ok(Compartment {
                    id: "limited".to_string(),
                    mask: Mask::Limits(Expr::And(
                        Box::new(Expr::Not(Box::new(Expr::Term(Limit {
                            axis: Axis::X,
                            op: Op::GreaterThan,
                            value: 10.0
                        })))),
                        Box::new(Expr::Term(Limit {
                            axis: Axis::Y,
                            op: Op::GreaterThan,
                            value: 3.0
                        }))
                    ))
                })
            )
        }

        #[test]
        fn limits_simple() {
            let src = "simple-example where x > 10.0";
            let res = lex_and_parse!(compartment_parser => src);
            assert_eq!(
                res,
                Ok(Compartment {
                    id: "simple-example".to_string(),
                    mask: Mask::Limits(Expr::Term(Limit {
                        axis: Axis::X,
                        op: Op::GreaterThan,
                        value: 10.0
                    })),
                })
            )
        }

        #[test]
        fn combination() {
            let src = "space combines not ((a) and b) or (c or d)";
            let res = lex_and_parse!(compartment_parser => src);
            assert_eq!(
                res,
                Ok(Compartment {
                    id: "space".to_string(),
                    mask: Mask::Combination(Expr::Or(
                        Box::new(Expr::Not(Box::new(Expr::And(
                            Box::new(Expr::Term("a".to_string(),)),
                            Box::new(Expr::Term("b".to_string(),)),
                        )),)),
                        Box::new(Expr::Or(
                            Box::new(Expr::Term("c".to_string(),)),
                            Box::new(Expr::Term("d".to_string(),)),
                        )),
                    ),),
                })
            )
        }

        #[test]
        fn combination_weird_whitespace() {
            let src = "space combines a";
            let res = lex_and_parse!(compartment_parser => src);
            assert_eq!(
                res,
                Ok(Compartment {
                    id: "space".to_string(),
                    mask: Mask::Combination(Expr::Term("a".to_string()))
                })
            )
        }

        #[test]
        fn section_single() {
            let src = "[ compartments ]

simple-example where x > 10.0";
            let res = lex_and_parse!(compartments_parser => src);
            assert_eq!(
                res,
                Ok(vec![Compartment {
                    id: "simple-example".to_string(),
                    mask: Mask::Limits(Expr::Term(Limit {
                        axis: Axis::X,
                        op: Op::GreaterThan,
                        value: 10.0
                    })),
                }])
            )
        }

        #[test]
        fn section_double() {
            let src = "[ compartments ]
simple-example where x > 10.0
simple-example where x > 10.0";
            let res = lex_and_parse!(compartments_parser => src);
            assert_eq!(
                res,
                Ok(vec![
                    Compartment {
                        id: "simple-example".to_string(),
                        mask: Mask::Limits(Expr::Term(Limit {
                            axis: Axis::X,
                            op: Op::GreaterThan,
                            value: 10.0
                        })),
                    };
                    2
                ])
            )
        }
    }

    mod segment {
        use super::*;

        fn base_segment() -> Segment {
            Segment {
                name: "3lyz".to_string(),
                tag: None,
                quantity: Quantity::Number(1),
                path: "structures/3lyz.gro".into(),
                compartment_ids: vec!["sphere".to_string()].into(),
                rules: vec!["rule".to_string()].into(),
            }
        }

        #[test]
        fn tag() {
            let src = r#"3lyz:tag 1 from "structures/3lyz.gro" in sphere satisfies rule"#;
            let res = lex_and_parse!(segment_parser => src);
            assert_eq!(
                res,
                Ok(Segment {
                    name: "3lyz".to_string(),
                    tag: Some("tag".to_string()),
                    ..base_segment()
                })
            )
        }

        #[test]
        fn no_tag() {
            let src = r#"3lyz 1 from "structures/3lyz.gro" in sphere satisfies rule"#;
            let res = lex_and_parse!(segment_parser => src);
            assert_eq!(
                res,
                Ok(Segment {
                    name: "3lyz".to_string(),
                    tag: None,
                    ..base_segment()
                })
            )
        }

        #[test]
        fn number() {
            let src = r#"3lyz 1000120 from "structures/3lyz.gro" in sphere satisfies rule"#;
            let res = lex_and_parse!(segment_parser => src);
            assert_eq!(
                res,
                Ok(Segment {
                    quantity: Quantity::Number(1000120),
                    ..base_segment()
                })
            )
        }

        #[test]
        fn concentration() {
            let src = r#"3lyz 0.05M from "structures/3lyz.gro" in sphere satisfies rule"#;
            let res = lex_and_parse!(segment_parser => src);
            assert_eq!(
                res,
                Ok(Segment {
                    quantity: Quantity::Concentration(0.05),
                    ..base_segment()
                })
            )
        }

        #[test]
        fn two_compartments() {
            let src = r#"3lyz 1 from "structures/3lyz.gro" in a,b satisfies rule"#;
            let res = lex_and_parse!(segment_parser => src);
            assert_eq!(
                res,
                Ok(Segment {
                    compartment_ids: vec!["a".to_string(), "b".to_string()].into(),
                    ..base_segment()
                })
            )
        }

        #[test]
        fn no_rules() {
            let src = r#"3lyz 1 from "structures/3lyz.gro" in sphere"#;
            let res = lex_and_parse!(segment_parser => src);
            assert_eq!(
                res,
                Ok(Segment {
                    rules: Default::default(),
                    ..base_segment()
                })
            )
        }

        #[test]
        fn two_rules() {
            let src = r#"3lyz 1 from "structures/3lyz.gro" in sphere satisfies rule1,rule2"#;
            let res = lex_and_parse!(segment_parser => src);
            assert_eq!(
                res,
                Ok(Segment {
                    rules: vec!["rule1".to_string(), "rule2".to_string()].into(),
                    ..base_segment()
                })
            )
        }

        #[test]
        fn section_normal() {
            let src = r#"[ segments ]
3lyz 1 from "structures/3lyz.gro" in sphere satisfies rule
3lyz 1 from "structures/3lyz.gro" in sphere satisfies rule"#;
            let res = lex_and_parse!(segments_parser => src);
            assert_eq!(res, Ok(vec![base_segment(); 2]));
        }

        #[test]
        fn section_weird_whitespace_lines() {
            let src = "[ segments ]
       \t
            3lyz 1 from \"structures/3lyz.gro\" in sphere satisfies rule     



    3lyz 1 from \"structures/3lyz.gro\" in sphere satisfies rule 
";
            let res = lex_and_parse!(segments_parser => src);
            assert_eq!(res, Ok(vec![base_segment(); 2]));
        }

        #[test]
        fn section_wrapped_line() {
            let src = r#"[ segments ]

3lyz 
    1 
    from "structures/3lyz.gro"
    in sphere 
    satisfies rule
3lyz 1 from 

    "structures/3lyz.gro"
in sphere satisfies rule  
"#;
            let res = lex_and_parse!(segments_parser => src);
            assert_eq!(res, Ok(vec![base_segment(); 2]));
        }
    }
}
