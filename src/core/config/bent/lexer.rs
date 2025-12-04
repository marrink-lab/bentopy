use std::str::FromStr;

use chumsky::{prelude::*, text::Char};

use crate::core::config::Point;

/// A shorter alias for this otherwise rather unwieldy type.
type E<'src> = extra::Err<Rich<'src, char>>;

// TODO: Move to chumsky::Spanned when that is released.
pub type Spanned<T> = (T, SimpleSpan);

#[derive(Debug, Clone, PartialEq)]
pub enum Number {
    Float(f64),
    Integer(u64),
}

impl Number {
    pub fn as_float(&self) -> f64 {
        match *self {
            Self::Float(f) => f,
            Self::Integer(i) => i as f64,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Token<'s> {
    // Values
    Ident(&'s str),
    Number(Number),
    Point(Point),
    Concentration(f64),
    String(&'s str),
    Section(&'s str),
    Placeholder(&'s str),

    // Interpunction
    ParenOpen,
    ParenClose,
    Comma,
    Colon,

    // Operators
    Not,
    And,
    Or,
    GreaterThan,
    LessThan,
}

impl std::fmt::Display for Token<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Token::Ident(ident) => write!(f, "word {ident:?}"),
            Token::Number(Number::Float(n)) => write!(f, "float {n}"),
            Token::Number(Number::Integer(n)) => write!(f, "integer {n}"),
            Token::Point([x, y, z]) => write!(f, "point ({x}, {y}, {z})"),
            Token::String(path) => write!(f, "string {path:?}"),
            Token::Concentration(conc) => write!(f, "concentration ({conc}M)"),
            Token::Section(header) => write!(f, "section header {header:?}"),
            Token::Placeholder(name) => write!(f, "placeholder {name:?}"),
            Token::ParenOpen => write!(f, "("),
            Token::ParenClose => write!(f, ")"),
            Token::Comma => write!(f, ","),
            Token::Colon => write!(f, ":"),
            Token::Not => write!(f, "not"),
            Token::And => write!(f, "and"),
            Token::Or => write!(f, "or"),
            Token::GreaterThan => write!(f, ">"),
            Token::LessThan => write!(f, "<"),
        }
    }
}

mod components {
    use super::*;

    // TODO: Consider removing this token and parsing it only where we expect a point. It could be
    // possible, actually, to have three number ids in a segment decl compartments list. Those
    // kinds of ids are technically legal, I guess. Although.. not really, they will always be
    // parsed as an ident. Hmm. Think about it for a while.
    pub(crate) fn point<'s>() -> impl Parser<'s, &'s str, Point, E<'s>> {
        let number = regex(r"[-+]?(\d*\.\d+|\d+\.?\d*)")
            .map(FromStr::from_str)
            .unwrapped()
            .labelled("number");
        number
            .separated_by(just(',').padded())
            .collect_exactly::<Point>()
    }

    pub(crate) fn identifier_builder<'s>(
        illegal: &'static [char],
    ) -> impl Parser<'s, &'s str, &'s str, E<'s>> {
        any()
            .filter(|c: &char| !c.is_whitespace() && !illegal.contains(c))
            .repeated()
            .at_least(1)
            .to_slice()
    }

    pub(crate) fn string<'s>() -> impl Parser<'s, &'s str, &'s str, E<'s>> {
        any()
            .filter(|c: &char| !c.is_newline() && *c != '"')
            .repeated()
            .at_least(1)
            .to_slice()
            .padded_by(just('"'))
    }

    pub(crate) fn identifier<'s>() -> impl Parser<'s, &'s str, &'s str, E<'s>> {
        identifier_builder(&['(', ')', '[', ']', ',', ':', '"', '<', '>'])
    }

    pub(crate) fn comment<'s>() -> impl Parser<'s, &'s str, &'s str, E<'s>> {
        choice((just('#'), just(';')))
            .padded()
            .ignore_then(any().filter(|c: &char| !c.is_newline()).repeated())
            .to_slice()
    }
}

pub fn lexer<'s>() -> impl Parser<'s, &'s str, Vec<Spanned<Token<'s>>>, extra::Err<Rich<'s, char>>>
{
    let operator = |s: &'static str| just(s).labelled(s);
    let integer = text::int(10)
        .map(FromStr::from_str)
        .unwrapped()
        // This is a bit cursed, but we need to make sure that an int is not followed by an
        // alphabetic character.
        .then_ignore(any().filter(|c: &char| c.is_alphabetic()).not())
        .boxed();
    let float = regex(r"[-+]?(\d*\.\d+)").map(FromStr::from_str).unwrapped();
    let number = choice((
        float.clone().map(Number::Float),
        integer.map(Number::Integer),
    ))
    .map(Token::Number);
    // TODO: Introduce mM and nM quantities as well.
    let concentration = float.then_ignore(just('M')).map(Token::Concentration);
    let section = components::identifier()
        .padded()
        .delimited_by(just('['), just(']'))
        .map(Token::Section);
    let placeholder = components::identifier()
        .padded()
        .delimited_by(just('<'), just('>'))
        .map(Token::Placeholder);
    let point = components::point().map(Token::Point).boxed();

    choice((
        // Unparseable
        placeholder,
        // Operators
        operator("not").to(Token::Not),
        operator("and").to(Token::And),
        operator("or").to(Token::Or),
        operator("<").to(Token::LessThan),
        operator(">").to(Token::GreaterThan),
        // Values
        point,
        concentration,
        number,
        section,
        components::string().map(Token::String),
        components::identifier().map(Token::Ident),
        // Separators
        just(',').to(Token::Comma),
        just(':').to(Token::Colon),
        just('(').to(Token::ParenOpen),
        just(')').to(Token::ParenClose),
    ))
    .map_with(|t, e| (t, e.span()))
    .padded()
    // Tokens usually follow each other directly, but there may be comments between them.
    .separated_by(choice((
        components::comment()
            .ignored()
            .padded()
            .repeated()
            .at_least(1),
        empty(),
    )))
    .allow_leading()
    .allow_trailing()
    .collect()
}
