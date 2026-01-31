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

    // There could be advantages to placing this in the parser, rather than the lexer. But for now
    // I think this is the best choice.
    /// A three-value point, values separated by commas, optionally delimited by parentheses.
    pub(crate) fn point<'s>() -> impl Parser<'s, &'s str, Point, E<'s>> {
        let number = regex(r"[-+]?(\d*\.\d+|\d+\.?\d*)")
            .map(FromStr::from_str)
            .unwrapped()
            .labelled("number");
        let point = number.separated_by(just(',').padded()).collect_exactly::<Point>();
        let parens_point = point.clone().padded().delimited_by(just('('), just(')'));
        choice((point, parens_point))
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
    let integer = any()
        .filter(move |c: &char| c.is_digit(10))
        .repeated()
        .at_least(1)
        .ignored()
        .or(just(char::digit_zero()).ignored())
        .to_slice()
        .map(FromStr::from_str)
        .unwrapped()
        .then_ignore(any().filter(|c: &char| c.is_alphabetic() || c == &'_').not())
        .boxed();
    let float = regex(r"[-+]?(\d*\.\d+)").map(FromStr::from_str).unwrapped().boxed();
    let number =
        choice((float.map(Number::Float), integer.map(Number::Integer))).map(Token::Number).boxed();
    // This is a bit silly, but if we want concentrations like 1M to work like 1.0M, we need a
    // leniant notion of a float that can also match integers. But just for concentrations.
    let leniant_float = regex(r"[-+]?(\d*\.?\d+)").map(f64::from_str).unwrapped().boxed();
    let unit = choice((
        just('M').to(1.0),
        just("mM").to(1e-3),
        just("uM").or(just("µM")).to(1e-6),
        just("nM").to(1e-9),
        just("pM").to(1e-12),
    ));
    let concentration =
        leniant_float.clone().then(unit).map(|(v, u)| Token::Concentration(v * u)).boxed();
    let section = components::identifier()
        .padded()
        .delimited_by(just('['), just(']'))
        .map(Token::Section)
        .boxed();
    let placeholder = components::identifier()
        .padded()
        .delimited_by(just('<'), just('>'))
        .map(Token::Placeholder)
        .boxed();
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
        components::comment().ignored().padded().repeated().at_least(1),
        empty(),
    )))
    .allow_leading()
    .allow_trailing()
    .collect()
}

#[cfg(test)]
mod test {
    use super::*;

    fn p<'s, T>(
        p: impl Parser<'s, &'s str, T, E<'s>>,
        s: &'s str,
    ) -> Result<T, Vec<Rich<'s, char>>> {
        p.parse(s).into_result()
    }

    mod point {
        use super::*;
        #[test]
        fn point_floats() {
            assert_eq!(p(components::point(), "1.0, 2.0, 3.0"), Ok([1.0, 2.0, 3.0]));
        }

        #[test]
        fn point_integers() {
            assert_eq!(p(components::point(), "1, 2, 3"), Ok([1.0, 2.0, 3.0]));
        }

        #[test]
        fn point_mixed() {
            assert_eq!(p(components::point(), "1.0, 2, 3"), Ok([1.0, 2.0, 3.0]));
        }

        #[test]
        fn point_whitespace() {
            assert_eq!(p(components::point(), "1.0,2,3"), Ok([1.0, 2.0, 3.0]));
            assert_eq!(p(components::point(), "1.0 ,2 ,3"), Ok([1.0, 2.0, 3.0]));
            assert_eq!(p(components::point(), "1.0 ,  2 ,  3"), Ok([1.0, 2.0, 3.0]));
        }

        #[test]
        fn point_with_parentheses() {
            assert_eq!(p(components::point(), "(1.0, 2.0, 3.0)"), Ok([1.0, 2.0, 3.0]));
        }

        #[test]
        fn point_with_parentheses_whitespace() {
            assert_eq!(p(components::point(), "( 1.0, 2.0, 3.0  )"), Ok([1.0, 2.0, 3.0]));
        }

        #[test]
        fn bad_points() {
            assert!(p(components::point(), "1.0, 2.0").is_err());
            assert!(p(components::point(), "1.0, (2.0").is_err());
            assert!(p(components::point(), "((1.0, 2.0, 3.0))").is_err());
            assert!(p(components::point(), "1.0, .2.0, 3.0").is_err());
        }
    }

    mod string {
        use super::*;
        #[test]
        fn string() {
            assert_eq!(p(components::string(), r#""abc""#), Ok("abc"));
        }

        #[test]
        fn bad_string_single_quotes() {
            assert!(p(components::string(), r#"'abc'"#).is_err());
        }

        #[test]
        fn bad_string_extra_quote() {
            assert!(p(components::string(), r#""ab"c""#).is_err());
        }
    }

    mod identifier {
        use super::*;

        #[test]
        fn identifier() {
            assert_eq!(p(components::identifier(), "abc"), Ok("abc"));
        }

        #[test]
        fn lysozyme() {
            assert_eq!(p(components::identifier(), "3lyz"), Ok("3lyz"));
        }

        #[test]
        fn zero_prefix() {
            assert_eq!(p(components::identifier(), "0abc"), Ok("0abc"));
        }

        #[test]
        fn zero_int_prefix() {
            assert_eq!(p(components::identifier(), "01abc"), Ok("01abc"));
        }

        #[test]
        fn zero_zero_prefix() {
            assert_eq!(p(components::identifier(), "00abc"), Ok("00abc"));
        }

        #[test]
        fn zero_zero_int_prefix() {
            assert_eq!(p(components::identifier(), "001abc"), Ok("001abc"));
        }

        #[test]
        fn zero_underscore_prefix() {
            assert_eq!(p(components::identifier(), "0_abc"), Ok("0_abc"));
        }

        #[test]
        fn zero_int_underscore_prefix() {
            assert_eq!(p(components::identifier(), "01_abc"), Ok("01_abc"));
        }

        #[test]
        fn zero_zero_underscore_prefix() {
            assert_eq!(p(components::identifier(), "00_abc"), Ok("00_abc"));
        }

        #[test]
        fn zero_zero_int_underscore_prefix() {
            assert_eq!(p(components::identifier(), "001_abc"), Ok("001_abc"));
        }

        #[test]
        fn trailing_numbers() {
            assert_eq!(p(components::identifier(), "abc0123"), Ok("abc0123"));
        }

        #[test]
        fn surrounding_numbers() {
            assert_eq!(p(components::identifier(), "012abc012"), Ok("012abc012"));
        }

        #[test]
        fn containing_double_quote() {
            assert!(p(components::identifier(), r#""ab"c""#).is_err());
        }
    }

    mod comment {
        use super::*;

        #[test]
        fn comment() {
            assert_eq!(p(components::comment(), "# abc"), Ok("# abc"));
            assert_eq!(p(components::comment(), "; abc"), Ok("; abc"));
            assert!(p(components::comment(), "// abc").is_err());
        }

        #[test]
        fn many() {
            assert_eq!(p(components::comment(), "# abc # a # b # c"), Ok("# abc # a # b # c"));
            assert_eq!(p(components::comment(), "; abc"), Ok("; abc"));
            assert!(p(components::comment(), "// abc").is_err());
        }

        #[test]
        fn weird_whitespace() {
            assert_eq!(p(components::comment(), "#abc"), Ok("#abc"));
            assert_eq!(
                p(components::comment(), "\t\t #       \t        abc\t\t"),
                Ok("\t\t #       \t        abc\t\t")
            );
        }
    }
}
