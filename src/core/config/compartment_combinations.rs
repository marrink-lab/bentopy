use std::str::FromStr;

use super::CompartmentID as Id;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum Token {
    Id(Id),
    Not,       // `!`
    Union,     // `union(`
    Intersect, // `intersect(`
    Close,     // `)`
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Expression {
    Id(Id),
    Not(Box<Self>),
    Union(Vec<Self>),
    Intersect(Vec<Self>),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParseExpressionError {
    UnexpectedEnd,
    UnexpectedClose,
    TrailingTokens,
}

impl std::error::Error for ParseExpressionError {}
impl std::fmt::Display for ParseExpressionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseExpressionError::UnexpectedEnd => {
                "encountered the end of the tokens, but expected more"
            }
            ParseExpressionError::UnexpectedClose => "encounterd an unexpected closing parenthesis",
            ParseExpressionError::TrailingTokens => {
                "some tokens remained after parsing the root expression"
            }
        }
        .fmt(f)
    }
}

fn consume_whitespace(input: &str) -> &str {
    input.trim_start()
}

/// Tokenize a single ID.
fn tokenize_id<'s>(mut input: &'s str, tokens: &mut Vec<Token>) -> &'s str {
    // Either, the end of the token is the end of `input`, whitespace, or a closing parenthesis.
    let end = input
        .find(|ch: char| ch.is_whitespace() || ch == ')')
        .unwrap_or(input.len());
    let id;
    (id, input) = input.split_at(end);
    // Make sure we're not accidentally emitting an empty id. This can happen if we're at the end
    // of the input string. That is usually wrong, but we'll catch that during the parse.
    if !id.is_empty() {
        tokens.push(Token::Id(id.to_string())); // Emit.
    }
    input
}

/// Tokenize the next operation, if there is one.
fn tokenize_operation<'s>(
    mut input: &'s str,
    tokens: &mut Vec<Token>,
) -> Result<Option<&'s str>, ParseExpressionError> {
    input = consume_whitespace(input);

    if input.starts_with('!') {
        tokens.push(Token::Not); // Emit.
        input = &input[1..];

        // We now expect another expression.
        return Ok(Some(tokenize_expr(input, tokens)?));
    }

    if input.starts_with("union(") {
        tokens.push(Token::Union); // Emit.
        input = &input[6..];

        // We now expect a list of operands.
        return Ok(Some(tokenize_operands(input, tokens)?));
    }

    if input.starts_with("intersect(") {
        tokens.push(Token::Intersect); // Emit.
        input = &input[10..];

        // We now expect a list of operands.
        return Ok(Some(tokenize_operands(input, tokens)?));
    }

    // Syntax errors.
    if input.starts_with(')') {
        return Err(ParseExpressionError::UnexpectedClose);
    }

    Ok(None)
}

/// Tokenize the operands of an operation until a [`Token::Close`] is found the same tree level.
fn tokenize_operands<'s>(
    mut input: &'s str,
    tokens: &mut Vec<Token>,
) -> Result<&'s str, ParseExpressionError> {
    while !input.is_empty() {
        // We now expect another expression.
        input = tokenize_expr(input, tokens)?;

        // If we have reached the end of this body, we are done.
        input = consume_whitespace(input);
        if input.starts_with(')') {
            tokens.push(Token::Close); // Emit.
            input = &input[1..];
            break;
        }
    }

    Ok(input)
}

/// Tokenize the next expression.
fn tokenize_expr<'s>(
    mut input: &'s str,
    tokens: &mut Vec<Token>,
) -> Result<&'s str, ParseExpressionError> {
    input = consume_whitespace(input);
    match tokenize_operation(input, tokens)? {
        Some(tail) => Ok(tail),
        None => Ok(tokenize_id(input, tokens)),
    }
}

/// Tokenize a compartment combinations definition input string.
fn tokenize(mut input: &str) -> Result<Vec<Token>, ParseExpressionError> {
    let mut tokens = Vec::new();
    while !input.is_empty() {
        input = tokenize_expr(input, &mut tokens)?;
    }
    Ok(tokens)
}

/// Returns a list of operands until the [`Token::Close`] if this level in the tree is encountered.
fn parse_operands_until_close(
    tokens: &mut &[Token],
) -> Result<Vec<Expression>, ParseExpressionError> {
    let mut operands = Vec::new();
    while let Some(item) = parse_next_item(tokens)? {
        operands.push(item);
    }
    Ok(operands)
}

/// Returns the next item or `None` when the [`Token::Close`] for this level in the three is hit.
fn parse_next_item(tokens: &mut &[Token]) -> Result<Option<Expression>, ParseExpressionError> {
    let token;
    (token, *tokens) = tokens
        .split_first()
        .ok_or(ParseExpressionError::UnexpectedEnd)?;
    let expr = match token {
        Token::Id(id) => Expression::Id(id.to_owned()),
        Token::Not => Expression::Not(Box::new(
            parse_next_item(tokens)?.ok_or(ParseExpressionError::UnexpectedClose)?,
        )),
        Token::Union => Expression::Union(parse_operands_until_close(tokens)?),
        Token::Intersect => Expression::Intersect(parse_operands_until_close(tokens)?),
        Token::Close => return Ok(None),
    };

    Ok(Some(expr))
}

fn parse(tokens: &mut &[Token]) -> Result<Expression, ParseExpressionError> {
    parse_next_item(tokens)?.ok_or(ParseExpressionError::UnexpectedClose)
}

impl FromStr for Expression {
    type Err = ParseExpressionError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let tokens = tokenize(s)?;
        let tokens = &mut tokens.as_slice();
        let expr = parse(tokens)?;
        if tokens.is_empty() {
            // No trailing tokens, as desired.
            Ok(expr)
        } else {
            Err(Self::Err::TrailingTokens)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokens() {
        let expected = Ok(vec![
            Token::Union,
            Token::Id("first".to_string()),
            Token::Intersect,
            Token::Id("second".to_string()),
            Token::Id("third".to_string()),
            Token::Close,
            Token::Not,
            Token::Union,
            Token::Id("fourth".to_string()),
            Token::Id("fifth".to_string()),
            Token::Close,
            Token::Close,
        ]);

        let clean = "union(first intersect(second third) !union(fourth fifth))";
        let extra_spaces =
            "    union(  first   intersect(  second   third)  !   union(  fourth  fifth  )  )  ";
        let messy = "     
   union(      first      
intersect(  

         second        third)  !        
        union(   fourth          

fifth          


  
) 
               )                
";

        assert_eq!(tokenize(clean), expected);
        assert_eq!(tokenize(extra_spaces), expected);
        assert_eq!(tokenize(messy), expected);
    }

    #[test]
    fn parsing() {
        let expected = Expression::Union(vec![
            Expression::Id("first".to_string()),
            Expression::Intersect(vec![
                Expression::Id("second".to_string()),
                Expression::Id("third".to_string()),
            ]),
            Expression::Not(Box::new(Expression::Union(vec![
                Expression::Id("fourth".to_string()),
                Expression::Id("fifth".to_string()),
            ]))),
        ]);

        let tokens = vec![
            Token::Union,
            Token::Id("first".to_string()),
            Token::Intersect,
            Token::Id("second".to_string()),
            Token::Id("third".to_string()),
            Token::Close,
            Token::Not,
            Token::Union,
            Token::Id("fourth".to_string()),
            Token::Id("fifth".to_string()),
            Token::Close,
            Token::Close,
        ];

        assert_eq!(parse(&mut tokens.as_slice()), Ok(expected));
    }

    #[test]
    fn parse_from_string() {
        let expected = Expression::Union(vec![
            Expression::Id("first".to_string()),
            Expression::Intersect(vec![
                Expression::Id("second".to_string()),
                Expression::Id("third".to_string()),
            ]),
            Expression::Not(Box::new(Expression::Union(vec![
                Expression::Id("fourth".to_string()),
                Expression::Id("fifth".to_string()),
            ]))),
        ]);

        let s = "union(first intersect(second third) !union(fourth fifth))";

        assert_eq!(Expression::from_str(s), Ok(expected));
    }

    #[test]
    fn parse_not() {
        let expected = Ok(Expression::Not(Box::new(Expression::Id("a".to_string()))));
        let unexpected_close = Err(ParseExpressionError::UnexpectedClose);

        assert_eq!(Expression::from_str("!a)"), unexpected_close);
        assert_eq!(Expression::from_str("!a "), expected);
        assert_eq!(Expression::from_str("!a"), expected);
    }
}
