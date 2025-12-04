use chumsky::Parser;
use chumsky::error::Rich;
use chumsky::input::{BorrowInput, Input};
use chumsky::span::SimpleSpan;

use lexer::{Spanned, Token};

use crate::core::config::Config;

mod lexer;
mod parser;
mod writer;

pub use writer::write;

pub(crate) fn make_input<'src>(
    eoi: SimpleSpan,
    toks: &'src [Spanned<Token<'src>>],
) -> impl BorrowInput<'src, Token = Token<'src>, Span = SimpleSpan> {
    toks.map(eoi, |(t, s)| (t, s))
}

mod report {
    use super::*;

    use ariadne::{Color, Fmt, Label, Report, ReportKind, Source};

    // Most of this is taken and modified from the mini_ml.rs example in the chumsky github page.
    // https://github.com/zesterer/chumsky/blob/0.11/examples/mini_ml.rs

    /// Report errors from a parsing step.
    pub(crate) fn error(path: &str, src: &str, err: &Rich<impl std::fmt::Display>) -> String {
        let msg = err.reason().to_string();
        let label = (
            err.found()
                .map(|c| c.to_string())
                .unwrap_or_else(|| "end of input".to_string()),
            *err.span(),
        );
        let extra_labels = err
            .contexts()
            .map(|(l, s)| (format!("while parsing this {l}"), *s));

        let cfg = ariadne::Config::new()
            .with_index_type(ariadne::IndexType::Char)
            .with_label_attach(ariadne::LabelAttach::Middle);
        Report::build(ReportKind::Error, (path, label.1.into_range()))
            .with_config(cfg)
            .with_message(&msg)
            .with_label(
                Label::new((path, label.1.into_range()))
                    .with_message(label.0)
                    .with_color(Color::Red),
            )
            .with_labels(extra_labels.into_iter().map(|label| {
                Label::new((path, label.1.into_range()))
                    .with_message(label.0)
                    .with_color(Color::Yellow)
            }))
            .finish()
            .eprint((path, Source::from(src)))
            .unwrap();
        msg
    }

    /// Report the result of a parsing step by showing the spans of source code associated with
    /// each item.
    pub(crate) fn result<T: std::fmt::Display>(path: &str, src: &str, items: &[Spanned<T>]) {
        let cfg = ariadne::Config::new()
            .with_compact(true)
            .with_index_type(ariadne::IndexType::Char)
            .with_label_attach(ariadne::LabelAttach::Middle);
        for (item, span) in items {
            Report::build(
                ReportKind::Custom("Result", Color::Green),
                (path, span.into_range()),
            )
            .with_config(cfg)
            .with_message(format!("found the following item: {item}"))
            .with_label(
                Label::new((path, span.into_range()))
                    .with_message(item.to_string())
                    .with_color(Color::Green),
            )
            .finish()
            .eprint((path, Source::from(src)))
            .unwrap();
        }
    }

    /// Report placeholders in an input file.
    pub(crate) fn placeholders<'s>(path: &str, src: &'s str, items: &[Spanned<Token<'s>>]) {
        let cfg = ariadne::Config::new()
            .with_compact(true)
            .with_index_type(ariadne::IndexType::Char)
            .with_label_attach(ariadne::LabelAttach::Middle);
        for (item, name, span) in items.iter().filter_map(|(item, span)| match item {
            Token::Placeholder(name) => Some((item, name, span)),
            _ => None,
        }) {
            Report::build(
                ReportKind::Custom("Placeholder", Color::Cyan),
                (path, span.into_range()),
            )
            .with_config(cfg)
            .with_message(format!(
                "replace {} with a valid value",
                format!("<{name}>").fg(Color::Cyan)
            ))
            .with_label(
                Label::new((path, span.into_range()))
                    .with_message(item.to_string())
                    .with_color(Color::Cyan),
            )
            .finish()
            .eprint((path, Source::from(src)))
            .unwrap();
        }
    }
}

pub fn parse_config(path: &str, src: &str) -> anyhow::Result<Config> {
    let tokens = lexer::lexer().parse(src).into_result().map_err(|errs| {
        // Display a nice report.
        let summary = report::error(path, src, &errs[0]);
        // Communicate the error upstream.
        anyhow::anyhow!("encountered an error while lexing {path:?}: {summary}")
    })?;

    if std::env::var("BENTOPY_SHOW_TOKENS").is_ok_and(|v| v.parse::<bool>().unwrap_or_default()) {
        eprintln!("Printing tokens.");
        report::result(path, src, &tokens);
    }

    // Report any placeholders in the stream of tokens.
    report::placeholders(path, src, &tokens);

    let tokens = make_input((0..src.len()).into(), &tokens);
    let config = parser::parser()
        .parse(tokens)
        .into_result()
        .map_err(|errs| {
            // Display a nice report.
            let summary = report::error(path, src, &errs[0]);
            // Communicate the error upstream.
            anyhow::anyhow!("encountered an error while parsing {path:?}: {summary}")
        })?;
    Ok(config)
}
