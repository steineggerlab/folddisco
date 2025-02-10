use std::fmt::{self, Display};
use std::cmp::Ordering;
use std::error;

/// A struct to define the context of an error message
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Context {
    /// When no context can be given
    None,
    /// When only a line (e.g. in a file) can be shown
    Show {
        /// The line to be shown to the user (e.g. filename)
        line: String,
    },
    /// When a full line is faulty and no special position can be annotated
    FullLine {
        /// The line number to recognise where the error is located
        linenumber: usize,
        /// The line to show the issue itself
        line: String,
    },
    /// When a special position can be annotated on a line.
    /// ```text
    ///      |
    /// 104  | ATOM      O  N   MET A   1      27.251  24.447   2.594  1.00 11.79           N
    ///      |        ^^^^
    ///        <-   -><-->
    /// ```
    /// The first space (annotated by `<-`, `->`) is the offset, in this case 7. The
    /// second space is the length, in this case 4.
    Line {
        /// The line number to recognise where the error is located.
        linenumber: usize,
        /// The line to show the issue itself.
        line: String,
        /// The offset of the special position to be annotated.
        offset: usize,
        /// The length of the special position to be annotated.
        length: usize,
    },
    /// To show multiple lines where an error occurred.
    Range {
        /// The linenumber of the first line
        start_linenumber: usize,
        /// The lines to show
        lines: Vec<String>,
        /// The possible offset of the first line, will be padded with spaces
        offset: usize,
    },
    /// To show multiple lines where an error occurred.
    RangeHighlights {
        /// The linenumber of the first line
        start_linenumber: usize,
        /// The lines to show
        lines: Vec<String>,
        /// Highlights defined by the line (relative to the set of lines given), start column in that line and length of highlight
        highlights: Vec<(usize, usize, usize)>,
    },
    /// To show multiple contexts
    Multiple {
        /// The contexts to show
        contexts: Vec<(Option<String>, Context)>,
    },
}

impl Context {
    /// Creates a new context when no context can be given
    pub const fn none() -> Context {
        Context::None
    }

    /// Creates a new context when only a line (eg filename) can be shown
    pub fn show(line: impl std::string::ToString) -> Context {
        Context::Show {
            line: line.to_string(),
        }
    }

    /// Creates a new context when a full line is faulty and no special position can be annotated
    pub fn full_line(linenumber: usize, line: impl std::string::ToString) -> Context {
        Context::FullLine {
            linenumber,
            line: line.to_string(),
        }
    }

    /// Creates a new context when a special position can be annotated on a line
    pub fn line(
        linenumber: usize,
        line: impl std::string::ToString,
        offset: usize,
        length: usize,
    ) -> Context {
        Context::Line {
            linenumber,
            line: line.to_string(),
            offset,
            length,
        }
    }

    /// Creates a new context to highlight a certain position
    #[allow(clippy::unwrap_used)]
    pub fn position(pos: &Position<'_>) -> Context {
        if pos.text.is_empty() {
            Context::Line {
                linenumber: pos.line,
                line: "".to_string(),
                offset: 0,
                length: 3,
            }
        } else {
            Context::Line {
                linenumber: pos.line,
                line: pos.text.lines().next().unwrap().to_string(),
                offset: 0,
                length: 3,
            }
        }
    }

    /// Creates a new context from a start and end point within a single file
    pub fn range(start: &Position<'_>, end: &Position<'_>) -> Context {
        if start.line == end.line {
            Context::Line {
                linenumber: start.line,
                line: start.text[..(end.column - start.column)].to_string(),
                offset: start.column,
                length: end.column - start.column,
            }
        } else {
            Context::Range {
                start_linenumber: start.line,
                lines: start
                    .text
                    .lines()
                    .take(end.line - start.line)
                    .map(ToString::to_string)
                    .collect::<Vec<String>>(),
                offset: start.column,
            }
        }
    }

    /// Display this context, with an optional note after the context.
    fn display(&self, f: &mut fmt::Formatter<'_>, note: Option<&str>) -> fmt::Result {
        let mut tail = true; // End with a tailing line ╵
        #[allow(
            clippy::cast_sign_loss,
            clippy::cast_precision_loss,
            clippy::cast_possible_truncation
        )]
        let get_margin = |n| ((n + 1) as f64).log10().max(1.0).ceil() as usize;
        let margin = match self {
            Context::None => 0,
            Context::Show { .. } => 2,
            Context::FullLine { linenumber: n, .. } => get_margin(*n),
            Context::Line { linenumber: n, .. } => get_margin(*n),
            Context::Range {
                start_linenumber: n,
                lines: l,
                ..
            } => get_margin(n + l.len()),
            Context::RangeHighlights {
                start_linenumber: n,
                lines: l,
                ..
            } => get_margin(n + l.len()),
            Context::Multiple { .. } => 0,
        };
        match self {
            Context::None => {
                return Ok(());
            }
            Context::Show { line } => {
                write!(f, "\n{:pad$} ╷\n{:pad$} │ {}", "", "", line, pad = margin)?
            }
            Context::FullLine { linenumber, line } => write!(
                f,
                "\n{:pad$} ╷\n{:<pad$} │ {}",
                "",
                linenumber,
                line,
                pad = margin
            )?,
            Context::Line {
                linenumber,
                line,
                offset,
                length,
            } => write!(
                f,
                "\n{:pad$} ╷\n{:<pad$} │ {}\n{:pad$} · {}{}",
                "",
                linenumber,
                line,
                "",
                " ".repeat(*offset),
                "─".repeat(*length),
                pad = margin
            )?,
            Context::Range {
                start_linenumber,
                lines,
                offset,
            } => {
                write!(f, "\n{:pad$} ╷", "", pad = margin)?;
                let mut number = *start_linenumber;
                write!(
                    f,
                    "\n{:<pad$} │ {}{}",
                    number,
                    " ".repeat(*offset),
                    lines[0],
                    pad = margin
                )?;
                for line in lines.iter().skip(1) {
                    number += 1;
                    write!(f, "\n{number:<margin$} │ {line}")?;
                }
            }
            Context::RangeHighlights {
                start_linenumber,
                lines,
                highlights,
            } => {
                write!(f, "\n{:pad$} ╷", "", pad = margin)?;
                let mut number = *start_linenumber;
                let mut highlights_peek = highlights.iter().peekable();
                #[allow(unused)]
                for (index, line) in lines.iter().enumerate() {
                    number += 1;
                    write!(f, "\n{number:<margin$} │ {line}")?;
                    let mut first = true;
                    let mut last_offset = 0;
                    while let Some(high) = highlights_peek.peek() {
                        if high.0 > index {
                            break;
                        }
                        if let Some(high) = highlights_peek.next() {
                            if first {
                                write!(f, "\n{:pad$} · ", "", pad = margin)?;
                                first = false;
                            }
                            if last_offset < high.1 {
                                write!(
                                    f,
                                    "{}{}",
                                    " ".repeat(high.1 - last_offset),
                                    "─".repeat(high.2)
                                )?;
                                last_offset = high.1 + high.2;
                            } else {
                                eprintln!("A highlight in a range error message is detected to overlap with a previous highlight, it is skipped.");
                                // Panicking on error gave the following very intense error message (in test code):
                                // `thread panicked while panicking. aborting. ... (exit code: 0xc000001d, STATUS_ILLEGAL_INSTRUCTION)`
                                // To prevent other people from panicking upon seeing this error message this error is not raised currently.
                            }
                        }
                    }
                }
            }
            Context::Multiple { contexts } => {
                for (note, context) in contexts {
                    context.display(f, note.as_deref())?;
                }
                tail = false;
            }
        }
        // Last line
        if let Some(note) = note {
            write!(f, "\n{:pad$} ╰{}", "", note, pad = margin)
        } else if tail {
            write!(f, "\n{:pad$} ╵", "", pad = margin)
        } else {
            Ok(())
        }
    }
}

impl fmt::Display for Context {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.display(f, None)
    }
}

#[derive(Debug, Eq, PartialEq, Copy, Clone)]
/// A position in a file for use in parsing/lexing
pub struct Position<'a> {
    /// The remaining text (as ref so no copies)
    pub text: &'a str,
    /// The current linenumber
    pub line: usize,
    /// The current column number
    pub column: usize,
}

#[derive(PartialEq, Eq, Debug, Copy, Clone, Default)]
/// The strictness to operate in, this defines at which [`ErrorLevel`] the program should stop execution upon finding an error.
pub enum StrictnessLevel {
    /// With `Strict` the program will always stop execution upon finding an error.
    Strict,
    /// With `Medium` the program will allow [`ErrorLevel::GeneralWarning`].
    #[default]
    Medium,
    /// With `Loose` the program will allow [`ErrorLevel::GeneralWarning`] and [`ErrorLevel::LooseWarning`].
    Loose,
}

impl Display for StrictnessLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                StrictnessLevel::Strict => "Strict",
                StrictnessLevel::Medium => "Medium",
                StrictnessLevel::Loose => "Loose",
            }
        )
    }
}

/// This indicates the level of the error, to handle it differently based on the level of the raised error.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum ErrorLevel {
    /// An error that breaks the execution of the program.
    BreakingError,
    /// An error that invalidates the output of the function generating the error. This concerns things like invalid
    /// characters, numeric literals etc.
    InvalidatingError,
    /// A warning that invalidates some strict invariants posed by the specification. These do not necessarily
    /// prevent the code from running, but will need to be checked.
    StrictWarning,
    /// A warning that invalidates some looser defined invariants. These are generally bad but sometimes occur
    /// due to other software packages not following the specifications to the letter.
    LooseWarning,
    /// A general warning.
    GeneralWarning,
}

impl ErrorLevel {
    /// Get the descriptor for this ErrorLevel (Error/Warning). This can be used to display to users to indicate
    /// the severity of the error.
    pub const fn descriptor(&self) -> &str {
        match self {
            ErrorLevel::BreakingError => "BreakingError",
            ErrorLevel::InvalidatingError => "InvalidatingError",
            ErrorLevel::StrictWarning => "StrictWarning",
            ErrorLevel::LooseWarning => "LooseWarning",
            ErrorLevel::GeneralWarning => "GeneralWarning",
        }
    }

    /// Tests if this errors is breaking with the given strictness level
    pub const fn fails(&self, level: StrictnessLevel) -> bool {
        match level {
            StrictnessLevel::Strict => true,
            StrictnessLevel::Medium => !matches!(self, ErrorLevel::GeneralWarning),
            StrictnessLevel::Loose => {
                !matches!(self, ErrorLevel::GeneralWarning | ErrorLevel::LooseWarning)
            }
        }
    }
}

impl fmt::Display for ErrorLevel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.descriptor())
    }
}


/// An error surfacing while handling a PDB
#[derive(PartialEq, Clone, Eq)]
pub struct PDBError {
    /// The level of the error, defining how it should be handled
    level: ErrorLevel,
    /// A short description of the error, generally used as title line
    short_description: String,
    /// A longer description of the error, presented below the context to give more information and helpful feedback
    long_description: String,
    /// The context, in the most general sense this produces output which leads the user to the right place in the code or file
    context: Context,
}

impl PDBError {
    /// Create a new PDBError
    ///
    /// ## Arguments
    /// * `level` - The level of the error, defining how it should be handled
    /// * `short_desc` - A short description of the error, generally used as title line
    /// * `long_desc` -  A longer description of the error, presented below the context to give more information and helpful feedback
    /// * `context` - The context, in the most general sense this produces output which leads the user to the right place in the code or file
    pub fn new(
        level: ErrorLevel,
        short_desc: impl std::string::ToString,
        long_descr: impl std::string::ToString,
        context: Context,
    ) -> PDBError {
        PDBError {
            level,
            short_description: short_desc.to_string(),
            long_description: long_descr.to_string(),
            context,
        }
    }

    /// The level of the error
    pub const fn level(&self) -> ErrorLevel {
        self.level
    }

    /// Tests if this errors is breaking with the given strictness level
    pub fn fails(&self, level: StrictnessLevel) -> bool {
        self.level.fails(level)
    }

    /// Gives the short description or title for this error
    pub fn short_description(&self) -> &str {
        &self.short_description
    }

    /// Gives the long description for this error
    pub fn long_description(&self) -> &str {
        &self.long_description
    }

    /// Gives the context for this error
    pub const fn context(&self) -> &Context {
        &self.context
    }
}

impl fmt::Debug for PDBError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}: {}{}\n{}\n",
            self.level, self.short_description, self.context, self.long_description
        )
    }
}

impl fmt::Display for PDBError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}: {}{}\n{}\n",
            self.level, self.short_description, self.context, self.long_description
        )
    }
}

impl error::Error for PDBError {}

impl PartialOrd for PDBError {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PDBError {
    fn cmp(&self, other: &Self) -> Ordering {
        self.level.cmp(&other.level)
    }
}

#[cfg(test)]
#[allow(clippy::print_stdout)]
mod tests {
    use super::*;
    use crate::Position;

    #[test]
    fn create_empty_error() {
        let a = PDBError::new(ErrorLevel::GeneralWarning, "test", "test", Context::none());
        println!("{a}");
        assert_eq!(format!("{a}"), "GeneralWarning: test\ntest\n");
        assert_eq!(a.level(), ErrorLevel::GeneralWarning);
        assert!(!a.fails(StrictnessLevel::Loose));
    }

    #[test]
    fn create_full_line_error() {
        let a = PDBError::new(
            ErrorLevel::StrictWarning,
            "test",
            "test",
            Context::full_line(1, "testing line"),
        );
        println!("{a}");
        assert_eq!(
            format!("{a}"),
            "StrictWarning: test\n  ╷\n1 │ testing line\n  ╵\ntest\n"
        );
        assert_eq!(a.level(), ErrorLevel::StrictWarning);
        assert!(a.fails(StrictnessLevel::Strict));
    }

    #[test]
    fn create_range_error() {
        let pos1 = Position {
            text: "hello world\nthis is a multiline\npiece of teXt",
            line: 1,
            column: 0,
        };
        let pos2 = Position {
            text: "",
            line: 4,
            column: 13,
        };
        let a = PDBError::new(
            ErrorLevel::LooseWarning,
            "test",
            "test error",
            Context::range(&pos1, &pos2),
        );
        println!("{a}");
        assert_eq!(format!("{a}"), "LooseWarning: test\n  ╷\n1 │ hello world\n2 │ this is a multiline\n3 │ piece of teXt\n  ╵\ntest error\n");
        assert_eq!(a.level(), ErrorLevel::LooseWarning);
        assert!(a.fails(StrictnessLevel::Strict));
        assert_eq!(pos2.text, "");
        assert_eq!(pos2.line, 4);
        assert_eq!(pos2.column, 13);
    }

    #[test]
    fn ordering_and_equality() {
        let a = PDBError::new(ErrorLevel::GeneralWarning, "test", "test", Context::none());
        let b = PDBError::new(ErrorLevel::LooseWarning, "test", "test", Context::none());
        let c = PDBError::new(ErrorLevel::LooseWarning, "test", "test", Context::none());
        let d = PDBError::new(ErrorLevel::BreakingError, "test", "test", Context::none());
        assert_ne!(a, b);
        assert_eq!(b, c);
        assert_ne!(c, d);
        assert!(a > b);
        assert!(c > d);
        assert!(c < a);
    }
}