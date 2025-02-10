pub mod lex_item;
pub mod error;
use lex_item::{DataBlock, DataItem, Item, Loop, SaveFrame, Single, Value};
use error::{Context, ErrorLevel, PDBError, Position};

/// Parse/lex a CIF file into CIF intermediate structure
pub fn lex_cif(text: &str) -> Result<DataBlock, PDBError> {
    parse_main(&mut Position {
        text,
        line: 1,
        column: 1,
    })
}

/// Parse a CIF file
fn parse_main(input: &mut Position<'_>) -> Result<DataBlock, PDBError> {
    trim_comments_and_whitespace(input);
    parse_data_block(input)
}

/// Parse a data block, the main item of a CIF file
fn parse_data_block(input: &mut Position<'_>) -> Result<DataBlock, PDBError> {
    if start_with(input, "data_").is_none() {
        return Err(PDBError::new(
            ErrorLevel::BreakingError,
            "Data Block not opened",
            "The data block should be opened with \"data_\".",
            Context::position(input),
        ));
        
    };
    let identifier = parse_identifier(input);
    let mut block = DataBlock {
        name: identifier.to_string(),
        items: Vec::new(),
    };
    loop {
        trim_comments_and_whitespace(input);
        if input.text.is_empty() {
            return Ok(block);
        }
        let item = parse_data_item_or_save_frame(input)?;
        block.items.push(item);
    }
}

/// Parse a main loop item, a data item or a save frame
fn parse_data_item_or_save_frame(input: &mut Position<'_>) -> Result<Item, PDBError> {
    let start = *input;
    if let Some(()) = start_with(input, "save_") {
        let mut frame = SaveFrame {
            name: parse_identifier(input).to_string(),
            items: Vec::new(),
        };
        while let Ok(item) = parse_data_item(input) {
            frame.items.push(item);
        }
        if let Some(()) = start_with(input, "save_") {
            Ok(Item::SaveFrame(frame))
        } else {
            Err(PDBError::new(
                ErrorLevel::BreakingError,
                "No matching \'save_\' found",
                "A save frame was instantiated but not closed (correctly)",
                Context::range(&start, input),
            ))
        }
    } else {
        let item = parse_data_item(input)?;
        Ok(Item::DataItem(item))
    }
}

/// Parse a data item, a loop or a single item
fn parse_data_item(input: &mut Position<'_>) -> Result<DataItem, PDBError> {
    let start = *input;
    trim_comments_and_whitespace(input);
    if let Some(()) = start_with(input, "loop_") {
        let mut loop_value = Loop {
            header: Vec::new(),
            data: Vec::new(),
        };
        let mut values = Vec::new();
        trim_comments_and_whitespace(input);

        while let Some(()) = start_with(input, "_") {
            let inner_name = parse_identifier(input);
            loop_value.header.push(inner_name.to_string());
            trim_comments_and_whitespace(input);
        }

        while let Ok(value) = parse_value(input) {
            values.push(value);
        }

        let columns = loop_value.header.len();
        if values.len() % columns == 0 {
            loop_value.data = Vec::with_capacity(values.len() / columns);
            let mut iter = values.into_iter().peekable();
            while iter.peek().is_some() {
                loop_value.data.push((&mut iter).take(columns).collect());
            }
        } else {
            return Err(PDBError::new(
                ErrorLevel::BreakingError,
                "Loop has incorrect number of data items",
                format!("A loop has to have a number of data items which is divisible by the amount of headers but here there are {} items left.", values.len() % columns),
                Context::range(&start, input),
            ));
        }

        Ok(DataItem::Loop(loop_value))
    } else if let Some(()) = start_with(input, "_") {
        let name = parse_identifier(input);
        if let Ok(value) = parse_value(input) {
            Ok(DataItem::Single(Single {
                name: name.to_string(),
                content: value,
            }))
        } else {
            Err(PDBError::new(
                ErrorLevel::BreakingError,
                "No valid Value",
                "A Data Item should contain a value or a loop.",
                Context::range(&start, input),
            ))
        }
    } else {
        Err(PDBError::new(
            ErrorLevel::BreakingError,
            "No valid DataItem",
            "A Data Item should be a tag with a value or a loop.",
            Context::position(input),
        ))
    }
}

/// Parse a value for a data item or inside a loop
fn parse_value(input: &mut Position<'_>) -> Result<Value, PDBError> {
    let start = *input;
    trim_comments_and_whitespace(input);
    if input.text.is_empty() {
        Err(PDBError::new(
            ErrorLevel::BreakingError,
            "Empty value",
            "No text left",
            Context::position(input),
        ))
        // The following are reserved words, and need to be checked with a branching position as otherwise it would consume the keyword if matched
    } else if start_with(&mut input.clone(), "data_").is_some()
        || start_with(&mut input.clone(), "global_").is_some()
        || start_with(&mut input.clone(), "loop_").is_some()
        || start_with(&mut input.clone(), "save_").is_some()
        || start_with(&mut input.clone(), "stop_").is_some()
    {
        Err(PDBError::new(
            ErrorLevel::BreakingError,
            "Use of reserved word",
            "\"data_\", \"global_\", \"loop_\", \"save_\" and \"stop_\" are reserved words.",
            Context::position(input),
        ))
    } else if input.text.starts_with('.') {
        // Technically there could be a number starting with a dot...
        let mut branch: Position<'_> = *input;
        if let Some(value) = parse_numeric(parse_identifier(&mut branch)) {
            input.text = branch.text;
            input.column = branch.column;
            Ok(value)
        } else {
            input.text = &input.text[1..];
            input.column += 1;
            Ok(Value::Inapplicable)
        }
    } else if input.text.starts_with('?') {
        input.text = &input.text[1..];
        input.column += 1;
        Ok(Value::Unknown)
    } else if input.text.starts_with('\'') {
        parse_enclosed(input, '\'').map(|text| Value::Text(text.to_string()))
    } else if input.text.starts_with('\"') {
        parse_enclosed(input, '\"').map(|text| Value::Text(text.to_string()))
    } else if input.text.starts_with(';') {
        parse_multiline_string(input).map(|text| Value::Text(text.to_string()))
    } else if is_ordinary(input.text.chars().next().unwrap()) {
        let text = parse_identifier(input);
        if let Some(value) = parse_numeric(text) {
            Ok(value)
        } else {
            Ok(Value::Text(text.to_string()))
        }
    } else {
        Err(PDBError::new(
            ErrorLevel::BreakingError,
            "Invalid value",
            "A value should be \'.\', \'?\', a string (possibly enclosed), numeric or a multiline string (starting with \';\'), but here is an invalid character.",
            Context::position(&start),
        ))
    }
}

/// Parse a numeric value from a string which is expected to be of non zero length and not containing whitespace
fn parse_numeric(text: &str) -> Option<Value> {
    let mut chars_to_remove = 0;
    let first_char = text.chars().next().unwrap();
    // Parse a possible sign
    let mut minus = false;
    if first_char == '-' {
        minus = true;
        chars_to_remove += 1;
    } else if first_char == '+' {
        chars_to_remove += 1;
    }

    // Parse the integer part
    let mut integer_set = false;
    let mut value = 0;
    for c in text.chars().skip(chars_to_remove) {
        if let Some(num) = c.to_digit(10) {
            integer_set = true;
            value *= 10;
            value += num;
            chars_to_remove += 1;
        } else {
            break;
        }
    }

    // Now take the decimal part
    let mut decimal_set = false;
    let mut decimal = 0.0;
    if text.len() > chars_to_remove && text.chars().nth(chars_to_remove).unwrap() == '.' {
        chars_to_remove += 1;
        let mut power: f32 = 1.0;
        for c in text.chars().skip(chars_to_remove) {
            if c.is_ascii_digit() {
                decimal_set = true;
                power *= 10.0;
                decimal += c.to_digit(10).unwrap() as f32 / power;
                chars_to_remove += 1;
            } else {
                break;
            }
        }
    }

    // Now take the exponent
    let mut exponent_set = false;
    let mut exponent = 0;
    if text.len() > chars_to_remove {
        let next_char = text.chars().nth(chars_to_remove).unwrap();
        if next_char == 'e' || next_char == 'E' {
            // Parse a possible sign
            chars_to_remove += 1;
            if text.len() == chars_to_remove {
                return None; // No number after the exponent
            }
            let exp_first_char = text.chars().nth(chars_to_remove).unwrap();
            let mut exp_minus = false;
            if exp_first_char == '-' {
                exp_minus = true;
                chars_to_remove += 1;
            } else if exp_first_char == '+' {
                chars_to_remove += 1;
            }

            // Parse the integer part
            let mut exp_value = 0;
            for c in text.chars().skip(chars_to_remove) {
                if let Some(num) = c.to_digit(10) {
                    exponent_set = true;
                    exp_value *= 10;
                    exp_value += num;
                    chars_to_remove += 1;
                } else {
                    break;
                }
            }
            #[allow(clippy::cast_possible_wrap)]
            if exp_minus {
                exponent = -(exp_value as i32);
            } else {
                exponent = exp_value as i32;
            }
        }
    }

    // Take the uncertainty
    let mut uncertainty_set = false;
    let mut uncertainty = 0;
    if text.len() > chars_to_remove && text.chars().nth(chars_to_remove).unwrap() == '(' {
        uncertainty_set = true;
        chars_to_remove += 1;
        for c in text.chars().skip(chars_to_remove) {
            if let Some(num) = c.to_digit(10) {
                uncertainty *= 10;
                uncertainty += num;
                chars_to_remove += 1;
            } else {
                break;
            }
        }
        if text.len() == chars_to_remove || text.chars().nth(chars_to_remove).unwrap() != ')' {
            return None;
        }
        chars_to_remove += 1;
    }

    if (!integer_set && !decimal_set) || text.len() != chars_to_remove {
        None
    } else {
        let mut number = value as f32 + decimal;
        if minus {
            number *= -1.0;
        }
        if exponent_set {
            number *= 10_f32.powi(exponent);
        }
        if uncertainty_set {
            Some(Value::NumericWithUncertainty(number, uncertainty))
        } else {
            Some(Value::Numeric(number))
        }
    }
}

/// Parse an identifier, basically all chars until the next whitespace
fn parse_identifier<'a>(input: &mut Position<'a>) -> &'a str {
    let mut chars_to_remove = 0;

    for c in input.text.chars() {
        if c.is_ascii_whitespace() {
            let identifier = &input.text[..chars_to_remove];
            input.text = &input.text[chars_to_remove..];
            input.column += chars_to_remove;
            return identifier;
        }
        chars_to_remove += 1;
    }

    let identifier = input.text;
    input.text = "";
    input.column += chars_to_remove;
    identifier
}

/// Check if the input starts with the given pattern, it is case insensitive by
/// lowercasing the input string, so the pattern should be lowercase otherwise
/// it can never match.
fn start_with(input: &mut Position<'_>, pattern: &str) -> Option<()> {
    if input.text.len() < pattern.len() {
        None
    } else {
        for (p, c) in pattern.chars().zip(input.text.chars()) {
            if p != c.to_ascii_lowercase() {
                return None;
            }
        }
        input.text = &input.text[pattern.len()..];
        input.column += pattern.len();
        Some(())
    }
}

/// Trim all allowed whitespace (including comments)
fn trim_comments_and_whitespace(input: &mut Position<'_>) {
    loop {
        trim_whitespace(input);
        if input.text.is_empty() {
            return;
        }
        if input.text.starts_with('#') {
            skip_to_eol(input);
        } else {
            return;
        }
    }
}

/// Parse a piece of text enclosed by a char, it assumes the first position also matches the char.
/// It will fail if it finds a newline in the text. SO it can be used for single or double quoted strings.
fn parse_enclosed<'a>(input: &mut Position<'a>, pat: char) -> Result<&'a str, PDBError> {
    let mut chars_to_remove = 1; //Assume the first position is 'pat'

    for c in input.text.chars().skip(1) {
        if c == pat {
            let trimmed = &input.text[1..chars_to_remove];
            input.text = &input.text[(chars_to_remove + 1)..];
            input.column += chars_to_remove + 1;
            return Ok(trimmed);
        } else if c == '\n' || c == '\r' {
            let mut end = *input;
            end.text = &input.text[(chars_to_remove + 1)..];
            end.column += chars_to_remove + 1;
            return Err(PDBError::new(
                ErrorLevel::BreakingError,
                "Invalid enclosing",
                format!("This element was enclosed by \'{pat}\' but the closing delimiter was not found."),
                Context::range(input, &end),
            ));
        }
        chars_to_remove += 1;
    }

    let mut end = *input;
    end.text = &input.text[chars_to_remove..];
    end.column += chars_to_remove;
    Err(PDBError::new(
        ErrorLevel::BreakingError,
        "Invalid enclosing",
        format!("This element was enclosed by \'{pat}\' but the closing delimiter was not found."),
        Context::range(input, &end),
    ))
}

/// Parse a multiline string <eol>; ...(text)... <eol>;, it assumes the first position is ';'
fn parse_multiline_string<'a>(input: &mut Position<'a>) -> Result<&'a str, PDBError> {
    let mut chars_to_remove = 1; //Assume the first position is ';'
    let mut eol = false;
    let mut iter = input.text.chars().skip(1).peekable();

    while let Some(c) = iter.next() {
        if eol && c == ';' {
            let trimmed = &input.text[1..chars_to_remove];
            input.text = &input.text[(chars_to_remove + 1)..];
            input.column += 1;
            return Ok(trimmed);
        } else if c == '\n' {
            if let Some('\r') = iter.peek() {
                chars_to_remove += 1;
                let _ = iter.next();
            }
            input.line += 1;
            input.column = 1;
            chars_to_remove += 1;
            eol = true;
        } else if c == '\r' {
            if let Some('\n') = iter.peek() {
                chars_to_remove += 1;
                let _ = iter.next();
            }
            input.line += 1;
            input.column = 1;
            chars_to_remove += 1;
            eol = true;
        } else {
            chars_to_remove += 1;
            input.column += 1;
            eol = false;
        }
    }

    let mut end = *input;
    end.text = &input.text[chars_to_remove..];
    end.column += chars_to_remove;
    Err(PDBError::new(
        ErrorLevel::BreakingError,
        "Multiline string not finished",
        "A multiline string has to be finished by \'<eol>;\'",
        Context::range(input, &end),
    ))
}

/// Skip forward until the next eol, \r\n and \n\r are but consumed in full
fn skip_to_eol(input: &mut Position<'_>) {
    let mut chars_to_remove = 1;
    let mut iter = input.text.chars().skip(1).peekable();

    while let Some(c) = iter.next() {
        if c == '\n' {
            if let Some('\r') = iter.peek() {
                chars_to_remove += 1
            }
            input.line += 1;
            input.column = 1;
            chars_to_remove += 1;
            input.text = &input.text[chars_to_remove..];
            return;
        } else if c == '\r' {
            if let Some('\n') = iter.peek() {
                chars_to_remove += 1
            }
            input.line += 1;
            input.column = 1;
            chars_to_remove += 1;
            input.text = &input.text[chars_to_remove..];
            return;
        }
        chars_to_remove += 1;
    }

    input.text = "";
    input.column += chars_to_remove;
}

/// Trim all whitespace (<space, \t, <eol>) from the start of the string
fn trim_whitespace(input: &mut Position<'_>) {
    let mut chars_to_remove = 0;
    let mut iter = input.text.chars().peekable();

    while let Some(c) = iter.next() {
        if c == ' ' || c == '\t' {
            input.column += 1;
            chars_to_remove += 1;
        } else if c == '\n' {
            if let Some('\r') = iter.peek() {
                chars_to_remove += 1;
                let _ = iter.next();
            }
            input.line += 1;
            input.column = 1;
            chars_to_remove += 1;
        } else if c == '\r' {
            if let Some('\n') = iter.peek() {
                chars_to_remove += 1;
                let _ = iter.next();
            }
            input.line += 1;
            input.column = 1;
            chars_to_remove += 1;
        } else {
            input.text = &input.text[chars_to_remove..];
            return;
        }
    }
    input.text = &input.text[chars_to_remove..];
}

/// Test if the character is an ordinary character, one which can start a line in a multiline string
const fn is_ordinary(c: char) -> bool {
    match c {
        '#' | '$' | '\'' | '\"' | '_' | '[' | ']' | ';' | ' ' | '\t' => false,
        _ => c.is_ascii_graphic(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! assert_numeric {
        ($res:expr, $exp:expr) => {
            if let Some(Value::Numeric(n)) = $res {
                if !close(n, $exp) {
                    panic!("assertion failed: {} is not close to {}", n, $exp);
                }
            } else {
                panic!("assertion failed: {:?} is Err", $res);
            }
        };
        ($res:expr, $exp:expr, $un:expr) => {
            if let Some(Value::NumericWithUncertainty(n, u)) = $res {
                if !close(n, $exp) {
                    panic!("assertion failed: {} is not close to {}", n, $exp);
                }
                if u != $un {
                    panic!("assertion failed: {} is not equal to {}", u, $un);
                }
            } else {
                panic!("assertion failed: {:?} is Err", $res);
            }
        };
    }

    #[test]
    fn trim_whitespace_only_spaces() {
        let mut pos = Position {
            text: "    a",
            line: 1,
            column: 1,
        };
        trim_whitespace(&mut pos);
        assert_eq!(pos.text, "a");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 5);
    }

    #[test]
    fn trim_whitespace_tabs_and_spaces() {
        let mut pos = Position {
            text: " \t \t a",
            line: 1,
            column: 1,
        };
        trim_whitespace(&mut pos);
        assert_eq!(pos.text, "a");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 6);
    }

    #[test]
    fn trim_whitespace_newlines() {
        let mut pos = Position {
            text: " \t \t \n \r \n\r \r\na",
            line: 1,
            column: 1,
        };
        trim_whitespace(&mut pos);
        assert_eq!(pos.text, "a");
        assert_eq!(pos.line, 5);
        assert_eq!(pos.column, 1);
    }

    #[test]
    fn skip_to_eol_test_n() {
        let mut pos = Position {
            text: "bla bla bla\n\na",
            line: 1,
            column: 1,
        };
        skip_to_eol(&mut pos);
        assert_eq!(pos.text, "\na");
        assert_eq!(pos.line, 2);
        assert_eq!(pos.column, 1);
    }

    #[test]
    fn skip_to_eol_test_r() {
        let mut pos = Position {
            text: "bla bla bla\r\ra",
            line: 1,
            column: 1,
        };
        skip_to_eol(&mut pos);
        assert_eq!(pos.text, "\ra");
        assert_eq!(pos.line, 2);
        assert_eq!(pos.column, 1);
    }

    #[test]
    fn skip_to_eol_test_nr() {
        let mut pos = Position {
            text: "bla bla bla\n\r\n\ra",
            line: 1,
            column: 1,
        };
        skip_to_eol(&mut pos);
        assert_eq!(pos.text, "\n\ra");
        assert_eq!(pos.line, 2);
        assert_eq!(pos.column, 1);
    }

    #[test]
    fn skip_to_eol_test_rn() {
        let mut pos = Position {
            text: "bla bla bla\r\n\r\na",
            line: 1,
            column: 1,
        };
        skip_to_eol(&mut pos);
        assert_eq!(pos.text, "\r\na");
        assert_eq!(pos.line, 2);
        assert_eq!(pos.column, 1);
    }

    #[test]
    fn skip_to_eol_test_end() {
        let mut pos = Position {
            text: "bla bla bla",
            line: 1,
            column: 1,
        };
        skip_to_eol(&mut pos);
        assert_eq!(pos.text, "");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 12);
    }

    #[test]
    fn trim_comments_and_whitespace_test() {
        let mut pos = Position {
            text: "  \n#comment\n  #comment\na",
            line: 1,
            column: 1,
        };
        trim_comments_and_whitespace(&mut pos);
        assert_eq!(pos.text, "a");
        assert_eq!(pos.line, 4);
        assert_eq!(pos.column, 1);
    }

    #[test]
    fn start_with_true() {
        let mut pos = Position {
            text: "BloCk_a",
            line: 1,
            column: 1,
        };
        let res = start_with(&mut pos, "block_");
        assert!(res.is_some());
        assert_eq!(pos.text, "a");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 7);
    }

    #[test]
    fn start_with_false() {
        let mut pos = Position {
            text: "loop_a",
            line: 1,
            column: 1,
        };
        let res = start_with(&mut pos, "block_");
        assert!(res.is_none());
        assert_eq!(pos.text, "loop_a");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 1);
    }

    #[test]
    fn parse_identifier_test0() {
        let mut pos = Position {
            text: "hello_world a",
            line: 1,
            column: 1,
        };
        let res = parse_identifier(&mut pos);
        assert_eq!(res, "hello_world");
        assert_eq!(pos.text, " a");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 12);
    }

    #[test]
    fn parse_identifier_test1() {
        let mut pos = Position {
            text: " a",
            line: 1,
            column: 1,
        };
        let res = parse_identifier(&mut pos);
        assert_eq!(res, "");
        assert_eq!(pos.text, " a");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 1);
    }

    #[test]
    fn parse_numeric_integer() {
        let res = parse_numeric("42");
        assert_numeric!(res, 42.0);
    }

    #[test]
    fn parse_numeric_float_no_decimal() {
        let res = parse_numeric("42.");
        assert_numeric!(res, 42.0);
    }

    #[test]
    fn parse_numeric_float_no_integer() {
        let res = parse_numeric(".42");
        assert_numeric!(res, 0.42);
    }

    #[test]
    fn parse_numeric_float_exp() {
        let res = parse_numeric("42e1");
        assert_numeric!(res, 420.0);
    }

    #[test]
    fn parse_numeric_float_no_decimal_exp() {
        let res = parse_numeric("42.e10");
        assert_numeric!(res, 42.0 * 10000000000.0);
    }

    #[test]
    fn parse_numeric_float_no_integer_exp() {
        let res = parse_numeric(".42e10");
        assert_numeric!(res, 0.420 * 10000000000.0);
    }

    #[test]
    fn parse_numeric_float_no_integer_positive_exp() {
        let res = parse_numeric(".42e+10");
        assert_numeric!(res, 0.420 * 10000000000.0);
    }

    #[test]
    fn parse_numeric_float_no_integer_negative_exp() {
        let res = parse_numeric(".42e-10");
        assert_numeric!(res, 0.420 / 10000000000.0);
    }

    #[test]
    fn parse_numeric_float_uncertainty() {
        let res = parse_numeric("42.0(9)");
        assert_numeric!(res, 42.0, 9);
    }

    #[test]
    fn parse_numeric_float_uncertainty_missing_bracket() {
        let res = parse_numeric("42.0(9");
        assert!(res.is_none());
    }

    #[test]
    fn parse_numeric_float_huge_uncertainty() {
        let res = parse_numeric("42.0(97845)");
        assert_numeric!(res, 42.0, 97845);
    }

    #[test]
    fn parse_numeric_missing_numbers0() {
        let res = parse_numeric(".");
        assert!(res.is_none());
    }

    #[test]
    fn parse_numeric_missing_numbers1() {
        let res = parse_numeric(".e");
        assert!(res.is_none());
    }

    #[test]
    fn parse_numeric_missing_numbers2() {
        let res = parse_numeric(".e42");
        assert!(res.is_none());
    }

    #[test]
    fn parse_enclosed_test() {
        let mut pos = Position {
            text: "\"hello world\"hello",
            line: 1,
            column: 1,
        };
        let res = parse_enclosed(&mut pos, '\"');
        assert_eq!(res, Ok("hello world"));
        assert_eq!(pos.text, "hello");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 14);
    }

    #[test]
    fn parse_value_inapplicable() {
        let mut pos = Position {
            text: ".hello hello",
            line: 1,
            column: 1,
        };
        let res = parse_value(&mut pos);
        assert_eq!(res, Ok(Value::Inapplicable));
        assert_eq!(pos.text, "hello hello");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 2);
    }

    #[test]
    fn parse_value_unknown() {
        let mut pos = Position {
            text: "?hello hello",
            line: 1,
            column: 1,
        };
        let res = parse_value(&mut pos);
        assert_eq!(res, Ok(Value::Unknown));
        assert_eq!(pos.text, "hello hello");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 2);
    }

    #[test]
    fn parse_char_string_simple() {
        let mut pos = Position {
            text: "hello hello",
            line: 1,
            column: 1,
        };
        let res = parse_value(&mut pos);
        assert_eq!(res, Ok(Value::Text("hello".to_string())));
        assert_eq!(pos.text, " hello");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 6);
    }

    #[test]
    fn parse_char_string_single_quoted() {
        let mut pos = Position {
            text: "\'hello hello\' hello",
            line: 1,
            column: 1,
        };
        let res = parse_value(&mut pos);
        assert_eq!(res, Ok(Value::Text("hello hello".to_string())));
        assert_eq!(pos.text, " hello");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 14);
    }

    #[test]
    fn parse_char_string_single_missing_quote() {
        let mut pos = Position {
            text: "\'hello hello",
            line: 1,
            column: 1,
        };
        let res = parse_value(&mut pos);
        assert!(res.is_err());
        assert_eq!(res.unwrap_err().short_description(), "Invalid enclosing")
    }

    #[test]
    fn parse_char_string_double_quoted() {
        let mut pos = Position {
            text: "\"hello hello\" hello",
            line: 1,
            column: 1,
        };
        let res = parse_value(&mut pos);
        assert_eq!(res, Ok(Value::Text("hello hello".to_string())));
        assert_eq!(pos.text, " hello");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 14);
    }

    #[test]
    fn parse_char_string_double_missing_quote() {
        let mut pos = Position {
            text: "\"hello hello\n",
            line: 1,
            column: 1,
        };
        let res = parse_value(&mut pos);
        assert!(res.is_err());
        assert_eq!(res.unwrap_err().short_description(), "Invalid enclosing")
    }

    #[test]
    fn parse_char_string_invalid_quoted() {
        let mut pos = Position {
            text: "\"hello\nhello\"hello",
            line: 1,
            column: 1,
        };
        let res = parse_value(&mut pos);
        assert!(res.is_err());
        assert_eq!(pos.text, "\"hello\nhello\"hello");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 1);
    }

    #[test]
    fn parse_value_numeric() {
        let mut pos = Position {
            text: "56.8 hello",
            line: 1,
            column: 1,
        };
        let res = parse_value(&mut pos);
        assert_eq!(res, Ok(Value::Numeric(56.8)));
        assert_eq!(pos.text, " hello");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 5);
    }

    #[test]
    fn parse_value_multiline_text() {
        let mut pos = Position {
            text: ";\n\tthis is a comment of considerable length\n; hello",
            line: 1,
            column: 1,
        };
        let res = parse_value(&mut pos);
        assert_eq!(
            res,
            Ok(Value::Text(
                "\n\tthis is a comment of considerable length\n".to_string()
            ))
        );
        assert_eq!(pos.text, " hello");
        assert_eq!(pos.line, 3);
        assert_eq!(pos.column, 2);
    }

    #[test]
    fn parse_value_multiline_text_with_semicolon() {
        let mut pos = Position {
            text: ";\n\tthis is a tricky comment; of considerable length!\n; hello",
            line: 1,
            column: 1,
        };
        let res = parse_value(&mut pos);
        assert_eq!(
            res,
            Ok(Value::Text(
                "\n\tthis is a tricky comment; of considerable length!\n".to_string()
            ))
        );
        assert_eq!(pos.text, " hello");
        assert_eq!(pos.line, 3);
        assert_eq!(pos.column, 2);
    }
    #[test]
    fn parse_value_multiline_text_with_newlines() {
        let mut pos = Position {
            text: ";\n\tthis is\na tricky\rcomment\n\rof considerable\r\nlength!\n; hello",
            line: 1,
            column: 1,
        };
        let res = parse_value(&mut pos);
        assert_eq!(
            res,
            Ok(Value::Text(
                "\n\tthis is\na tricky\rcomment\n\rof considerable\r\nlength!\n".to_string()
            ))
        );
        assert_eq!(pos.text, " hello");
        assert_eq!(pos.line, 7);
        assert_eq!(pos.column, 2);
    }

    #[test]
    fn parse_value_multiline_text_missing_end() {
        let mut pos = Position {
            text: ";\n\tthis is\na tricky\rcomment\n\rof considerable\r\nlength!\n",
            line: 1,
            column: 1,
        };
        let res = parse_value(&mut pos);
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().short_description(),
            "Multiline string not finished"
        )
    }

    #[test]
    fn classify_char_test() {
        assert!(is_ordinary('a'));
        assert!(is_ordinary('!'));
        assert!(is_ordinary('h'));
        assert!(is_ordinary('%'));
        assert!(is_ordinary('~'));
        assert!(!is_ordinary(' '));
        assert!(!is_ordinary('\t'));
        assert!(!is_ordinary(';'));
        assert!(!is_ordinary('#'));
        assert!(!is_ordinary('\''));
        assert!(!is_ordinary('\"'));
        assert!(!is_ordinary('$'));
        assert!(!is_ordinary('_'));
        assert!(!is_ordinary('['));
        assert!(!is_ordinary(']'));
    }

    #[test]
    fn parse_data_single_item_numeric() {
        let mut pos = Position {
            text: "_tag\n42.3",
            line: 1,
            column: 1,
        };
        let res = parse_data_item(&mut pos);
        assert_eq!(
            res,
            Ok(DataItem::Single(Single {
                name: "tag".to_string(),
                content: Value::Numeric(42.3)
            }))
        );
        assert_eq!(pos.text, "");
        assert_eq!(pos.line, 2);
        assert_eq!(pos.column, 5);
    }

    #[test]
    fn parse_data_single_item_numeric_no_leading_zero() {
        let mut pos = Position {
            text: "_tag\n+.16",
            line: 1,
            column: 1,
        };
        let res = parse_data_item(&mut pos);
        assert_eq!(
            res,
            Ok(DataItem::Single(Single {
                name: "tag".to_string(),
                content: Value::Numeric(0.16)
            }))
        );
        assert_eq!(pos.text, "");
        assert_eq!(pos.line, 2);
        assert_eq!(pos.column, 5);
    }

    #[test]
    fn parse_data_single_item_string() {
        let mut pos = Position {
            text: "_tag\t\"of course I would\"",
            line: 1,
            column: 1,
        };
        let res = parse_data_item(&mut pos);
        assert_eq!(
            res,
            Ok(DataItem::Single(Single {
                name: "tag".to_string(),
                content: Value::Text("of course I would".to_string())
            }))
        );
        assert_eq!(pos.text, "");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 25);
    }

    #[test]
    fn parse_data_single_missing_value() {
        let mut pos = Position {
            text: "_tag\t",
            line: 1,
            column: 1,
        };
        let res = parse_data_item(&mut pos);
        assert!(res.is_err());
        assert_eq!(res.unwrap_err().short_description(), "No valid Value")
    }

    #[test]
    fn parse_data_single_item_multiline_string() {
        let mut pos = Position {
            text: "_long__tag\n;\tOf course I would\nAlso on multiple lines ;-)\n;",
            line: 1,
            column: 1,
        };
        let res = parse_data_item(&mut pos);
        assert_eq!(
            res,
            Ok(DataItem::Single(Single {
                name: "long__tag".to_string(),
                content: Value::Text(
                    "\tOf course I would\nAlso on multiple lines ;-)\n".to_string()
                )
            }))
        );
        assert_eq!(pos.text, "");
        assert_eq!(pos.line, 4);
        assert_eq!(pos.column, 2);
    }

    #[test]
    fn parse_data_item_loop() {
        let mut pos = Position {
            text: "loop_\n\t_first\n\t_second\n\t_last\n#Some comment because I need to put that in here as well!\n. 23.2 ?\nHello 25.9 ?\nHey 30.3 N",
            line: 1,
            column: 1,
        };
        let res = parse_data_item(&mut pos);
        assert_eq!(
            res,
            Ok(DataItem::Loop(Loop {
                header: vec![
                    "first".to_string(),
                    "second".to_string(),
                    "last".to_string()
                ],
                data: vec![
                    vec![Value::Inapplicable, Value::Numeric(23.2), Value::Unknown,],
                    vec![
                        Value::Text("Hello".to_string()),
                        Value::Numeric(25.9),
                        Value::Unknown,
                    ],
                    vec![
                        Value::Text("Hey".to_string()),
                        Value::Numeric(30.3),
                        Value::Text("N".to_string())
                    ]
                ]
            }))
        );
        assert_eq!(pos.text, "");
        assert_eq!(pos.line, 8);
        assert_eq!(pos.column, 11);
    }

    #[test]
    fn parse_invalid_data_item_loop() {
        let mut pos = Position {
            text: "loop_\n\t_first\n\t_second\n\t_last\n. 23.2 #?\nHello 25.9 ?\nHey 30.3 N",
            line: 1,
            column: 1,
        };
        let res = parse_data_item(&mut pos);
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().short_description(),
            "Loop has incorrect number of data items"
        )
    }

    #[test]
    fn parse_data_item_or_save_frame_data_item() {
        let mut pos = Position {
            text: "_data ?",
            line: 1,
            column: 1,
        };
        let res = parse_data_item_or_save_frame(&mut pos);
        assert_eq!(
            res,
            Ok(Item::DataItem(DataItem::Single(Single {
                name: "data".to_string(),
                content: Value::Unknown
            })))
        );
        assert_eq!(pos.text, "");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 8);
    }

    #[test]
    fn parse_data_item_or_save_frame_save_frame() {
        let mut pos = Position {
            text: "save_something_to_save _data . save_",
            line: 1,
            column: 1,
        };
        let res = parse_data_item_or_save_frame(&mut pos);
        assert_eq!(
            res,
            Ok(Item::SaveFrame(SaveFrame {
                name: "something_to_save".to_string(),
                items: vec![DataItem::Single(Single {
                    name: "data".to_string(),
                    content: Value::Inapplicable
                })]
            }))
        );
        assert_eq!(pos.text, "");
        assert_eq!(pos.line, 1);
        assert_eq!(pos.column, 37);
    }

    #[test]
    fn parse_invalid_save_frame() {
        let mut pos = Position {
            text: "save_something_to_save _data . safe_",
            line: 1,
            column: 1,
        };
        let res = parse_data_item_or_save_frame(&mut pos);
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().short_description(),
            "No matching \'save_\' found"
        )
    }

    #[test]
    fn parse_data_block_test() {
        let mut pos = Position {
            text: "data_1UBQ\n#\n_entry.id   1UBQ",
            line: 1,
            column: 1,
        };
        let res = parse_data_block(&mut pos);
        assert_eq!(
            res,
            Ok(DataBlock {
                name: "1UBQ".to_string(),
                items: vec![Item::DataItem(DataItem::Single(Single {
                    name: "entry.id".to_string(),
                    content: Value::Text("1UBQ".to_string())
                }))]
            })
        );
        assert_eq!(pos.text, "");
        assert_eq!(pos.line, 3);
        assert_eq!(pos.column, 17);
    }

    #[test]
    fn parse_invalid_data_block_test() {
        let mut pos = Position {
            text: "1UBQ\n#\n_entry.id   1UBQ",
            line: 1,
            column: 1,
        };
        let res = parse_data_block(&mut pos);
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().short_description(),
            "Data Block not opened"
        )
    }

    fn close(a: f32, b: f32) -> bool {
        let dif = a / b;
        (1.0 - dif) > -0.000000000000001 && (dif - 1.0) < 0.000000000000001
    }
}
