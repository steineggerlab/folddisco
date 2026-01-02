// File: log.rs
// Created: 2023-12-05 21:05:10
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved

// Colored string for log
pub const INFO: &str = "\x1b[1;32m[INFO]\x1b[0m";
pub const FAIL: &str = "\x1b[1;31m[FAIL]\x1b[0m";
pub const WARN: &str = "\x1b[1;33m[WARN]\x1b[0m";
pub const DONE: &str = "\x1b[1;34m[DONE]\x1b[0m";

pub fn log_msg(prefix: &str, msg: &str) -> String { format!("{} {}", prefix, msg) }
pub fn print_log_msg(prefix: &str, msg: &str) { eprintln!("{}", log_msg(prefix, msg)); }

// Macro for measuring time for a function
// measure_time gets expression and optional verbose flag
// Usage: measure_time!(expr) or measure_time!(expr, verbose)
#[macro_export]
macro_rules! measure_time {
    ($x:expr) => {{
        let start = std::time::Instant::now();
        let result = $x;
        let duration = start.elapsed();
        // Print only the name of the function
        let mut fn_name = stringify!($x);
        fn_name = fn_name.split("(").collect::<Vec<&str>>()[0];
        // Don't print quotation marks around the function name
        eprintln!("\x1b[1;32m[INFO]\x1b[0m {}: {:?}", fn_name, duration);
        result
    }};
    ($x:expr, $verbose:expr) => {{
        if $verbose {
            let start = std::time::Instant::now();
            let result = $x;
            let duration = start.elapsed();
            // Print only the name of the function
            let mut fn_name = stringify!($x);
            fn_name = fn_name.split("(").collect::<Vec<&str>>()[0];
            // Don't print quotation marks around the function name
            eprintln!("\x1b[1;32m[INFO]\x1b[0m {}: {:?}", fn_name, duration);
            result
        } else {
            $x
        }
    }}
}

#[cfg(test)]
mod tests {
    use super::*;
    // Test log with colored prefix
    #[test]
    fn test_colored_log() {
        let msg = "Hello, world!";
        let info = log_msg(INFO, msg);
        let fail = log_msg(FAIL, msg);
        let warn = log_msg(WARN, msg);
        let done = log_msg(DONE, msg);

        assert_eq!(info, "\x1b[1;32m[INFO]\x1b[0m Hello, world!");
        assert_eq!(fail, "\x1b[1;31m[FAIL]\x1b[0m Hello, world!");
        assert_eq!(warn, "\x1b[1;33m[WARN]\x1b[0m Hello, world!");
        assert_eq!(done, "\x1b[1;34m[DONE]\x1b[0m Hello, world!");
        println!("{}", info);
        println!("{}", fail);
        println!("{}", warn);
        println!("{}", done);
    }
}
