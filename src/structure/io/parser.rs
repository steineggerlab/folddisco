use crate::structure::atom::Atom;

pub fn parse_line(line: &String) -> Result<Atom, &str> {
    // Check line length is enough to index b_factor
    match line.len() {
        67..=80 => {}
        _ => return Err("Line length is not 80 characters long"),
    }
    // Parse line
    let x = line[30..38].trim().parse::<f32>();
    let y = line[38..46].trim().parse::<f32>();
    let z = line[46..54].trim().parse::<f32>();
    let atom_name = parse_atom(&line[12..16]);
    let atom_serial = line[6..11].trim().parse::<u64>();
    let res_name = parse_residue(&line[17..20]);
    let res_serial = line[22..26].trim().parse::<u64>();
    let chain = line[21..22].as_bytes()[0];
    let b_factor = line[60..66].trim().parse::<f32>();
    // let occupancy = &line[54..60]; // NOT USING occupancy yet

    // Check if all the parsing was successful
    match (
        x,
        y,
        z,
        atom_name,
        atom_serial,
        res_name,
        res_serial,
        b_factor,
    ) {
        (
            Ok(x),
            Ok(y),
            Ok(z),
            Ok(atom_name),
            Ok(atom_serial),
            Ok(res_name),
            Ok(res_serial),
            Ok(b_factor),
        ) => Ok(Atom::new(
            x,
            y,
            z,
            atom_name,
            atom_serial,
            res_name,
            res_serial,
            chain,
            b_factor,
        )),
        _ => Err("Error parsing line"),
    }
}

pub fn parse_atom(name: &str) -> Result<[u8; 4], &str> {
    let bytes = name.as_bytes();
    // Check atom name is 4 ASCII characters
    match bytes.len() {
        4 => Ok([bytes[0], bytes[1], bytes[2], bytes[3]]),
        _ => Err("Atom name is not 4 characters long"),
    }
}

pub fn parse_residue(name: &str) -> Result<[u8; 3], &str> {
    let bytes = name.as_bytes();
    // Check residue name is 3 ASCII characters
    match bytes.len() {
        3 => Ok([bytes[0], bytes[1], bytes[2]]),
        _ => Err("Residue name is not 3 characters long"),
    }
}

#[cfg(test)]
mod parser_tests {
    use super::*;

    #[test]
    fn test_parse_atom() {
        let atom_name = "CA  ";
        let atom_name_bytes = parse_atom(atom_name).unwrap();
        assert_eq!(atom_name_bytes, [67, 65, 32, 32]);
    }
    #[test]
    fn test_parse_atom_fail() {
        let atom_name = "CA";
        let atom_name_bytes = parse_atom(atom_name);
        assert!(atom_name_bytes.is_err());
        let atom_name = "CA   ";
        let atom_name_bytes = parse_atom(atom_name);
        assert!(atom_name_bytes.is_err());
    }

    #[test]
    fn test_parse_residue() {
        let res_name = "ALA";
        let res_name_bytes = parse_residue(res_name);
        assert!(res_name_bytes.is_ok());
        assert_eq!(res_name_bytes.unwrap(), [65, 76, 65]);
    }
    #[test]
    fn test_parse_residue_fail() {
        let res_name = "ALAN";
        let res_name_bytes = parse_residue(res_name);
        assert!(res_name_bytes.is_err());
        let res_name = "AL";
        let res_name_bytes = parse_residue(res_name);
        assert!(res_name_bytes.is_err());
    }

    #[test]
    fn test_parse_line_success() {
        let line =
            "ATOM      1  N   ALA A 340      -2.311   2.993 -33.448  1.00  6.00           N  "
                .to_string();
        let atom = parse_line(&line).unwrap();
        assert_eq!(atom.atom_name, [32, 78, 32, 32]); // N
        assert_eq!(atom.res_name, [65, 76, 65]); // ALA
        assert_eq!(atom.chain, 65); // A
        assert_eq!(atom.atom_serial, 1); // 1
        assert_eq!(atom.res_serial, 340); // 340
        assert_eq!(atom.x, -2.311);
        assert_eq!(atom.y, 2.993);
        assert_eq!(atom.z, -33.448);
        assert_eq!(atom.b_factor, 6.00);
    }

    #[test]
    fn test_parse_line_fail_length() {
        // Short line
        let line = "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1".to_string();
        let atom = parse_line(&line);
        assert!(atom.is_err());
        // Long line
        let line = "ATOM      1  N   ALA A   1               10.000  10.000  10.000  1.00  0.00           N  1".to_string();
        let atom = parse_line(&line);
        assert!(atom.is_err());
    }

    #[test]
    fn test_parse_line_float() {
        // Error in X
        let line =
            "ATOM      1  N   ALA A   1      1A.000  10.000  10.000  1.00  0.00           N  "
                .to_string();
        let atom = parse_line(&line);
        assert!(atom.is_err());
        // Error in Y
        let line =
            "ATOM      1  N   ALA A   1      10.000  1A.000  10.000  1.00  0.00           N  "
                .to_string();
        let atom = parse_line(&line);
        assert!(atom.is_err());
        // Error in Z
        let line =
            "ATOM      1  N   ALA A   1      10.000  10.000  1A.000  1.00  0.00           N  "
                .to_string();
        let atom = parse_line(&line);
        assert!(atom.is_err());
        // Error in B-factor
        let line =
            "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1A.00  0.00           N  "
                .to_string();
        let atom = parse_line(&line);
        assert!(atom.is_err());
    }

    #[test]
    fn test_parse_line_fail_serial() {
        // Error in residue serial
        let line =
            "ATOM      1  N   ALA A 3A0      -2.311   2.993 -33.448  1.00  6.00           N  "
                .to_string();
        let atom = parse_line(&line);
        assert!(atom.is_err());
        // Error in atom serial
        let line =
            "ATOM      K  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N  "
                .to_string();
        let atom = parse_line(&line);
        assert!(atom.is_err());
    }
}
