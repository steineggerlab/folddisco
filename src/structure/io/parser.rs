use super::super::atom::*;

pub fn parse_line(line: &String) -> Atom {
    // ??? Error, None type handling ???
    let x = line[30..38].trim().parse::<f32>().unwrap();
    let y = line[38..46].trim().parse::<f32>().unwrap();
    let z= line[46..54].trim().parse::<f32>().unwrap();
    let atom_name  = parse_atom(&line[12..16]);
    let atom_serial = line[6..11].trim().parse::<u64>().unwrap();
    let res_name  = parse_resi(&line[17..20]) ;
    let res_serial = line[22..26].trim().parse::<u64>().unwrap();
    let chain = line[21..22].as_bytes()[0];
    let b_factor = line[60..66].trim().parse::<f32>().unwrap();
    let occupancy = &line[54..60] ;    

    Atom::build(x, y, z, atom_name, atom_serial, 
                    res_name, res_serial, chain, b_factor)
}

pub fn parse_atom(name:&str) -> [u8;4] {
    let bytes = name.as_bytes();
    [bytes[0], bytes[1], bytes[2],bytes[3]] 

}
pub fn parse_resi(name:&str) -> [u8;3] {
    let bytes = name.as_bytes();
    [bytes[0], bytes[1], bytes[2]]
}
