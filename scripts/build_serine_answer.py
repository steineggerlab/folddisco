#!/usr/bin/env python3
"""Build a serine-protease answer set for a folddisco benchmark index.

The zinc-finger answer set (data/zinc_answer.tsv) is a UniProt keyword export.
This one is derived from sequence instead, so it needs no network access and,
more importantly, it is defined independently of anything folddisco computes --
a geometry-derived answer set would make the benchmark circular.

A protein is an answer if its sequence matches BOTH PROSITE chymotrypsin-family
active-site patterns, with the His site before the Ser site:

  PS00134  TRYPSIN_HIS  [LIVM]-[ST]-A-[STAG]-H-C
  PS00135  TRYPSIN_SER  [DNSTAGC]-[GSTAPIMVQH]-x(2)-G-[DE]-S-G-[GS]-[SAPHV]-
                        [LIVMFYWH]-[LIVMFYSTANQH]

Caveats, both of which make benchmark precision CONSERVATIVE rather than optimistic:
  * Only the S1 (chymotrypsin) clan is covered. A geometric Ser-His-Asp triad
    query legitimately also hits S8/S9/S28 peptidases and alpha/beta hydrolases;
    those count as false positives here. Config-to-config comparisons are unaffected.
  * Divergent family members, and catalytic domains split across an AlphaFold
    fragment boundary, are missed.
  * The output is specific to the index it was built from. Regenerate it for any
    other index.

Regenerate (about 35 s for the 23,391-structure human proteome):

  cargo run --release --example seqdump -- index/h_sapiens > h_sapiens_seq.tsv
  python3 scripts/build_serine_answer.py h_sapiens_seq.tsv > data/serine_answer.tsv

Output columns: uniprot_accession, family_tag, structure_name, his_site, ser_site,
sequence_length. Only column 0 is read by `folddisco benchmark --afdb-to-uniprot`;
the rest are there so a reviewer can audit any row by hand.
"""
import re
import sys

HIS = re.compile(r'[LIVM][ST]A[STAG]HC')
SER = re.compile(r'[DNSTAGC][GSTAPIMVQH]..G[DE]SG[GS][SAPHV][LIVMFYWH][LIVMFYSTANQH]')


def main(seq_tsv):
    rows = []
    for line in open(seq_tsv):
        name, _, seq = line.rstrip('\n').partition('\t')
        if not seq:
            continue
        his, ser = HIS.search(seq), SER.search(seq)
        # the chymotrypsin fold orders the triad His < Asp < Ser in sequence
        if his and ser and his.start() < ser.start():
            acc = name.split('-')[1] if name.count('-') >= 2 else name
            rows.append((acc, 'S1_serine_protease', name,
                         f'his@{his.start() + 1}', f'ser@{ser.start() + 1}', str(len(seq))))
    seen = set()
    for row in sorted(rows):
        if row[0] in seen:
            continue
        seen.add(row[0])
        print('\t'.join(row))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
