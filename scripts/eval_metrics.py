#!/usr/bin/env python3
"""Ranking metrics for a folddisco result file.

`folddisco benchmark` reports the two metrics this project decides on - Sens@1FP
(with `--fp 1`) and F1 over the whole list (without it) - and this script
reproduces both from the same inputs, plus the depths in between that neither
call gives you.

It replicates `benchmark --afdb-to-uniprot` exactly: an identifier is the
basename with its structure extension removed, then the second dash-separated
field (`AF-P17538-F1-model_v4.pdb` -> `P17538`), and duplicates are dropped
keeping first occurrence, so the AlphaFold fragments of one protein collapse to
a single accession. Ranking order is the order of lines in the result file.

Reported, in this order:
  sens@1fp       tp@1fp / answer_len - the author's headline metric. Fragile: it
                 is decided by the rank of one false positive, and its deltas
                 keep their sign in only about half of 5 % answer-dropout
                 replicates. Report it; do not decide on it alone.
  f1             set F1 over the whole list, the decision metric. Has no rank
                 dependence at all, so a reordering that promotes real hits into
                 the top can leave it unchanged - pair it with tp@10fp/tp@100fp.
  precision, recall, tp/fp/fn over the whole list
  tp@kfp         true positives above the k-th false positive, k in 1..500
  tp/fp/prec@topN at fixed depths N in 100..2000
  ap             average precision: for every rank i holding a true positive,
                 accumulate precision@i, then divide by the size of the FULL
                 answer set. Kept because it is rank-aware and stable, but it is
                 not a metric this project uses.

Absolute values depend on the answer set as much as on the search - the same
ranked lists score F1 0.93 against one zinc-finger set and 0.59 against another,
with opposite signs between configurations - so always name the answer set with
the number.

Usage:
  python3 scripts/eval_metrics.py <result.tsv> <answer.tsv> <index.lookup>

Reproduce a documented row (matched 4-residue zinc query, --sensitive, F1 0.9641):

  IDX=index/h_sapiens_folddisco
  folddisco query -i $IDX -p query/1G2F.pdb -q F207,F212,F225,F229 -t 12 \
    --covered-node 3 --max-node 4 --rmsd 1.0 --per-structure --sensitive > result.tsv
  python3 scripts/eval_metrics.py result.tsv <zinc answers>.tsv $IDX.lookup
"""
import sys, json

def parse(x):
    x = x.split('/')[-1]
    for ext in ('.pdb','.cif','.fcz','.ent'):
        if x.endswith(ext): x = x[:-4]; break
        if x.endswith(ext+'.gz'): x = x[:-7]; break
    s = x.split('-')
    return s[1] if len(s) >= 2 else x

def read_col(path, col=0):
    out, seen = [], set()
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line: continue
            v = parse(line.split('\t')[col])
            if v not in seen:
                seen.add(v); out.append(v)
    return out

def tp_at_kfp(result, answer, k):
    tp = fp = 0
    for t in result:
        if t in answer: tp += 1
        else:
            fp += 1
            if fp >= k: break
    return tp

def main(res, ans, lookup):
    result = read_col(res)
    answer = set(read_col(ans))
    allids = set(read_col(lookup, 1))
    n_ans = len(answer)
    # the two metrics of record first
    tp_1fp = tp_at_kfp(result, answer, 1)
    tp_all = sum(1 for t in result if t in answer)
    prec = tp_all/len(result) if result else 0.0
    rec = tp_all/n_ans if n_ans else 0.0
    o = {'sens@1fp': round(tp_1fp/n_ans, 4) if n_ans else None,
         'f1': round(2*prec*rec/(prec+rec), 4) if prec+rec else 0.0,
         'hits_raw': sum(1 for _ in open(res)), 'result_len': len(result),
         'answer_len': n_ans, 'total_ids': len(allids)}
    # TP@kFP
    for k in (1,2,5,10,20,50,100,200,500):
        o[f'tp@{k}fp'] = tp_at_kfp(result, answer, k)
    # precision@N / FP@N
    for n in (100,200,500,1000,2000):
        head = result[:n]
        tp = sum(1 for t in head if t in answer)
        o[f'tp@top{n}'] = tp
        o[f'fp@top{n}'] = len(head)-tp
        o[f'prec@top{n}'] = round(tp/len(head), 4) if head else None
    o['tp_all'] = tp_all; o['fp_all'] = len(result)-tp_all
    o['fn_all'] = n_ans - tp_all
    o['precision'] = round(prec,4); o['recall'] = round(rec,4)
    # average precision over the ranked list (AP wrt full answer set)
    tp = 0; ap = 0.0
    for i,t in enumerate(result,1):
        if t in answer:
            tp += 1; ap += tp/i
    o['ap'] = round(ap/n_ans, 4)
    print(json.dumps(o))

if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
