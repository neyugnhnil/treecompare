#!/usr/bin/env python3
"""
write_mrbayes_nexus.py
create a MrBayes-ready NEXUS written to --out-nexus
"""

import argparse, os, sys, re
from typing import List, Tuple
from fasta import read_fasta_aligned, FastaError


def parse_rate_flags(model: str) -> tuple[bool, bool]:
    up = (model or "").upper()
    return ("+G" in up or "+R" in up), ("+I" in up)


def map_dna_model_to_mrbayes(model: str) -> list[str]:
    base = (model or "").upper().split("+")[0]
    use_gamma, use_invar = parse_rate_flags(model)
    if base == "JC69":
        cmds = ["lset nst=1 rates=equal;", "prset statefreqpr=fixed(equal);"]
    elif base == "F81":
        cmds = ["lset nst=1 rates=equal;"]
    elif base in ("K80", "K2P"):
        cmds = ["lset nst=2 rates=equal;", "prset statefreqpr=fixed(equal);"]
    elif base == "HKY":
        cmds = ["lset nst=2 rates=equal;"]
    else:
        cmds = ["lset nst=6 rates=equal;"]
    if use_gamma and use_invar:
        cmds[0] = cmds[0].replace("rates=equal", "rates=invgamma")
    elif use_gamma:
        cmds[0] = cmds[0].replace("rates=equal", "rates=gamma")
    elif use_invar:
        cmds[0] = cmds[0].replace("rates=equal", "rates=propinv")
    return cmds


def map_aa_model_to_mrbayes(model: str) -> list[str]:
    base = (model or "").upper().split("+")[0]
    use_gamma, use_invar = parse_rate_flags(model)
    token = {
        "JTT": "jones",
        "WAG": "wag",
        "LG": "lg",
        "DAYHOFF": "dayhoff",
    }.get(base, "lg")
    cmds = ["lset rates=equal;", f"prset aamodelpr=fixed({token});"]
    if use_gamma and use_invar:
        cmds[0] = "lset rates=invgamma;"
    elif use_gamma:
        cmds[0] = "lset rates=gamma;"
    elif use_invar:
        cmds[0] = "lset rates=propinv;"
    return cmds


MIN_BRANCH = "0.000001"


def _clip_negative_branches(newick: str) -> str:
    return re.sub(r":-([0-9\\.eE+-]+)", f":{MIN_BRANCH}", newick)


def write_nexus(out_path: str, records: List[Tuple[str, str]],
                datatype: str, model: str, starter: str, gap_symbol: str):
    ntax = len(records)
    nchar = len(records[0][1])
    datatype_token = "dna" if datatype.lower() == "dna" else "protein"
    model_lines = map_dna_model_to_mrbayes(model) if datatype_token == "dna" else map_aa_model_to_mrbayes(model)

    start_clause = "mcmcp starttree=random;"
    trees_block = ""
    if isinstance(starter, str) and starter.lower() != "random" and os.path.exists(starter):
        with open(starter, "r", encoding="utf-8") as fh:
            newick = fh.read().strip()
        if newick:
            newick = _clip_negative_branches(newick)
            trees_block = (
                "begin trees;\n"
                "  tree usertree = " + newick.rstrip(" ;") + ";\n"
                "end;\n\n"
            )
            start_clause = "mcmcp starttree=user;"

    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("#NEXUS\n\n")
        out.write("begin data;\n")
        out.write(f"  dimensions ntax={ntax} nchar={nchar};\n")
        out.write(f"  format datatype={datatype_token} gap={gap_symbol} missing=?;\n")
        out.write("  matrix\n")
        for name, seq in records:
            seq_norm = seq.replace(" ", "").replace("\t", "").replace("\r", "").replace("\n", "")
            seq_norm = seq_norm.upper()
            out.write(f"  {name}    {seq_norm}\n")
        out.write("  ;\n")
        out.write("end;\n\n")
        if trees_block:
            out.write(trees_block)
        out.write("begin mrbayes;\n")
        out.write("  set autoclose=yes nowarn=yes;\n")
        for line in model_lines:
            out.write("  " + line + "\n")
        out.write("  " + start_clause + "\n")
        out.write("end;\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--datatype", required=True, choices=["dna", "aa"])
    ap.add_argument("--model", required=True)
    ap.add_argument("--starter", required=True)
    ap.add_argument("--out-nexus", required=True)
    ap.add_argument("--gap", required=True)
    args = ap.parse_args()

    try:
        recs = read_fasta_aligned(args.fasta, min_distinct=3)
    except FastaError as e:
        sys.stderr.write(f"[error] {e}\n")
        return 2

    write_nexus(args.out_nexus, recs, args.datatype, args.model, args.starter, args.gap)
    return 0

if __name__ == "__main__":
    sys.exit(main())
