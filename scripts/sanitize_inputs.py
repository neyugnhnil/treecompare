#!/usr/bin/env python3
import argparse, os, sys
from typing import Dict
from fasta import read_fasta


def sanitize_label(label: str, used: Dict[str, int]) -> str:
    base = "".join(ch if (ch.isalnum() or ch == "_") else "_" for ch in label)
    if not base:
        base = "TAXON"
    idx = used.get(base, 0)
    used[base] = idx + 1
    if idx == 0:
        return base
    return f"{base}_{idx}"

def write_fasta(records, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            fh.write(f"{seq}\n")

def sanitize_metadata(meta_path, out_path, mapping):
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    if not meta_path or not os.path.exists(meta_path):
        with open(out_path, "w", encoding="utf-8") as fh:
            fh.write("")
        return
    with open(meta_path, "r", encoding="utf-8") as src:
        lines = [ln.rstrip("\n") for ln in src]
    if not lines:
        open(out_path, "w", encoding="utf-8").close()
        return
    header = lines[0]
    body = []
    for ln in lines[1:]:
        if not ln:
            continue
        cols = ln.split("\t")
        token = mapping.get(cols[0])
        if not token:
            continue
        cols[0] = token
        body.append("\t".join(cols))
    with open(out_path, "w", encoding="utf-8") as out:
        out.write(header + "\n")
        for row in body:
            out.write(row + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gap", default="-", help="gap symbol in FASTA")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--map-out", required=True)
    ap.add_argument("--meta", default=None)
    ap.add_argument("--out-meta", required=True)
    args = ap.parse_args()

    try:
        records = read_fasta(args.fasta)
    except Exception as e:
        sys.stderr.write(f"[error] could not read FASTA: {e}\n")
        return 2

    used = {}
    sanitized_records = []
    name_map = {}
    gap_lower = (args.gap or "-").lower()[0]
    gap_upper = gap_lower.upper()

    for header, seq in records:
        token = sanitize_label(header, used)
        name_map[header] = token
        seq_upper = seq.upper().replace(gap_lower, "-").replace(gap_upper, "-")
        sanitized_records.append((token, seq_upper))

    write_fasta(sanitized_records, args.out_fasta)

    os.makedirs(os.path.dirname(os.path.abspath(args.map_out)) or ".", exist_ok=True)
    with open(args.map_out, "w", encoding="utf-8") as fh:
        for original, token in name_map.items():
            fh.write(f"{token}\t{original}\n")

    sanitize_metadata(args.meta, args.out_meta, name_map)
    return 0

if __name__ == "__main__":
    sys.exit(main())
