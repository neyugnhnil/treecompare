# treecompare/scripts/fasta.py
from typing import List, Tuple

class FastaError(Exception):
    pass

def read_fasta(path: str, *, strip_whitespace: bool = True) -> List[Tuple[str, str]]:
    """
    standard fasta reader that returns list of (header, sequence). raises FastaError on structural issues.
    """
    try:
        headers: List[str] = []
        seqs: List[str] = []
        name = None
        buf: List[str] = []

        with open(path, "r", encoding="utf-8") as fh:
            for i, line in enumerate(fh, 1):
                s = line.rstrip("\n")
                if not s:
                    continue
                if s.startswith(">"):
                    if name is not None:
                        seq = "".join(buf)
                        if strip_whitespace:
                            seq = seq.replace(" ", "").replace("\t", "")
                        if not seq:
                            raise FastaError(f"Empty sequence for header '{name}'")
                        headers.append(name)
                        seqs.append(seq)
                        buf = []
                    name = s[1:].strip()
                    if not name:
                        raise FastaError(f"Blank header at line {i}")
                else:
                    if name is None:
                        raise FastaError(f"Sequence line before any header at line {i}")
                    buf.append(s.strip() if strip_whitespace else s)
        if name is not None:
            seq = "".join(buf)
            if strip_whitespace:
                seq = seq.replace(" ", "").replace("\t", "")
            if not seq:
                raise FastaError(f"Empty sequence for header '{name}'")
            headers.append(name)
            seqs.append(seq)

        if not headers:
            raise FastaError("No sequences found in FASTA.")
        return list(zip(headers, seqs))
    except FileNotFoundError:
        raise FastaError(f"FASTA not found: {path}")

def read_fasta_aligned(path: str, *, min_distinct: int = 3) -> List[Tuple[str, str]]:
    """
    does read_fasta, also checks:
      - all sequences same length
      - at least `min_distinct` distinct sequences, default 3.
    """
    recs = read_fasta(path)
    lengths = {len(seq) for _, seq in recs}
    if len(lengths) != 1:
        raise FastaError("Sequences must be aligned (varying lengths found).")
    distinct = {seq for _, seq in recs}
    if len(distinct) < min_distinct:
        raise FastaError(f"Need at least {min_distinct} distinct sequences.")
    return recs
