# workflow/common.py
# shared utilities/constants for Snakemake rules

SANITIZED_FASTA = "resources/alignment.sanitized.fasta"
SANITIZED_META = "resources/metafile.sanitized.tsv"
HEADER_MAP = "resources/header_map.tsv"
QC_REPORT = "results/qc/alignment_qc.txt"


def tree_names(cfg):
    trees = cfg.get("trees", [])
    return [t["name"] for t in trees]


def tree_targets(cfg):
    names = tree_names(cfg)
    nwk = [f"results/{n}/tree.nwk" for n in names]
    stats = [f"results/{n}/stats.json" for n in names]
    mantel = [f"results/{n}/mantel.tsv" for n in names]
    return nwk + stats + mantel


def tree_starter_inputs(cfg):
    mapping = {}
    seen = set()
    for tree in cfg.get("trees", []):
        if not isinstance(tree, dict):
            continue
        name = tree.get("name")
        if not name or name in seen:
            continue
        seen.add(name)
        starter = tree.get("starter")
        if isinstance(starter, str):
            starter_norm = starter.strip()
            if starter_norm and starter_norm.lower() != "random":
                mapping[name] = f"results/{starter_norm}/tree.nwk"
    return mapping
