#!/usr/bin/env python3
"""
dispatcher.py
snakemake should have passed --config, --tree-name, --fasta, --out-tree and --out-stats
- reads a single tree spec (by --tree-name) from config.yaml
- dispatches to the appropriate backend tools with needed specs
- ensures outputs land at output paths
- creates placeholders if trees/stats could not be generated to not break downstream methods
"""

import argparse, json
import os, subprocess, sys
import time
import yaml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def run(cmd):
    '''run a command as a subprocess'''
    subprocess.run(cmd, check=True)

def load_tree_spec(cfg_path, tree_name):
    '''reads config.yaml, finds the tree entry matching --tree-name, 
    and returns both the overall config dict and that tree's spec bloc'''
    with open(cfg_path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}
    trees = [t for t in cfg.get("trees", []) if isinstance(t, dict)]
    for t in trees:
        if t.get("name") == tree_name:
            return cfg, t
    sys.exit(f"[error] tree '{tree_name}' not found in {cfg_path}")

def starter_path_or_random(starter):
    '''converts a starter name to results/<starter>/tree.nwk, unless it’s “random”.'''
    if isinstance(starter, str) and starter.lower().strip() == "random":
        return "random"
    # starter references a previously built tree's Newick file
    return os.path.join("results", starter, "tree.nwk")

def maybe_write_minimal_stats(path, payload):
    """dumps payload if it can't find path."""
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)

def main():

    ############## prep ############################################
    # parse CLI args #
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--tree-name", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out-tree", required=True)
    ap.add_argument("--out-stats", required=True)
    ap.add_argument("--seed", type=int, default=None)
    args = ap.parse_args()

    # extract and prep specs 
    cfg, spec = load_tree_spec(args.config, args.tree_name)

    data_block = cfg.get("data", {})
    datatype = (data_block.get("datatype") or "").lower()
    method = spec["method"]
    model = spec.get("model")
    starter = spec.get("starter", "random")
    burninfrac = spec.get("burninfrac")
    starter_arg = starter_path_or_random(starter)

    # does out dir exist 
    os.makedirs(os.path.dirname(args.out_tree), exist_ok=True)

    # base stats (will be used as placeholder) 
    t0 = time.time()
    base_stats = {
        "tree": args.tree_name,
        "method": method,
        "datatype": datatype,
        "model": model,
        "starter": starter if isinstance(starter, str) else None,
    }

    ############### branch on method ################################
    try:
        if method in ("UPGMA", "NJ"):
            run([
                "Rscript", os.path.join(SCRIPT_DIR, "build_dist_tree.R"),
                "--fasta", args.fasta,
                "--datatype", datatype,
                "--model", model,
                "--method", method,
                "--out-tree", args.out_tree,
                "--out-stats", args.out_stats,
            ])
            #seed doesn't matter, distance matrix is deterministic
        elif method == "MP":
            cmd = [
                "Rscript", os.path.join(SCRIPT_DIR, "build_mp.R"),
                "--fasta", args.fasta,
                "--datatype", datatype,
                "--starter", starter_arg,
                "--out-tree", args.out_tree,
                "--out-stats", args.out_stats,
            ]
            if args.seed is not None:
                cmd += ["--seed", str(args.seed)]
            run(cmd)
        elif method == "ML":
            cmd = [
                "bash", os.path.join(SCRIPT_DIR, "run_iqtree.sh"),
                "--fasta", args.fasta,
                "--datatype", datatype,
                "--model", model,
                "--starter", starter_arg,
                "--out-tree", args.out_tree,
                "--out-stats", args.out_stats,
            ]
            if args.seed is not None:
                cmd += ["--seed", str(args.seed)]
            run(cmd)
        elif method == "BI":
            result_dir = os.path.dirname(args.out_tree)
            nexus_path = os.path.join(result_dir, "input.nex")
            # write nexus
            run([
                sys.executable, os.path.join(SCRIPT_DIR, "write_mrbayes_nexus.py"),
                "--fasta", args.fasta,
                "--datatype", datatype,
                "--model", model,
                "--starter", starter_arg,
                "--gap", data_block.get("gap_symbol") or "-",
                "--out-nexus", nexus_path,
            ])
            # use nexus file for mrbayes
            mrbayes_cmd = [
                "bash", os.path.join(SCRIPT_DIR, "run_mrbayes.sh"),
                "--nexus", nexus_path,
                "--burninfrac", str(burninfrac if burninfrac is not None else 0.1),
                "--out-tree", args.out_tree,
                "--out-stats", args.out_stats,
            ]
            if args.seed is not None:
                mrbayes_cmd += ["--seed", str(args.seed)]
            run(mrbayes_cmd)
        else:
            sys.exit(f"[error] unsupported method: {method}") # hard exit if unsupported method

    # if any subprocess raises calledprocesserror
    except subprocess.CalledProcessError as e:
        sys.exit(e.returncode)
    # success

    # calculate elapsed time after method finishes.
    elapsed = time.time() - t0

    # if stats.json not created by method scripts, create placeholder.
    base_stats.update({
        "elapsed_sec": round(elapsed, 3),
        "out_tree": args.out_tree,
    })
    maybe_write_minimal_stats(args.out_stats, base_stats)

    return 0


if __name__ == "__main__":
    sys.exit(main())
