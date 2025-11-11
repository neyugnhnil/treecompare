#!/usr/bin/env python3
"""
validate_config.py 
gatekeeper before snakemake pipeline starts.
runs structural/semantic checks on user input, normalizes values.

- uses scripts/fasta.py to read & validate the aligned FASTA
- hard errors causes exit
- soft issues collected as warnings, printed once, user can continue or abort
- fills defaults and coerces values so downstream steps can assume sanity for the most part

"""

import argparse, os, sys, yaml
from collections import Counter
from fasta import read_fasta_aligned, FastaError

DNA_ALLOWED = set("ACGTRYMKSWHBVDN")
AA_ALLOWED = set("ARNDCEQGHILKMFPSTWYVX")

# dicts for allowed params and default fallback values
ALLOWED_METHODS = {"UPGMA", "NJ", "MP", "ML", "BI"}
ALLOWED_KEYS_BY_METHOD = {
    "UPGMA": {"name", "method", "model"},
    "NJ": {"name", "method", "model"},
    "MP": {"name", "method", "starter"},
    "ML": {"name", "method", "model", "starter"},
    "BI": {"name", "method", "model", "starter", "burninfrac"},
}
DEFAULTS = {
    "UPGMA": {"model": {"dna": "JC69", "aa": "JTT"}},
    "NJ": {"model": {"dna": "JC69", "aa": "JTT"}},
    "MP": {"starter": "random"},
    "ML": {"model": {"dna": "GTR+F+R4", "aa": "LG+F+R4"}, "starter": "random"},
    "BI": {"model": {"dna": "GTR+G", "aa": "LG"}, "starter": "random", "burninfrac": 0.1},
}

# function to coerce yes/no input
def as_bool(x, default=False):
    '''coerces user booleans, returning (value, warning)
    warning exists only if user input is invalid'''
    if isinstance(x, bool):
        return x, None
    
    if isinstance(x, str):
        t = x.strip().lower()
        if t in {"true", "t", "yes", "y", "1"}:  return True, None
        if t in {"false", "f", "no", "n", "0"}:  return False, None
    
    return default, f"value '{x}' invalid; using default {default}"

# get default value function
def default_for(datatype, method, key):
    '''fetches datatype-specific defaults from DEFAULTS dict'''
    v = DEFAULTS.get(method, {}).get(key)
    
    if isinstance(v, dict):
        return v.get(datatype)
   
    return v

# substitution model validator, returns default model if invalid
def normalized_model(model, datatype, method):
    if not isinstance(model, str):
        return None, None
    cleaned = model.strip().upper()
    warning = None
    if method in {"UPGMA", "NJ"}:
        allowed = {"JC69", "F81", "K80", "HKY", "TN93", "GTR"} if datatype == "dna" else {"JTT", "LG", "WAG", "DAYHOFF"}
        if cleaned not in allowed:
            warning = f"model '{model}' invalid for {method}; using default"
            cleaned = default_for(datatype, method, "model")
    elif method == "ML":
        allowed = {"GTR+F+R4", "GTR+G", "GTR+I+G", "LG+F+R4", "LG+G", "LG+I+G"}
        if cleaned not in allowed:
            warning = f"model '{model}' invalid for {method}; using default"
            cleaned = default_for(datatype, method, "model")
    elif method == "BI":
        allowed = {"GTR+G", "GTR+I+G", "LG", "LG+G"}
        if cleaned not in allowed:
            warning = f"model '{model}' invalid for {method}; using default"
            cleaned = default_for(datatype, method, "model")
    return cleaned, warning

# metafile validation
def validate_metafile(path, fasta_headers_set, warnings):
    """soft-validate metafile for mantel use."""
    
    if not isinstance(path, str) or not path:
        warnings.append("mantel enabled but no metafile provided; mantel tests will not be performed.")
        return False
    
    if not os.path.exists(path):
        warnings.append(f"mantel enabled but metafile not found: {path}; mantel tests will not be performed.")
        return False
    
    try:
        with open(path, "r", encoding="utf-8") as fh:
            lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
        if not lines:
            warnings.append(f"metafile is empty: {path}; mantel tests will not be performed.")
            return False
       
        header = lines[0].split("\t")
        if len(header) < 2:
            warnings.append(f"metafile must be TSV with header and ≥2 columns: {path}; mantel tests will not be performed.")
            return False
        
        first_col = [ln.split("\t")[0] for ln in lines[1:]]
        missing = [v for v in first_col if v not in fasta_headers_set]
        if missing:
            warnings.append(f"metafile first column has values not in FASTA headers; mantel tests will not be performed.")
            return False
        
        return True # everything checks out

    except Exception as e:
        warnings.append(f"could not parse metafile ({e}); mantel tests will not be performed.")
        return False

# sequence validation
def normalize_gap_symbol(value, warnings):
    default = "-"
    if isinstance(value, str):
        trimmed = value.strip()
        if len(trimmed) == 1:
            return trimmed
    warnings.append(f"[data.gap_symbol] invalid value '{value}'; using default '-'")
    return default

def ensure_allowed_characters(records, datatype, gap_symbol):
    gap = gap_symbol.upper()
    allowed = DNA_ALLOWED if datatype == "dna" else AA_ALLOWED
    allowed_with_gap = allowed | {gap}
    for header, seq in records:
        seq_up = seq.upper()
        invalid = sorted({ch for ch in seq_up if ch not in allowed_with_gap})
        if invalid:
            allowed_str = "".join(sorted(allowed_with_gap))
            sys.stderr.write(
                f"[error] sequence '{header}' contains invalid characters for {datatype.upper()} data: {', '.join(invalid)} "
                f"(allowed: {allowed_str})\n"
            )
            sys.exit(2)

# function that writes sanitized config if necessary once user accepts warnings
def write_sanitized_config(cfg, config_path):
    project_dir = os.path.dirname(os.path.abspath(config_path))
    resources_dir = os.path.join(project_dir, "resources")
    os.makedirs(resources_dir, exist_ok=True)
    sanitized_path = os.path.join(resources_dir, "config.sanitized.yaml")
    with open(sanitized_path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    rel_path = os.path.relpath(sanitized_path, project_dir)
    print(f"\n[sanitize] wrote sanitized config to {rel_path}\n")


def main():
    ###### CLI args parsing #############################################
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="path to config.yaml")
    ap.add_argument("--yes", action="store_true", help="auto-continue on soft warnings")
    args = ap.parse_args()

    warnings = []

    ########## load YAML config ###########################################
    try:
        with open(args.config, "r", encoding="utf-8") as fh:
            cfg = yaml.safe_load(fh) or {}
    except Exception as e:
        sys.stderr.write(f"[error] could not read config: {e}\n")
        sys.exit(2)

    ###### seed ##########################################################
    # attempts to coerce runtime.seed from config into positive int
    # invalid seed -> fallback to None

    seed_raw = cfg.get("runtime", {}).get("seed", None)
    seed = None # this means no seed enforced

    if seed_raw is None:
        seed = None  
    else:
        try:
            seed = int(seed_raw)
            if seed <= 0:
                raise ValueError
        except Exception:
            warnings.append(f"[runtime.seed] '{seed_raw}' invalid; no seed enforced")
            seed = None

    ################## "data" block checks #####################################
    data = cfg.get("data", {})

    # fasta path is present?
    fasta_path = data.get("fasta")
    if not isinstance(fasta_path, str) or not fasta_path:
        sys.stderr.write("[error] data.fasta must be a path to an aligned FASTA\n")
        sys.exit(2)

    # read aligned FASTA, checks equal length, 3 distinct taxa 
    gap_symbol = normalize_gap_symbol(data.get("gap_symbol", "-"), warnings)
    data["gap_symbol"] = gap_symbol

    try:
        recs = read_fasta_aligned(fasta_path, min_distinct=3)
    except FastaError as e:
        sys.stderr.write(f"[error] invalid FASTA: {e}\n")
        sys.exit(2)
    fasta_headers = [h for (h, _) in recs]
    fasta_headers_set = set(fasta_headers)
    if len(fasta_headers) < 4:
        warnings.append("FASTA has fewer than 4 taxa. BI (MrBayes) requires 4. BI trees will be skipped.")
        bi_allowed = False
    else:
        bi_allowed = True

    # data type hard check (dna/aa) 
    datatype = data.get("datatype")
    if not isinstance(datatype, str) or datatype.lower() not in {"dna", "aa"}:
        sys.stderr.write("[error] data.datatype must be 'dna' or 'aa'\n")
        sys.exit(2)
    datatype = datatype.lower()
    ensure_allowed_characters(recs, datatype, gap_symbol)

    # post-alignment qc yes/no validation 
    pa_qc = data.get("post_alignment_qc", False)
    val, note = as_bool(pa_qc, False)
    if note:
        warnings.append(f"[data.post_alignment_qc] {note}")
    data["post_alignment_qc"] = val  # normalize in "data" dict for later

    # mantel soft checks 
    mantel = cfg.get("mantel", {})
    mant_enabled, note = as_bool(mantel.get("enabled", False), False)
    if note:
        warnings.append(f"[mantel.enabled] {note}")
    method = mantel.get("method", "spearman")
    if method not in {"spearman", "pearson"}:
        warnings.append(f"[mantel.method] '{method}' invalid; using 'spearman'")
        method = "spearman"
    perms_raw = mantel.get("permutations", 1000) # permutations must be positive int
    try:
        perms = int(perms_raw)
        if perms <= 0:
            raise ValueError
    except Exception:
        warnings.append(f"[mantel.permutations] '{perms_raw}' invalid; using 1000")
        perms = 1000
    if mant_enabled:
        validate_metafile(data.get("metafile"), fasta_headers_set, warnings) 
        # failures -> warning that mantel won't be performed

    ################## "trees" block checks ###################################
    trees = cfg.get("trees", []) # pulls list of tree blocks as list of dicts
    
    # tree names 
    # keep only dict entries with a non-empty name; warn on skips
    valid_entries = []
    for t in trees:
        if not isinstance(t, dict):
            warnings.append("[trees] non-dict entry skipped"); continue
        nm = t.get("name")
        if not isinstance(nm, str) or not nm.strip():
            warnings.append("[trees] entry without a valid 'name' skipped"); continue
        valid_entries.append(t)
    # duplicate names not allowed
    names = [t["name"] for t in valid_entries]
    dups = [n for n, c in Counter(names).items() if c > 1]
    if dups:
        sys.stderr.write(f"[error] duplicate tree names are not allowed: {', '.join(dups)}\n"); sys.exit(2)
    # only trees with non-duplicate names remain in trees dict

    # tree params
    validated, seen = [], set()
    for t in valid_entries:
        name = t["name"]
        
        # invalid method names -> skip
        method_name = t.get("method")
        if method_name not in ALLOWED_METHODS:
            warnings.append(f"[{name}] method '{method_name}' invalid; this tree will be skipped"); continue
        if method_name == "BI" and not bi_allowed:
            warnings.append(f"[{name}] BI requires ≥4 taxa; skipping this tree because only {len(fasta_headers)} were found.")
            continue

        # start saving tree in vt dict if name and method are valid so far
        vt = {"name": name, "method": method_name}

        # irrelevant keys ignored
        allowed = ALLOWED_KEYS_BY_METHOD[method_name]
        extras = set(t.keys()) - allowed
        if extras:
            warnings.append(f"[{name}] ignoring keys not applicable to {method_name}: {sorted(extras)}")

        # missing model name (when applicable) goes to default
        if "model" in allowed:
            model = t.get("model")
            if not isinstance(model, str) or not model.strip():
                model = default_for(datatype, method_name, "model")
                warnings.append(f"[{name}] model missing; using default '{model}' for {datatype.upper()}")
            else:
                normalized, warn = normalized_model(model, datatype, method_name)
                if warn:
                    warnings.append(f"[{name}] {warn} '{default_for(datatype, method_name, 'model')}' for {datatype.upper()}")
                    model = normalized
                else:
                    model = normalized
            vt["model"] = model

        # starter is "random" or a previously seen tree, default to random if invalid
        if "starter" in allowed:
            raw = t.get("starter", default_for(datatype, method_name, "starter")) # default if doesnt exist
            st = raw.strip() if isinstance(raw, str) else ""
            if st.lower() == "random":
                vt["starter"] = "random"
            elif st in seen:
                vt["starter"] = st
            else:
                if st:
                    warnings.append(f"[{name}] starter '{st}' not found among previously defined trees; using 'random'")
                vt["starter"] = "random"

        # burninfrac (BI only), must be in (0,1)
        if "burninfrac" in allowed:
            raw_bf = t.get("burninfrac", DEFAULTS["BI"]["burninfrac"])
            try:
                bf = float(raw_bf)
                if not (0.0 < bf < 1.0):
                    raise ValueError
            except Exception:
                bf = DEFAULTS["BI"]["burninfrac"]
                warnings.append(f"[{name}] burninfrac '{raw_bf}' invalid; using {bf}")
            vt["burninfrac"] = bf

        validated.append(vt); seen.add(name)
    # hard exit if no trees remain
    if not validated:
        sys.stderr.write("[error] no tree requests survived validation!\n"); sys.exit(2)
    
    cfg["trees"] = validated

    ######### end: show warnings once (quiet if none) ###################################
    if warnings:
        print("\n!!!validation warnings!!!\n")
        for w in warnings:
            print(f"- {w}")
        if args.yes:
            print("--yes provided; continuing despite warnings.")
        else:
            try:
                resp = input("continue anyways? [y/n]: ").strip().lower()
            except EOFError:
                resp = ""
            if resp not in {"y", "yes"}:
                print("aborting at user request.")
                sys.exit(1)

    write_sanitized_config(cfg, args.config)

    # success
    return 0

if __name__ == "__main__":
    sys.exit(main())
