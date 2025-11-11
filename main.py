#!/usr/bin/env python3
import argparse, os, sys, subprocess

def main():

    # arg parsing
    ap = argparse.ArgumentParser()
    ap.add_argument("project_dir", help="project directory containing config.yaml")
    ap.add_argument("--snakemake", default="snakemake", help="snakemake executable")
    ap.add_argument("--snakefile", default="Snakefile", help="path to Snakefile")
    ap.add_argument("--yes", action="store_true", help="auto-continue when validate_config.py reports warnings")
    args, snakemake_extra = ap.parse_known_args() # everything unknown lands in snakemake_extra

    # define repo root, project directory, and config file
    repo_root = os.path.dirname(os.path.abspath(__file__))
    proj = os.path.abspath(args.project_dir)
    cfg = os.path.join(proj, "config.yaml")
    if not os.path.exists(cfg):
        sys.exit(f"[error] config.yaml not found at: {cfg}")

    # validate_config.py, outcomes: hard exit, soft warnings + prompt, or nothing
    validate_script = os.path.join(repo_root, "scripts", "validate_config.py")
    validate_cmd = [sys.executable, validate_script, "--config", cfg]
    if args.yes:
        validate_cmd.append("--yes")
    rc = subprocess.run(validate_cmd, cwd=proj).returncode
    if rc != 0:
        sys.exit(rc)

    # prefer sanitized config in /resources if exists
    sanitized_rel = os.path.join("resources", "config.sanitized.yaml")
    sanitized_abs = os.path.join(proj, sanitized_rel)
    configfile_arg = sanitized_rel if os.path.exists(sanitized_abs) else "config.yaml"

    # run snakemake DAG
    default_cores = str(max(1, (os.cpu_count() or 1)))
    snakemake_cmd = [
        args.snakemake,
        "-s", os.path.join(os.path.dirname(__file__), args.snakefile),
        "--directory", proj,
        "--configfile", configfile_arg,
        "--printshellcmds",
        "--use-conda",
        "-j", default_cores,
    ]

    # native snakemake options
    if snakemake_extra:
        extra = snakemake_extra
        if extra[0] == "--": # some people do this
            extra = extra[1:]
        snakemake_cmd.extend(extra)

    env = os.environ.copy()
    env["TREECOMPARE_CONFIG"] = os.path.abspath(os.path.join(proj, configfile_arg))
    rc = subprocess.run(snakemake_cmd, cwd=proj, env=env).returncode
    sys.exit(rc)

if __name__ == "__main__":
    main()
