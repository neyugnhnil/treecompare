# Snakefile
# sets up global context

import os

CONFIG_PATH = os.environ.get("TREECOMPARE_CONFIG", "config.yaml")

configfile: CONFIG_PATH
PROJECT_CONFIG = CONFIG_PATH
from workflow.common import tree_names, tree_targets, QC_REPORT

REPO_DIR = workflow.basedir
SCRIPTS_DIR = os.path.join(REPO_DIR, "scripts")
ENV_YAML = os.path.join(REPO_DIR, "env.yaml")
QC_ENABLED = bool(config.get("data", {}).get("post_alignment_qc", False))
QC_TARGETS = [QC_REPORT] if QC_ENABLED else []

# sets final target of workflow (tree files, stats files, mantel files, html report)
rule all:
    input:
        tree_targets(config)
        + QC_TARGETS
        + ["results/report/index.html"]


# rule-sets
include: "workflow/sanitize.rules"
include: "workflow/qc.rules"
include: "workflow/trees.rules"
include: "workflow/mantel.rules"
include: "workflow/report.rules"
