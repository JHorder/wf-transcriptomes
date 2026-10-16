"""Check that the DGE gene expression prefilter parameters reach the DE process."""

import json
from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[4]
FILTER_PARAMS = ("de_min_count", "de_min_samples")


def test_de_gene_filter_params_default_to_disabled():
    """Schema and config declare non-negative integers defaulting to 0 (off)."""
    schema = json.loads((REPO_ROOT / "nextflow_schema.json").read_text())
    properties = schema["definitions"]["de_analysis_options"]["properties"]
    config = (REPO_ROOT / "nextflow.config").read_text()

    for name in FILTER_PARAMS:
        assert properties[name]["type"] == "integer"
        assert properties[name]["default"] == 0
        assert properties[name]["minimum"] == 0
        assert re.search(rf"^\s*{name}\s*=\s*0\s*$", config, re.MULTILINE)


def test_de_gene_filter_params_passed_to_de_process():
    """The DE process passes both parameters to supeRglue de_analysis."""
    process = (REPO_ROOT / "subworkflows" / "differential_expression.nf").read_text()
    r_spec = (
        REPO_ROOT / "bin" / "workflow_glue_r" / "R" / "de_analysis.R"
    ).read_text()

    for name in FILTER_PARAMS:
        assert f"--{name} ${{params.{name}}}" in process
        assert f'flag = "--{name}"' in r_spec
