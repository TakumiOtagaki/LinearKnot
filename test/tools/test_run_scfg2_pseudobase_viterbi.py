#!/usr/bin/env python3
"""Parser/aggregation tests for the PseudoBase++ Viterbi runner."""

from __future__ import annotations

import os
import sys
import subprocess
import tempfile
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_scfg2_pseudobase_viterbi import (  # noqa: E402
    aggregate,
    evaluate_case,
    load_resumable_rows,
    parse_prediction_output,
    run_captured_command,
    select_case_ids,
    status_summary,
    validate_case,
    write_rows,
)


SYNTHETIC = """\
SUMMARY\tok=1\treachable=1\tlength=8\tbeam_width=20\tanchor_cap=2\thotspot_width=4\tpasses=1\tmax_results=2\trecords=2\tgroups=2\tproposal_seconds=0.1\tviterbi_seconds=0.2\ttotal_seconds=0.3\terror=
RECORD_HEADER\tindex\treachable\tscaffold_pairs\theuristic_score\tpass\torigin\tenergy\tseconds\tscaffold_G\tgenerated_Gprime\tunion_R
RECORD\t0\t1\t0\t0\t0\t0\t-1\t0.1\t........\t(....)..\t(....)..
RECORD\t1\t1\t1\t2\t0\t1\t-2\t0.1\t(......)\t.(....).\t((....))
GROUP_HEADER\tindex\trecord_count\tminimum_energy\tmaximum_energy\tenergy_spread\tselected_record\tselected
GROUP\t0\t1\t-1\t-1\t0\t0\t0
GROUP\t1\t1\t-2\t-2\t0\t1\t1
"""


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("SCFG2_BUILD_DIR", root / "build")).resolve()
    sequence = "UUAUAUGGGUAUUAUUCUCA"
    arguments = [
        sequence,
        str(root / "params/rna_Turner04.par"),
        "64", "2", "4", "1", "4", "none",
    ]
    sequential_output = subprocess.run(
        [str(build / "scfg2_hotspot_viterbi_tool"), *arguments],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout
    parallel_output = subprocess.run(
        [
            sys.executable,
            str(root / "test/tools/run_scfg2_parallel_hotspot_viterbi.py"),
            "--proposal-tool",
            str(build / "scfg2_hotspot_proposal_tool"),
            "--worker-tool",
            str(build / "scfg2_fixed_scaffold_viterbi_tool"),
            "--candidate-jobs", "2",
            *arguments,
        ],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout
    sequential = parse_prediction_output(sequential_output)
    parallel = parse_prediction_output(parallel_output)
    record_fields = (
        "reachable", "scaffold_pairs", "heuristic_score", "pass", "origin",
        "energy", "scaffold_G", "generated_Gprime", "union_R",
    )
    assert len(sequential["records"]) == len(parallel["records"])
    for index in sequential["records"]:
        assert {
            field: sequential["records"][index][field]
            for field in record_fields
        } == {
            field: parallel["records"][index][field]
            for field in record_fields
        }
    group_fields = (
        "record_count", "minimum_energy", "maximum_energy", "energy_spread",
        "selected_record", "selected",
    )
    assert sequential["summary"]["case3_truncated_records"] == "0"
    assert parallel["summary"]["case3_truncated_records"] == "0"
    assert len(sequential["groups"]) == len(parallel["groups"])
    for index in sequential["groups"]:
        assert {
            field: sequential["groups"][index][field]
            for field in group_fields
        } == {
            field: parallel["groups"][index][field]
            for field in group_fields
        }

    with tempfile.TemporaryDirectory() as tmpdir:
        fake_linearfold = Path(tmpdir) / "linearfold_predict"
        fake_linearfold.write_text(
            "#!/usr/bin/env python3\n"
            "print('((................))')\n"
            "print('-1.25')\n",
            encoding="utf-8",
        )
        fake_linearfold.chmod(0o755)
        augmented_output = subprocess.run(
            [
                sys.executable,
                str(root / "test/tools/run_scfg2_parallel_hotspot_viterbi.py"),
                "--proposal-tool",
                str(build / "scfg2_hotspot_proposal_tool"),
                "--worker-tool",
                str(build / "scfg2_fixed_scaffold_viterbi_tool"),
                "--candidate-jobs", "2",
                "--linearfold-bin", str(fake_linearfold),
                "--linearfold-beam-width", "100",
                *arguments,
            ],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        ).stdout
        augmented = parse_prediction_output(augmented_output)
        assert len(augmented["records"]) == len(parallel["records"]) + 1
        added = augmented["records"][len(parallel["records"])]
        assert added["proposal_source"] == "linearfold_mfe"
        assert added["scaffold_G"] == "((................))"
        assert augmented["summary"]["linearfold_scaffold_added"] == "1"

    failed_parallel = subprocess.run(
        [
            sys.executable,
            str(root / "test/tools/run_scfg2_parallel_hotspot_viterbi.py"),
            "--proposal-tool",
            str(build / "scfg2_hotspot_proposal_tool"),
            "--worker-tool", "/bin/false",
            "--candidate-jobs", "2",
            *arguments,
        ],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
    )
    assert failed_parallel.returncode != 0
    assert "SUMMARY\tok=0" in failed_parallel.stdout

    dataset = {"PKB1": {}, "PKB2": {}}
    assert select_case_ids(dataset, None, None, True) == ["PKB1", "PKB2"]
    with tempfile.TemporaryDirectory() as tmpdir:
        ids_path = Path(tmpdir) / "ids.txt"
        ids_path.write_text("# frozen order\nPKB2\nPKB1\n", encoding="utf-8")
        assert select_case_ids(
            dataset, None, ids_path, False) == ["PKB2", "PKB1"]

        result_path = Path(tmpdir) / "result.tsv"
        provenance = {
            "dataset_sha256": "dataset-hash",
            "tool_sha256": "tool-hash",
        }
        write_rows(result_path, [{
            **provenance,
            "pkb_short": "PKB1",
            "status": "ok",
            "all_f1": 1.0,
        }])
        resumed = load_resumable_rows(
            result_path, provenance, ["PKB1", "PKB2"])
        assert list(resumed) == ["PKB1"]
        try:
            load_resumable_rows(
                result_path,
                {**provenance, "tool_sha256": "different"},
                ["PKB1"],
            )
        except ValueError as error:
            assert "provenance mismatch" in str(error)
        else:
            raise AssertionError("mismatched provenance was resumed")

        command_result = run_captured_command(
            [sys.executable, "-c", "print('native output')"],
            2.0,
            Path(tmpdir) / "stdout.txt",
            Path(tmpdir) / "stderr.txt",
            Path(tmpdir) / "rss.txt",
        )
        assert command_result.exit_status == 0
        assert not command_result.timed_out
        assert command_result.stdout == "native output\n"
        assert command_result.peak_rss_kb is not None

        timeout_result = run_captured_command(
            [sys.executable, "-c", "import time; time.sleep(2)"],
            0.05,
            Path(tmpdir) / "timeout.stdout.txt",
            Path(tmpdir) / "timeout.stderr.txt",
            Path(tmpdir) / "timeout.rss.txt",
        )
        assert timeout_result.timed_out
        assert timeout_result.exit_status != 0

        heterogeneous = Path(tmpdir) / "heterogeneous.tsv"
        write_rows(heterogeneous, [
            {"pkb_short": "PKB1", "status": "ok", "all_f1": 1.0},
            {"pkb_short": "PKB2", "status": "timeout", "error": "limit"},
        ])
        status = status_summary([
            {"status": "ok"}, {"status": "timeout"}])
        assert status["completed_cases"] == 1
        assert status["timeout_cases"] == 1
        assert status["failed_cases"] == 1

    invalid = {
        "pkb_short": "ambiguous",
        "pipeline_ready": "true",
        "sequence": "ACYU",
        "length": "4",
        "dotbracket": "....",
    }
    try:
        validate_case(invalid)
    except ValueError as error:
        assert "unsupported residues" in str(error)
    else:
        raise AssertionError("ambiguous residue was accepted")

    parsed = parse_prediction_output(SYNTHETIC)
    assert parsed["selected_record"]["index"] == "1"
    row = {
        "pkb_short": "synthetic",
        "length": "8",
        "rna_type": "test",
        "dotbracket": "((....))",
    }
    evaluated = evaluate_case(row, parsed)
    assert evaluated["all_f1"] == 1.0
    assert evaluated["best_completed_f1"] == 1.0
    assert evaluated["selection_f1_regret"] == 0.0
    assert evaluated["exact"] == 1
    assert evaluated["nonempty_scaffold_records"] == 1
    assert evaluated["any_nonempty_scaffold_subset_of_reference"] == 1
    assert evaluated["best_compatible_scaffold_pair_recall"] == 0.5
    assert evaluated["case3_truncated_records"] == ""
    summary = aggregate([evaluated])
    assert summary["macro_all_f1"] == 1.0
    assert summary["micro_all_f1"] == 1.0
    assert summary["exact_rate"] == 1.0
    print("scfg2_pseudobase_viterbi_runner=ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
