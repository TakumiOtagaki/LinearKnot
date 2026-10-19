#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parents[1]
EXPORT_SCRIPT = SCRIPT_DIR / "export_legacy_replay_artifact.py"
PARSE_STDOUT_SCRIPT = SCRIPT_DIR / "parse_cparty_stdout.sh"


W_RE = re.compile(
    r"^\[MFE_TRACE_W\] "
    r"j=(?P<j>-?\d+) "
    r"m1=(?P<m1>-?\d+) "
    r"m2=(?P<m2>-?\d+) "
    r"m3=(?P<m3>-?\d+) "
    r"best_k_m2=(?P<best_k_m2>-?\d+) "
    r"best_k_m3=(?P<best_k_m3>-?\d+) "
    r"chosen=(?P<chosen>-?\d+) "
    r"Wj=(?P<Wj>-?\d+)$"
)

WMB_RE = re.compile(
    r"^\[legacy-trace-WMB\] "
    r"i=(?P<i>-?\d+) "
    r"j=(?P<j>-?\d+) "
    r"m2=(?P<m2>-?\d+) "
    r"mWMBP=(?P<mWMBP>-?\d+) "
    r"chosen=(?P<chosen>-?\d+) "
    r"best_l=(?P<best_l>-?\d+) "
    r"best_Bp_lj=(?P<best_Bp_lj>-?\d+) "
    r"best_be=(?P<best_be>-?\d+) "
    r"best_wmbp=(?P<best_wmbp>-?\d+) "
    r"best_wi=(?P<best_wi>-?\d+)$"
)

WMBP_RE = re.compile(
    r"^\[MFE_TRACE_WMBP\] "
    r"i=(?P<i>-?\d+) "
    r"j=(?P<j>-?\d+) "
    r"m1=(?P<m1>-?\d+) "
    r"m2=(?P<m2>-?\d+) "
    r"m3=(?P<m3>-?\d+) "
    r"m4=(?P<m4>-?\d+) "
    r"total=(?P<total>-?\d+) "
    r"chosen=(?P<chosen>-?\d+)$"
)

FWD_WI_RE = re.compile(
    r"^\[LEGACY_FWD_WI\] "
    r"i=(?P<i>-?\d+) "
    r"j=(?P<j>-?\d+) "
    r"m1=(?P<m1>-?\d+) "
    r"m2=(?P<m2>-?\d+) "
    r"m3=(?P<m3>-?\d+) "
    r"m4=(?P<m4>-?\d+) "
    r"m5=(?P<m5>-?\d+) "
    r"best=(?P<best>-?\d+)$"
)

FWD_WIP_RE = re.compile(
    r"^\[LEGACY_FWD_WIP\] "
    r"i=(?P<i>-?\d+) "
    r"j=(?P<j>-?\d+) "
    r"m1=(?P<m1>-?\d+) "
    r"m2=(?P<m2>-?\d+) "
    r"m3=(?P<m3>-?\d+) "
    r"m4=(?P<m4>-?\d+) "
    r"m5=(?P<m5>-?\d+) "
    r"m6=(?P<m6>-?\d+) "
    r"m7=(?P<m7>-?\d+) "
    r"best=(?P<best>-?\d+)$"
)

BT_FREE_RE = re.compile(
    r"^\[LEGACY_BT_FREE\] "
    r"j=(?P<j>-?\d+) "
    r"best_row=(?P<best_row>-?\d+) "
    r"best_i=(?P<best_i>-?\d+) "
    r"best_score=(?P<best_score>-?\d+)$"
)

BT_BE_RE = re.compile(
    r"^\[LEGACY_BT_BE\] "
    r"i=(?P<i>-?\d+) "
    r"j=(?P<j>-?\d+) "
    r"best_row=(?P<best_row>-?\d+) "
    r"best_score=(?P<best_score>-?\d+)$"
)

BT_BE_DETAIL_RE = re.compile(
    r"^\[LEGACY_BT_BE_DETAIL\] "
    r"i=(?P<detail_i>-?\d+) "
    r"j=(?P<detail_j>-?\d+) "
    r"ip=(?P<detail_ip>-?\d+) "
    r"jp=(?P<detail_jp>-?\d+) "
    r"best_row=(?P<detail_best_row>-?\d+) "
    r"best_l=(?P<detail_best_l>-?\d+) "
    r"best_score=(?P<detail_best_score>-?\d+)$"
)


@dataclass(frozen=True)
class Probe:
    kind: str
    i: int
    j: int

    def key(self) -> str:
        if self.kind in {"W", "BT_FREE"}:
            return f"{self.kind}:{self.j}"
        return f"{self.kind}:{self.i},{self.j}"


PROBE_KIND_TO_RE = {
    "W": W_RE,
    "WMB": WMB_RE,
    "WMBP": WMBP_RE,
    "FWD_WI": FWD_WI_RE,
    "FWD_WIP": FWD_WIP_RE,
    "BT_FREE": BT_FREE_RE,
    "BT_BE": BT_BE_RE,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare legacy replay on the same public restricted structure with "
            "and without CPARTY_LEGACY_REPLAY_PROVENANCE_G, and summarize the "
            "first traced forward/backtrack differences."
        )
    )
    parser.add_argument("legacy_bin")
    parser.add_argument("param_file")
    parser.add_argument("seq")
    parser.add_argument("input_g")
    parser.add_argument(
        "--restricted",
        help=(
            "Public restricted structure to replay. Defaults to the exported "
            "printed_round_projection from export_legacy_replay_artifact.py."
        ),
    )
    parser.add_argument(
        "--probe",
        action="append",
        default=[],
        help=(
            "Explicit trace probe. Forms: W:70, WMB:3,70, WMBP:3,70, "
            "FWD_WI:19,24, FWD_WIP:19,24, BT_FREE:70, BT_BE:18,28."
        ),
    )
    parser.add_argument(
        "--no-auto",
        action="store_true",
        help="Disable automatic W/BT_FREE root probes and WMB child follow-up.",
    )
    return parser.parse_args()


def parse_cparty_stdout(stdout_text: str) -> dict[str, str]:
    run = subprocess.run(
        [str(PARSE_STDOUT_SCRIPT)],
        input=stdout_text,
        text=True,
        capture_output=True,
        check=True,
    )
    fields = run.stdout.rstrip("\n").split("\t")
    if len(fields) != 6:
        raise RuntimeError(f"unexpected parse_cparty_stdout output: {run.stdout!r}")
    return {
        "seq": fields[0],
        "restricted": fields[1],
        "mfe_structure": fields[2],
        "mfe_energy": fields[3],
        "pf_structure": fields[4],
        "pf_energy": fields[5],
    }


def export_legacy_replay_artifact(legacy_bin: str, param_file: str, seq: str, input_g: str) -> dict[str, str]:
    run = subprocess.run(
        [str(EXPORT_SCRIPT), legacy_bin, param_file, seq, input_g, "--format", "kv"],
        text=True,
        capture_output=True,
        check=True,
    )
    fields: dict[str, str] = {}
    for line in run.stdout.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        fields[key] = value
    return fields


def parse_probe(text: str) -> Probe:
    if ":" not in text:
        raise ValueError(f"invalid probe '{text}': expected KIND:... form")
    kind, payload = text.split(":", 1)
    kind = kind.strip().upper()
    if kind not in PROBE_KIND_TO_RE and kind != "BT_FREE":
        raise ValueError(f"unsupported probe kind '{kind}'")
    if kind in {"W", "BT_FREE"}:
        return Probe(kind, -1, int(payload))
    if "," not in payload:
        raise ValueError(f"invalid probe '{text}': expected i,j payload")
    i_text, j_text = payload.split(",", 1)
    return Probe(kind, int(i_text), int(j_text))


def parse_pair_map(structure: str) -> dict[int, int]:
    stack: list[int] = []
    pair_map: dict[int, int] = {}
    for idx, ch in enumerate(structure, start=1):
        if ch == "(":
            stack.append(idx)
        elif ch == ")":
            if not stack:
                raise ValueError(f"unbalanced structure at position {idx}: {structure}")
            left = stack.pop()
            pair_map[left] = idx
            pair_map[idx] = left
        elif ch != ".":
            raise ValueError(f"unsupported structure glyph '{ch}' in restricted structure: {structure}")
    if stack:
        raise ValueError(f"unbalanced structure, unmatched '(' at positions {stack}: {structure}")
    return pair_map


def build_trace_env(probe: Probe) -> dict[str, str]:
    if probe.kind == "W":
        return {
            "CPARTY_MFE_TRACE_W": "1",
            "CPARTY_MFE_TRACE_W_J": str(probe.j),
        }
    if probe.kind == "WMB":
        return {"CPARTY_LEGACY_TRACE_WMB": f"{probe.i},{probe.j}"}
    if probe.kind == "WMBP":
        return {"CPARTY_MFE_TRACE_WMBP": f"{probe.i},{probe.j}"}
    if probe.kind == "FWD_WI":
        return {"CPARTY_LEGACY_FWD_WI": f"{probe.i},{probe.j}"}
    if probe.kind == "FWD_WIP":
        return {"CPARTY_LEGACY_FWD_WIP": f"{probe.i},{probe.j}"}
    if probe.kind == "BT_FREE":
        return {"CPARTY_LEGACY_BT_FREE_J": str(probe.j)}
    if probe.kind == "BT_BE":
        return {"CPARTY_LEGACY_BT_BE": f"{probe.i},{probe.j}"}
    raise ValueError(f"unsupported probe kind '{probe.kind}'")


def run_legacy(
    legacy_bin: str,
    param_file: str,
    seq: str,
    restricted: str,
    provenance_g: str | None,
    extra_env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    if provenance_g is not None:
        env["CPARTY_LEGACY_REPLAY_PROVENANCE_G"] = provenance_g
    if extra_env:
        env.update(extra_env)
    return subprocess.run(
        [legacy_bin, "-d2", "-r", restricted, "-P", param_file, seq],
        text=True,
        capture_output=True,
        env=env,
        check=True,
    )


def parse_probe_fields(probe: Probe, stderr_text: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    regex = PROBE_KIND_TO_RE.get(probe.kind)
    if regex is not None:
        for line in stderr_text.splitlines():
            match = regex.match(line.strip())
            if match:
                fields.update(match.groupdict())
    if probe.kind == "BT_BE":
        for line in stderr_text.splitlines():
            match = BT_BE_DETAIL_RE.match(line.strip())
            if match:
                fields.update(match.groupdict())
    return fields


def changed_keys(left: dict[str, str], right: dict[str, str]) -> list[str]:
    keys = sorted(set(left) | set(right))
    return [key for key in keys if left.get(key) != right.get(key)]


def format_fields(fields: dict[str, str]) -> str:
    if not fields:
        return "<none>"
    return ", ".join(f"{key}={fields[key]}" for key in sorted(fields))


def auto_followup_probes(
    root_probe: Probe,
    off_fields: dict[str, str],
    on_fields: dict[str, str],
    restricted_pair_map: dict[int, int],
) -> list[Probe]:
    probes: list[Probe] = []
    if root_probe.kind == "W":
        for key in ("best_k_m3", "best_k_m2"):
            off_value = int(off_fields.get(key, "-1"))
            on_value = int(on_fields.get(key, "-1"))
            if off_value > 0:
                probes.append(Probe("WMB", off_value, root_probe.j))
            if on_value > 0:
                probes.append(Probe("WMB", on_value, root_probe.j))
    elif root_probe.kind == "WMB":
        for fields in (off_fields, on_fields):
            best_l = int(fields.get("best_l", "-1"))
            best_bp_lj = int(fields.get("best_Bp_lj", "-1"))
            outer_left = restricted_pair_map.get(root_probe.j, -1)
            inner_left = restricted_pair_map.get(best_bp_lj, -1)
            if best_l > 0:
                probes.append(Probe("WMBP", root_probe.i, best_l))
            if outer_left > 0 and inner_left > 0:
                probes.append(Probe("BT_BE", outer_left, inner_left))
    else:
        return []
    deduped: list[Probe] = []
    seen: set[str] = set()
    for probe in probes:
        probe_key = probe.key()
        if probe_key in seen:
            continue
        seen.add(probe_key)
        deduped.append(probe)
    return deduped


def main() -> int:
    args = parse_args()

    artifact = export_legacy_replay_artifact(args.legacy_bin, args.param_file, args.seq, args.input_g)
    restricted = args.restricted or artifact["printed_round_projection"]
    restricted_pair_map = parse_pair_map(restricted)

    print(f"seq={args.seq}")
    print(f"input_g={args.input_g}")
    print(f"restricted={restricted}")
    print(f"legacy_reported_mfe_r={artifact['legacy_reported_mfe_r']}")
    print(f"legacy_reported_mfe_energy={artifact['legacy_reported_mfe_energy']}")
    print(f"printed_round_projection={artifact['printed_round_projection']}")
    print(f"printed_square_projection={artifact['printed_square_projection']}")
    print(f"emitted_lane_g={artifact['emitted_lane_g']}")
    print(f"emitted_lane_gprime={artifact['emitted_lane_gprime']}")

    off_run = run_legacy(args.legacy_bin, args.param_file, args.seq, restricted, None)
    on_run = run_legacy(args.legacy_bin, args.param_file, args.seq, restricted, args.input_g)
    off_stdout = parse_cparty_stdout(off_run.stdout)
    on_stdout = parse_cparty_stdout(on_run.stdout)
    print(f"replay_off_mfe_energy={off_stdout['mfe_energy']}")
    print(f"replay_on_mfe_energy={on_stdout['mfe_energy']}")
    print(f"replay_energy_delta={float(on_stdout['mfe_energy']) - float(off_stdout['mfe_energy']):.2f}")

    probes: list[Probe] = []
    if not args.no_auto:
        probes.extend([Probe("W", -1, len(args.seq)), Probe("BT_FREE", -1, len(args.seq))])
    probes.extend(parse_probe(text) for text in args.probe)

    ordered: list[Probe] = []
    seen: set[str] = set()
    while probes:
        probe = probes.pop(0)
        probe_key = probe.key()
        if probe_key in seen:
            continue
        seen.add(probe_key)
        ordered.append(probe)

        off_probe_run = run_legacy(
            args.legacy_bin, args.param_file, args.seq, restricted, None, build_trace_env(probe)
        )
        on_probe_run = run_legacy(
            args.legacy_bin, args.param_file, args.seq, restricted, args.input_g, build_trace_env(probe)
        )
        off_fields = parse_probe_fields(probe, off_probe_run.stderr)
        on_fields = parse_probe_fields(probe, on_probe_run.stderr)
        diff_keys = changed_keys(off_fields, on_fields)

        print(f"probe={probe.key()}")
        print(f"probe_off={format_fields(off_fields)}")
        print(f"probe_on={format_fields(on_fields)}")
        print(f"probe_changed_keys={','.join(diff_keys) if diff_keys else '<none>'}")

        if not args.no_auto:
            probes.extend(auto_followup_probes(probe, off_fields, on_fields, restricted_pair_map))

    print(f"probe_count={len(ordered)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
