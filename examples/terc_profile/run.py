"""Recompute the Terc profile with the installed LinearKnot command."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin-dir", type=Path, help="directory containing linearknot; default: PATH")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "build/examples/terc_profile/output")
    args = parser.parse_args()
    executable = args.bin_dir.resolve() / "linearknot" if args.bin_dir else shutil.which("linearknot")
    if not executable or not Path(executable).is_file():
        parser.error("install linearknot or provide --bin-dir build")
    executable = Path(executable).resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    arguments = ["profile", "-i", str(EXAMPLE / "input.txt"), "--beam", "200",
                 "--ranking", "outer", "--outer-pilot-exponent", "2/3",
                 "--max-hairpin", "200", "--max-case3-depth", "200"]
    environment = dict(os.environ)
    environment.pop("CPARTY_FIXED_ENERGY_PARAM_FILE", None)
    result = subprocess.run([str(executable), *arguments], cwd=output, env=environment,
                            capture_output=True, text=True, check=True, timeout=180)
    # This example has one record; omit the CLI's trailing blank separator.
    (output / "loop_profile.txt").write_text(result.stdout.rstrip("\n") + "\n")
    if result.stderr:
        print(result.stderr, end="")
    display = ["linearknot", *arguments]
    display[display.index(str(EXAMPLE / "input.txt"))] = "examples/terc_profile/input.txt"
    manifest = {
        "repository_version": (ROOT / "VERSION").read_text().strip(),
        "input_sha256": hashlib.sha256((EXAMPLE / "input.txt").read_bytes()).hexdigest(),
        "output_sha256": hashlib.sha256((output / "loop_profile.txt").read_bytes()).hexdigest(),
        "binary_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "parameters": "bundled Turner2004; default completed-WMB weak closure",
        "scope": "full 388-nt RNA, beam=200, outer ranking, pilot exponent=2/3, local caps=200",
        "command": display,
        "stdout": "loop_profile.txt",
        "stdout_normalization": "one terminal newline for the single record",
    }
    (output / "run.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("example_outputs=" + str(output))


if __name__ == "__main__":
    main()
