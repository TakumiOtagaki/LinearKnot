"""Generate the gallery's raw outputs with the installed LinearKnot subcommands."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "examples/hairpin_exterior/input.txt"
DEFAULT_OUTPUT = ROOT / "build/examples/hairpin_exterior/output"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin-dir", type=Path, help="directory containing linearknot; default: PATH")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    common = [
        "-i", str(INPUT), "--beam", "100", "--ranking", "outer",
        "--outer-pilot-exponent", "2/3", "--max-hairpin", "30",
        "--max-case3-depth", "30",
    ]
    jobs = [
        ("fold", ["--bpp", str(output / "bpp.tsv")], "prediction.txt"),
        ("profile", [], "loop_profile.txt"),
        ("access", ["--access-len", "1,3,5", "--probability"], "raccess_probability.txt"),
        ("access", ["--access-len", "1,3,5"], "raccess_energy.txt"),
        ("access", ["--access-len", "1,3,5", "--probability", "--by-loop"], "raccess_byloop.txt"),
        ("sample", ["--count", "20", "--seed", "42"], "samples.txt"),
        ("sample", ["--count", "1000", "--seed", "42", "--unique"], "sample_counts.txt"),
        ("subopt", ["--top-k", "8", "--stats"], "subopt.txt"),
        ("mea", ["--gamma", "1"], "mea.txt"),
        ("centroid", ["--gamma", "1"], "centroid.txt"),
    ]
    environment = dict(os.environ)
    environment.pop("CPARTY_FIXED_ENERGY_PARAM_FILE", None)
    executable = (args.bin_dir.resolve() / "linearknot") if args.bin_dir else shutil.which("linearknot")
    if not executable or not Path(executable).is_file():
        parser.error("cannot find linearknot; install it or provide --bin-dir build")
    executable = Path(executable).resolve()
    version = (ROOT / "VERSION").read_text().strip()
    reported = subprocess.run([str(executable), "--version"], capture_output=True,
                              text=True, check=True).stdout.strip()
    if reported != "LinearKnot " + version:
        parser.error(f"expected LinearKnot {version}, found {reported}; rebuild the executable")
    commands = []
    hashes = {}
    # Run from the output directory to avoid relying on checkout-relative assets.
    for name, extra, filename in jobs:
        result = subprocess.run(
            [str(executable), name, *common, *extra], cwd=output, env=environment,
            text=True, capture_output=True, timeout=120, check=True,
        )
        (output / filename).write_text(result.stdout.rstrip("\n") + "\n")
        if name == "subopt":
            (output / "subopt_stats.txt").write_text(result.stderr)
        elif result.stderr:
            print(result.stderr, end="")
        hashes["linearknot"] = hashlib.sha256(executable.read_bytes()).hexdigest()
        # Keep the checked-in manifest free of machine-specific absolute paths.
        display = ["linearknot", name, *common, *extra]
        display[display.index(str(INPUT))] = "examples/hairpin_exterior/input.txt"
        if "--bpp" in display:
            display[display.index("--bpp") + 1] = "bpp.tsv"
        commands.append({"argv": display, "stdout": filename})
    manifest = {
        "repository_version": version,
        "input_sha256": hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        "parameters": "bundled Turner2004 (external parameter-file override removed)",
        "scope": "beam=100, outer ranking, pilot exponent=2/3, local caps=30",
        "binary_sha256": hashes,
        "commands": commands,
        "stdout_normalization": "one terminal newline for the single record",
    }
    (output / "run.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("example_outputs=" + str(output))


if __name__ == "__main__":
    main()
