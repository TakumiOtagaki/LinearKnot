"""Exercise every installed CLI without relying on the source checkout."""

import argparse
import os
from pathlib import Path
import subprocess
import sysconfig
import tempfile


COMMANDS = ("linearknot",)
SUBCOMMANDS = ("fold", "access", "profile", "sample", "subopt", "mea", "centroid")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin-dir", type=Path)
    args = parser.parse_args()
    bin_dir = (args.bin_dir or Path(sysconfig.get_path("scripts"))).resolve()
    environment = dict(os.environ)
    environment.pop("CPARTY_FIXED_ENERGY_PARAM_FILE", None)
    common = ["--sequence", "GGGGGAACCAACCC", "--scaffold", "((.....))....."]
    probes = [["--help"], ["--version"]]
    probes.extend([[name, "--help"] for name in SUBCOMMANDS])
    probes.extend([
        ["fold", *common],
        ["fold", *common, "--exact"],
        ["profile", *common],
        ["access", *common, "--access-len", "1,3"],
        ["sample", *common, "--count", "2"],
        ["subopt", *common, "--top-k", "2"],
        ["mea", *common],
        ["centroid", *common],
    ])
    executable = bin_dir / "linearknot"
    if not executable.is_file():
        raise RuntimeError("missing installed command: " + str(executable))
    with tempfile.TemporaryDirectory(prefix="linearknot-installed-") as cwd:
        for arguments in probes:
            result = subprocess.run(
                [str(executable), *arguments], cwd=cwd, env=environment,
                capture_output=True, text=True, timeout=60,
            )
            if result.returncode:
                raise RuntimeError(
                    "{} failed ({}):\n{}{}".format(
                        "linearknot " + " ".join(arguments), result.returncode, result.stdout, result.stderr,
                    )
                )
            if not (result.stdout.strip() or result.stderr.strip()):
                raise RuntimeError("linearknot produced no output")
            if arguments[0] == "fold" and "--help" not in arguments:
                if "(([[[..))..]]]" not in result.stdout:
                    raise RuntimeError("hairpin/exterior pseudoknot was not predicted: " + result.stdout)
                if "ensemble free energy:" not in result.stdout:
                    raise RuntimeError("bundled-parameter prediction changed: " + result.stdout)
    print("installed_cli_smoke=ok (1 command, {} subcommands, {} probes)".format(len(SUBCOMMANDS), len(probes)))


if __name__ == "__main__":
    main()
