#!/usr/bin/env python3
"""Shared pin and recorded build identity for external legacy observations."""

import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(checkout, *args):
    return subprocess.check_output(["git", "-C", str(checkout), *args])


def configuration():
    config = json.loads((HERE / "config.json").read_text())
    if (config["schema_version"] != 1 or
            config["wmb_closure_policy"] != "legacy-compatible-wmb" or
            not config["patches"] or len(set(config["patches"])) != len(config["patches"])):
        raise ValueError("invalid legacy oracle configuration")
    for name in config["patches"]:
        if Path(name).name != name or not (HERE / "patches" / name).is_file():
            raise ValueError(f"invalid oracle patch: {name}")
    return config


def record(checkout, destination):
    config = configuration()
    head = git(checkout, "rev-parse", "HEAD").decode().strip()
    if head != config["upstream_commit"]:
        raise ValueError(f"unexpected legacy oracle commit: {head}")
    sources = git(checkout, "diff", "--name-only", "HEAD").decode().splitlines()
    objects = ["build/CMakeFiles/CParty.dir/src/" + name + ".cc.o" for name in
               ("part_func", "pseudo_loop", "sparse_tree", "SHAPE", "dot_plot", "mea", "centroid")]
    result = dict(config,
                  config_sha256=digest(HERE / "config.json"),
                  patches=[{"name": name, "sha256": digest(HERE / "patches" / name)}
                           for name in config["patches"]],
                  patched_diff_sha256=hashlib.sha256(git(checkout, "diff", "--binary", "HEAD")).hexdigest(),
                  source_sha256={name: digest(checkout / name) for name in sources},
                  artifact_sha256={name: digest(checkout / name) for name in
                                   ["build/CParty", "build/libRNA.a", *objects]})
    destination.write_text(json.dumps(result, indent=2) + "\n")
    return result


def verify_checkout(checkout):
    config = configuration()
    manifest_path = checkout.parent / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for key in ("schema_version", "upstream_url", "upstream_commit", "wmb_closure_policy"):
        if manifest[key] != config[key]:
            raise ValueError(f"legacy oracle manifest disagrees with config: {key}")
    if (manifest["config_sha256"] != digest(HERE / "config.json") or
            manifest["patches"] != [{"name": name, "sha256": digest(HERE / "patches" / name)}
                                    for name in config["patches"]]):
        raise ValueError("legacy oracle patches changed; rebuild the oracle")
    head = git(checkout, "rev-parse", "HEAD").decode().strip()
    if head != config["upstream_commit"]:
        raise ValueError(f"unexpected legacy oracle commit: {head}")
    actual_diff = hashlib.sha256(git(checkout, "diff", "--binary", "HEAD")).hexdigest()
    if actual_diff != manifest["patched_diff_sha256"]:
        raise ValueError("legacy oracle checkout changed after its recorded build")
    for group in ("source_sha256", "artifact_sha256"):
        for name, expected in manifest[group].items():
            if digest(checkout / name) != expected:
                raise ValueError(f"legacy oracle source/build changed: {name}")
    return manifest


if __name__ == "__main__":
    if sys.argv[1:] == ["config"]:
        config = configuration()
        print("\n".join([config["upstream_url"], config["upstream_commit"], *config["patches"]]))
    elif len(sys.argv) == 4 and sys.argv[1] == "record":
        record(Path(sys.argv[2]), Path(sys.argv[3]))
    else:
        raise SystemExit("usage: oracle_config.py config | record CHECKOUT MANIFEST")
