#!/usr/bin/env python3
"""Use only installed LinearKnot targets and run installed CLIs outside the repo."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def run(args, cwd, env):
    result = subprocess.run(args, cwd=cwd, env=env, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(f"{args} failed:\n{result.stdout}\n{result.stderr}")
    return result.stdout


def main():
    root = Path(__file__).resolve().parents[2]
    build = Path(sys.argv[1] if len(sys.argv) > 1 else
                 os.environ.get("SCFG2_BUILD_DIR", root / "build")).resolve()
    version = (root / "VERSION").read_text().strip()
    env = os.environ.copy()
    env.pop("CPARTY_FIXED_ENERGY_PARAM_FILE", None)
    with tempfile.TemporaryDirectory(prefix="scfg2-install-") as directory:
        work = Path(directory)
        prefix = work / "prefix"
        run(["cmake", "--install", str(build), "--prefix", str(prefix)], work, env)
        for original, installed in [
                ("LICENSE", "LICENSE"),
                ("THIRD_PARTY_NOTICES.md", "THIRD_PARTY_NOTICES.md"),
                ("src/ViennaRNA/COPYING", "licenses/ViennaRNA/COPYING"),
                ("src/ViennaRNA/AUTHORS", "licenses/ViennaRNA/AUTHORS")]:
            notice = prefix / "share/LinearKnot" / installed
            if not notice.is_file() or notice.read_bytes() != (root / original).read_bytes():
                raise RuntimeError("installed license notice is missing or changed: " + installed)
        if not (prefix / "lib/libRNA.a").is_file():
            raise RuntimeError("installed Vienna dependency is missing")
        source = work / "consumer"
        source.mkdir()
        (source / "CMakeLists.txt").write_text('''
cmake_minimum_required(VERSION 3.8)
project(InstalledLinearKnotConsumer LANGUAGES CXX)
find_package(LinearKnot @VERSION@ EXACT CONFIG REQUIRED)
foreach(component RNA scfg2_exact scfg2_beam scfg2_hotspot_proposal scfg2_hotspot_viterbi)
  if(NOT TARGET LinearKnot::${component})
    message(FATAL_ERROR "Missing installed target: ${component}")
  endif()
endforeach()
add_executable(consumer main.cc)
target_link_libraries(consumer PRIVATE LinearKnot::scfg2_beam)
'''.replace("@VERSION@", version))
        (source / "main.cc").write_text('''
#include "scfg2/api/beam.hh"
#include "scfg2/api/subopt.hh"
#include "scfg2/api/sampling.hh"
#include "scfg2/api/profiles.hh"
#include <cmath>
#include <iostream>
#include <string>
int main() {
  auto exact = cparty::scfg2::api::inside("cuucgg", "(....)");
  auto beam = cparty::scfg2::api::beam_outside("CTTCGG", "(....)");
  auto subopt = cparty::scfg2::api::exact_subopt("CUUCGG", "(....)");
  if (!exact.ok || !beam.ok || !subopt.ok || subopt.structures.size() != 1 ||
      std::abs(subopt.structures.front().energy - 3.7) > 1e-10 ||
      std::abs(exact.ensemble_free_energy - 3.7) > 1e-10 ||
      std::abs(beam.ensemble_free_energy - 3.7) > 1e-10) return 1;
  std::size_t sampled = 0, enumerated = 0;
  auto samples = cparty::scfg2::api::beam_sample_each("CUUCGG", "(....)", 100, 7,
      [&](auto &&s) { ++sampled; return false; });
  auto ranks = cparty::scfg2::api::exact_subopt_each("CUUCGG", "(....)",
      [&](auto &&s) { ++enumerated; return false; });
  if (!samples.ok || samples.samples_emitted != 1 || sampled != 1 || !samples.stopped_by_consumer ||
      !ranks.ok || !ranks.reachable || ranks.structures_emitted != 1 || enumerated != 1 ||
      ranks.stop_reason != cparty::scfg2::api::SuboptStopReason::Consumer) return 2;
  namespace api = cparty::scfg2::api;
  const auto profiles = api::beam_structural_profiles("AAAAAA", "......", {1, 3}, {}, true);
  const auto exact_profiles = api::exact_structural_profiles("AAAAAA", "......", {1, 3},
      cparty::scfg2::WmbClosurePolicy::RequireWeakClosure, true);
  if (!profiles.ok || !exact_profiles.ok || profiles.accessibility_by_loop.size() != 2 ||
      profiles.accessibility_by_loop != exact_profiles.accessibility_by_loop ||
      profiles.accessibility_by_loop[1][1][4] != 1.0 ||
      std::string(api::accessibility_loop_context_names[4]) != "Exterior") return 3;
  std::cout << "installed API ok\\n";
}
''')
        consumer_build = work / "consumer-build"
        run(["cmake", "-S", str(source), "-B", str(consumer_build),
             f"-DCMAKE_PREFIX_PATH={prefix}"], work, env)
        run(["cmake", "--build", str(consumer_build), "-j8"], work, env)
        run([str(consumer_build / "consumer")], work, env)

        invalid_parameters = [work / "missing.par", work / "empty.par"]
        invalid_parameters[1].write_text("")
        for index, content in enumerate([
                "x", "\n", "arbitrary parameter text\n",
                "## RNAfold parameter file v2.0\n",
                "## RNAfold parameter file v2.0\n#END\n",
                "## RNAfold parameter file v2.0\n# unknown\n1 2 3\n"]):
            path = work / f"invalid-{index}.par"
            path.write_text(content)
            invalid_parameters.append(path)
        common = ["--sequence", "cuucgg", "--scaffold", "(....)"]
        for tool in ["fold", "profile", "access",
                     "sample", "subopt", "mea", "centroid"]:
            extra = ["--access-len", "1,3", "--by-loop"] if tool == "access" else []
            binary = [str(prefix / "bin/linearknot"), tool]
            run([*binary, *common, *extra], work, env)
            bad = subprocess.run([*binary, "--sequence", "CUUCGG", "--scaffold",
                                  "[....]", *extra], cwd=work, env=env,
                                 text=True, capture_output=True)
            if not bad.returncode or "round brackets" not in bad.stderr:
                raise RuntimeError(f"{tool} accepted square scaffold: {bad}")
            for path in [*invalid_parameters, ""]:
                for mode in [[], ["--exact"]]:
                    bad = subprocess.run([*binary, *common, *extra, *mode],
                                         cwd=work, text=True, capture_output=True,
                                         env=dict(env, CPARTY_FIXED_ENERGY_PARAM_FILE=str(path)))
                    if (not bad.returncode or
                            "CPARTY_FIXED_ENERGY_PARAM_FILE" not in bad.stderr):
                        raise RuntimeError(f"{tool} ignored invalid parameter file {path!r}: {bad}")
            run([*binary, *common, *extra], work, dict(env,
                CPARTY_FIXED_ENERGY_PARAM_FILE=str(root / "params/rna_Turner04.par")))

        binary = [str(prefix / "bin/linearknot"), "fold"]
        direct = run([*binary, *common, "--exact"], work, env)
        record = work / "input.txt"
        record.write_text(">test\ncuucgg\n(....)\n")
        file_input = run([*binary, "-i", str(record), "--exact"], work, env)
        if ("ensemble free energy: 3.7000 kcal/mol" not in direct or
                "ensemble free energy: 3.7000 kcal/mol" not in file_input):
            raise RuntimeError("file and direct sequence spelling disagree")
    print("scfg2_installed_api=ok")


if __name__ == "__main__":
    main()
