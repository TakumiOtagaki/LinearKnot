# LinearKnot

LinearKnot computes RNA pseudoknot ensembles conditioned on a noncrossing
input scaffold. It provides partition functions, base-pair probabilities,
minimum-energy prediction, structural profiles, accessibility, sampling,
ranked suboptimal structures, and MEA/gamma-centroid decoding.

The Python distribution installs the native `linearknot` command, with
`fold`, `access`, `profile`, `sample`, `subopt`, `mea`, and `centroid` subcommands.
Use **LinearKnot** for the project name; never abbreviate it to **LinKnot**
or **linknot**.

After a version is published to PyPI, install it with one of:

```bash
pip install linearknot
uv add linearknot
uv tool install linearknot
```

Run a prediction directly, or through `uv run` in a uv project:

```bash
linearknot fold --sequence GGGGGAACCAACCC --scaffold '((.....)).....'
uv run linearknot fold --sequence GGGGGAACCAACCC --scaffold '((.....)).....'
linearknot access --help
```

This example predicts `(([[[..))..]]]`: the new pairs `(3,14)`, `(4,13)`,
and `(5,12)` connect the hairpin to the exterior, crossing the scaffold pairs
`(1,9)` and `(2,8)` (positions are one-based).

Scaffolds use balanced round brackets and must be noncrossing. The default
beam mode uses outer ranking; `--exact` selects unpruned computation. Energy
parameters are bundled, so no external ViennaRNA installation or parameter
file is needed.

Source inputs, saved outputs, and reproducible PNG/SVG figures for BPP,
loop profiles, accessibility, and suboptimal structures, along with sample
draws and MEA/centroid predictions, are in the repository's
[example gallery](https://github.com/TakumiOtagaki/LinearKnot/tree/main/examples).

Python 3.9 or newer is required by the installer. The package provides a native
CLI program; it does not provide an importable folding API. To use the tools
from Python, call them with `subprocess.run` and pass arguments as a list.

The package workflow defines separate Linux and macOS wheel builds for
x86_64 and arm64. When a compatible wheel is unavailable, installing the
source distribution requires a C/C++17 compiler; CMake and Ninja are supplied
by the build backend when needed. Windows builds are not currently supported.

To install from a local checkout:

```bash
pip install .
uv add /path/to/LinearKnot
uv tool install /path/to/LinearKnot
```

The distribution name is `linearknot`. See the
[repository](https://github.com/TakumiOtagaki/LinearKnot) for source,
documentation, and release status.

The distribution carries the project's GPLv3 license and the bundled
ViennaRNA's separate COPYING/AUTHORS notices; see
[third-party notices](https://github.com/TakumiOtagaki/LinearKnot/blob/main/THIRD_PARTY_NOTICES.md).
