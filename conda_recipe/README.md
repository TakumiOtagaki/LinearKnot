# LinearKnot Conda package

The package is named `linearknot`; it installs `linearknot` with seven
subcommands: `fold`, `access`, `profile`, `sample`, `subopt`, `mea`, and
`centroid`, plus the C++ libraries and public headers.

## Build locally

Use a clean checkout/export so local build directories and research data are
not copied by `source: path: ..`. A C++17 compiler is required. The version
comes from the repository's `VERSION` file, shared with CMake.

```bash
conda install -c conda-forge conda-build
conda build ./conda_recipe --no-anaconda-upload
conda install --use-local linearknot
linearknot fold --sequence GGGGGAACCAACCC --scaffold '((.....)).....'
```

The recipe builds with `BUILD_TESTING=OFF` and `-j8`. Package tests exercise
all seven installed subcommands and both beam/exact folding modes, using the
bundled parameter tables without a source-tree-relative parameter path.
The bundled ViennaRNA C sources use C99 explicitly. The recipe disables
LTO so installed static libraries contain ordinary object code that does
not require the package builder's compiler-specific LTO plugin. Native
builds keep their existing LTO default (`LINEARKNOT_ENABLE_LTO=ON`).

## Release workflow

`.github/workflows/CondaUpload.yml` builds and tests Linux/macOS x86_64 and
arm64 packages, and saves them as workflow artifacts. It runs on a published
release or manually with `workflow_dispatch`. Release tags must match
`VERSION` (`v0.1.1` or `0.1.1`, for example).

For Anaconda upload, configure this repository's `ANACONDA_USER` Actions
variable and `ANACONDA_API_TOKEN` secret for that account. When the account
variable is absent, packages are built and saved without publishing. Manual
workflow runs never upload. The upstream CParty account is not hardcoded,
and existing packages are not force-overwritten.

To upload an artifact manually after downloading it:

```bash
anaconda upload --user YOUR_ACCOUNT /path/to/linearknot-VERSION-BUILD.conda
```

Metadata uses the existing repository `LICENSE` (GPL-3.0); the former recipe's
MIT label did not match that file. This cleanup does not change the license.

The workflow follows the official [runner labels](https://docs.github.com/en/actions/reference/runners/github-hosted-runners)
and [setup-miniconda activation instructions](https://github.com/conda-incubator/setup-miniconda).
Version extraction uses conda-build's [load_file_regex](https://docs.conda.io/projects/conda-build/en/latest/resources/define-metadata.html#loading-data-from-other-files).
