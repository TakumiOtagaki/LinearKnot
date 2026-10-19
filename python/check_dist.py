"""Check distribution metadata, source scope, and installed command payloads."""

import argparse
from email.parser import BytesParser
from pathlib import Path, PurePosixPath
import tarfile
import zipfile

from smoke_test import COMMANDS


def metadata(data):
    result = BytesParser().parsebytes(data)
    if not result["Name"] or not result["Version"]:
        raise RuntimeError("missing distribution name/version")
    if result["Name"] != "linearknot":
        raise RuntimeError("distribution name must be linearknot")
    return result


def check_source(path, expected_version):
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        roots = {PurePosixPath(m.name).parts[0] for m in members}
        if len(roots) != 1:
            raise RuntimeError("source archive must have one root")
        root = roots.pop()
        names = set()
        for member in members:
            relative = PurePosixPath(member.name).relative_to(root)
            if ".." in relative.parts or member.issym() or member.islnk():
                raise RuntimeError("unexpected archive path/link: " + member.name)
            if not relative.parts:
                continue
            if relative.parts[0] not in {
                "CMakeLists.txt", "VERSION", "LICENSE", "THIRD_PARTY_NOTICES.md", "pyproject.toml",
                "cmake", "src", "python", "PKG-INFO",
            }:
                raise RuntimeError("unexpected source payload: " + str(relative))
            names.add(str(relative))
        required = {
            "CMakeLists.txt", "VERSION", "LICENSE", "THIRD_PARTY_NOTICES.md", "pyproject.toml",
            "cmake/LinearKnotConfig.cmake", "cmake/BundledViennaRNA.cmake",
            "cmake/vienna_config.h.in",
            "src/scfg2/api/beam.cc", "src/ViennaRNA/params/default.c",
            "python/README.md", "python/smoke_test.py", "python/check_dist.py",
            "src/ViennaRNA/COPYING", "src/ViennaRNA/AUTHORS",
            "PKG-INFO",
        }
        if not required <= names:
            raise RuntimeError("missing source files: " + repr(sorted(required - names)))
        package = metadata(archive.extractfile(root + "/PKG-INFO").read())
        source_version = archive.extractfile(root + "/VERSION").read().decode().strip()
        if package["Version"] != source_version:
            raise RuntimeError("source metadata and VERSION disagree")
        if expected_version and source_version != expected_version:
            raise RuntimeError("source version differs from checkout")
    print("sdist=ok {} ({} files)".format(path.name, len(names)))


def check_wheel(path, expected_version):
    with zipfile.ZipFile(path) as archive:
        names = [info.filename for info in archive.infolist() if not info.is_dir()]
        metadata_files = [name for name in names if name.endswith(".dist-info/METADATA")]
        if len(metadata_files) != 1:
            raise RuntimeError("expected one wheel METADATA")
        package = metadata(archive.read(metadata_files[0]))
        if expected_version and package["Version"] != expected_version:
            raise RuntimeError("wheel version differs from checkout")
        wheel = metadata_files[0].removesuffix("METADATA") + "WHEEL"
        wheel_metadata = BytesParser().parsebytes(archive.read(wheel))
        if wheel_metadata["Root-Is-Purelib"] != "false":
            raise RuntimeError("native commands must have a platform wheel")
        tags = wheel_metadata.get_all("Tag", [])
        if not tags or any(not tag.startswith("py3-none-") or tag.endswith("-any") for tag in tags):
            raise RuntimeError("unexpected native wheel tags: " + repr(tags))
        if any(not tag.split("-", 2)[2].startswith(("manylinux_", "musllinux_", "macosx_")) for tag in tags):
            raise RuntimeError("release wheels require repaired Linux or macOS tags: " + repr(tags))
        scripts = [name for name in names if ".data/scripts/" in name]
        if {PurePosixPath(name).name for name in scripts} != set(COMMANDS):
            raise RuntimeError("wheel must contain exactly the linearknot native command")
        for name in scripts:
            if archive.read(name)[:4] not in {
                b"\x7fELF", b"\xcf\xfa\xed\xfe", b"\xfe\xed\xfa\xcf",
                b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca",
            }:
                raise RuntimeError("command is not a native executable: " + name)
            if not ((archive.getinfo(name).external_attr >> 16) & 0o111):
                raise RuntimeError("command is not executable: " + name)
        license_files = {name.split(".dist-info/licenses/", 1)[1]
                         for name in names if ".dist-info/licenses/" in name}
        required_licenses = {"LICENSE", "THIRD_PARTY_NOTICES.md",
                             "src/ViennaRNA/COPYING", "src/ViennaRNA/AUTHORS"}
        if not required_licenses <= license_files:
            raise RuntimeError("missing distribution notices: " + repr(sorted(required_licenses - license_files)))
        if any(name.startswith(("include/", "lib/", "paper/", "documents/", "src/")) for name in names):
            raise RuntimeError("unexpected development/research files in wheel")
    print("wheel=ok {} ({} commands)".format(path.name, len(scripts)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("distributions", type=Path, nargs="+")
    parser.add_argument("--version-file", type=Path)
    args = parser.parse_args()
    expected_version = args.version_file.read_text().strip() if args.version_file else None
    for path in args.distributions:
        if path.name.endswith(".tar.gz"):
            check_source(path, expected_version)
        elif path.suffix == ".whl":
            check_wheel(path, expected_version)
        else:
            raise RuntimeError("unsupported distribution: " + str(path))


if __name__ == "__main__":
    main()
