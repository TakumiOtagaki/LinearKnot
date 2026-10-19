"""Plot the gallery's saved CLI output; no folding library is imported."""

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
GALLERY = ROOT / "examples/hairpin_exterior"
SCAFFOLD_COLOR = "#0f766e"
ADDED_COLOR = "#4f46e5"
CONTEXTS = ["Stem", "Hairpin", "Bulge", "Internal", "Multiloop", "Exterior", "PseudoLoop"]
CONTEXT_COLORS = ["#0f766e", "#f59e0b", "#f97316", "#ef4444", "#a855f7", "#cbd5e1", "#4f46e5"]


def pair_map(structure):
    stacks = {opening: [] for opening in "([{<"}
    closing = {")": "(", "]": "[", "}": "{", ">": "<"}
    pairs = set()
    for position, symbol in enumerate(structure, 1):
        if symbol in stacks:
            stacks[symbol].append(position)
        elif symbol in closing:
            pairs.add((stacks[closing[symbol]].pop(), position))
        elif symbol != ".":
            raise ValueError("unsupported structure symbol: " + symbol)
    if any(stacks.values()):
        raise ValueError("unbalanced structure")
    return pairs


def read_access(path):
    rows = {}
    metadata = {}
    for line in path.read_text().splitlines():
        if line.startswith("# ") and "=" in line:
            key, value = line[2:].split("=", 1)
            metadata[key] = value
        elif line and line[0].isdigit():
            fields = line.split("\t")
            if len(fields) == 1:
                continue
            start = int(fields[0]) + 1  # CLI starts are zero-based.
            for item in fields[1].split(";"):
                if item:
                    length, value = item.split(",")
                    rows[(start, int(length))] = float(value)
    return rows, metadata


def read_loop_access(path):
    rows = {}
    codes = ("E", "H", "B", "I", "M", "P")
    for line in path.read_text().splitlines():
        if not line[:1].isdigit():
            continue
        fields = line.split("\t")
        if len(fields) == 1:
            continue
        assert len(fields) == 3 and fields[2].startswith("byloop:")
        start = int(fields[0]) + 1
        for item in fields[2][7:].rstrip(";").split(";"):
            length, *values = item.split(",")
            assert tuple(value.split("=")[0] for value in values) == codes
            rows[(start, int(length))] = {
                code: float(value.split("=")[1]) for code, value in zip(codes, values)
            }
    return rows


def arcs(ax, sequence, pairs, scaffold, title, labels=True):
    n = len(sequence)
    ax.axhline(0, color="#cbd5e1", lw=1)
    for i, j in sorted(pairs):
        color = SCAFFOLD_COLOR if (i, j) in scaffold else ADDED_COLOR
        ax.add_patch(Arc(((i + j) / 2, 0), j - i, (j - i) * 0.55,
                         theta1=0, theta2=180, color=color, lw=2))
    for i, base in enumerate(sequence, 1):
        ax.text(i, -0.38, base, ha="center", va="top", family="monospace", fontsize=10)
        if labels:
            ax.text(i, -0.9, str(i), ha="center", va="top", fontsize=8, color="#64748b")
    ax.set(xlim=(0.3, n + 0.7), ylim=(-1.6, n * 0.3 + 0.3))
    ax.set_title(title, loc="left", pad=12)
    ax.axis("off")


def save(fig, directory, name):
    fig.savefig(directory / (name + ".png"), dpi=180, bbox_inches="tight")
    svg = directory / (name + ".svg")
    fig.savefig(svg, bbox_inches="tight", metadata={"Date": None})
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=GALLERY / "output")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "build/examples/hairpin_exterior/figures")
    args = parser.parse_args()
    data = args.data_dir
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    _, sequence, scaffold = GALLERY.joinpath("input.txt").read_text().splitlines()
    n = len(sequence)
    fixed_pairs = pair_map(scaffold)
    prediction = data.joinpath("prediction.txt").read_text().splitlines()[3]
    best, energy = prediction.split(" ", 1)
    assert len(best) == n and fixed_pairs <= pair_map(best)
    plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans", "svg.hashsalt": "linearknot-gallery"})

    probabilities = np.zeros((n, n))
    for line in data.joinpath("bpp.tsv").read_text().splitlines():
        if not line or line.startswith(">"):
            continue
        i, j, raw = line.split()
        i, j, p = int(i), int(j), float(raw)
        assert 1 <= i < j <= n and math.isfinite(p) and 0 <= p <= 1 + 1e-9
        probabilities[i-1, j-1] = probabilities[j-1, i-1] = p
    paired = probabilities.sum(axis=1)
    assert np.all(paired <= 1 + 1e-9)
    for i, j in fixed_pairs:
        assert abs(probabilities[i-1, j-1] - 1) < 1e-9

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1, 1.1]})
    arcs(axes[0], sequence, pair_map(best), fixed_pairs, "MFE extension  " + energy + " kcal/mol")
    axes[0].text(0.5, -0.05, "scaffold G: teal   added pairs: indigo\n" + best,
                 transform=axes[0].transAxes, ha="center", family="monospace", fontsize=10)
    heat = axes[1].imshow(probabilities, origin="upper", vmin=0, vmax=1, cmap="Blues",
                          extent=(0.5, n+0.5, n+0.5, 0.5), interpolation="nearest")
    axes[1].set(title="Union base-pair probabilities", xlabel="Position j (1-based)", ylabel="Position i (1-based)")
    axes[1].set_xticks(range(1, n+1, 2))
    axes[1].set_yticks(range(1, n+1, 2))
    fig.colorbar(heat, ax=axes[1], label="P(pair | scaffold G)", fraction=0.046)
    fig.tight_layout()
    save(fig, output, "bpp")

    profiles = {}
    for line in data.joinpath("loop_profile.txt").read_text().splitlines():
        fields = line.split()
        if fields and fields[0] in CONTEXTS:
            profiles[fields[0]] = np.array([float(value) for value in fields[1:]])
    assert set(profiles) == set(CONTEXTS)
    profile = np.array([profiles[context] for context in CONTEXTS])
    assert profile.shape == (7, n) and np.all(np.isfinite(profile))
    assert np.all(profile >= 0) and np.all(profile <= 1 + 1e-9)
    np.testing.assert_allclose(profile.sum(axis=0), 1, atol=1e-9)
    np.testing.assert_allclose(profiles["Stem"], paired, atol=1e-9)
    x = np.arange(1, n+1)
    fig, ax = plt.subplots(figsize=(10, 3.8))
    bottom = np.zeros(n)
    for context, color in zip(CONTEXTS, CONTEXT_COLORS):
        ax.bar(x, profiles[context], bottom=bottom, width=0.88, color=color, label=context)
        bottom += profiles[context]
    ax.set(xlim=(0.3, n+0.7), ylim=(0, 1), ylabel="Context probability", xlabel="Position (1-based) / base",
           title="linearknot profile: structural contexts, including PseudoLoop")
    ax.set_xticks(x, [str(i) + "\n" + base for i, base in zip(x, sequence)])
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=4, frameon=False)
    fig.tight_layout()
    save(fig, output, "loop_profile")

    access, _ = read_access(data / "raccess_probability.txt")
    opening, metadata = read_access(data / "raccess_energy.txt")
    assert access.keys() == opening.keys()
    rt = float(metadata["folding_rt"])
    for (start, length), p in access.items():
        assert 1 <= start <= n-length+1 and math.isfinite(p) and 0 <= p <= 1
        if length == 1:
            assert abs(p - (1 - paired[start-1])) < 1e-9
        if p == 0:
            assert opening[(start, length)] == math.inf
        else:
            assert abs(opening[(start, length)] + rt * math.log(p)) < 1e-9
    fig, ax = plt.subplots(figsize=(10, 3.8))
    for length, color, marker in zip((1, 3, 5), ("#0f766e", "#4f46e5", "#f59e0b"), ("o", "s", "^")):
        starts = list(range(1, n-length+2))
        assert all((start, length) in access for start in starts)
        ax.plot(starts, [access[(start, length)] for start in starts], color=color, marker=marker,
                ms=5, label="L=" + str(length))
    ax.set(xlim=(0.5, n+0.5), ylim=(-0.025, 1.05), ylabel="P(entire window unpaired | G)",
           xlabel="Window start (1-based) / base", title="linearknot access: joint window accessibility")
    ax.set_xticks(x, [str(i) + "\n" + base for i, base in zip(x, sequence)])
    ax.legend(frameon=False, loc="upper right")
    ax.grid(axis="y", alpha=0.15)
    fig.tight_layout()
    save(fig, output, "raccess")

    loop_access = read_loop_access(data / "raccess_byloop.txt")
    grouped_total, _ = read_access(data / "raccess_byloop.txt")
    assert loop_access.keys() == access.keys() == grouped_total.keys()
    for key, row in loop_access.items():
        assert all(math.isfinite(p) and 0 <= p <= 1 for p in row.values())
        assert abs(sum(row.values()) - access[key]) < 1e-9
        assert abs(grouped_total[key] - access[key]) < 1e-9
        if key[1] == 1:
            for code, context in zip("HBIMEP", CONTEXTS[1:]):
                assert abs(row[code] - profiles[context][key[0]-1]) < 1e-9
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    for ax, length in zip(axes, (3, 5)):
        starts = np.arange(1, n-length+2)
        bottom = np.zeros(len(starts))
        for code, context, color in zip("HBIMEP", CONTEXTS[1:], CONTEXT_COLORS[1:]):
            values = np.array([loop_access[(int(start), length)][code] for start in starts])
            ax.bar(starts, values, bottom=bottom, width=0.88, color=color, label=context)
            bottom += values
        ax.set(xlim=(0.5, n+0.5), ylim=(0, max(float(bottom.max()), 1e-4)*1.25),
               ylabel="Joint probability", title="L=" + str(length))
        ax.grid(axis="y", alpha=0.15)
    axes[0].legend(frameon=False, ncol=3, loc="upper right")
    axes[-1].set_xticks(x, [str(i) + "\n" + base for i, base in zip(x, sequence)])
    axes[-1].set_xlabel("Window start (1-based) / base")
    fig.suptitle("linearknot access: window accessibility by loop context")
    fig.tight_layout()
    save(fig, output, "raccess_byloop")

    lines = data.joinpath("subopt.txt").read_text().splitlines()[3:]
    ranked = []
    for line in lines:
        if line:
            structure, value = line.split(" ", 1)
            e = float(value.strip("()"))
            assert len(structure) == n and fixed_pairs <= pair_map(structure) and math.isfinite(e)
            ranked.append((structure, e))
    assert len(ranked) == 8 and len({structure for structure, _ in ranked}) == 8
    assert all(left[1] <= right[1] + 1e-9 for left, right in zip(ranked, ranked[1:]))
    assert ranked[0][0] == best
    fig, axes = plt.subplots(3, 2, figsize=(10, 8))
    for rank, (ax, (structure, e)) in enumerate(zip(axes.flat, ranked[:6]), 1):
        arcs(ax, sequence, pair_map(structure), fixed_pairs, "#{}   {:.2f} kcal/mol".format(rank, e), labels=False)
        ax.text(0.5, 0.02, structure, transform=ax.transAxes, ha="center", family="monospace")
    fig.suptitle("linearknot subopt: first six of eight ranked structures\nscaffold G: teal   added pairs: indigo", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, output, "subopt")
    print("example_cross_checks=ok; figures=" + str(output))


if __name__ == "__main__":
    main()
