"""Plot the full Terc profile and a zoom from one saved full-length run."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "linearknot-gallery-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = Path(__file__).resolve().parent
CONTEXTS = ("Stem", "Hairpin", "Bulge", "Internal", "Multiloop", "Exterior", "PseudoLoop")
# Match the colors in examples/plot.py and the other README figures.
COLORS = ("#0f766e", "#f59e0b", "#f97316", "#ef4444", "#a855f7", "#cbd5e1", "#4f46e5")
REGION = (40, 160)


def read_profile(path, name, scaffold):
    lines = [line.split() for line in path.read_text().splitlines() if line.strip()]
    if lines[0] != [name] or len(lines) != 8:
        raise ValueError("expected one profile record matching input.txt")
    rows = {fields[0]: np.array([float(value) for value in fields[1:]]) for fields in lines[1:]}
    if set(rows) != set(CONTEXTS):
        raise ValueError("expected all seven structural contexts")
    profile = np.array([rows[context] for context in CONTEXTS])
    if profile.shape != (7, len(scaffold)) or not np.all(np.isfinite(profile)):
        raise ValueError("profile length or numeric values do not match the input")
    if np.any(profile < 0) or np.any(profile > 1):
        raise ValueError("context probability is outside [0, 1]")
    np.testing.assert_allclose(profile.sum(axis=0), 1, atol=1e-9, rtol=0)
    fixed = np.array([base != "." for base in scaffold])
    np.testing.assert_allclose(profile[0, fixed], 1, atol=1e-9, rtol=0)
    return profile


def draw_profile(ax, profile, start, end):
    positions = np.arange(start, end + 1)
    edges = np.arange(start - 0.5, end + 1)
    bottom = np.zeros(len(positions))
    for row, color in zip(profile, COLORS):
        values = row[start - 1:end]
        ax.stairs(bottom + values, edges, baseline=bottom, fill=True,
                  color=color, linewidth=0, antialiased=False)
        bottom += values
    ax.set(xlim=(start - 0.5, end + 0.5), ylim=(0, 1), yticks=(0, 0.5, 1))
    ax.set_yticklabels(("0", "0.5", "1"))
    ax.axhline(0.5, color="white", alpha=0.45, lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#94a3b8")
        ax.spines[side].set_linewidth(0.7)
    ax.tick_params(axis="both", colors="#475569", length=3)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=EXAMPLE / "output")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "build/examples/terc_profile/figures")
    args = parser.parse_args()
    name, sequence, scaffold = (EXAMPLE / "input.txt").read_text().splitlines()
    if len(sequence) != len(scaffold):
        raise ValueError("sequence/scaffold length mismatch")
    profile_path = args.data_dir / "loop_profile.txt"
    profile = read_profile(profile_path, name, scaffold)
    manifest = json.loads((args.data_dir / "run.json").read_text())
    for path, field in ((EXAMPLE / "input.txt", "input_sha256"), (profile_path, "output_sha256")):
        if hashlib.sha256(path.read_bytes()).hexdigest() != manifest[field]:
            raise ValueError("input/output no longer matches run.json: " + str(path))
    start, end = REGION
    peak_index = start - 1 + int(np.argmax(profile[-1, start - 1:end]))
    peak = float(profile[-1, peak_index])
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "text.color": "#0f172a", "svg.hashsalt": "linearknot-terc-profile"})
    fig = plt.figure(figsize=(11.6, 6.0), facecolor="white")
    grid = fig.add_gridspec(2, 1, height_ratios=(1, 2.05), hspace=0.75,
                           left=0.075, right=0.98, top=0.79, bottom=0.18)
    overview, zoom = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
    fig.text(0.075, 0.945, "Structural context along rat telomerase RNA", fontsize=18, weight="bold")
    fig.text(0.075, 0.895, "388 nt · seven probabilities at each base · conditioned on the fixed scaffold",
             fontsize=11, color="#64748b")
    draw_profile(overview, profile, 1, len(sequence))
    overview.set_xticks((1, 50, 100, 150, 200, 250, 300, 350, len(sequence)))
    overview.set_title("A   Full RNA", loc="left", fontsize=12, weight="bold", pad=11)
    overview.set_title("Box: positions 40–160, enlarged below", loc="right", fontsize=10,
                       color="#64748b", pad=11)
    overview.add_patch(Rectangle((start - 0.5, -0.025), end - start + 1, 1.05,
                                 fill=False, edgecolor="#334155", lw=1.3,
                                 linestyle=(0, (3, 2)), clip_on=False))
    draw_profile(zoom, profile, start, end)
    zoom.set_xticks(np.arange(start, end + 1, 20))
    zoom.set_title("B   Pseudoknot region", loc="left", fontsize=12, weight="bold", pad=14)
    zoom.set_title(f"PseudoLoop peak: {peak:.0%} at base {peak_index + 1}", loc="right",
                   fontsize=11, color=COLORS[-1], pad=14)
    zoom.set_xlabel("Position in the full RNA (nt, 1-based)", labelpad=10)
    fig.text(0.018, 0.47, "Context probability", rotation=90, va="center", fontsize=11)
    handles = [Patch(facecolor=color, label=context) for context, color in zip(CONTEXTS, COLORS)]
    legend = fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.53, 0.035),
                        ncol=7, frameon=False, handlelength=1.1, handletextpad=0.5,
                        columnspacing=1.2, fontsize=10.5)
    legend.get_texts()[-1].set_color(COLORS[-1])
    legend.get_texts()[-1].set_weight("bold")
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / "loop_profile.png", dpi=180, facecolor="white")
    svg = output / "loop_profile.svg"
    fig.savefig(svg, metadata={"Date": None}, facecolor="white")
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
    plt.close(fig)
    print(json.dumps({"profile_shape": list(profile.shape),
                      "normalization_max_error": float(np.max(np.abs(profile.sum(axis=0) - 1))),
                      "zoom_region": [start, end], "pseudoloop_peak_position": peak_index + 1,
                      "pseudoloop_peak_probability": peak, "figure_directory": str(output)}))


if __name__ == "__main__":
    main()
