import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# BIN1 variant position
VARIANT_POS = 127609280
INTERVAL_START = 127347136

# Style constants
COLOR_REF = "#888888"
COLOR_ALT = "#E24B4A"
COLOR_POS_DIFF = "#FF6B8A"   # pink — risk allele increases signal
COLOR_NEG_DIFF = "#4A90D9"   # blue  — risk allele decreases signal
DARK_BG = "#0f0f0f"
DARK_PAPER = "#1a1a1a"
GRID_COLOR = "#2a2a2a"
FONT_COLOR = "#cccccc"

# Tracks to display: (track index, label, assay)
RNA_TRACKS = [
    (0, "Brain (polyA+, +)", "RNA-seq"),
    (2, "Brain (polyA+, −)", "RNA-seq"),
    (4, "Frontal Cortex GTEx (polyA+)", "RNA-seq"),
]
DNASE_TRACKS = [
    (0, "Brain", "DNase-seq"),
    (1, "Frontal Cortex", "DNase-seq"),
]


def _get_x_and_slices(zoom_bp: int = 25000, bin_size: int = 32):
    """Return genomic x-axis coordinates and array slice for the zoom window."""
    variant_bin = VARIANT_POS - INTERVAL_START
    half = zoom_bp
    start_bin = variant_bin - half
    end_bin = variant_bin + half
    x_full = np.arange(INTERVAL_START + start_bin, INTERVAL_START + end_bin)
    n_bins = (end_bin - start_bin) // bin_size
    trimmed = n_bins * bin_size
    x = x_full[:trimmed].reshape(n_bins, bin_size).mean(axis=1).astype(int)
    return x, start_bin, end_bin, n_bins, bin_size


def _downsample(arr, start_bin, end_bin, n_bins, bin_size):
    """Slice and downsample a 1-D track array."""
    sliced = arr[start_bin:end_bin]
    trimmed = n_bins * bin_size
    return sliced[:trimmed].reshape(n_bins, bin_size).mean(axis=1)


def _apply_dark_theme(fig):
    fig.update_layout(
        paper_bgcolor=DARK_PAPER,
        plot_bgcolor=DARK_BG,
        font=dict(color=FONT_COLOR, family="monospace"),
        legend=dict(bgcolor=DARK_BG, bordercolor=GRID_COLOR, borderwidth=1),
    )
    fig.update_xaxes(gridcolor=GRID_COLOR, zerolinecolor=GRID_COLOR)
    fig.update_yaxes(gridcolor=GRID_COLOR, zerolinecolor=GRID_COLOR)
    return fig


def _add_bin1_annotation(fig, x_pos, label="BIN1 gene"):
    """Arrow annotation pointing at the BIN1 variant, placed in paper y-space."""
    fig.add_annotation(
        x=x_pos,
        y=0.97,
        xref="x",
        yref="paper",
        text=f"▼ {label}",
        showarrow=False,
        font=dict(color="#ffaa00", size=11, family="monospace"),
        bgcolor="#1a1a1a",
        bordercolor="#ffaa00",
        borderwidth=1,
        borderpad=4,
    )


def plot_ref_alt(variant_output, zoom_bp: int = 25000, bin_size: int = 32):
    """
    Overlaid REF (grey dashed) vs ALT (red) tracks for RNA-seq and DNase-seq.
    REF is dashed so it stays visible even when signals are nearly identical.
    Returns a Plotly figure.
    """
    all_tracks = RNA_TRACKS + DNASE_TRACKS
    n_rows = len(all_tracks)
    titles = [f"{label} — {assay}" for _, label, assay in all_tracks]

    fig = make_subplots(
        rows=n_rows,
        cols=1,
        shared_xaxes=True,
        subplot_titles=titles,
        vertical_spacing=0.04,
    )

    x, start_bin, end_bin, n_bins, bin_size_ = _get_x_and_slices(zoom_bp, bin_size)

    def add_track_pair(row, ref_arr, alt_arr, track_idx, show_legend):
        ref_ds = _downsample(ref_arr[:, track_idx], start_bin, end_bin, n_bins, bin_size_)
        alt_ds = _downsample(alt_arr[:, track_idx], start_bin, end_bin, n_bins, bin_size_)
        # REF: dashed grey — visible even when nearly identical to ALT
        fig.add_trace(go.Scatter(
            x=x, y=ref_ds, name="REF allele (C)",
            line=dict(color=COLOR_REF, width=1.5, dash="dash"),
            showlegend=show_legend, legendgroup="ref",
        ), row=row, col=1)
        # ALT: solid red on top
        fig.add_trace(go.Scatter(
            x=x, y=alt_ds, name="ALT allele (T) — Alzheimer's risk",
            line=dict(color=COLOR_ALT, width=1.5),
            showlegend=show_legend, legendgroup="alt",
        ), row=row, col=1)

    ref_rna = variant_output.reference.rna_seq.values
    alt_rna = variant_output.alternate.rna_seq.values
    ref_dnase = variant_output.reference.dnase.values
    alt_dnase = variant_output.alternate.dnase.values

    for row_i, (track_idx, _, assay) in enumerate(all_tracks, start=1):
        show_legend = row_i == 1
        if assay == "RNA-seq":
            add_track_pair(row_i, ref_rna, alt_rna, track_idx, show_legend)
        else:
            add_track_pair(row_i, ref_dnase, alt_dnase, track_idx, show_legend)

    for row_i in range(1, n_rows + 1):
        fig.add_vline(
            x=VARIANT_POS, line_dash="dash", line_color="#ffaa00",
            line_width=1, row=row_i, col=1,
        )

    _add_bin1_annotation(fig, VARIANT_POS + 3000)

    fig.update_layout(
        title=dict(
            text=(
                "BIN1 rs6733839 (chr2:127,609,280 C→T) — Predicted Gene Expression: REF vs ALT"
                "<br><sup>Grey dashed = reference allele (C) &nbsp;|&nbsp; "
                "Red = Alzheimer's risk allele (T) &nbsp;|&nbsp; "
                "Yellow line = variant position</sup>"
            ),
            font=dict(size=15),
        ),
        height=220 * n_rows,
        margin=dict(l=60, r=30, t=100, b=40),
    )
    _apply_dark_theme(fig)
    return fig


def plot_diff(variant_output, zoom_bp: int = 25000, bin_size: int = 32):
    """
    ALT − REF difference tracks for RNA-seq and DNase-seq.
    Pink = risk allele increases signal, blue = decreases signal.
    Returns a Plotly figure.
    """
    all_tracks = RNA_TRACKS + DNASE_TRACKS
    n_rows = len(all_tracks)
    titles = [f"Δ {label} — {assay}" for _, label, assay in all_tracks]

    fig = make_subplots(
        rows=n_rows,
        cols=1,
        shared_xaxes=True,
        subplot_titles=titles,
        vertical_spacing=0.04,
    )

    x, start_bin, end_bin, n_bins, bin_size_ = _get_x_and_slices(zoom_bp, bin_size)

    def add_diff_track(row, ref_arr, alt_arr, track_idx):
        ref_ds = _downsample(ref_arr[:, track_idx], start_bin, end_bin, n_bins, bin_size_)
        alt_ds = _downsample(alt_arr[:, track_idx], start_bin, end_bin, n_bins, bin_size_)
        diff = alt_ds - ref_ds

        pos = np.where(diff >= 0, diff, 0)
        neg = np.where(diff < 0, diff, 0)

        fig.add_trace(go.Bar(
            x=x, y=pos,
            name="Risk allele increases signal",
            marker_color=COLOR_POS_DIFF,
            showlegend=(row == 1), legendgroup="pos",
        ), row=row, col=1)
        fig.add_trace(go.Bar(
            x=x, y=neg,
            name="Risk allele decreases signal",
            marker_color=COLOR_NEG_DIFF,
            showlegend=(row == 1), legendgroup="neg",
        ), row=row, col=1)

    ref_rna = variant_output.reference.rna_seq.values
    alt_rna = variant_output.alternate.rna_seq.values
    ref_dnase = variant_output.reference.dnase.values
    alt_dnase = variant_output.alternate.dnase.values

    for row_i, (track_idx, _, assay) in enumerate(all_tracks, start=1):
        if assay == "RNA-seq":
            add_diff_track(row_i, ref_rna, alt_rna, track_idx)
        else:
            add_diff_track(row_i, ref_dnase, alt_dnase, track_idx)

    for row_i in range(1, n_rows + 1):
        fig.add_hline(y=0, line_color=GRID_COLOR, line_width=1, row=row_i, col=1)
        fig.add_vline(
            x=VARIANT_POS, line_dash="dash", line_color="#ffaa00",
            line_width=1, row=row_i, col=1,
        )

    _add_bin1_annotation(fig, VARIANT_POS + 3000)

    fig.update_layout(
        title=dict(
            text=(
                "BIN1 rs6733839 — How the Alzheimer's Risk Allele Changes Brain Gene Expression (ALT − REF)"
                "<br><sup>Pink = risk allele (T) increases signal &nbsp;|&nbsp; "
                "Blue = risk allele suppresses signal &nbsp;|&nbsp; "
                "Yellow line = variant position (chr2:127,609,280)</sup>"
            ),
            font=dict(size=15),
        ),
        barmode="overlay",
        bargap=0,
        height=220 * n_rows,
        margin=dict(l=60, r=30, t=100, b=40),
    )
    _apply_dark_theme(fig)
    return fig


if __name__ == "__main__":
    from api_call import run_bin1_variant

    print("Fetching variant predictions...")
    output = run_bin1_variant()

    print("Generating REF vs ALT plot...")
    fig_compare = plot_ref_alt(output)
    fig_compare.write_html("bin1_ref_alt.html")
    print("  Saved: bin1_ref_alt.html")

    print("Generating difference plot...")
    fig_diff = plot_diff(output)
    fig_diff.write_html("bin1_diff.html")
    print("  Saved: bin1_diff.html")
