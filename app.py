import streamlit as st
from api_call import run_bin1_variant
from viz import plot_ref_alt, plot_diff

st.set_page_config(
    page_title="BIN1 rs6733839 — AlphaGenome Explorer",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Dark theme override ──────────────────────────────────────────────────────
st.markdown("""
<style>
    html, body, [data-testid="stAppViewContainer"], [data-testid="stApp"] {
        background-color: #0f0f0f;
        color: #cccccc;
    }
    [data-testid="stSidebar"] {
        background-color: #1a1a1a;
    }
    h1, h2, h3, h4 { color: #eeeeee; }
    hr { border-color: #2a2a2a; }
    .variant-badge {
        display: inline-block;
        background: #E24B4A22;
        border: 1px solid #E24B4A;
        border-radius: 6px;
        padding: 4px 12px;
        font-family: monospace;
        color: #E24B4A;
        font-size: 1rem;
        margin-bottom: 0.5rem;
    }
    .metric-box {
        background: #1a1a1a;
        border: 1px solid #2a2a2a;
        border-radius: 8px;
        padding: 12px 18px;
        margin: 4px 0;
    }
</style>
""", unsafe_allow_html=True)

# ── Sidebar ──────────────────────────────────────────────────────────────────
with st.sidebar:
    st.title("Controls")
    st.markdown("---")

    zoom_kb = st.slider(
        "Zoom window (kb each side)",
        min_value=5,
        max_value=100,
        value=25,
        step=5,
    )
    zoom_bp = zoom_kb * 1000

    bin_size = st.select_slider(
        "Bin size (bp)",
        options=[16, 32, 64, 128, 256],
        value=32,
    )

    st.markdown("---")
    st.markdown("**Variant**")
    st.markdown('<div class="variant-badge">chr2:127609280 C→T</div>', unsafe_allow_html=True)
    st.markdown("**Gene:** BIN1  \n**rsID:** rs6733839  \n**Disease:** Alzheimer's (GWAS)")
    st.markdown("---")
    st.markdown("**Tissues**")
    st.markdown("- Brain (UBERON:0000955)  \n- Frontal Cortex (UBERON:0001870)")
    st.markdown("**Outputs**")
    st.markdown("- RNA-seq (polyA+, total, GTEx)  \n- DNase-seq (chromatin accessibility)")

# ── Header ───────────────────────────────────────────────────────────────────
st.title("BIN1 rs6733839 — Alzheimer's Risk Variant")
st.markdown(
    "Predicted effect of the **rs6733839** risk SNP on gene expression and chromatin "
    "accessibility in brain tissue, generated with the "
    "[AlphaGenome](https://deepmind.google/technologies/alphagenome/) API."
)

col1, col2, col3 = st.columns(3)
with col1:
    st.markdown('<div class="metric-box">📍 <b>Position</b><br>chr2:127,609,280</div>', unsafe_allow_html=True)
with col2:
    st.markdown('<div class="metric-box">🧬 <b>Alleles</b><br>REF: C &nbsp;|&nbsp; ALT: T</div>', unsafe_allow_html=True)
with col3:
    st.markdown('<div class="metric-box">🧠 <b>Context</b><br>Brain / Frontal Cortex</div>', unsafe_allow_html=True)

st.markdown("---")

# ── Data loading (cached) ────────────────────────────────────────────────────
@st.cache_resource(show_spinner="Querying AlphaGenome API...")
def load_variant_output():
    return run_bin1_variant()

variant_output = load_variant_output()

# ── Tabs ─────────────────────────────────────────────────────────────────────
tab1, tab2 = st.tabs(["REF vs ALT", "Difference (ALT − REF)"])

with tab1:
    st.markdown(
        "**Grey** = reference allele (C) &nbsp;|&nbsp; **Red** = alternate allele (T)  \n"
        "Yellow dashed line marks the variant position."
    )
    fig_compare = plot_ref_alt(variant_output, zoom_bp=zoom_bp, bin_size=bin_size)
    st.plotly_chart(fig_compare, use_container_width=True)

with tab2:
    st.markdown(
        "**Red bars** = ALT increases signal over REF &nbsp;|&nbsp; "
        "**Blue bars** = ALT decreases signal  \n"
        "Yellow dashed line marks the variant position."
    )
    fig_diff = plot_diff(variant_output, zoom_bp=zoom_bp, bin_size=bin_size)
    st.plotly_chart(fig_diff, use_container_width=True)
