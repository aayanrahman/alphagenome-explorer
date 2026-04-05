# BIN1 rs6733839 AlphaGenome Explorer

A Streamlit web app that visualizes how the **BIN1 rs6733839** Alzheimer's risk variant disrupts gene expression and chromatin accessibility in brain tissue, using the [AlphaGenome](https://deepmind.google/technologies/alphagenome/) API from Google DeepMind.

---

## What it does

The BIN1 rs6733839 SNP (chr2:127,609,280 C→T) is one of the strongest non-APOE genetic risk factors for late-onset Alzheimer's disease, identified through GWAS. This app queries AlphaGenome's deep learning model to predict how the ALT allele (T) changes:

- **RNA-seq signal** — gene expression in brain and frontal cortex
- **DNase-seq signal** — chromatin accessibility (open chromatin) in the same tissues

Two views are provided:
1. **REF vs ALT** — overlaid tracks showing the reference (C) and alternate (T) allele predictions side-by-side
2. **Difference (ALT − REF)** — explicit difference tracks highlighting where the risk allele increases or decreases signal

---

## Stack

| Layer | Tool |
|---|---|
| Genomic predictions | AlphaGenome SDK (Google DeepMind) |
| Visualization | Plotly |
| Web interface | Streamlit |
| Environment | python-dotenv |

---

## Project structure

```
.
├── api_call.py   # AlphaGenome API call — predict_variant for BIN1 rs6733839
├── viz.py        # Plotly visualization functions (REF vs ALT + diff tracks)
├── app.py        # Streamlit web app
├── .env          # API key (never committed)
└── README.md
```

---

## Setup

### 1. Clone and enter the repo

```bash
git clone <repo-url>
cd alphagenome-explorer
```

### 2. Create the conda environment

```bash
conda create -n alphagenome-env python=3.11
conda activate alphagenome-env
```

### 3. Install dependencies

```bash
pip install alphagenome python-dotenv plotly streamlit
```

### 4. Add your API key

Create a `.env` file in the project root:

```
ALPHAGENOME_API_KEY=your_key_here
```

> Get an API key from the [AlphaGenome API access page](https://deepmind.google/technologies/alphagenome/).

---

## Running

### Web app

```bash
conda activate alphagenome-env
streamlit run app.py
```

Then open `http://localhost:8501` in your browser.

### API call only (CLI)

```bash
python api_call.py
```

Prints the REF and ALT output tensor shapes to confirm the API is working.

### Generate static HTML plots

```bash
python viz.py
```

Writes `bin1_ref_alt.html` and `bin1_diff.html` to the project root.

---

## How it works

### `api_call.py`

- Loads the AlphaGenome API key from `.env`
- Creates a `DnaClient` via `dna_client.create(api_key=...)`
- Defines `run_bin1_variant()` which calls `client.predict_variant()` with:
  - **Interval:** chr2:127,347,136–127,871,424 (500 kb window centered on the variant)
  - **Variant:** chr2:127,609,280 C→T
  - **Outputs:** `RNA_SEQ`, `DNASE`
  - **Ontology terms:** `UBERON:0000955` (brain), `UBERON:0001870` (frontal cortex)
    - Note: `CL:0000129` (microglia) is not supported by the current API — frontal cortex was used as the closest available brain region
- Returns a `VariantOutput` object with `.reference` and `.alternate` fields, each containing `rna_seq` and `dnase` `TrackData` objects

Output shapes: `(524288, 5)` for RNA-seq, `(524288, 2)` for DNase-seq (524,288 bp at 1 bp resolution).

### `viz.py`

- `plot_ref_alt(variant_output, zoom_bp, bin_size)` — overlaid REF (grey `#888`) vs ALT (red `#E24B4A`) tracks across 5 subplots
- `plot_diff(variant_output, zoom_bp, bin_size)` — ALT − REF difference bar tracks; positive diff in red, negative in blue
- Data is zoomed to ±`zoom_bp` around the variant and downsampled to `bin_size` bp resolution for display
- All plots use Plotly with a dark theme; a yellow dashed line marks the variant position on every subplot

**Tracks shown:**
| Row | Track | Assay |
|---|---|---|
| 1 | Brain polyA+ RNA-seq (+) | RNA-seq |
| 2 | Brain polyA+ RNA-seq (−) | RNA-seq |
| 3 | Frontal Cortex GTEx polyA+ | RNA-seq |
| 4 | Brain | DNase-seq |
| 5 | Frontal Cortex | DNase-seq |

### `app.py`

- Streamlit app with a dark CSS override
- Sidebar controls for zoom window (5–100 kb) and bin size (16–256 bp)
- `@st.cache_resource` caches the API call so re-renders don't re-query the model
- Two tabs: **REF vs ALT** and **Difference (ALT − REF)**, each rendering the corresponding Plotly figure

---

## Notes

- The AlphaGenome model takes a 500 kb sequence window — the interval is centered on the variant position
- The `VariantOutput` object holds the full 524,288-position prediction; zooming and downsampling happen only at render time in `viz.py`
- `.env` is gitignored — never commit your API key
