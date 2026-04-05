# AlphaGenome Demo Project

## What this is
A web demo visualizing how the BIN1 rs6733839 Alzheimer's risk variant 
disrupts gene expression in brain tissue, using the AlphaGenome API.

## Stack
- Python + AlphaGenome SDK for API calls
- Plotly for visualizations  
- Streamlit for the web interface
- .env file for API keys (never commit this)

## Key files
- `api_call.py` — AlphaGenome predict_variant logic
- `viz.py` — Plotly visualization functions
- `app.py` — Streamlit web app

## Style
- Dark theme throughout
- REF tracks in grey (#888), ALT tracks in red (#E24B4A)
- Clean minimal UI, no clutter
- All plots use Plotly (not matplotlib)