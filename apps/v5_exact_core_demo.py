import streamlit as st
import pandas as pd
from pathlib import Path

st.set_page_config(page_title="Kira v5 Exact-Core Demo", layout="wide")
st.title("Kira v5 Exact-Core Demo")
st.markdown("**Narrow, reproducible demo of the v5 selectivity exact core**")

demo_dir = Path("data/processed")
summary_path = demo_dir / "selectivity_v5_exact_core_summary.json"
csv_path = demo_dir / "selectivity_v5_exact_core_rows.csv"

if summary_path.exists():
    summary = pd.read_json(summary_path, typ="series")
    st.subheader("Summary Statistics")
    st.json(summary.to_dict())
else:
    st.warning("Summary JSON not found. Run the exact_core pipeline first.")

if csv_path.exists():
    df = pd.read_csv(csv_path)
    st.subheader("Trainable Exact-Core Rows")
    st.dataframe(df.head(20), use_container_width=True)

    st.subheader("Label Status")
    st.bar_chart(df["label_status"].value_counts())

    trainable = df[df["is_trainable_exact_core"] == 1]
    st.subheader("Trainable Rows by Pair")
    pair_counts = trainable.groupby("pair_id").size()
    st.bar_chart(pair_counts)

    st.subheader("Positive vs Negative by Pair")
    pos_neg = trainable.groupby("pair_id")["is_selective"].value_counts().unstack(fill_value=0)
    st.dataframe(pos_neg)

    st.caption(f"Total rows: {len(df)} | Trainable: {len(trainable)}")
else:
    st.warning("CSV not found. Run the exact_core pipeline first.")
