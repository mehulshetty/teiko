#!/usr/bin/env python3
"""Interactive Streamlit dashboard for clinical trial data analysis."""

import os
import streamlit as st
import pandas as pd

from analysis import get_data_overview, get_statistical_analysis, get_subset_analysis

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(SCRIPT_DIR, "loblaw.db")

st.set_page_config(
    page_title="Clinical Trial Analysis",
    layout="wide",
)

st.title("Loblaw Bio Dashboard")
st.caption("Immune cell population analysis for Loblaw Bio clinical trial")

# Check that the database exists
if not os.path.exists(DB_PATH):
    st.error(
        f"Database not found at `{DB_PATH}`. "
        "Please run `python load_data.py` first to create the database."
    )
    st.stop()

tab1, tab2, tab3 = st.tabs([
    "Data Overview",
    "Statistical Analysis",
    "Subset Analysis",
])

# ---------------------------------------------------------------------------
# Tab 1: Data Overview (Part 2)
# ---------------------------------------------------------------------------
with tab1:
    st.header("Cell Population Relative Frequencies")
    st.markdown(
        "For each sample, the total cell count is the sum across all five populations. "
        "The percentage shows each population's relative frequency within that sample."
    )

    @st.cache_data
    def load_summary():
        return get_data_overview(DB_PATH)

    summary_df = load_summary()

    # Filters
    col_filter1, col_filter2 = st.columns(2)
    with col_filter1:
        pop_filter = st.multiselect(
            "Filter by population",
            options=sorted(summary_df["population"].unique()),
            default=sorted(summary_df["population"].unique()),
        )
    with col_filter2:
        sample_search = st.text_input("Search by sample ID", placeholder="e.g. sample00001")

    filtered = summary_df
    if pop_filter:
        filtered = filtered[filtered["population"].isin(pop_filter)]
    if sample_search:
        filtered = filtered[filtered["sample"].str.contains(sample_search, case=False)]

    st.dataframe(
        filtered,
        width='stretch',
        height=500,
        column_config={
            "sample": "Sample ID",
            "total_count": st.column_config.NumberColumn("Total Count", format="%d"),
            "population": "Population",
            "count": st.column_config.NumberColumn("Count", format="%d"),
            "percentage": st.column_config.NumberColumn("Percentage (%)", format="%.2f"),
        },
    )
    st.caption(f"Showing {len(filtered):,} of {len(summary_df):,} rows")

# ---------------------------------------------------------------------------
# Tab 2: Statistical Analysis (Part 3)
# ---------------------------------------------------------------------------
with tab2:
    st.header("Responders vs Non-Responders")
    st.markdown(
        "Comparison of cell population relative frequencies for **melanoma** patients "
        "receiving **miraclib** (PBMC samples only). "
        "A Mann-Whitney U test (non-parametric) is used to assess significance at α = 0.05."
    )

    @st.cache_data
    def load_stats():
        df, fig, stats_df = get_statistical_analysis(DB_PATH)
        return df, fig, stats_df

    stats_data, stats_fig, stats_results = load_stats()

    st.pyplot(stats_fig)

    st.subheader("Statistical Test Results")
    display_stats = stats_results.copy()
    display_stats.columns = ["Population", "U Statistic", "p-value", "Significant (p<0.05)"]
    st.dataframe(
        display_stats,
        width='stretch',
        hide_index=True,
        column_config={
            "U Statistic": st.column_config.NumberColumn(format="%.1f"),
            "p-value": st.column_config.NumberColumn(format="%.6f"),
        },
    )

    sig_pops = stats_results[stats_results["significant"]]["population"].tolist()
    if sig_pops:
        st.success(
            f"Populations with significant differences (p < 0.05): **{', '.join(sig_pops)}**"
        )
    else:
        st.info("No populations showed a statistically significant difference at α = 0.05.")

# ---------------------------------------------------------------------------
# Tab 3: Subset Analysis (Part 4)
# ---------------------------------------------------------------------------
with tab3:
    st.header("Baseline Subset Analysis")
    st.markdown(
        "Melanoma patients treated with **miraclib**, **PBMC** samples at "
        "**baseline** (time_from_treatment_start = 0)."
    )

    @st.cache_data
    def load_subset():
        return get_subset_analysis(DB_PATH)

    subset = load_subset()

    col1, col2, col3 = st.columns(3)

    with col1:
        st.subheader("Samples per Project")
        proj = subset["project_counts"]
        for _, row in proj.iterrows():
            st.metric(label=row["project"], value=row["sample_count"])
        st.bar_chart(proj.set_index("project"), y="sample_count", color="#4a90d9")

    with col2:
        st.subheader("Response Distribution")
        resp = subset["response_counts"]
        for _, row in resp.iterrows():
            label = "Responders" if row["response"] == "yes" else "Non-Responders"
            st.metric(label=label, value=row["subject_count"])
        chart_resp = resp.copy()
        chart_resp["response"] = chart_resp["response"].map({"yes": "Responders", "no": "Non-Responders"})
        st.bar_chart(chart_resp.set_index("response"), y="subject_count", color="#2ecc71")

    with col3:
        st.subheader("Sex Distribution")
        sex = subset["sex_counts"]
        for _, row in sex.iterrows():
            label = "Male" if row["sex"] == "M" else "Female"
            st.metric(label=label, value=row["subject_count"])
        chart_sex = sex.copy()
        chart_sex["sex"] = chart_sex["sex"].map({"M": "Male", "F": "Female"})
        st.bar_chart(chart_sex.set_index("sex"), y="subject_count", color="#e67e22")
