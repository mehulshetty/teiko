#!/usr/bin/env python3
"""Analysis functions for clinical trial data (Parts 2-4)."""

import sqlite3
import os

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(SCRIPT_DIR, "loblaw.db")

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def get_data_overview(db_path=DB_PATH):
    """Relative frequency of each cell type in each sample.

    Returns a DataFrame with columns:
        sample, total_count, population, count, percentage
    """
    conn = sqlite3.connect(db_path)
    query = """
        SELECT
            sample,
            SUM(count) OVER (PARTITION BY sample) AS total_count,
            population,
            count,
            ROUND(count * 100.0 / SUM(count) OVER (PARTITION BY sample), 2) AS percentage
        FROM cell_counts
        ORDER BY sample, population
    """
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df


def get_statistical_analysis(db_path=DB_PATH):
    """Compare responders vs non-responders for melanoma+miraclib+PBMC.

    Returns:
        df: DataFrame with columns (sample, subject, response, time_from_treatment_start,
            population, count, percentage)
        fig: matplotlib Figure with boxplots
        stats_df: DataFrame with Mann-Whitney U test results per population
    """
    conn = sqlite3.connect(db_path)
    query = """
        SELECT
            cc.sample,
            sub.subject,
            sub.response,
            sa.time_from_treatment_start,
            cc.population,
            cc.count,
            ROUND(cc.count * 100.0 / SUM(cc.count) OVER (PARTITION BY cc.sample), 2) AS percentage
        FROM cell_counts cc
        JOIN samples sa ON cc.sample = sa.sample
        JOIN subjects sub ON sa.subject = sub.subject
        WHERE sub.condition = 'melanoma'
          AND sub.treatment = 'miraclib'
          AND sa.sample_type = 'PBMC'
        ORDER BY cc.sample, cc.population
    """
    df = pd.read_sql_query(query, conn)
    conn.close()

    # Build boxplots and run statistical tests
    fig, axes = plt.subplots(1, 5, figsize=(22, 5))
    fig.suptitle(
        "Cell Population Frequencies: Responders vs Non-Responders\n"
        "(Melanoma, Miraclib, PBMC)",
        fontsize=14,
        y=1.02,
    )

    stats_rows = []
    for i, pop in enumerate(POPULATIONS):
        pop_data = df[df["population"] == pop]
        responders = pop_data[pop_data["response"] == "yes"]["percentage"]
        non_responders = pop_data[pop_data["response"] == "no"]["percentage"]

        stat, pval = mannwhitneyu(responders, non_responders, alternative="two-sided")
        significant = pval < 0.05

        stats_rows.append(
            {
                "population": pop,
                "u_statistic": stat,
                "p_value": round(pval, 6),
                "significant": significant,
            }
        )

        sns.boxplot(
            data=pop_data,
            x="response",
            y="percentage",
            hue="response",
            order=["yes", "no"],
            hue_order=["yes", "no"],
            ax=axes[i],
            palette={"yes": "#2ecc71", "no": "#e74c3c"},
            legend=False,
        )
        title = f"{pop}\n(p={pval:.4f})"
        if significant:
            title += " *"
        axes[i].set_title(title, fontsize=11)
        axes[i].set_xlabel("Response")
        axes[i].set_ylabel("Relative Frequency (%)")

    fig.tight_layout()

    stats_df = pd.DataFrame(stats_rows)
    return df, fig, stats_df


def get_subset_analysis(db_path=DB_PATH):
    """Baseline melanoma PBMC samples treated with miraclib.

    Returns a dict with three DataFrames:
        project_counts: samples per project
        response_counts: subjects per response status
        sex_counts: subjects per sex
    """
    conn = sqlite3.connect(db_path)

    base_where = """
        WHERE sub.condition = 'melanoma'
          AND sa.sample_type = 'PBMC'
          AND sa.time_from_treatment_start = 0
          AND sub.treatment = 'miraclib'
    """

    # Samples per project
    project_counts = pd.read_sql_query(
        f"""
        SELECT sub.project, COUNT(*) AS sample_count
        FROM samples sa
        JOIN subjects sub ON sa.subject = sub.subject
        {base_where}
        GROUP BY sub.project
        ORDER BY sub.project
        """,
        conn,
    )

    # Responder vs non-responder subjects
    response_counts = pd.read_sql_query(
        f"""
        SELECT sub.response, COUNT(DISTINCT sa.subject) AS subject_count
        FROM samples sa
        JOIN subjects sub ON sa.subject = sub.subject
        {base_where}
        GROUP BY sub.response
        ORDER BY sub.response
        """,
        conn,
    )

    # Male vs female subjects
    sex_counts = pd.read_sql_query(
        f"""
        SELECT sub.sex, COUNT(DISTINCT sa.subject) AS subject_count
        FROM samples sa
        JOIN subjects sub ON sa.subject = sub.subject
        {base_where}
        GROUP BY sub.sex
        ORDER BY sub.sex
        """,
        conn,
    )

    conn.close()

    return {
        "project_counts": project_counts,
        "response_counts": response_counts,
        "sex_counts": sex_counts,
    }
