# Loblaw Bio — Clinical Trial Analysis

Analysis of immune cell populations for Bob Loblaw's clinical trial, exploring how the drug candidate **miraclib** affects five immune cell types (B cells, CD8 T cells, CD4 T cells, NK cells, and monocytes) across melanoma, carcinoma, and healthy patients.

## Getting Started

From a GitHub Codespace (or any terminal):

```bash
pip install -r requirements.txt
python load_data.py
streamlit run dashboard.py
```

That's it. `load_data.py` creates the SQLite database (`loblaw.db`) from the raw CSV, and the dashboard reads from it. If you're in Codespaces, Streamlit will open on port 8501 and then you should get a prompt to open it in your browser.

## Database Schema

The raw CSV (`data/cell-count.csv`) is a flat file where every row is a sample, with subject-level info (age, sex, treatment, etc.) repeated across that subject's three timepoint rows, and the five cell populations spread across columns. The database normalizes this into three tables:

```text
subjects (1 row per patient)
  ├── subject     TEXT  PK
  ├── project     TEXT
  ├── condition   TEXT       
  ├── age         INTEGER
  ├── sex         TEXT       
  ├── treatment   TEXT       
  └── response    TEXT       

samples (1 row per biological sample — 3 per subject at days 0, 7, 14)
  ├── sample                    TEXT  PK
  ├── subject                   TEXT  FK 
  ├── sample_type               TEXT       
  └── time_from_treatment_start INTEGER

cell_counts (1 row per population per sample — 5 per sample)
  ├── id          INTEGER  PK
  ├── sample      TEXT     FK 
  ├── population  TEXT          
  └── count       INTEGER
```

**Why this design?**

- **No repeated data.** In the CSV, a subject's demographics get copied across all their sample rows. Here, that info lives in `subjects` once. If Bob updates a patient's metadata, it's one row, not three (or more).

- **Cell counts are rows, not columns.** The CSV has `b_cell`, `cd8_t_cell`, etc. as separate columns. The database "unpivots" them so each population is its own row in `cell_counts`. This means queries like "give me the percentage of each population" are just a `GROUP BY` — no need to name every column. More importantly, if the lab starts tracking a sixth cell type tomorrow, no schema change is needed. It's just more rows.

- **Scales naturally.** With hundreds of projects and thousands of samples, this structure stays clean. Subject-level queries (demographics, treatment arms) hit the small `subjects` table. Sample-level queries filter `samples` without touching cell data at all. The heavy `cell_counts` table is indexed on `sample` and `population`, so joins and filters stay fast. You could add new tables (e.g., `projects` with site/PI info, or `assays` for different measurement types) without touching existing ones.

- **Flexible for new analytics.** Want to add gene expression data? New assay results? Longitudinal metadata? Each becomes its own table joining back to `subjects` or `samples`. The normalized structure means new analysis types don't require reshaping what's already there.

## Code Structure

```text
├── data/
│   └── cell-count.csv        # Raw trial data
├── load_data.py               # Creates loblaw.db 
├── analysis.py                # All the computation 
├── dashboard.py               # Streamlit UI
└── requirements.txt
```

The split is intentional:

- **`load_data.py`** is standalone and uses only the Python standard library (`csv`, `sqlite3`). No pandas, no dependencies. It does one thing: build the database. You can re-run it anytime to reset to a clean state.

- **`analysis.py`** holds three functions — `get_data_overview()`, `get_statistical_analysis()`, and `get_subset_analysis()` — each corresponding to a part of the analysis. They all take a database path, run SQL, and return DataFrames (and a matplotlib figure for the stats part). This makes them easy to test, reuse, or call from a notebook if you'd rather work that way.

- **`dashboard.py`** is purely presentation. It imports the analysis functions, caches their results with `@st.cache_data`, and lays everything out in three tabs. No SQL or computation happens here — if you wanted to swap Streamlit for something else, the analysis code wouldn't change at all.
