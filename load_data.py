#!/usr/bin/env python3
"""Load cell-count.csv into a normalized SQLite database."""

import csv
import sqlite3
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(SCRIPT_DIR, "loblaw.db")
CSV_PATH = os.path.join(SCRIPT_DIR, "data", "cell-count.csv")

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def create_schema(cursor):
    """Create the normalized database schema."""
    cursor.executescript("""
        CREATE TABLE IF NOT EXISTS subjects (
            subject TEXT PRIMARY KEY,
            project TEXT NOT NULL,
            condition TEXT NOT NULL,
            age INTEGER NOT NULL,
            sex TEXT NOT NULL,
            treatment TEXT NOT NULL,
            response TEXT
        );

        CREATE TABLE IF NOT EXISTS samples (
            sample TEXT PRIMARY KEY,
            subject TEXT NOT NULL,
            sample_type TEXT NOT NULL,
            time_from_treatment_start INTEGER NOT NULL,
            FOREIGN KEY (subject) REFERENCES subjects(subject)
        );

        CREATE TABLE IF NOT EXISTS cell_counts (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sample TEXT NOT NULL,
            population TEXT NOT NULL,
            count INTEGER NOT NULL,
            FOREIGN KEY (sample) REFERENCES samples(sample)
        );

        CREATE INDEX IF NOT EXISTS idx_samples_subject ON samples(subject);
        CREATE INDEX IF NOT EXISTS idx_cell_counts_sample ON cell_counts(sample);
        CREATE INDEX IF NOT EXISTS idx_cell_counts_population ON cell_counts(population);
    """)


def load_csv(cursor, csv_path):
    """Read CSV and insert into normalized tables."""
    seen_subjects = set()

    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            subj = row["subject"]

            # Insert subject only once
            if subj not in seen_subjects:
                cursor.execute(
                    "INSERT INTO subjects (subject, project, condition, age, sex, treatment, response) "
                    "VALUES (?, ?, ?, ?, ?, ?, ?)",
                    (
                        subj,
                        row["project"],
                        row["condition"],
                        int(row["age"]),
                        row["sex"],
                        row["treatment"],
                        row["response"] if row["response"] else None,
                    ),
                )
                seen_subjects.add(subj)

            # Insert sample
            cursor.execute(
                "INSERT INTO samples (sample, subject, sample_type, time_from_treatment_start) "
                "VALUES (?, ?, ?, ?)",
                (
                    row["sample"],
                    subj,
                    row["sample_type"],
                    int(row["time_from_treatment_start"]),
                ),
            )

            # Unpivot cell counts into rows
            cell_count_rows = [
                (row["sample"], pop, int(row[pop])) for pop in POPULATIONS
            ]
            cursor.executemany(
                "INSERT INTO cell_counts (sample, population, count) VALUES (?, ?, ?)",
                cell_count_rows,
            )


def main():
    # Remove existing DB for re-runs
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)

    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON")
    cursor = conn.cursor()

    create_schema(cursor)
    load_csv(cursor, CSV_PATH)

    conn.commit()

    conn.close()
    print(f"\nDatabase created at {DB_PATH}")


if __name__ == "__main__":
    main()
