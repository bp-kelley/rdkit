#!/usr/bin/env python
"""Build an xlsx (summary + per-molecule + per-rotor + generation-vs-ranking) from
the `fragConfGenBench bank` stats.tsv.  Usage: bank_stats_xlsx.py <outdir>
Reads <outdir>/stats.tsv, writes <outdir>/bank_stats.xlsx."""
import sys, os
import pandas as pd
import numpy as np

outdir = sys.argv[1] if len(sys.argv) > 1 else "./bank_results"
tsv = os.path.join(outdir, "stats.tsv")
df = pd.read_csv(tsv, sep="\t")

# -1 sentinels -> NaN for rank columns (means "not present in ensemble")
for c in ["best_rank", "rank_lt1", "rank_lt2"]:
    df.loc[df[c] < 0, c] = np.nan

methods = list(df["method"].unique())


def rotor_bin(r, maxbin=14):
    return min(int(r), maxbin)


df["rot_bin"] = df["n_rotors"].apply(rotor_bin)


def overall(g):
    n = len(g)
    hit1 = g["hit_lt1"].sum()
    hit2 = g["hit_lt2"].sum()
    # of the <1A hits, where does the good pose rank by energy?
    hits = g[g["hit_lt1"] == 1]
    return pd.Series({
        "n_mols": n,
        "pct_lt1": 100.0 * hit1 / n if n else np.nan,
        "pct_lt2": 100.0 * hit2 / n if n else np.nan,
        "mean_gen_ms": g["gen_ms"].mean(),
        "median_gen_ms": g["gen_ms"].median(),
        "mean_n_confs": g["n_confs"].mean(),
        "median_n_confs": g["n_confs"].median(),
        "median_best_rms": g["best_rms"].median(),
        # ranking diagnostics (only over mols that DO contain a <1A pose)
        "median_rank_lt1": hits["rank_lt1"].median(),
        "pct_lt1_is_rank0": 100.0 * (hits["rank_lt1"] == 0).sum() / len(hits) if len(hits) else np.nan,
        "pct_lt1_in_top10": 100.0 * (hits["rank_lt1"] < 10).sum() / len(hits) if len(hits) else np.nan,
    })


summary = df.groupby("method").apply(overall).reset_index()

# per-rotor-bin accuracy/timing per method
def per_rotor(metric, aggfn="mean"):
    p = df.pivot_table(index="rot_bin", columns="method", values=metric, aggfunc=aggfn)
    return p

rot_pct1 = df.pivot_table(index="rot_bin", columns="method", values="hit_lt1",
                          aggfunc="mean").mul(100.0)
rot_pct2 = df.pivot_table(index="rot_bin", columns="method", values="hit_lt2",
                          aggfunc="mean").mul(100.0)
rot_nconf = per_rotor("n_confs")
rot_ms = per_rotor("gen_ms")
rot_n = df.pivot_table(index="rot_bin", columns="method", values="name",
                       aggfunc="count")

# generation vs ranking decomposition (per rotor bin) --
# split the accuracy gap into "can't generate a <1A pose" vs "generate but bury it".
gr_rows = []
for method in methods:
    g = df[df["method"] == method]
    for rb, gg in g.groupby("rot_bin"):
        n = len(gg)
        gen = gg["hit_lt1"].sum()               # a <1A pose EXISTS in the ensemble
        hits = gg[gg["hit_lt1"] == 1]
        top10 = (hits["rank_lt1"] < 10).sum()   # ...and is energy-ranked in top 10
        gr_rows.append({
            "method": method, "rot_bin": rb, "n_mols": n,
            "pct_generate_lt1": 100.0 * gen / n if n else np.nan,
            "pct_rank_top10_of_gen": 100.0 * top10 / gen if gen else np.nan,
            "median_rank_lt1": hits["rank_lt1"].median(),
        })
genrank = pd.DataFrame(gr_rows).sort_values(["method", "rot_bin"])

# wide per-molecule table
key = "name"
pivot_cols = ["n_confs", "gen_ms", "best_rms", "rank_lt1", "hit_lt1", "hit_lt2"]
wide = df.pivot_table(index=["name", "n_rotors"], columns="method",
                      values=pivot_cols, aggfunc="first")
wide.columns = [f"{m}_{c}" for c, m in wide.columns]
wide = wide.reset_index()
# gap flags: molecules with no sub-1A pose (generation failures worth inspecting)

xlsx = os.path.join(outdir, "bank_stats.xlsx")
with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
    summary.to_excel(xw, sheet_name="summary", index=False)
    rot_pct1.to_excel(xw, sheet_name="pct_lt1_by_rotor")
    rot_pct2.to_excel(xw, sheet_name="pct_lt2_by_rotor")
    rot_nconf.to_excel(xw, sheet_name="nconfs_by_rotor")
    rot_ms.to_excel(xw, sheet_name="gen_ms_by_rotor")
    rot_n.to_excel(xw, sheet_name="nmols_by_rotor")
    genrank.to_excel(xw, sheet_name="generation_vs_ranking", index=False)
    wide.sort_values("n_rotors").to_excel(xw, sheet_name="per_molecule", index=False)
    df.to_excel(xw, sheet_name="raw", index=False)

# also echo the key numbers to stdout for the terminal analysis
pd.set_option("display.width", 200, "display.max_columns", 40)
print("=== SUMMARY (per method) ===")
print(summary.to_string(index=False))
print("\n=== %<1A by rotor bin ===")
print(rot_pct1.round(1).to_string())
print("\n=== mean n_confs by rotor bin ===")
print(rot_nconf.round(1).to_string())
print("\n=== generation vs ranking (pct_generate_lt1 = a <1A pose exists; "
      "pct_rank_top10_of_gen = of those, how many rank in energy top-10) ===")
print(genrank.round(1).to_string(index=False))
print(f"\nwrote {xlsx}")
