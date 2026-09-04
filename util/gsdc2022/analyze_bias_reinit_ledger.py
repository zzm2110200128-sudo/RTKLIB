#!/usr/bin/env python3
"""E7c udbias_ppp() 模糊度重建账本统计。"""

import argparse
import csv
import math
import statistics
from collections import Counter
from pathlib import Path

STEADY_S = 300.0


def pct(values, p):
    if not values:
        return None
    values = sorted(values)
    x = (len(values) - 1) * p
    lo, hi = math.floor(x), math.ceil(x)
    return values[lo] if lo == hi else values[lo] * (hi - x) + values[hi] * (x - lo)


def fmt(value):
    return "" if value is None else f"{value:.3f}"


def load(path):
    with path.open(newline="", encoding="utf-8", errors="replace") as stream:
        reader = csv.DictReader(stream)
        if len(reader.fieldnames or []) != 27:
            raise ValueError(f"{path}: 预期 27 列，实际 {len(reader.fieldnames or [])}")
        rows = list(reader)
    if rows:
        t0 = float(rows[0]["tow"])
        for row in rows:
            row["_rel"] = float(row["tow"]) - t0
    return rows


def validate(rows, label):
    bad_p = sum(row["did_init"] == "1" and abs(float(row["p_post"]) - 3600.0) > 1e-6
                for row in rows)
    bad_action = sum((row["did_init"] == "1") != row["action"].startswith(("init_", "reinit_"))
                     for row in rows)
    bad_missing = sum(row["action"] == "no_bias" and
                      (float(row["Lc"]) != 0.0 or float(row["Pc"]) != 0.0) for row in rows)
    if bad_p or bad_action or bad_missing:
        raise ValueError(f"{label}: p_post={bad_p}, action={bad_action}, no_bias_LcPc={bad_missing}")


def summarize(label, rows):
    print(f"\n=== {label} ===")
    for title, pool in (("ALL", rows), (">=300s", [r for r in rows if r["_rel"] >= STEADY_S])):
        actions = Counter(r["action"] for r in pool)
        init = [r for r in pool if r["did_init"] == "1"]
        sources = Counter(r["slip_source"] for r in init)
        continued = [r for r in init if float(r["x_pre"]) != 0.0]
        deltas = [abs(float(r["dx_init"])) for r in continued]
        print(f"\n-- {title}: rows={len(pool)} init={len(init)} --")
        print("actions:", " ".join(f"{k}={v}" for k, v in actions.most_common()))
        print("init sources:", " ".join(f"{k}={v}" for k, v in sources.most_common()))
        print(f"prior_nonzero={len(continued)} |new_bias-old_bias| P50/P95="
              f"{fmt(pct(deltas,.5))}/{fmt(pct(deltas,.95))} m; "
              f"cleared={sum(r['cleared']=='1' for r in pool)}; "
              f"phase-code-jump epochs={len({r['tow'] for r in pool if float(r['jump_corr']) != 0.0})}")
        print(f"{'sys':>4} {'source':<16} {'init':>7} {'prior!=0':>9} {'delta P50':>10} {'delta P95':>10}")
        for sysc in ("G", "E"):
            for source in ("raw_bit0", "raw_half_only", "detector_only", "clean"):
                group = [r for r in init if r["sys"] == sysc and r["slip_source"] == source]
                if not group:
                    continue
                old = [r for r in group if float(r["x_pre"]) != 0.0]
                ds = [abs(float(r["dx_init"])) for r in old]
                print(f"{sysc:>4} {source:<16} {len(group):>7} {len(old):>9} "
                      f"{fmt(pct(ds,.5)):>10} {fmt(pct(ds,.95)):>10}")


def parse_case(text):
    parts = text.split("|", 1)
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("格式：LABEL|bias.csv")
    return parts[0], Path(parts[1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", type=parse_case, required=True)
    args = parser.parse_args()
    for label, path in args.case:
        rows = load(path)
        validate(rows, label)
        summarize(label, rows)


if __name__ == "__main__":
    main()
