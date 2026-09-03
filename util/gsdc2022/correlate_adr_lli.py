#!/usr/bin/env python3
"""逐观测匹配 Android ADR state 与转换后的 RINEX 3 LLI。"""

import argparse
import math
from collections import Counter, defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

from analyze_android_adr import decode_state, read_raw_rows
from analyze_rinex_lli import read_observation_types


SPEED_OF_LIGHT_MPS = 299_792_458.0
GPS_UTC_OFFSET_SECONDS = 18
CONSTELLATION_PREFIX = {1: "G", 2: "S", 3: "R", 5: "C", 6: "E", 7: "I"}


def android_satellite_id(constellation_type: int, svid: int):
    if constellation_type == 4:
        qzss_prn = svid - 192 if svid >= 193 else svid
        return f"J{qzss_prn:02d}"
    prefix = CONSTELLATION_PREFIX.get(constellation_type)
    return f"{prefix}{svid:02d}" if prefix else None


def rinex_epoch_utc_ms(line: str):
    fields = line.split()
    year, month, day, hour, minute = map(int, fields[1:6])
    second = float(fields[6])
    gpst = datetime(year, month, day, hour, minute, tzinfo=timezone.utc) + timedelta(
        seconds=second
    )
    return round((gpst - timedelta(seconds=GPS_UTC_OFFSET_SECONDS)).timestamp() * 1000)


def read_rinex_phases(path: Path):
    """建立 (UTC毫秒,卫星) 到相位观测列表的索引。"""
    observation_types = read_observation_types(path)
    phases = defaultdict(list)
    current_time_ms = None
    in_body = False

    with path.open(encoding="ascii", errors="ignore") as stream:
        for line in stream:
            if not in_body:
                label = line[60:].strip() if len(line) >= 60 else ""
                if label == "END OF HEADER":
                    in_body = True
                continue
            if line.startswith(">"):
                current_time_ms = rinex_epoch_utc_ms(line)
                continue
            satellite = line[:3]
            if current_time_ms is None or len(satellite) != 3 or not satellite[1:].isdigit():
                continue

            satellite_phases = []
            malformed = False
            for index, observation_type in enumerate(observation_types.get(satellite[0], [])):
                if not observation_type.startswith("L"):
                    continue
                field = line[3 + 16 * index : 3 + 16 * (index + 1)]
                value_text = field[:14].strip()
                if not value_text:
                    continue
                lli = field[14] if len(field) >= 15 and field[14].strip() else "blank"
                try:
                    value = float(value_text)
                except ValueError:
                    # 少量极端异常观测会溢出RINEX的14字符数值栏，使整行错位。
                    malformed = True
                    break
                satellite_phases.append((observation_type, value, lli))
            if not malformed:
                phases[(current_time_ms, satellite)].extend(satellite_phases)
    return phases


def correlate(
    log_path: Path,
    rinex_path: Path,
    tolerance_cycles: float,
    time_tolerance_ms: int,
):
    phases = read_rinex_phases(rinex_path)
    pair_counts = Counter()
    type_counts = Counter()
    examples = []
    raw_with_adr = 0
    matched = 0

    for row in read_raw_rows(log_path):
        state_text = row.get("AccumulatedDeltaRangeState", "").strip()
        adr_text = row.get("AccumulatedDeltaRangeMeters", "").strip()
        frequency_text = row.get("CarrierFrequencyHz", "").strip()
        if not state_text or not adr_text or not frequency_text:
            continue
        raw_with_adr += 1

        satellite = android_satellite_id(
            int(row["ConstellationType"]), int(row["Svid"])
        )
        if satellite is None:
            continue
        raw_time_ms = int(row["utcTimeMillis"])
        candidates = []
        for offset_ms in range(-time_tolerance_ms, time_tolerance_ms + 1):
            candidates.extend(phases.get((raw_time_ms + offset_ms, satellite), []))
        if not candidates:
            continue

        frequency_hz = float(frequency_text)
        adr_cycles = float(adr_text) / (SPEED_OF_LIGHT_MPS / frequency_hz)
        observation_type, phase_cycles, lli = min(
            candidates, key=lambda item: abs(item[1] - adr_cycles)
        )
        difference = abs(phase_cycles - adr_cycles)
        if not math.isfinite(difference) or difference > tolerance_cycles:
            continue

        state = int(state_text)
        matched += 1
        pair_counts[(state, lli)] += 1
        type_counts[(observation_type, state, lli)] += 1
        if len(examples) < 8:
            examples.append(
                (raw_time_ms, satellite, observation_type, adr_cycles, phase_cycles, difference, state, lli)
            )

    return raw_with_adr, matched, pair_counts, type_counts, examples


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", required=True, type=Path, help="Android gnss_log.txt")
    parser.add_argument("--rinex", required=True, type=Path, help="对应 gnss_rinex.*o")
    parser.add_argument(
        "--tolerance-cycles", type=float, default=0.1, help="ADR与RINEX相位最大差，默认0.1周"
    )
    parser.add_argument(
        "--time-tolerance-ms",
        type=int,
        default=2,
        help="Raw与RINEX纪元的最大时间差，默认前后2毫秒",
    )
    args = parser.parse_args()

    raw_count, matched, pair_counts, type_counts, examples = correlate(
        args.log, args.rinex, args.tolerance_cycles, args.time_tolerance_ms
    )
    print(f"含ADR的Raw记录: {raw_count}")
    print(f"与RINEX相位精确匹配: {matched} ({matched / raw_count:.2%})")
    print("\nADR state -> RINEX LLI：")
    states = sorted({state for state, lli in pair_counts})
    llis = sorted({lli for state, lli in pair_counts})
    for state in states:
        total = sum(pair_counts[(state, lli)] for lli in llis)
        mapping = ", ".join(
            f"LLI={lli}: {pair_counts[(state, lli)]} ({pair_counts[(state, lli)] / total:.2%})"
            for lli in llis
            if pair_counts[(state, lli)]
        )
        print(f"  {state:>3} {decode_state(state):<55} n={total:<7} {mapping}")

    print("\n匹配示例：")
    for time_ms, satellite, obs_type, adr, phase, difference, state, lli in examples:
        print(
            f"  {time_ms} {satellite} {obs_type}: ADR={adr:.3f}周 "
            f"RINEX={phase:.3f}周 差={difference:.6f}周 state={state} LLI={lli}"
        )


if __name__ == "__main__":
    main()
