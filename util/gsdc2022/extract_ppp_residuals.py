#!/usr/bin/env python3
"""把 RTKLIB PPP .pos.stat 中的 $SAT 记录提取为结构化 CSV。

字段语义（核对 src/ppp.c）：
- 无电离层组合（IFLC）模式下 NF=1：每星每历元只输出 1 行，
  frequency_slot 列的取值恒为打印的 j+1=1，表示唯一的 IF 组合槽位，
  不是 L2/L5 之类的具体频点，无法据此做 L1/L5 分频统计；
- pseudorange_residual_m 与 carrier_phase_residual_m 分别来自
  ppp_res() 中对 IF 组合伪距 Pc（code=1）和 IF 组合相位 Lc（code=0）
  计算的残差 resp[j]/resc[j]（j=0）。ppp_res 每历元会先做验前残差、
  滤波后再做验后残差，数组被多次覆盖，故文件中保存的是最后一次（验后）
  通过的值：code 列≈码噪声/多路径（m 级），phase 列被滤波几乎完全吸收，
  通常 ≤0.001 m 甚至精确为 0，不能代表真实相位噪声；
- vsat=1 表示该星相位残差成功入列（相位可用），与伪距可用性无关；
- cn0_dbhz 为 L1 槽位的 C/N0（dB-Hz）。
"""

import argparse
import csv
from pathlib import Path


FIELD_NAMES = [
    "gps_week",
    "gps_tow_s",
    "gps_time_s",
    "satellite",
    "system",
    "frequency_slot",
    "azimuth_deg",
    "elevation_deg",
    "pseudorange_residual_m",
    "carrier_phase_residual_m",
    "vsat",
    "cn0_dbhz",
    "ambiguity_fix",
    "slip",
    "lock_count",
    "outage_count",
    "slip_count",
    "reject_count",
    "phase_bias_m",
    "phase_bias_variance_m2",
    "interchannel_bias",
]


def parse_sat_line(line: str, path: Path, line_number: int):
    """按 src/ppp.c 中 $SAT 的实际输出顺序解析一行。"""
    fields = line.rstrip().split(",")
    if len(fields) != 20 or fields[0] != "$SAT":
        raise ValueError(f"{path}:{line_number} 不是预期的 PPP $SAT 格式：{line.rstrip()}")

    gps_week = int(fields[1])
    gps_tow = float(fields[2])
    satellite = fields[3]
    return {
        "gps_week": gps_week,
        "gps_tow_s": gps_tow,
        "gps_time_s": gps_week * 604800.0 + gps_tow,
        "satellite": satellite,
        "system": satellite[0],
        "frequency_slot": int(fields[4]),
        "azimuth_deg": float(fields[5]),
        "elevation_deg": float(fields[6]),
        "pseudorange_residual_m": float(fields[7]),
        "carrier_phase_residual_m": float(fields[8]),
        "vsat": int(fields[9]),
        "cn0_dbhz": float(fields[10]),
        "ambiguity_fix": int(fields[11]),
        "slip": int(fields[12]),
        "lock_count": int(fields[13]),
        "outage_count": int(fields[14]),
        "slip_count": int(fields[15]),
        "reject_count": int(fields[16]),
        "phase_bias_m": float(fields[17]),
        "phase_bias_variance_m2": float(fields[18]),
        "interchannel_bias": float(fields[19]),
    }


def read_sat_records(path: Path):
    records = []
    with path.open(encoding="utf-8", errors="ignore") as stream:
        for line_number, line in enumerate(stream, start=1):
            if line.startswith("$SAT,"):
                records.append(parse_sat_line(line, path, line_number))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stat", required=True, type=Path, help="输入 .pos.stat 文件")
    parser.add_argument("--output", required=True, type=Path, help="输出 CSV 文件")
    args = parser.parse_args()

    records = read_sat_records(args.stat)
    if not records:
        raise SystemExit(f"没有找到 $SAT 记录：{args.stat}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELD_NAMES)
        writer.writeheader()
        writer.writerows(records)

    valid_count = sum(record["vsat"] == 1 for record in records)
    invalid_count = len(records) - valid_count
    slip_count = sum(record["slip"] != 0 for record in records)
    satellites = {record["satellite"] for record in records}

    print(f"输入文件: {args.stat}")
    print(f"输出文件: {args.output}")
    print(f"$SAT 总记录数: {len(records)}")
    print(f"卫星数: {len(satellites)}")
    print(f"vsat=1: {valid_count} ({valid_count / len(records):.2%})")
    print(f"vsat=0: {invalid_count} ({invalid_count / len(records):.2%})")
    print(f"slip!=0: {slip_count} ({slip_count / len(records):.2%})")


if __name__ == "__main__":
    main()
