#!/usr/bin/env python3
"""Descriptive tables for a two-sample SV VCF (stdlib only, read-only on the input).

Usage:
    33_sv_describe.py --vcf CALLS.vcf[.gz] --out-prefix PREFIX [--autosomes-only]

Writes PREFIX.type_size.tsv, PREFIX.per_chrom.tsv, PREFIX.per_sample.tsv,
PREFIX.sharing.tsv and PREFIX.size_log10.tsv. Existing output files are never
overwritten. Nothing is executed from the input; it is only parsed as text.
"""
import argparse, collections, gzip, math, os, re, sys

BINS = [(0, 50, "<50"), (50, 100, "50-99"), (100, 500, "100-499"), (500, 5000, "500-4999"),
        (5000, 100000, "5000-99999"), (100000, float("inf"), ">=100000")]
BIN_ORDER = [b[2] for b in BINS] + ["no length"]
AUTOSOME = re.compile(r"^NC_0308(0[89]|[12][0-9]|3[0-6])\.1$")  # NC_030808.1 .. NC_030836.1


def size_class(n):
    for lo, hi, name in BINS:
        if lo <= n < hi:
            return name
    raise ValueError(n)


def parse_info(s):
    d = {}
    for item in s.split(";"):
        if not item:
            continue
        k, _, v = item.partition("=")
        d[k] = v
    return d


def sv_length(info, pos, ref, alt):
    """Absolute SV length, or None when the record has no length (e.g. breakends)."""
    t = info.get("SVTYPE", "")
    if t == "BND":
        return None
    v = info.get("SVLEN")
    if v:
        try:
            return abs(int(float(v.split(",")[0])))
        except ValueError:
            pass
    if "END" in info:
        try:
            return abs(int(info["END"]) - pos)
        except ValueError:
            pass
    if alt and not alt.startswith("<") and "[" not in alt and "]" not in alt:
        return abs(len(alt) - len(ref))
    return None


def gt_state(gt):
    """ref / het / hom / missing from a GT string (phased or not)."""
    alleles = re.split(r"[/|]", gt.split(":")[0])
    if any(a in (".", "") for a in alleles):
        return "missing"
    try:
        nz = [int(a) != 0 for a in alleles]
    except ValueError:
        return "missing"
    if not any(nz):
        return "ref"
    return "hom" if all(nz) else "het"


def _open(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


def describe(path, autosomes_only=False):
    r = {"n_records": 0, "n_malformed": 0, "type_counts": collections.Counter(),
         "type_size": collections.Counter(), "per_chrom": collections.defaultdict(collections.Counter),
         "per_sample": {}, "sharing": collections.Counter(), "log10": collections.Counter(),
         "samples": []}
    samples = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                r["samples"] = samples
                r["per_sample"] = {s: collections.Counter() for s in samples}
                continue
            if len(f) < 10:
                r["n_malformed"] += 1
                continue
            chrom = f[0]
            if autosomes_only and not AUTOSOME.match(chrom):
                continue
            try:
                pos = int(f[1])
            except ValueError:
                r["n_malformed"] += 1
                continue
            info = parse_info(f[7])
            t = info.get("SVTYPE", "UNK")
            n = sv_length(info, pos, f[3], f[4])
            r["n_records"] += 1
            r["type_counts"][t] += 1
            r["per_chrom"][chrom][t] += 1
            r["type_size"][(t, "no length" if n is None else size_class(n))] += 1
            if n is not None and n > 0:
                r["log10"][(t, round(math.floor(math.log10(n) * 10) / 10, 1))] += 1
            fmt = f[8].split(":")
            gi = fmt.index("GT") if "GT" in fmt else None
            states = []
            for s, col in zip(samples, f[9:]):
                st = gt_state(col.split(":")[gi]) if gi is not None and len(col.split(":")) > gi else "missing"
                states.append(st)
                r["per_sample"][s][st] += 1
                if st in ("het", "hom"):
                    r["per_sample"][s][t] += 1
            if len(samples) == 2:
                a, b = (st in ("het", "hom") for st in states)
                key = "shared" if a and b else "A_only" if a else "B_only" if b else "neither"
                r["sharing"][(t, key)] += 1
    return r


def _w(path, rows):
    # 'x' mode: refuse to overwrite an existing file
    with open(path, "x") as fh:
        for row in rows:
            fh.write("\t".join(str(c) for c in row) + "\n")


def write_tables(r, prefix):
    types = sorted(r["type_counts"])
    _w(prefix + ".type_size.tsv", [["type"] + BIN_ORDER + ["total"]] + [
        [t] + [r["type_size"][(t, b)] for b in BIN_ORDER] + [r["type_counts"][t]] for t in types])
    _w(prefix + ".per_chrom.tsv", [["chrom"] + types + ["total"]] + [
        [c] + [r["per_chrom"][c][t] for t in types] + [sum(r["per_chrom"][c].values())]
        for c in sorted(r["per_chrom"])])
    _w(prefix + ".per_sample.tsv", [["sample", "het", "hom", "ref", "missing"] + types] + [
        [s, r["per_sample"][s]["het"], r["per_sample"][s]["hom"], r["per_sample"][s]["ref"],
         r["per_sample"][s]["missing"]] + [r["per_sample"][s][t] for t in types] for s in r["samples"]])
    keys = ["shared", "A_only", "B_only", "neither"]
    _w(prefix + ".sharing.tsv", [["type"] + keys] + [[t] + [r["sharing"][(t, k)] for k in keys] for t in types])
    bins = sorted({b for (_, b) in r["log10"]})
    _w(prefix + ".size_log10.tsv", [["log10_size_lower"] + types] + [
        [b] + [r["log10"][(t, b)] for t in types] for b in bins])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--autosomes-only", action="store_true")
    a = ap.parse_args(argv)
    try:
        r = describe(a.vcf, a.autosomes_only)
        write_tables(r, a.out_prefix)
    except (OSError, EOFError) as e:  # unreadable/corrupt input, or an output file already exists
        print("error: %s" % e, file=sys.stderr)
        return 2
    print("records=%d malformed=%d samples=%s (A=%s, B=%s)" % (
        r["n_records"], r["n_malformed"], ",".join(r["samples"]),
        *(r["samples"] + ["", ""])[:2]), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
