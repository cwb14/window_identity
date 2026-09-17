"""ks_summary.py: near-zero Ks must be kept, not silently dropped.

KaKs_Calculator 3.0 prints Ks as 'NA' whenever the estimate is < 1e-6 (base.cpp, parseOutput).
It is a display convention, not a failure flag: no method assigns Ks its NA sentinel, and with
zero synonymous differences YN00 falls back to Jukes-Cantor, which returns exactly 0. Dropping
every 'NA' removed 51% of PaA/PiA gene pairs, all near-identical, and biased the median upward.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import ks_summary as ksum  # noqa: E402


def row(name, ks, syn_subs, ka="0.01", subs="3"):
    r = dict.fromkeys(ksum.KAKS_COLUMNS, "NA")
    r.update({"Sequence": name, "Method": "YN", "Ka": ka, "Ks": ks,
              "Substitutions": subs, "Syn-Subs": syn_subs, "Length": "900"})
    return r


def write_table(path, rows, header=True):
    with open(path, "w") as fh:
        if header:
            fh.write("\t".join(ksum.KAKS_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in ksum.KAKS_COLUMNS) + "\n")


def test_numeric_ks_is_returned_as_is():
    assert ksum.parse_ks(row("a-b", "0.0870", "28")) == 0.087


def test_na_ks_without_synonymous_substitutions_is_zero():
    # Identical CDS: Ka, Ks and Syn-Subs all NA, 0 substitutions.
    assert ksum.parse_ks(row("a-b", "NA", "NA", ka="NA", subs="0")) == 0.0
    # Nonsynonymous-only differences: Ka > 0, but no synonymous change.
    assert ksum.parse_ks(row("a-b", "NA", "NA", ka="0.00127", subs="3")) == 0.0


def test_na_ks_is_zero_whatever_syn_subs_says():
    # 'NA' only ever means estimate < 1e-6; Syn-Subs is irrelevant to that.
    assert ksum.parse_ks(row("a-b", "NA", "5.18556e-05")) == 0.0
    assert ksum.parse_ks(row("a-b", "NA", "3.2")) == 0.0


def test_unparseable_ks_is_not_silently_zero():
    assert ksum.parse_ks(row("a-b", "", "3")) is None
    assert ksum.parse_ks(row("a-b", "garbage", "3")) is None


def test_read_ks_keeps_zeros_and_drops_unparseable_and_saturation(tmp_path):
    p = tmp_path / "A.B.kaks.tsv"
    write_table(p, [row("a1-b1", "0.1", "20"),
                    row("a2-b2", "NA", "NA", ka="NA", subs="0"),   # identical -> 0
                    row("a3-b3", "garbage", "3"),                  # unparseable
                    row("a4-b4", "2.5", "400")])                   # saturated
    assert ksum.read_ks(str(p), max_ks=2.0, verbose=False) == [0.1, 0.0]


def test_headerless_table_is_read_with_zeros(tmp_path):
    p = tmp_path / "A.B.kaks.tsv"
    write_table(p, [row("a1-b1", "NA", "NA", ka="NA", subs="0"),
                    row("a2-b2", "0.2", "30")], header=False)
    assert ksum.read_ks(str(p), max_ks=2.0, verbose=False) == [0.0, 0.2]


def test_iter_ks_yields_pair_names(tmp_path):
    p = tmp_path / "A.B.kaks.tsv"
    write_table(p, [row("a1-b1", "0.3", "40"), row("a2-b2", "NA", "7"),
                    row("a3-b3", "garbage", "7")])
    assert list(ksum.iter_ks(str(p))) == [("a1-b1", 0.3), ("a2-b2", 0.0), ("a3-b3", None)]
