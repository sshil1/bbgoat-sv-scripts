"""Unit tests for scripts/33_sv_describe.py (stdlib only). Run: python3 -I -m unittest discover -s tests"""
import gzip, importlib.util, io, os, tempfile, unittest

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("svd", os.path.join(HERE, "..", "scripts", "33_sv_describe.py"))
svd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(svd)

HDR = ("##fileformat=VCFv4.2\n"
       '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="t">\n'
       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n")


def rec(chrom, pos, svtype, svlen, gta, gtb, ref="N", alt="<DEL>", extra=""):
    info = "SVTYPE=%s" % svtype + (";SVLEN=%s" % svlen if svlen is not None else "") + extra
    return "\t".join([chrom, str(pos), ".", ref, alt, ".", "PASS", info, "GT", gta, gtb]) + "\n"


def write(text, gz=False):
    fd, path = tempfile.mkstemp(suffix=".vcf.gz" if gz else ".vcf")
    os.close(fd)
    with (gzip.open if gz else open)(path, "wt") as fh:
        fh.write(text)
    return path


class SizeClass(unittest.TestCase):
    def test_edges(self):
        self.assertEqual(svd.size_class(50), "50-99")
        self.assertEqual(svd.size_class(99), "50-99")
        self.assertEqual(svd.size_class(100), "100-499")
        self.assertEqual(svd.size_class(4999), "500-4999")
        self.assertEqual(svd.size_class(5000), "5000-99999")
        self.assertEqual(svd.size_class(100000), ">=100000")
        self.assertEqual(svd.size_class(49), "<50")


class SvLen(unittest.TestCase):
    def test_svlen_negative_and_list(self):
        self.assertEqual(svd.sv_length({"SVTYPE": "DEL", "SVLEN": "-120"}, 100, "N", "<DEL>"), 120)
        self.assertEqual(svd.sv_length({"SVTYPE": "INS", "SVLEN": "300,300"}, 100, "N", "<INS>"), 300)

    def test_fallback_end(self):
        self.assertEqual(svd.sv_length({"SVTYPE": "DEL", "END": "1099"}, 100, "N", "<DEL>"), 999)

    def test_fallback_alleles(self):
        self.assertEqual(svd.sv_length({"SVTYPE": "INS"}, 100, "A", "A" + "C" * 60), 60)

    def test_breakend_has_no_length(self):
        self.assertIsNone(svd.sv_length({"SVTYPE": "BND"}, 100, "N", "N[chr2:5["))


class Genotype(unittest.TestCase):
    def test_states(self):
        self.assertEqual(svd.gt_state("0/0"), "ref")
        self.assertEqual(svd.gt_state("0|1"), "het")
        self.assertEqual(svd.gt_state("1/1"), "hom")
        self.assertEqual(svd.gt_state("./."), "missing")
        self.assertEqual(svd.gt_state("."), "missing")
        self.assertEqual(svd.gt_state("1/2"), "hom")
        self.assertEqual(svd.gt_state("0/."), "missing")


class Describe(unittest.TestCase):
    def setUp(self):
        text = HDR
        text += rec("NC_030808.1", 100, "DEL", -60, "0/1", "0/1")          # shared, 50-99
        text += rec("NC_030808.1", 900, "INS", 250, "1/1", "./.", alt="<INS>")  # A only
        text += rec("NC_030809.1", 50, "DEL", -7000, "./.", "0/1")          # B only
        text += rec("NC_030809.1", 80, "INV", 2000, "0/1", "0/1", alt="<INV>")
        text += rec("NW_017189516.1", 5, "DEL", -80, "0/1", "0/0")          # X scaffold
        text += rec("NC_030810.1", 10, "BND", None, "0/1", "0/1", alt="N[x:1[")
        self.path = write(text)

    def tearDown(self):
        os.unlink(self.path)

    def test_counts_all_contigs(self):
        r = svd.describe(self.path, autosomes_only=False)
        self.assertEqual(r["n_records"], 6)
        self.assertEqual(r["type_counts"]["DEL"], 3)
        self.assertEqual(r["type_counts"]["INS"], 1)

    def test_autosomes_only_drops_scaffold(self):
        r = svd.describe(self.path, autosomes_only=True)
        self.assertEqual(r["n_records"], 5)
        self.assertEqual(r["type_counts"]["DEL"], 2)

    def test_type_size_matrix(self):
        r = svd.describe(self.path, autosomes_only=True)
        self.assertEqual(r["type_size"][("DEL", "50-99")], 1)
        self.assertEqual(r["type_size"][("DEL", "5000-99999")], 1)
        self.assertEqual(r["type_size"][("INS", "100-499")], 1)
        self.assertEqual(r["type_size"][("BND", "no length")], 1)

    def test_sample_states(self):
        r = svd.describe(self.path, autosomes_only=True)
        a = r["per_sample"]["A"]
        self.assertEqual((a["het"], a["hom"], a["missing"]), (3, 1, 1))
        b = r["per_sample"]["B"]
        self.assertEqual(b["missing"], 1)

    def test_sharing(self):
        r = svd.describe(self.path, autosomes_only=True)
        s = r["sharing"]
        self.assertEqual(s[("DEL", "shared")], 1)
        self.assertEqual(s[("INS", "A_only")], 1)
        self.assertEqual(s[("DEL", "B_only")], 1)

    def test_per_chrom(self):
        r = svd.describe(self.path, autosomes_only=True)
        self.assertEqual(r["per_chrom"]["NC_030808.1"]["DEL"], 1)
        self.assertEqual(r["per_chrom"]["NC_030809.1"]["INV"], 1)

    def test_gz_input(self):
        with open(self.path) as fh:
            p = write(fh.read(), gz=True)
        try:
            self.assertEqual(svd.describe(p, autosomes_only=True)["n_records"], 5)
        finally:
            os.unlink(p)

    def test_malformed_line_is_counted_not_fatal(self):
        p = write(HDR + "chr\tonly\tthree\n" + rec("NC_030808.1", 1, "DEL", -60, "0/1", "0/1"))
        try:
            r = svd.describe(p, autosomes_only=False)
            self.assertEqual(r["n_records"], 1)
            self.assertEqual(r["n_malformed"], 1)
        finally:
            os.unlink(p)


class Output(unittest.TestCase):
    def test_write_tables_refuses_to_overwrite(self):
        d = tempfile.mkdtemp()
        p = write(HDR + rec("NC_030808.1", 1, "DEL", -60, "0/1", "0/1"))
        try:
            r = svd.describe(p, autosomes_only=False)
            prefix = os.path.join(d, "out")
            svd.write_tables(r, prefix)
            self.assertTrue(os.path.exists(prefix + ".type_size.tsv"))
            with self.assertRaises(FileExistsError):
                svd.write_tables(r, prefix)
        finally:
            os.unlink(p)


if __name__ == "__main__":
    unittest.main()
