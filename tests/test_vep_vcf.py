import pytest
import tests.config as config
from vepvcf.vep_vcf import VEP_VCF, parse_filter_string, eval_filter


class TestParseFilterString(object):
    def test_basic(self):
        assert parse_filter_string("#AF < 0.1") == ('allele_data["AF"] < 0.1', ["AF"])

    def test_multiple(self):
        assert parse_filter_string("#foo == 1 and #bar == 2") == (
            'allele_data["foo"] == 1 and allele_data["bar"] == 2',
            ["foo", "bar"],
        )

    def test_replacement(self):
        assert parse_filter_string("@lof") == (
            '(allele_data["LoF"] == "HC" or allele_data["Consequence"] == "start_lost")',
            ["LoF", "Consequence"],
        )

    def test_plus(self):
        assert parse_filter_string("@stop_gained+") == (
            'allele_data["Consequence"] in ["splice_acceptor_variant","splice_donor_variant","stop_gained","transcript_ablation"]',
            ["Consequence"],
        )

    def test_file(self, filter_file):
        assert parse_filter_string("#Consequence in ~" + filter_file) == (
            'allele_data["Consequence"] in TILDE_FILE_DATA["'
            + config.TEST_DIR
            + '/data/filter_file.txt"]',
            ["Consequence"],
        )

    def test_file_fail(self):
        with pytest.raises(IOError, match=r"No such file or directory"):
            parse_filter_string("~foobarfoobar")


class TestEvalFilter(object):
    def test_basic(self):
        f = parse_filter_string('#foo == "bar"')[0]

        assert eval_filter(f, {"foo": "bar"}) == True
        assert eval_filter(f, {"foo": "barb"}) == False

    def test_file(self, filter_file):
        f = parse_filter_string("#Consequence in ~" + filter_file)[0]

        assert eval_filter(f, {"Consequence": "foo"}) == False
        assert eval_filter(f, {"Consequence": "intron_variant"}) == True


class TestVEPVCF(object):
    def test_init(self, vcf_file):
        v = VEP_VCF(vcf_file=vcf_file)
        assert v.__class__.__name__ == "VEP_VCF"
        assert v.vcf.__class__.__name__ == "VCF"
        assert v.vep_info_field == "CSQ"
        assert len(v.vep_cols) == 31
        assert len(v.info_fields) == 42

    def test_failed_init(self):
        with pytest.raises(IOError, match=r"Error opening"):
            v = VEP_VCF(vcf_file="foo")

    def test_vep_info_field(self, info_field_vcf_file):
        v = VEP_VCF(vcf_file=info_field_vcf_file, vep_info_field="FOO")
        assert len(v.vep_cols) == 31

    def test_no_vep_info_field(self, vcf_file):
        v = VEP_VCF(vcf_file=vcf_file, vep_info_field="FOO")
        assert len(v.vep_cols) == 0

    def test_parse_info_data(self, test_vep_vcf, ma_var):
        # no info_fields given, gets all
        parsed_data = test_vep_vcf.parse_info_data(ma_var, [])
        assert sorted(parsed_data.keys()) == [1, 2]
        assert len(parsed_data[1].keys()) == 18

        # info_fields given, only returns those we asked for
        assert test_vep_vcf.parse_info_data(
            ma_var, ["AC", "EX_TARGET", "GRCH37_POS"]
        ) == {
            1: {"AC": 4, "EX_TARGET": True, "GRCH37_POS": 97544543},
            2: {"AC": 49, "EX_TARGET": True, "GRCH37_POS": 97544543},
        }

    def test_parse_vep_data(self, test_vep_vcf, ma_var):
        vep_data = test_vep_vcf.parse_vep_data(ma_var)

        assert sorted(vep_data.keys()) == [1]
        assert len(vep_data[1]) == 1
        assert len(vep_data[1][0].keys()) == 31
        assert vep_data[1][0]["BIOTYPE"] == "protein_coding"

    def test_parse_vep_data_no_vep(self, vcf_file, ma_var):
        vcf = VEP_VCF(vcf_file, vep_info_field="foo")
        vep_data = vcf.parse_vep_data(ma_var)
        assert vep_data == {}

    def test_parse_vcf_data(self, test_vep_vcf, ma_var):
        vcf_data = test_vep_vcf.parse_vcf_data(ma_var, 0)

        assert vcf_data == {
            "ALT": "A",
            "CHROM": "1",
            "FILTER": None,
            "ID": "rs114096998",
            "POS": 97078987,
            "QUAL": 100.0,
            "REF": "G",
        }

    # 	def test_no_allele_num(self)

    def test_get_all_allele_data(self, test_vep_vcf, ma_var):
        ad = test_vep_vcf.get_all_allele_data(ma_var, ["AC"])

        assert sorted(ad.keys()) == [1, 2]
        assert len(ad[1]) == 1
        assert len(ad[1][0].keys()) == 39
        assert ad[1][0]["BIOTYPE"] == "protein_coding"
        assert ad[1][0]["AC"] == 4
        assert ad[1][0]["POS"] == 97078987

    def test_get_all_allele_data_no_vep(self, vcf_file, ma_var):
        vcf = VEP_VCF(vcf_file, vep_info_field="foo")
        ad = vcf.get_all_allele_data(ma_var, ["AC"])

        assert sorted(ad.keys()) == [1, 2]
        assert len(ad[1]) == 1
        assert len(ad[1][0].keys()) == 8
        assert "BIOTYPE" not in ad[1][0]
        assert ad[1][0]["AC"] == 4
        assert ad[1][0]["POS"] == 97078987
