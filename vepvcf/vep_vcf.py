#!/usr/bin/env python3

import re
from cyvcf2 import VCF

# used as shortcuts in filter strings
FILTER_REPLACEMENTS = {
    "@lof": '(#LoF == "HC" or #Consequence == "start_lost")',
    "@pathogenic": '#clinvar_CLNSIG[-9:] == "athogenic"',
}

# used to generate shortcuts based on VEP SO term consequences
# spit out from VEP module with the following perl one-liner:
# perl -I ~/software/ensembl-vep/ -MBio::EnsEMBL::Variation::Utils::Constants -e 'for my $c(sort {$a->{rank} <=> $b->{rank}} values(%Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES)) { print "\"".$c->{SO_term}."\": ".$c->{rank}.",\n"}' > vep_consequence_ranks.txt
VEP_SO_TERMS = {
    "transcript_ablation": 1,
    "splice_acceptor_variant": 3,
    "splice_donor_variant": 3,
    "stop_gained": 4,
    "frameshift_variant": 5,
    "stop_lost": 6,
    "start_lost": 7,
    "transcript_amplification": 8,
    "inframe_insertion": 10,
    "inframe_deletion": 11,
    "protein_altering_variant": 12,
    "missense_variant": 12,
    "splice_region_variant": 13,
    "incomplete_terminal_codon_variant": 14,
    "stop_retained_variant": 15,
    "synonymous_variant": 15,
    "start_retained_variant": 15,
    "coding_sequence_variant": 16,
    "mature_miRNA_variant": 17,
    "5_prime_UTR_variant": 18,
    "3_prime_UTR_variant": 19,
    "non_coding_transcript_exon_variant": 20,
    "intron_variant": 21,
    "NMD_transcript_variant": 22,
    "non_coding_transcript_variant": 23,
    "upstream_gene_variant": 24,
    "downstream_gene_variant": 25,
    "TFBS_ablation": 26,
    "TFBS_amplification": 28,
    "TF_binding_site_variant": 30,
    "regulatory_region_ablation": 31,
    "regulatory_region_amplification": 33,
    "regulatory_region_variant": 36,
    "feature_elongation": 36,
    "feature_truncation": 37,
    "intergenic_variant": 38,
    "sequence_variant": 39,
}

TILDE_FILE_DATA = {}

# generate consequence term based shortcuts
# e.g. @missense_variant+ ===> everything as or more "severe" than missense_variant in the VEP_SO_TERMS list
for t in VEP_SO_TERMS.keys():
    FILTER_REPLACEMENTS["@" + t + "+"] = (
        "#Consequence in ["
        + ",".join(
            list(
                map(
                    lambda x: '"' + x + '"',
                    filter(
                        lambda x: VEP_SO_TERMS[x] <= VEP_SO_TERMS[t],
                        sorted(VEP_SO_TERMS.keys()),
                    ),
                )
            )
        )
        + "]"
    )


def parse_filter_string(filter_string):
    for k, v in FILTER_REPLACEMENTS.items():
        filter_string = filter_string.replace(k, v)
    filter_fields = list(
        map(lambda x: x.replace("#", ""), re.findall(r"#[\w\-\_]+", filter_string))
    )
    filter_eval = re.sub(r"#([\w\-\_]+)", r'allele_data["\1"]', filter_string)

    m = re.search(r"\~(\S+)", filter_string)
    if m:
        global TILDE_FILE_DATA
        for group in m.groups():
            try:
                with open(group, "r") as f:
                    TILDE_FILE_DATA[group] = list(map(lambda x: x.rstrip(), f))
            except:
                raise

        filter_eval = re.sub(r"\~(\S+)", r'TILDE_FILE_DATA["\1"]', filter_eval)

    return (filter_eval, filter_fields)


def eval_filter(filter_eval, allele_data):
    filter_pass = False
    try:
        if eval(filter_eval):
            filter_pass = True
    except:
        pass
    return filter_pass


class VEP_VCF(object):
    def __init__(self, vcf_file, vep_info_field="CSQ"):  # , vep_info_field='CSQ'):
        self.vep_info_field = vep_info_field
        self.vcf = VCF(vcf_file)
        self._vep_cols = None
        self._info_fields = None

    @property
    def vep_cols(self):
        if self._vep_cols is None:
            if self.vep_info_field in self.vcf:
                matches = re.search(
                    'Format: (.+?)"', self.vcf[self.vep_info_field]["Description"]
                )
                if matches:
                    self._vep_cols = matches.groups()[0].split("|")
            else:
                self._vep_cols = []

        return self._vep_cols

    @property
    def info_fields(self):
        if self._info_fields is None:
            info_fields = {}

            for header in self.vcf.header_iter():
                if header.type == "INFO":
                    header_dict = header.info()
                    info_fields[header_dict["ID"]] = header_dict

            self._info_fields = info_fields

        return self._info_fields

    def parse_vep_data(self, var):
        # init data structure to store by allele
        # we're going to key on the ALT allele index but start at 1
        # this way it's a more obvious match up to the genotypes where 0=REF, 1=ALT1, 2=ALT2 etc.
        data = {}
        alleles = {var.REF: 0}
        allele_lengths = set([len(var.REF)])
        allele_prefixes = set([var.REF[0]])
        for i, a in enumerate(var.ALT):
            alleles[a] = i + 1
            allele_lengths.add(len(a))
            allele_prefixes.add(a[0])

        # are any of the alleles longer than 1bp?
        trimmed_alleles = {}
        if filter(lambda x: x > 1, allele_lengths) and len(allele_prefixes) == 1:
            for a in alleles:
                trimmed_alleles["-" if len(a) == 1 else a[1:]] = alleles[a]
        else:
            trimmed_alleles.update(alleles)

        vep_data_string = var.INFO.get(self.vep_info_field)
        if vep_data_string:
            for chunk in vep_data_string.split(","):
                vep_data = {}
                for c, v in zip(self.vep_cols, chunk.split("|")):
                    vep_data[c] = v

                # figure out which allele this is
                allele_num = None
                if "ALLELE_NUM" in vep_data:
                    allele_num = vep_data["ALLELE_NUM"]
                elif len(var.ALT) == 1:
                    allele_num = 1
                elif "Allele" in vep_data:
                    a = vep_data["Allele"]

                    if (a in alleles and a in trimmed_alleles) or (
                        a not in alleles and a not in trimmed_alleles
                    ):
                        raise Exception(
                            "Unable to identify allele from VEP data - you can avoid this issue by running VEP with --allele_number"
                        )

                    allele_num = alleles.get(a, trimmed_alleles[a])
                else:
                    raise Exception(
                        "Unable to identify allele from VEP data - you can avoid this issue by running VEP with --allele_number"
                    )

                allele_num = int(allele_num)
                if allele_num not in data:
                    data[allele_num] = []
                data[allele_num].append(vep_data)

        return data

    def parse_info_data(self, var, info_fields=[]):
        data = {}
        for i, a in enumerate(var.ALT):
            data[i + 1] = {}

        if len(info_fields) == 0:
            info_fields = self.info_fields.keys()

        for info_field in info_fields:
            values = var.INFO.get(info_field)

            if values is not None:
                if type(values) is not tuple:
                    values = tuple([values])

                if len(values) == 1:
                    for k in data:
                        data[k][info_field] = values[0]
                elif len(values) == len(var.ALT):
                    for i, v in enumerate(values):
                        data[i + 1][info_field] = v
                else:
                    raise Exception(
                        "Number of values for INFO field "
                        + info_field
                        + " does not match number of ALTs"
                    )

        return data

    def parse_vcf_data(self, var, alt_index):
        vcf_data = {}
        for f in ("CHROM", "POS", "ID", "REF", "QUAL", "FILTER"):
            vcf_data[f] = getattr(var, f)
        vcf_data["ALT"] = var.ALT[alt_index]

        return vcf_data

    def get_all_allele_data(self, var, filter_fields=[]):
        info_data = self.parse_info_data(var=var, info_fields=filter_fields)

        vep_data = self.parse_vep_data(var=var)

        all_allele_data = {}

        for allele_i in info_data:
            all_allele_data[allele_i] = []

            for allele_vep_data in vep_data.get(allele_i, [{}]):
                allele_vep_data.update(info_data[allele_i])
                allele_vep_data.update(self.parse_vcf_data(var, allele_i - 1))

                all_allele_data[allele_i].append(allele_vep_data)

        return all_allele_data


if __name__ == "__main__":
    main()
