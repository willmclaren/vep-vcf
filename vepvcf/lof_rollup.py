#!/usr/bin/env python3

import argparse
import re
from vepvcf.vep_vcf import VEP_VCF, parse_filter_string, eval_filter
from vepvcf.utils import smart_open, smart_write


def main():
    parser = argparse.ArgumentParser(
        description="Generate LOF rollup file from VEP-annotated VCF file"
    )
    parser.add_argument("vcf", metavar="VCF", nargs="+", help="VCF file")
    parser.add_argument(
        "-i",
        "--info_field",
        default="CSQ",
        help="VCF INFO field name containing VEP annotation (default: %(default)s)",
    )
    parser.add_argument(
        "-f",
        "--filter",
        default='#PICK == "1" and #AF <= 0.05 and @lof',
        help='Filter string eval\'d to select variants, use #[field_name] to denote VEP fields (default: "%(default)s")',
    )
    parser.add_argument("-r", "--region", help="Limit to region [chr:start-end] of VCF")
    parser.add_argument("-o", "--output_file", help="Output file name")
    parser.add_argument(
        "-g", "--group_file", help="Create EPACTS group file with this name"
    )
    parser.add_argument(
        "-q", "--freq_file", help="Create frequency output file with this name"
    )
    parser.add_argument(
        "-p",
        "--pops_file",
        help="Additionally use population assignments in this file for frequency calculations",
    )
    parser.add_argument(
        "-s",
        "--symbol_only",
        action="store_true",
        help="By default genes are shown as [gene_symbol]([ensembl_id]); use this flag to show only [gene_symbol]",
    )
    args = parser.parse_args()

    (filter_eval, filter_fields) = parse_filter_string(args.filter)

    lof_table = {}
    epacts_groups = {}
    samples = []
    gene_AN = {}

    for vcf_file in args.vcf:
        try:
            vep_vcf = VEP_VCF(vcf_file=vcf_file, vep_info_field=args.info_field)
        except:
            raise

        for s in vep_vcf.vcf.samples:
            if s not in samples:
                samples.append(s)

        gene_data = {}
        line_count = 0
        prev_chrom = None

        if args.region:
            vcf_iter = vep_vcf.vcf(args.region)
        else:
            vcf_iter = vep_vcf.vcf

        for var in vcf_iter:

            # process data at each chromosome boundary to save memory
            if prev_chrom is not None and var.CHROM != prev_chrom:
                process_gene_data(
                    gene_data,
                    vep_vcf.vcf.samples,
                    lof_table,
                    epacts_groups,
                    args.group_file,
                )
            prev_chrom = var.CHROM

            var_gene_data = get_var_gene_data(
                vep_vcf=vep_vcf,
                var=var,
                filter_eval=filter_eval,
                filter_fields=filter_fields,
                symbol_only=args.symbol_only,
            )

            # now merge var_gene_data into gene_data
            for gene in var_gene_data:
                if gene not in gene_data:
                    gene_data[gene] = []

                gene_data[gene] += var_gene_data[gene]

                if "AN" in var.INFO:
                    if gene not in gene_AN or var.INFO["AN"] > gene_AN[gene]:
                        gene_AN[gene] = var.INFO["AN"]

        # process remaining data
        process_gene_data(
            gene_data, vep_vcf.vcf.samples, lof_table, epacts_groups, args.group_file
        )

    # dump table
    genes = sorted(lof_table.keys())

    with smart_open(args.output_file) as fh:
        smart_write(fh, "#\t" + "\t".join(genes) + "\n")

        for s in samples:
            smart_write(
                fh,
                s + "\t" + "\t".join(map(lambda x: str(lof_table[x][s]), genes)) + "\n",
            )

        fh.close()

    if args.group_file:
        with smart_open(args.group_file) as fh:
            for k, v in epacts_groups.items():
                smart_write(fh, "\t".join([k] + v) + "\n")
            fh.close()

    if args.freq_file:
        sample_pops = {}
        all_pops = set([])
        all_pops.add("AAA")

        for s in samples:
            sample_pops[s] = set(["AAA"])
            all_pops.add("AAA")

        if args.pops_file:
            with open(args.pops_file, "r") as fh:

                for l in fh:
                    data = l.rstrip().split("\t")

                    s = data.pop(0)
                    if s not in samples:
                        continue

                    for pop in data:
                        sample_pops[s].add(pop)
                        all_pops.add(pop)

        with smart_open(args.freq_file) as fh:

            headers = ["#Gene"]
            for pop in sorted(all_pops):
                if pop == "AAA":
                    pop = "ALL"

                headers += ["AC_" + pop, "AN_" + pop, "AF_" + pop]

            smart_write(fh, "\t".join(headers) + "\n")

            for g in genes:
                by_pop = {}

                for s in samples:
                    v = lof_table[g][s]

                    for pop in sorted(sample_pops[s]):
                        if pop not in by_pop:
                            by_pop[pop] = {"ac": 0, "an": 0}
                        by_pop[pop]["an"] += 2
                        by_pop[pop]["ac"] += v

                # an = gene_AN.get(gene, count)

                smart_write(fh, g)

                for pop in sorted(all_pops):
                    ac = by_pop[pop]["ac"]
                    an = by_pop[pop]["an"]
                    af = float(ac) / float(an)

                    smart_write(fh, "\t" + "\t".join(list(map(str, [ac, an, af]))))

                smart_write(fh, "\n")

            fh.close()


def process_gene_data(gene_data, samples, lof_table, epacts_groups, group_file):
    # process gene data
    for gene in gene_data.keys():
        this_gene_data = gene_data[gene]
        lof_table[gene] = get_sample_lof_counts(this_gene_data, samples)

        if group_file:
            epacts_groups[gene] = create_epacts_group(this_gene_data)

    gene_data = {}


def get_sample_lof_counts(gene_data, samples):
    max_counts = {}
    for s in samples:
        max_counts[s] = 0

    for var_data in gene_data:
        var = var_data["var"]

        alleles = [var.REF] + var.ALT
        required_alt = alleles[var_data["allele_i"]]

        for sample, bases in zip(samples, var.gt_bases):
            count = 0

            for a in re.split(r"\||\/|\\", bases):
                if a == required_alt:
                    count += 1

            if count > max_counts[sample]:
                max_counts[sample] = count

    return max_counts


def create_epacts_group(gene_data):
    group = []

    for var_data in gene_data:
        var = var_data["var"]
        alleles = [var.REF] + var.ALT

        # [CHROM]:[POS]_[REF]/[ALT]
        group.append(
            "{}:{}_{}/{}".format(
                var.CHROM, var.POS, var.REF, alleles[var_data["allele_i"]]
            )
        )

    return group


def get_var_gene_data(vep_vcf, var, filter_eval, filter_fields, symbol_only=False):
    var_gene_data = {}

    all_allele_data = vep_vcf.get_all_allele_data(var, filter_fields)

    for allele_i in all_allele_data:

        # track if we've added this variant allele for this gene already
        gene_varallele_added = {}

        for allele_data in all_allele_data[allele_i]:
            gene = get_gene_key(allele_data, symbol_only)
            if not gene or gene in gene_varallele_added:
                continue

            if eval_filter(filter_eval, allele_data):
                if gene not in var_gene_data:
                    var_gene_data[gene] = []

                var_gene_data[gene].append({"var": var, "allele_i": allele_i})

                gene_varallele_added[gene] = True

    return var_gene_data


def get_gene_key(allele_data, symbol_only=False):
    if symbol_only:
        return allele_data.get("SYMBOL", None)

    gene_key = allele_data.get("Gene", None)

    if "SYMBOL" in allele_data:
        gene_key = allele_data["SYMBOL"] + "(" + gene_key + ")"

    return gene_key


if __name__ == "__main__":
    main()
