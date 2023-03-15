#!/usr/bin/env python3

import argparse
from vepvcf.vep_vcf import VEP_VCF, parse_filter_string, eval_filter
from vepvcf.utils import smart_open, smart_write


def main():
	parser = argparse.ArgumentParser(description="Filter VCF file, optionally with VEP-annotations")
	parser.add_argument('vcf', metavar='VCF', help="VCF file")
	parser.add_argument('-i', '--info_field', default="CSQ", help="VCF INFO field name containing VEP annotation (default: %(default)s)")
	parser.add_argument('-f', '--filter', default='1', help="Filter string eval'd to select variants, use #[field_name] to denote VEP or INFO fields (default: \"%(default)s\")")
	parser.add_argument('-r', '--region', help="Limit to region [chr:start-end] of VCF")
	parser.add_argument('-o', '--output_file', help="Output file name")
	parser.add_argument('-c', '--count_only', action="store_true", help="Output a count of matched entries only")
	parser.add_argument('-n', '--no_header', action="store_true", help="Don't write the VCF header")
	args = parser.parse_args()

	(filter_eval, filter_fields) = parse_filter_string(args.filter)

	try:
		vep_vcf = VEP_VCF(vcf_file=args.vcf, vep_info_field=args.info_field)
	except:
		raise

	if args.region:
		vcf_iter = vep_vcf.vcf(args.region)
	else:
		vcf_iter = vep_vcf.vcf

	with smart_open(args.output_file) as out_fh:
		if not(args.count_only or args.no_header):
			smart_write(out_fh, vep_vcf.vcf.raw_header)

		passed_count = 0

		for var in vcf_iter:
			allele_data = vep_vcf.get_all_allele_data(
				var=var,
				filter_fields=filter_fields
			)

			passed = []

			for allele_i in sorted(allele_data.keys()):
				passed += list(filter(lambda x: eval_filter(filter_eval, x), allele_data[allele_i]))

			if len(passed) > 0:
				passed_count += 1

				if not args.count_only:
					smart_write(out_fh, str(var))

	if args.count_only:
		smart_write(out_fh, str(passed_count) + "\n")



if __name__ == '__main__':
	main()
	