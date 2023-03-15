import pytest
import sys
import tempfile
import os
import csv
import binascii

import vepvcf.lof_rollup as lof_rollup

class TestHelpers(object):
	def test_get_gene_key(self):
		assert lof_rollup.get_gene_key({'foo': 'bar'}) == None
		assert lof_rollup.get_gene_key({'Gene': 'foo'}) == 'foo'
		assert lof_rollup.get_gene_key({'Gene': 'foo', 'SYMBOL': 'bar'}) == 'bar(foo)'
		assert lof_rollup.get_gene_key({'Gene': 'foo', 'SYMBOL': 'bar'}, symbol_only=True) == 'bar'

	def test_create_epacts_group(self, ma_var):
		assert lof_rollup.create_epacts_group([{'var': ma_var, 'allele_i': 1}]) == ['1:97078987_G/A']
		assert lof_rollup.create_epacts_group([{'var': ma_var, 'allele_i': 2}]) == ['1:97078987_G/T']


class TestGetVarGeneData(object):
	def test_pass_filter(self, test_vep_vcf, ma_var):
		assert lof_rollup.get_var_gene_data(test_vep_vcf, ma_var, 'True', []) == {'DPYD(ENSG00000188641)': [{'allele_i': 1, 'var': ma_var}]}

	def test_fail_filter(self, test_vep_vcf, ma_var):
		assert lof_rollup.get_var_gene_data(test_vep_vcf, ma_var, 'False', []) == {}

	def test_symbol_only(self, test_vep_vcf, ma_var):
		assert lof_rollup.get_var_gene_data(test_vep_vcf, ma_var, 'True', [], symbol_only=True) == {'DPYD': [{'allele_i': 1, 'var': ma_var}]}


class TestGetSampleLofCounts(object):
	def test_basic(self, test_vep_vcf):
		var = next(test_vep_vcf.vcf('1:97078076-97078076'))
		gd = lof_rollup.get_var_gene_data(test_vep_vcf, var, 'True', [])

		counts = lof_rollup.get_sample_lof_counts(list(gd.values())[0], test_vep_vcf.vcf.samples)
		assert len(list(filter(lambda x: counts[x] != 0, counts.keys()))) == 15
		assert max(counts.values()) == 1

	def test_homs(self, test_vep_vcf):
		var = next(test_vep_vcf.vcf('1:97078196-97078196'))
		gd = lof_rollup.get_var_gene_data(test_vep_vcf, var, 'True', [])

		counts = lof_rollup.get_sample_lof_counts(list(gd.values())[0], test_vep_vcf.vcf.samples)
		assert len(list(filter(lambda x: counts[x] == 2, counts.keys()))) == 33
		assert max(counts.values()) == 2


class TestMain(object):
	def test_lof(self, vcf_file):
		[tmp, fn] = tempfile.mkstemp()

		sys.argv = ['test', '-f', '#AF <= 0.05 and @lof', '-o', fn, vcf_file]
		lof_rollup.main()

		with open(fn, 'r') as fh:
			reader = csv.reader(fh, delimiter='\t')

			assert next(reader) == ['#', 'DPYD(ENSG00000188641)']

			data = {}
			for row in reader:
				data[row[0]] = row[1]

			assert sorted(list(filter(lambda x: data[x] != '0', data))) == ['HG00185', 'HG00349', 'HG00357']

	def test_groups(self, vcf_file):
		[tmp, fn] = tempfile.mkstemp()

		sys.argv = ['test', '-f', '#AF <= 0.05 and @lof', '-o', fn, '-g', fn + '_groups', vcf_file]
		lof_rollup.main()

		with open(fn+'_groups', 'r') as fh:
			reader = csv.reader(fh, delimiter='\t')
			assert next(reader) == [
				'DPYD(ENSG00000188641)',
				'1:97382464_T/C',
				'1:97450058_C/T',
				'1:97828127_G/A',
				'1:97883353_G/A'
			]

	def test_freqs(self, vcf_file):
		[tmp, fn] = tempfile.mkstemp()

		sys.argv = ['test', '-f', '#AF <= 0.05 and @lof', '-o', fn, '-q', fn + '_freqs', vcf_file]
		lof_rollup.main()

		with open(fn+'_freqs', 'r') as fh:
			reader = csv.reader(fh, delimiter='\t')
			assert next(reader) == ['#Gene', 'AC_ALL', 'AN_ALL', 'AF_ALL']
			assert next(reader) == [
				'DPYD(ENSG00000188641)',
				'3', '382', str(3.0 / 382)
			]

	def test_pops(self, vcf_file, pops_file):
		[tmp, fn] = tempfile.mkstemp()

		sys.argv = ['test', '-f', '#AF <= 0.05 and @lof', '-o', fn, '-q', fn + '_freqs', '-p', pops_file, vcf_file]
		lof_rollup.main()

		with open(fn+'_freqs', 'r') as fh:
			reader = csv.reader(fh, delimiter='\t')
			assert next(reader) == [
				'#Gene',
				'AC_ALL', 'AN_ALL', 'AF_ALL',
				'AC_EAS', 'AN_EAS', 'AF_EAS',
				'AC_EUR', 'AN_EUR', 'AF_EUR'
			]
			assert next(reader) == [
				'DPYD(ENSG00000188641)',
				'3', '382', str(3.0 / 382),
				'0', '12', str(0.0 / 12),
				'3', '370', str(3.0 / 370)
			]

	def test_lof_gzip_out(self, vcf_file):
		[tmp, fn] = tempfile.mkstemp(suffix=".txt.gz")

		sys.argv = ['test', '-f', '#AF <= 0.05 and @lof', '-o', fn, vcf_file]
		lof_rollup.main()

		assert os.path.isfile(fn)
		with open(fn, 'rb') as test_f:
			assert binascii.hexlify(test_f.read(2)) == b'1f8b'

	def test_symbol_only(self, vcf_file):
		[tmp, fn] = tempfile.mkstemp()

		sys.argv = ['test', '--symbol_only', '-f', '#AF <= 0.05 and @lof', '-o', fn, vcf_file]
		lof_rollup.main()

		with open(fn, 'r') as fh:
			reader = csv.reader(fh, delimiter='\t')
			assert next(reader) == ['#', 'DPYD']
			