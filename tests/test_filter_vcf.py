import pytest
import sys
import tempfile
import os
import binascii

import vepvcf.filter_vcf as filter_vcf

class TestFilterVcf(object):
	def test_basic_filter(self, capsys, vcf_file, test_region):
		sys.argv = ['test', '-r', test_region, '-c', '-f', '#AF > 0.5', vcf_file]
		filter_vcf.main()

		assert capsys.readouterr().out.rstrip() == "0"

		sys.argv = ['test', '-r', test_region, '-c', '-f', '#AF < 0.5', vcf_file]
		filter_vcf.main()
		
		assert capsys.readouterr().out.rstrip() == "1"

	def test_count(self, capsys, vcf_file):
		sys.argv = ['test', '-c', vcf_file]
		filter_vcf.main()

		assert capsys.readouterr().out.rstrip() == "25354"

	def test_filter_file(self, capsys, vcf_file, filter_file):
		sys.argv = ['test', '-c', vcf_file, '-f', '#Consequence in ~' + filter_file]
		filter_vcf.main()
		assert capsys.readouterr().out.rstrip() == "23279"

		sys.argv = ['test', '-c', vcf_file, '-f', '#Consequence not in ~' + filter_file]
		filter_vcf.main()
		assert capsys.readouterr().out.rstrip() == "2075"

	def test_region(self, capsys, vcf_file, test_region):
		sys.argv = ['test', '-c', '-r', test_region, vcf_file]
		filter_vcf.main()

		assert capsys.readouterr().out.rstrip() == "1"

	def test_info_field(self, capsys, info_field_vcf_file):
		sys.argv = ['test', '-i', 'FOO', '-f', '#Allele == "C"', '-c', info_field_vcf_file]
		filter_vcf.main()

		assert capsys.readouterr().out.rstrip() == "2"

	def test_no_header(self, capsys, info_field_vcf_file):
		sys.argv = ['test', '-i', 'FOO', '-n', info_field_vcf_file]
		filter_vcf.main()

		assert len(capsys.readouterr().out.split('\n')) == 6

	def test_output_file(self, info_field_vcf_file):
		[tmp, fn] = tempfile.mkstemp(suffix=".vcf")
		sys.argv = ['test', '-i', 'FOO', '-o', fn, info_field_vcf_file]
		filter_vcf.main()

		assert os.path.isfile(fn)
		assert os.path.getsize(fn) == 13401

	def test_output_file_gzip(self, info_field_vcf_file):
		[tmp, fn] = tempfile.mkstemp(suffix=".vcf.gz")
		sys.argv = ['test', '-i', 'FOO', '-o', fn, info_field_vcf_file]
		filter_vcf.main()

		assert os.path.isfile(fn)
		with open(fn, 'rb') as test_f:
			assert binascii.hexlify(test_f.read(2)) == b'1f8b'



