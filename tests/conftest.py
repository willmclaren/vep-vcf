import pytest
import tests.config as config
from vepvcf.vep_vcf import VEP_VCF, parse_filter_string

@pytest.fixture(scope="module")
def vcf_file():
	return config.TEST_DIR + '/data/test.vcf.gz'

@pytest.fixture(scope="module")
def test_vep_vcf(vcf_file):
	return VEP_VCF(vcf_file)

@pytest.fixture(scope="module")
def info_field_vcf_file():
	return config.TEST_DIR + '/data/test_info_field.vcf.gz'

@pytest.fixture(scope="module")
def test_region():
	return '1:97078987-97078987'

@pytest.fixture(scope="module")
def ma_var(test_vep_vcf, test_region):
	return next(test_vep_vcf.vcf(test_region))

@pytest.fixture(scope="module")
def pops_file():
	return config.TEST_DIR + '/data/continental_pops.txt'

@pytest.fixture(scope="module")
def filter_file():
	return config.TEST_DIR + '/data/filter_file.txt'