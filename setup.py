from setuptools import find_packages, setup

setup(
    name="vepvcf",
    version="0.1.0",
    packages=find_packages(),
    include_package_data=False,
    entry_points={
        "console_scripts": [
            "filter_vcf = vepvcf.filter_vcf:main",
            "lof_rollup = vepvcf.lof_rollup:main",
        ]
    },
)
