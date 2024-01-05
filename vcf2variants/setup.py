from setuptools import setup

setup(name='vcf2variants',
    version='1.2',
    description='A python tool to convert vcf to varda variant files.',
    url='http://github.com/varda/varda2_preprocessing',
    author='Mark Santcroos',
    author_email='m.a.santcroos@lumc.nl',
    license='MIT',
    packages=['vcf2variants'],
    zip_safe=False,
    entry_points = {
        'console_scripts': ['vcf2variants=vcf2variants:main'],
    },
    install_requires=[
        'natsort',
        'vcfphasesets',
    ],
    extras_require={
        'testing': [
            'pytest',
        ],
    },
)
