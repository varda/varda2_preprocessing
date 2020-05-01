from setuptools import setup

setup(name='gvcf2coverage',
    version='0.1',
    description='A python tool to extract coverage from gvcf files.',
    url='http://github.com/varda/varda2_preprocessing',
    author='Mark Santcroos',
    author_email='m.a.santcroos@lumc.nl',
    license='MIT',
    packages=['gvcf2coverage'],
    zip_safe=False,
    entry_points = {
        'console_scripts': ['gvcf2coverage=gvcf2coverage'],
    },
    install_requires=[
        'cyvcf2'
    ],
)
