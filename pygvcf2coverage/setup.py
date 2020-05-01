from setuptools import setup

setup(name='pygvcf2coverage',
    version='0.2',
    description='A python tool to extract coverage from gvcf files.',
    url='http://github.com/varda/varda2_preprocessing',
    author='Mark Santcroos',
    author_email='m.a.santcroos@lumc.nl',
    license='MIT',
    packages=['pygvcf2coverage'],
    zip_safe=False,
    entry_points = {
        'console_scripts': ['pygvcf2coverage=pygvcf2coverage:main'],
    },
    install_requires=[
        'cyvcf2'
    ],
)
