from setuptools import setup
import sys

if sys.version_info[:2] != (2, 7):
    sys.exit('Sorry, only Python 2.7 is supported')

setup(
    name="TETyper",
    version="1.0.0",
    url="https://github.com/aesheppard/TETyper",
    description="transposable element typing pipeline",
    license="GPL'",
    keywords="TE typing NGS",
    classifiers=[
        'Development Status :: 5 - Alpha',
        'License :: OSI Approved :: GPL',
        'Programming Language :: Python :: 2.7'
        ],
    install_requires=['biopython', 'pysam', 'pyvcf'],
    scripts=['TETyper.py'])
