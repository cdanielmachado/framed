from setuptools import setup, find_packages

setup(
    name = 'framed',
    version = '0.2.1',
    package_dir = {'':'src'},
    packages = find_packages('src'),
    install_requires = [], # e.g. ["numpy>=1.3.0"]
    author = 'Daniel Machado',
    author_email = 'cdanielmachado@gmail.com',
    description = 'FRAMED - FRAmework for Metabolic Engineering and Design',
    license = 'Apache License Version 2.0',
    keywords = 'biology metabolism bioinformatics',
    url = 'https://github.com/cdanielmachado/framed',
    long_description = open('README.rst').read(),
    classifiers = [
        'Development Status :: 4 - Beta',
        'Topic :: Utilities',
        'Programming Language :: Python :: 2.5',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: Apache Software License',
    ],
)
