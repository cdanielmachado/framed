from setuptools import setup, find_packages

setup(
    name = 'framed',
    version = '0.3.5',
    package_dir = {'':'src'},
    packages = find_packages('src'),
    install_requires = ["numpy",
                        "scipy",
                        "sympy",
                        "python-libsbml",
                        "matplotlib",
                        "pandas"],
    scripts = ["scripts/framed-smetana-pipeline"],
    author = 'Daniel Machado',
    author_email = 'cdanielmachado@gmail.com',
    description = 'framed - metabolic modeling for python',
    license = 'Apache License Version 2.0',
    keywords = 'metabolism modeling',
    url = 'https://github.com/cdanielmachado/framed',
    long_description = open('README.rst').read(),
    classifiers = [
        'Development Status :: 4 - Beta',
        'Topic :: Utilities',
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: Apache Software License',
    ],
)
