from pathlib    import Path
import setuptools 
import sys
# modified from nextstrain/augur repo
min_version = (3, 6)

if sys.version_info < min_version:
    error = """
Python {0} or above is required.

Make sure you have an up-to-date pip installed.  
""".format('.'.join(str(n) for n in min_version)), sys.exit(error)

base_dir = Path(__file__).parent.resolve()
version_file = base_dir / "perobarosa/__version__.py"
readme_file = base_dir / "README.md"

# Eval the version file to get __version__; avoids importing our own package
with version_file.open() as f:
    exec(f.read())

with readme_file.open(encoding = "utf-8") as f:
    long_description = f.read()

setuptools.setup(
    name = "peroba",
    version = __version__,
    author = "Leonardo de Oliveira Martins",
    author_email = "Leonardo.de-Oliveira-Martins@quadram.ac.uk",
    description = "Phylogenetic analysis pipeline for viral epigenomics at the Quadram Institute Biosciences",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    keywords = "phylogenetics, COVID19",
    url = "https://github.com/quadram-institute-bioscience/peroba",
    project_urls = {
        "Source": "https://github.com/quadram-institute-bioscience/peroba",
    },
    packages = setuptools.find_packages(),
    include_package_data=True,
    package_data = {'peroba': ['data/*','data/report/*']},
    data_files = [("", ["LICENSE"])],
    python_requires = '>={}'.format('.'.join(str(n) for n in min_version)),
    license='GPLv3+',
    install_requires=[
           'biopython >= 1.70',
           'ete3',
           #'pastml',
           'numpy',
#           'matplotlib',
           'pandas',
           #'seaborn',
           #'basemap',
           #'geopandas',
#           'treeswift',
           'xxhash',
           #'pandas_profiling',
           'scikit-learn'
       ],
    classifiers = [
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        # Python 3 only
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6"
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    # Install a "peroba" program which calls peroba.__main__.main()
    #   https://setuptools.readthedocs.io/en/latest/setuptools.html#automatic-script-creation
    entry_points = {
        "console_scripts": [ "peroba = perobarosa.peroba:main" ]
    }
)
