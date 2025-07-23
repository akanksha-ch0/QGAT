from setuptools import setup, find_packages

setup(
    name="QGAT",
    version="0.1",
    author="Akanksha Choudhary",
    description="QGAT: QTL and Gene Annotation Tool",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "QGAT": [
            "QTLdb/*.bed.gz",      # Include all .bed QTLdb files
            "QTLdb/*.gff",      # Include all .gff QTLdb files
            "QTLdb/*.gff.gz",   # ✅ NEW: Include compressed .gff.gz files too
            "input/*",
            "output/*",
        ],
    },
    install_requires=[
        "pandas",
        "argparse",
        "matplotlib",  # Required for plotting
    ],
    entry_points={
        "console_scripts": [
            "QGAT=QGAT.main:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
