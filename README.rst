pairLiftOver
============
pairLiftOver is a Python package that converts the two-dimensional genomic coordinates
of chromatin contact pairs between assemblies. It supports 4 input data formats:
`4DN pairs <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_,
`allValidPairs <https://nservant.github.io/HiC-Pro/RESULTS.html>`_, `cool <https://open2c.github.io/cooler/>`_,
and `hic <https://github.com/aidenlab/juicer/wiki/Data>`_. The default output of pairLiftOver is a
sorted pairs file in the standard 4DN pairs format, containing seven columns: “readID”,
“chr1”, “pos1”, “chr2”, “pos2”, “strand1”, and “strand2”. However, you can also choose to
output a matrix file in cool or hic format by setting the parameter “--output-format”.

Installation
============
pairLiftOver is developed and tested on UNIX-like operating system, and following packages
are required:

- python 3.6+
- cooler 0.8.6
- pairtools 0.3.0
- pyliftover
- hic-straw

We recommend using `conda <https://conda.io/miniconda.html>`_ to manage these packages. After
you have installed conda on your system, execute the commands below::

    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge
    $ conda create -n pairliftover python=3.7 pairtools=0.3.0 cooler=0.8.6 pyliftover
    $ conda activate pairliftover
    $ pip install pairLiftOver hic-straw

Usage
=====
Open a terminal, type ``pairLiftOver -h`` for help information.

Citation
========
Wang, X., Luan, Y., Yue, F. pairLiftOver: a Python package to convert two-dimensional genomic coordinates between assemblies. bioRxiv.

