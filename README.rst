.. image:: https://readthedocs.org/projects/aspring/badge/?version=latest
    :alt: ReadTheDocs
    :target: https://aspring.readthedocs.io/en/stable/
.. image:: https://img.shields.io/coveralls/github/PhyloSofS-Team/aspring/main.svg
    :alt: Coveralls
    :target: https://coveralls.io/r/PhyloSofS-Team/aspring
.. image:: https://img.shields.io/pypi/v/aspring.svg
    :alt: PyPI-Server
    :target: https://pypi.org/project/aspring/
.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/

|

=======
aspring
=======


    Alternatively Spliced Pseudo Repeat IN-Gene


**ASPRING** is a computational tool for detecting Alternative Splicing Repetitive 
Units (ASRUs) on a gene. It analyzes the outputs of ThorAxe, which provides 
Multiple Sequence Alignments of exonic regions and generates information about 
their use on alternative isoforms. You can run aspring through the 
command line to find duplication events of such exonic regions. 
This tool will provide information on the duplicated regions through a 
couple of tables.

Requirements
============

**ASPRING** is a Python package that can be installed from PyPI using the pip 
package manager. Its pipeline uses **R** and the **HH-suite3**. In particular, 
`Rscript` from **R** and `hhmake` and `hhalign` from **HH-suite3** must be in 
the `PATH`. 

You can see the installation instructions for both in the following links:

- R: https://www.r-project.org/
- HH-suite3: https://github.com/soedinglab/hh-suite

If you have miniconda_ installed, you can use it to install both **HH-suite3** 
and **R**. For example:

.. code-block:: shell

    conda install -c conda-forge -c bioconda hhsuite
    conda install -c conda-forge r-base=4.2.2

Also, you will need to know the path where the **HH-suite3** scripts are 
installed, as **ASPRING** needs to access `reformat.pl`. If you 
installed **HH-suite3** using the base environment miniconda_ on a Linux, 
then the path will be something like `~/miniconda3/scripts`.

The **ASPRING** package ships with a renv_ environment 
(located at `src/aspring/R_scripts`) that will be automatically executed, 
so you do not need to worry about installing the **R** dependencies.

Outputs
=======

For a given `gene` (Ensembl Gene ID), ASPRING returns:


- `{gene}_ASRUs_table.csv`
- `{gene}_instances_table.csv`
- `{gene}_duplication_pairs.csv`
- `{gene}_eventsDup_withCols.txt`
- `DupRaw/{gene}` folder containing the `s-exon_A.s-exon_B.hhr` files (HMM-HMM alignments) 



.. _pyscaffold-notes:

Note
====

This project has been set up using PyScaffold 4.4. For details and usage
information on PyScaffold see https://pyscaffold.org/.


.. _miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _renv: https://rstudio.github.io/renv/articles/renv.html