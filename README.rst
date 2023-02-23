.. image:: https://readthedocs.org/projects/aspring/badge/?version=latest
    :alt: ReadTheDocs
    :target: https://aspring.readthedocs.io/en/stable/
.. image:: https://img.shields.io/coveralls/github/PhyloSofS-Team/aspring/main.svg
    :alt: Coveralls
    :target: https://coveralls.io/r/PhyloSofS-Team/aspring
.. image:: https://img.shields.io/pypi/v/aspring.svg
    :alt: PyPI-Server
    :target: https://pypi.org/project/aspring/
.. image:: https://img.shields.io/docker/v/diegozea/aspring?label=docker
    :alt: Docker
    :target: https://hub.docker.com/r/diegozea/aspring
.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/


=================
🌼 **ASPRING** 🌼
=================


*Alternatively Spliced Pseudo Repeat IN-Gene*


**ASPRING** is a computational tool for detecting Alternative Splicing Repetitive Units
(ASRUs) on a gene. It analyzes the outputs of ThorAxe, which provides Multiple Sequence
Alignments of exonic regions and generates information about their use on alternative
isoforms. You can run aspring through the command line to find duplication events of such
exonic regions. This tool will provide information on the duplicated regions through a
couple of tables.

.. note:: 

    **ASPRING** requires **ThorAxe outputs** for a single query gene to run. If you don't
    have ThorAxe outputs, you can visit the `ThorAxe documentation`_ to learn how to install
    and run ThorAxe on your data. ThorAxe can also be run using the `Ases web server`_,
    which provides a user-friendly interface for running ThorAxe online. Note that you can
    skip the PhyloSofS step of Ases to obtain results more quickly for use with ``aspring``.
    Once you have ThorAxe outputs, you can use ``aspring`` to identify ASRUs for your query
    gene.



Requirements
============

**ASPRING** is a Python package that can be installed from PyPI using the pip package
manager. Its pipeline uses **R** and the **HH-suite3**. In particular, ``Rscript`` from **R**
and ``hhmake`` and ``hhalign`` from **HH-suite3** must be in the ``PATH``. 

You can see the installation instructions for both in the following links:

- R: https://www.r-project.org/
- HH-suite3: https://github.com/soedinglab/hh-suite

If you have miniconda_ installed, you can use it to install both **HH-suite3** and **R**.
For example:

.. code-block:: bash

    conda install -c conda-forge -c bioconda hhsuite conda install -c conda-forge
    r-base=4.2.2

Also, you will need to know the path where the **HH-suite3** scripts are installed, as
**ASPRING** needs to access ``reformat.pl``. If you installed **HH-suite3** using the base
environment miniconda_ on a Linux, then the path will be something like
`~/miniconda3/scripts`.

The **ASPRING** package ships with a renv_ environment (located at ``src/aspring/R_scripts``)
that will be automatically executed, so you do not need to worry about installing the **R**
dependencies.

Nomenclature
============

.. image:: https://raw.githubusercontent.com/PhyloSofS-Team/aspring/main/docs/_static/nomenclature_example.png
  :alt: ASPRING nomenclature explained using a ThorAxe Evolutionary Splicing Graph (ESG). The image shows an ASRU composed by two ASPRs, one of them composed by multiple s-exons.

The figure shows an example of an **Alternative Splicing Repetitive Unit** (**ASRU**)
composed by two **Alternatively Spliced Pseudo Repeats** (**ASPRs**), one of them composed
by multiple **s-exons**.

The nodes are the s-exons. The opaque red boxes are the ASPRs, and the transparent red box
is the ASRU.

The **ASPRs** are repetitive units identified by ASPRING that consist of one or more s-exons
alternatively spliced in different isoforms. Note that ASPRs are called **instances** on the
output tables.

How to use aspring
==================

``aspring`` is a Python-based command-line tool that helps identify Alternative Splicing
Repetitive Units (ASRUs) from Thoraxe outputs for a single query gene. The tool executes
several steps that involve converting data, creating HMM profiles, aligning profiles,
parsing and filtering alignments, and generating ASRUs and Alternatively Spliced Pseudo
Repeats (ASPRs) tables for the query gene.

Here is how to run the script after installing the package:

1. Open your terminal.
2. Run the script with the following command:

   ::

       aspring --gene GENE_NAME --path_data PATH_TO_THORAXE_OUTPUTS --path_hhsuite_scripts PATH_TO_HHSUITE_SCRIPTS

   To query a specific gene, replace ``GENE_NAME`` with the corresponding Ensemble
   ID. If ThorAxe outputs are in the current working directory, ``--path_data``
   parameter can be avoided. In case the ``reformat.pl`` script from the HH-suite3
   is in the path indicated by the ``HHSUITE_SCRIPTS`` environment variable, 
   ``--path_hhsuite_scripts`` parameter can be omitted.
   Otherwise, replace ``PATH_TO_HHSUITE_SCRIPTS`` with the
   path to the ``reformat.pl`` script directory. Replace ``PATH_TO_THORAXE_OUTPUTS``
   with the path to the ThorAxe outputs directory for the gene of  interest.

   Optional arguments are available to customize the behavior of the script. Run the command
   ``aspring --help`` to see the full list of options.

The script will execute several steps and generate output files containing ASRU and
ASPR tables for the query gene.


Docker
------

To ease the use and installation of ASPRING, we have created a Docker image that
you can easily download and run. The aspring Docker image is available on
`Docker Hub`_. To run the ``aspring`` tool using the Docker image, you must have 
Docker installed on your system. You can download and install **Docker** from the 
`official website`_. Once Docker is installed, you can run ``aspring`` using the 
following command:

.. code-block:: bash

  sudo docker run --mount type=bind,source=$(pwd),target=/data diegozea/aspring aspring --gene GENE_NAME

In this command, we use the ``docker run`` command to run ``aspring``. We are
mounting the current working directory using the ``--mount`` option, which is
necessary for providing access to the data files required by ``aspring``. The
``--mount`` option takes two parameters: ``type`` and ``source``. ``type`` specifies 
the type of mount to use. In this case, we use a ``bind`` mount, which allows us 
to mount a directory from the host system to the container; that is a
**requirement** to enable ``aspring`` to see the input files and to let it save
the output files in your filesystem. ``source`` specifies the source directory to
mount. In this case, we use ``$(pwd)`` to select the current working directory as
the source. We are also specifying the ``target`` directory as ``/data`` in the
container. This means that the files from the current working directory on the
host system will be available in the ``/data`` directory in the container.

The aspring tool requires **R** and the **HH-suite3**, which are already
installed in the Docker image. Therefore, there is no need to specify
``--path_hhsuite_scripts`` or ``--path_data``; the last one is set to ``/data`` by
default.

Finally, we specify the ``--gene`` option with ``GENE_NAME`` to run aspring on that gene.


Pipeline
--------

ASPRING is a tool for detecting Alternative Splicing Repetitive Units (ASRUs) on a gene. The
pipeline consists of nine steps, each of which can be executed separately, but it is
recommended to run the main script ``aspring`` to execute the entire pipeline. Only steps 1,
2, and 3 require **HH-suite3** and step 6 requires **R**. You can use the ``-h`` argument to
show the arguments for each step.

The pipeline steps are:

1. ``step_01_preprocess``: Reformat s-exons fasta files to a2m.
2. ``step_02_hmm_maker``: Generates a Hidden Markov Model (HMM) profile for each s-exon.
3. ``step_03_hmm_aligner``: HMM-HMM alignment of all the s-exons combinations.
4. ``step_04_gettable``: Parses the alignment files and creates a table.
5. ``step_05_filter``: Filter the table to keep gene duplication pairs based on identity,
   coverage, p-value and number of species in the MSAs.
6. ``step_06_stats``: Generates statistics on the filtered duplicated regions.
7. ``step_07_reformat``: Reformat the previous outputs to add the information about the
   duplicated regions.
8. ``step_08_ASRUs``: Identifies the Alternative Splicing Repetitive Units (ASRUs) on the
   gene.
9. ``step_09_clean``: Removes the intermediate files generated during the pipeline.

Note that the main script ``aspring`` **runs the entire pipeline** automatically. However,
the user can also execute the scripts of each pipeline step individually for more control
over the pipeline.


Outputs
=======

For a given ``gene`` (Ensembl Gene ID), ASPRING returns:

- ``{gene}_ASRUs_table.csv``
- ``{gene}_instances_table.csv``
- ``{gene}_duplication_pairs.csv``
- ``{gene}_eventsDup_withCols.txt``
- ``DupRaw/{gene}`` folder containing the ``s-exon_A.s-exon_B.hhr`` files (HMM-HMM alignments) 

{gene}_ASRUs_table.csv
----------------------

This table provides information on the Alternatively Spliced Repeat Units (ASRUs) detected
for the given ``gene``. Each row corresponds to a distinct ASRU and provides the following
information:

- ``gene``: The Ensembl Gene ID for the given gene.
- ``ASRU``: The set of duplicated s-exons, a.k.a Alternatively Spliced Pseudo Repeats (ASPRs)
  that belong to the ASRU.
- ``Nbinstances``: The number of Alternatively Spliced Pseudo Repeats of the ASRU that were
  found in the exonic regions of the gene.
- ``max``: The length of the longest ASPR instance of the ASRU, in residues.
- ``min``: The length of the shortest ASPR instance of the ASRU, in residues.
- ``moy``: The mean length of the instances of the ASRU, in amino acid residues.
- ``median``: The median length of the instances of the ASRU, in residues.
- ``std``: The standard deviation of the lengths of the instances of the ASRU, in amino acid
  residues.
- ``eventsRank``: The rank/position of the alternative splicing events involving the ASRU in
  the ``ases.csv`` output table from ThorAxe — from the most to the least conserved/frequent.

{gene}_instances_table.csv
--------------------------

This table provides information on the instances of ASRUs detected for the given ``gene``.
Each row corresponds to a distinct instance and provides the following information:

- ``instance``: The sequence of the ASPR instance, in the form of a string of amino acid
  residues.
- ``size``: The length of the ASPR instance, in amino acid residues.
- ``NbSex``: The number of exonic regions where the ASPR instance was detected.
- ``ASRU``: The set of homologous/duplicated s-exons that belong to the ASRU to which the ASPR
  instance belongs.
- ``gene``: The Ensembl Gene ID for the given gene.
 
{gene}_duplication_pairs.csv
----------------------------

This table provides information on the pairs of exonic regions that were involved in the
duplication events. Each row corresponds to a distinct pair of s-exons and provides the
following information:

- ``S_exon_Q``: The identifier of the first s-exon.
- ``S_exon_T``: The identifier of the second s-exon.
- ``Gene``: The Ensembl Gene ID for the given gene.
- ``Prob``: The probability score of the alignment of the exonic region pair.
- ``E-value``: The E-value associated with the alignment of the exonic region pair.
- ``P-value``: The P-value associated with the alignment of the exonic region pair.
- ``Score``: The alignment score of the alignment of the exonic region pair.
- ``Cols_Q``: The alignment columns corresponding to the first s-exon, in the format
  "start-end".
- ``Cols_T``: The alignment columns corresponding to the second s-exon, in the format
  "start-end".
- ``Length_Q``: The length of the first s-exon, in amino acid residues.
- ``Length_T``: The length of the second s-exon, in amino acid residues.
- ``Identities``: The percentage of identical residues in the alignment of the exonic region
  pair.
- ``IdCons``: The percentage of conserved residues in the alignment of the exonic region pair.
- ``Similarity``: The fraction of similar residues in the alignment of the exonic region pair.
- ``NoSpecies_Q``: The number of species in which the first s-exon is conserved.
- ``NoSpecies_T``: The number of species in which the second s-exon is conserved.

{gene}_eventsDup_withCols.txt
-----------------------------

This table provides detailed information on the alternative splicing events in with the
ASRUs are involved. Each row corresponds to a distinct event and provides the following
information:

- ``gene``: The Ensembl Gene ID for the given gene.
- ``sexA``: The index of the first s-exon in the ASRU.
- ``sexB``: The index of the second s-exon in the ASRU.
- ``rank``: The rank of the alternative splicing event, as ordered in the ThorAxe output table
  from the most to the least conserved/frequent.
- ``type``: The type of the alternative splicing events, e.g "alternative".
- ``statusA``: The status of the path with the first s-exon, which can be alternative or
  canonical.
- ``statusB``: The status of the path with the first s-exon, which can be alternative or
  canonical.
- ``lePathA``: Number of s-exons in the path with the first s-exon.
- ``lePathB``: Number of s-exons in the path with the second s-exon.
- ``exclu``: A boolean indicating whether the event involves mutually exclusive s-exons.
- ``pval``: The P-value associated with the alignment of the exonic region pair.
- ``ncols``: The number of columns in the alignment.
- ``leA``: The length of the first s-exon, in amino acid residues.
- ``leB``: The length of the second s-exon, in amino acid residues.
- ``typePair``: The type of the alternative splicing event.
- ``ColA``: The alignment columns corresponding to the first s-exon, in the format
  "start-end".
- ``ColB``: The alignment columns corresponding to the second s-exon, in the format
  "start-end".


.. _pyscaffold-notes:

Note
====

This project has been set up using PyScaffold 4.4. For details and usage information on
PyScaffold see https://pyscaffold.org/.


.. _miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _renv: https://rstudio.github.io/renv/articles/renv.html
.. _ThorAxe documentation: https://phylosofs-team.github.io/thoraxe/
.. _Ases web server: http://www.lcqb.upmc.fr/Ases
.. _Docker Hub: https://hub.docker.com/r/diegozea/aspring
.. _official website: https://www.docker.com/
