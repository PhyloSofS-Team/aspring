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

.. code-block:: shell

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

.. _nomenclature-example:

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

.. code-block:: shell
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

.. image::
    :target: nomenclature-example
    :alt: ASPRING nomenclature explained using a ThorAxe Evolutionary Splicing Graph (ESG). The image shows an ASRU composed by two ASPRs, one of them composed by multiple s-exons.

    iVBORw0KGgoAAAANSUhEUgAAATIAAAKECAMAAACzY0YCAAAABGdBTUEAALGPC/xhBQAAAAFzUkdC
    AK7OHOkAAAHWaVRYdFhNTDpjb20uYWRvYmUueG1wAAAAAAA8eDp4bXBtZXRhIHhtbG5zOng9ImFk
    b2JlOm5zOm1ldGEvIiB4OnhtcHRrPSJYTVAgQ29yZSA2LjAuMCI+CiAgIDxyZGY6UkRGIHhtbG5z
    OnJkZj0iaHR0cDovL3d3dy53My5vcmcvMTk5OS8wMi8yMi1yZGYtc3ludGF4LW5zIyI+CiAgICAg
    IDxyZGY6RGVzY3JpcHRpb24gcmRmOmFib3V0PSIiCiAgICAgICAgICAgIHhtbG5zOmV4aWY9Imh0
    dHA6Ly9ucy5hZG9iZS5jb20vZXhpZi8xLjAvIj4KICAgICAgICAgPGV4aWY6UGl4ZWxZRGltZW5z
    aW9uPjY0NDwvZXhpZjpQaXhlbFlEaW1lbnNpb24+CiAgICAgICAgIDxleGlmOlBpeGVsWERpbWVu
    c2lvbj4zMDY8L2V4aWY6UGl4ZWxYRGltZW5zaW9uPgogICAgICAgICA8ZXhpZjpVc2VyQ29tbWVu
    dD5TY3JlZW5zaG90PC9leGlmOlVzZXJDb21tZW50PgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4K
    ICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KJca5fwAAAAlwSFlzAAAWJQAAFiUBSVIk8AAAAwBQ
    TFRF/83KbgFIcAFJQQFUAAAARAJWXahw/////wAAcAJKocRNagFFQABS///7jYxe//zb/+Xj/4aA
    /vi5/5+a/xsW/0xF///+dQ9O/8nH//b1n8NLWaduz4qd//70+Pn5W6hvy6hPazl5aCJS+/z7jTFj
    eRNS/f79fh5WgLyR//vS/vnGWKZtbBdPcbSCcQlKUhdi1urc/vnA//rL/8zG/sbFEBAQXCZs/vzk
    9PX1/Mm///3rrlyACAgIqFV6u6TC352r2JWllz9sTA9dVyBnrcxl7vDvw67JtJq74NXjum6LMDAw
    z73UPgBQeriJRwlZk5OT4KOtiYlbW1tb+ce0ZDByQUFB56iy6ejp5OLkZzFY0Kxa7LO4giVaTExM
    5/DSeRdSoUx0ODg45LiHZmVlg1mPckR/mMil3e3ha7B819fX+L/B6ryU4ezH87q9yePS9MOqvnWP
    78CfxeHM5fHo3Oi96q+1xn6W3rR5ubm50ObWy8vL/1RM/5SNj8OefVGJq46zm0dw08LXkjlosbGx
    ICAg2LFuiStft9J1FxcXybbO7PPb/25omnikioqKb29vZa12YKpzyoOZo4Or0tLSvr6+2Mvcp8hY
    o86xrNO43Nzd1K9ltWaGxsTEkm2dwNiOwnmSwNZ+8ri7dnZ2pZpzS6JzVqVrqqqqoqKjyN2bmsFe
    oWB9iEdr/wgIKCgot9nBVFRUsmKDR0dH9vrvdCNWimKVvt3I8PXi0eKogYGB29DfmJiY/7m01bei
    6b+zgLZkLIWKtc5eIoCEaBJNJCQk8Pbp5uZzU6Z6rKB7m5ydfHx8sM+Bya2W0N+Q+PGS2+a0/19X
    9e6BnceI9PTHhjNIb7WPpchx/PWlvdJrYa2E7++rPDVpm0xAfy9f0OO4q3mJvtqwocnLsGCC3bqp
    on9SqWGAfbS3+fnXkFp0sGqGjLtX1OjNkkVtuqWJrtCbeYCeVJ2gbEJp/gIE1h8qksDD8VA3vaOg
    O2iB50MtWH947quUMl16/ywl+4555H5Fm0w+xJV8Cw5E6QAAQhBJREFUeNrsnF9MWlkex7ncQyiR
    GJvZKU7CGqkgLxN8oJC6ovyLIGgYbBEYwSqjPGjTFB+kjdMAWodIak2w9aFV4zom003/xHHGZGpb
    kib74DbdbB8naZpJ6kwyTcY0O5smm8k87J5zL38uigjXSRU5XxDuOScxl09+93e+53cPcPhYRYqD
    EWBkGBlGhpFhZFgYGUaGkWFkGBkWRoaRYWQYGUaGhZFhZBgZRoaRYWRYGBlGhpFhZBgZFkaGkWFk
    GBlGhoWRYWQYGUaGkWEdBmSnZ1vM+dXS3oqRMYmZ/GrjHg9J+DRGlpF7gSck95J/HSPLqG1AyNtL
    pKQdI9sFGUliZIUiI+HFySONAxIEjaRa6B0dYWQ5kdmMDkcDzxZqD6sbeKRB7VAbSJ6toQEeGW0k
    RrYTmWHAZDJNGRwtp92mNZsjHDPFZtTChqmWcGzMtGbAyHYgIx2mWCjktzVMtZkGHOrw7ExoZjZs
    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