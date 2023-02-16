Pipeline steps
==============

.. note::
    
    Note that the main script ``aspring`` **runs the entire pipeline** automatically.


step_01_preprocess
------------------

.. argparse::
   :ref: aspring.step_01_preprocess.get_arg_parser
   :prog: step_01_preprocess

step_02_hmm_maker
-----------------

.. argparse::
   :ref: aspring.step_02_hmm_maker.get_arg_parser
   :prog: step_02_hmm_maker

step_03_hmm_aligner
-------------------

.. argparse::
    :ref: aspring.step_03_hmm_aligner.get_arg_parser
    :prog: step_03_hmm_aligner

step_04_gettable
----------------

.. argparse::
    :ref: aspring.step_04_gettable.get_arg_parser
    :prog: step_04_gettable

step_05_filter
--------------

.. argparse::
    :ref: aspring.step_05_filter.get_arg_parser
    :prog: step_05_filter

step_06_stats
-------------

.. argparse::
    :ref: aspring.step_06_stats.get_arg_parser
    :prog: step_06_stats

step_07_reformat
----------------

.. argparse::
    :ref: aspring.step_07_reformat.get_arg_parser
    :prog: step_07_reformat

step_08_ASRUs
-------------

.. argparse::
    :ref: aspring.step_08_ASRUs.get_arg_parser
    :prog: step_08_ASRUs

step_09_clean
-------------

.. argparse::
    :ref: aspring.step_09_clean.get_arg_parser
    :prog: step_09_clean