import pytest
import os
import pandas as pd
import shutil

from aspring.main import run_pipeline


@pytest.fixture(scope='module')
def clean_up(request):
    filename = request.module.__file__
    test_dir = os.path.dirname(filename)
    path_data = os.path.join(test_dir, 'data')
    path_gene = os.path.join(path_data, 'data')
    path_dupraw = os.path.join(path_data, 'DupRaw')
    yield path_gene, path_dupraw
    os.path.isdir(path_gene) and shutil.rmtree(path_gene)
    os.path.isdir(path_dupraw) and shutil.rmtree(path_dupraw)


@pytest.fixture(scope='module')
def setup_and_run_pipeline(request, clean_up):
    # example : python3 ASPRIN/scripts/withoutDB/main.py
    # --path_data ~/lab_bench --path_hhsuite ~/hhsuite
    # --gene ENSG00000107643 --len 5 --id 100 --norealign 1 --glo_loc 1
    # --mact 0 --id_pair 50 --idCons_pair 50 --pval 0.001 --nbSpe 2 --cov 0.8
    filename = request.module.__file__
    test_dir = os.path.dirname(filename)
    path_data = os.path.join(test_dir, 'data')
    path_gene = os.path.join(path_data, 'data', 'ENSG00000007866')
    path_ref = os.path.join(path_data, 'ENSG00000007866', 'ASPRING_reference')

    # Local
    path_diego = '/home/diego/bin/miniconda3/scripts'
    if os.path.isdir(path_diego):
        path_hhsuite_scripts = path_diego
    else:
        # GitHUb action CI
        CONDA = os.environ.get('CONDA')
        path_hhsuite_scripts = os.path.join(CONDA, "/envs/test/scripts")
    # show the content of the directory
    print(f"path_hhsuite_scripts : {path_hhsuite_scripts} : ",
          os.listdir(path_hhsuite_scripts))

    gene = 'ENSG00000007866'
    msa_len = 5
    msa_id_threshold = 100
    re_align = 1
    glo_loc = 1
    mact = 0
    id_pair = 50
    idCons_pair = 50
    pval = 0.001
    nbSpe = 2
    cov = 0.8
    run_pipeline(gene, path_data, path_hhsuite_scripts, msa_len,
                 msa_id_threshold, re_align, glo_loc, mact, id_pair,
                 idCons_pair, pval, nbSpe, cov)

    # return the paths to the created files
    return {
        'path_data': path_data,
        'path_gene': path_gene,
        'path_ref': path_ref
    }


def test_instances(setup_and_run_pipeline):
    path_gene = setup_and_run_pipeline['path_gene']
    path_ref = setup_and_run_pipeline['path_ref']

    filename = "ENSG00000007866_instances_table.csv"
    out = os.path.join(path_gene, filename)
    ref = os.path.join(path_ref, filename)

    columns_to_compare = ['instance', 'size', 'NbSex', 'gene']

    df_out = pd.read_csv(out, index_col=False)
    df_ref = pd.read_csv(ref, index_col=False)

    df_out_sorted = df_out[columns_to_compare].sort_values(
        by=['instance', 'gene', 'size', 'NbSex'])
    df_ref_sorted = df_ref[columns_to_compare].sort_values(
        by=['instance', 'gene', 'size', 'NbSex'])

    assert (df_out_sorted == df_ref_sorted).all().all()


def test_events(setup_and_run_pipeline):
    path_gene = setup_and_run_pipeline['path_gene']
    path_ref = setup_and_run_pipeline['path_ref']

    filename = "ENSG00000007866_eventsDup_withCols.txt"
    out = os.path.join(path_gene, filename)
    ref = os.path.join(path_ref, filename)

    with open(out, 'r') as f_out, open(ref, 'r') as f_ref:
        string_out = f_out.read()
        string_ref = f_ref.read()

    assert string_out == string_ref


def test_duplications(setup_and_run_pipeline):
    path_gene = setup_and_run_pipeline['path_gene']
    path_ref = setup_and_run_pipeline['path_ref']

    filename = "ENSG00000007866_duplication_pairs.csv"
    out = os.path.join(path_gene, filename)
    ref = os.path.join(path_ref, filename)

    df_out = pd.read_csv(out)
    df_ref = pd.read_csv(ref)

    assert df_out.shape == df_ref.shape


def test_asrus(setup_and_run_pipeline):
    path_gene = setup_and_run_pipeline['path_gene']
    path_ref = setup_and_run_pipeline['path_ref']

    filename = "ENSG00000007866_ASRUs_table.csv"
    out = os.path.join(path_gene, filename)
    ref = os.path.join(path_ref, filename)

    df_out = pd.read_csv(out, index_col=False)
    df_ref = pd.read_csv(ref, index_col=False)

    columns_to_compare = [
        'gene', 'Nbinstances', 'max', 'min', 'moy', 'median', 'std',
        'eventsRank'
    ]

    df_out_sorted = df_out[columns_to_compare].sort_values(by=['gene'])
    df_ref_sorted = df_ref[columns_to_compare].sort_values(by=['gene'])

    assert (df_out[columns_to_compare] == df_ref[columns_to_compare]
            ).all().sum() == len(columns_to_compare)