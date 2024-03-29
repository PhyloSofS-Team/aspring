import warnings
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
    filename = request.module.__file__
    test_dir = os.path.dirname(filename)
    path_data = os.path.join(test_dir, 'data')
    path_gene = os.path.join(path_data, 'data', 'ENSG00000007866')
    path_ref = os.path.join(path_data, 'ENSG00000007866', 'ASPRING_reference')

    path_hhsuite_scripts = None
    possible_paths = [
        # Local
        '/home/diego/bin/miniconda3/scripts',
        # GitHUb action CI
        '/usr/share/miniconda/envs/test/scripts',
        '/usr/local/miniconda/envs/test/scripts']
    for path in possible_paths:
        if os.path.isdir(path):
            path_hhsuite_scripts = path
            break
    if path_hhsuite_scripts is None:
        raise Exception(
            'path_hhsuite_scripts is None: please set it manually to test aspring in your machine')

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
    # show the files in the data folder
    print("path_gene :", os.listdir(path_gene))

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
        by=['instance', 'gene', 'size', 'NbSex']).reset_index(drop=True)
    df_ref_sorted = df_ref[columns_to_compare].sort_values(
        by=['instance', 'gene', 'size', 'NbSex']).reset_index(drop=True)

    assert (df_out_sorted == df_ref_sorted).all().all()


def test_events(setup_and_run_pipeline):
    path_gene = setup_and_run_pipeline['path_gene']
    path_ref = setup_and_run_pipeline['path_ref']

    filename = "ENSG00000007866_eventsDup_withCols.txt"
    out = os.path.join(path_gene, filename)
    ref = os.path.join(path_ref, filename)

    columns_to_compare = ['gene', 'rank', 'type', 'lePathA', 'lePathB', 'exclu',
                          'ncols', 'leA', 'leB', 'typePair', 'ColA', 'ColB']

    df_out = pd.read_csv(out, index_col=False)
    df_ref = pd.read_csv(ref, index_col=False)

    df_out_sorted = df_out[columns_to_compare]
    df_ref_sorted = df_ref[columns_to_compare]

    assert (df_out_sorted == df_ref_sorted).all().all()


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


def is_pdb(path):
    """Check if the file at the given path is a PDB file."""
    pdb_keywords = ['HEADER', 'OBSLTE', 'TITLE', 'SPLT', 'CAVEAT', 'COMPND', 'SOURCE',
                    'KEYWDS', 'EXPDTA', 'NUMMDL', 'MDLTYP', 'AUTHOR', 'REVDAT', 'SPRSDE',
                    'JRNL', 'REMARK', 'DBREF', 'DBREF1', 'DBREF2', 'SEQADV', 'SEQRES',
                    'MODRES', 'HET', 'FORMUL', 'HETNAM', 'HETSYN', 'HELIX', 'SHEET',
                    'SSBOND', 'LINK', 'CISPEP', 'SITE', 'CRYST1', 'MTRIX1', 'MTRIX2',
                    'MTRIX3', 'ORIGX1', 'ORIGX2', 'ORIGX3', 'SCALE1', 'SCALE2', 'SCALE3',
                    'MODEL', 'ATOM', 'ANISOU', 'TER', 'HETATM', 'ENDMDL', 'CONECT',
                    'MASTER', 'END']
    with open(path, 'r') as f:
        for line in f:
            if not line.startswith(tuple(pdb_keywords)):
                print(f'The line that failed is: {line}')
                return False
    return True


def is_pml(path):
    """Check if the file at the given path is a PyMol script."""
    # list of the PyMol commands that are used in the write_pml function
    pml_keywords = ['load', 'bg_color', 'color', 'select', 'join', 'deselect']
    with open(path, 'r') as f:
        for line in f:
            if not line.startswith(tuple(pml_keywords)):
                return False
    return True


def test_structures_folder(setup_and_run_pipeline):
    path_gene = setup_and_run_pipeline['path_gene']
    structures_folder = os.path.join(path_gene, 'structures')
    assert os.path.isdir(structures_folder)

    files = os.listdir(structures_folder)
    if files:
        for file in files:
            # test that the file has the right extension
            assert file.endswith('.pdb') or file.endswith('.pml')
            # test that the file is not empty
            assert os.path.getsize(os.path.join(structures_folder, file)) > 0
            # test that the file has the right name
            if file.endswith('.pdb'):
                assert file.startswith('AF-')
                assert file.endswith('-F1-model_v4.pdb')
                # test that the file is a proper PDB file
                assert is_pdb(os.path.join(structures_folder, file))
            if file.endswith('.pml'):
                assert '_AF-' in file
                assert file.endswith('-F1-model_v4.pml')
                # test that the file is a proper PyMol script
                assert is_pml(os.path.join(structures_folder, file))
    else:
        # warn the user that no structure was generated withouth failing the test
        with pytest.warns(UserWarning):
            warnings.warn('No structures were downloaded for this gene.')
