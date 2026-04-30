import os
import shutil
import subprocess

import pytest

SMILES_TO_HELM = 'sketcher_app/smiles_to_helm'

DB_ENV_VAR = 'SCHRODINGER_CUSTOM_MONOMER_DB_PATH'

SKETCHER_SOURCE_DIR = os.environ['SKETCHER_SOURCE_DIR']


@pytest.fixture(autouse=True)
def tmp_db(tmp_path):
    """
    Run each test case with its own temporary monomer database, both to avoid
    picking up a user's custom database, and because some tests have the side
    effect of writing to the database.
    """
    db_path = str(tmp_path / "test.db")
    old_db = os.environ.get(DB_ENV_VAR)
    try:
        os.environ[DB_ENV_VAR] = db_path
        yield
    finally:
        if old_db is None:
            del os.environ[DB_ENV_VAR]
        else:
            os.environ[DB_ENV_VAR] = old_db


@pytest.fixture()
def sqlite_db(tmp_path):
    """
    Make a copy of the sqlite test db and return its path. This is to prevent
    modifications to the original.
    """
    src = f'{SKETCHER_SOURCE_DIR}/test/testfiles/X1.db'
    dest = tmp_path / 'X1.db'
    shutil.copy(src, dest)
    return str(dest)


@pytest.fixture()
def json_db():
    return f'{SKETCHER_SOURCE_DIR}/test/testfiles/X1.json'


def smiles_to_helm(smiles: str, db: str | None = None) -> str:
    """
    Helper function to run smiles_to_helm with the given input and optional
    monomer database file.

    :param smiles: stdin to feed into smiles_to_helm
    :param db: custom monomer databases relative to this test module
    :return: stdout from smiles_to_helm
    """
    cmd = [SMILES_TO_HELM]
    if db:
        cmd += ['--db', db]
    p = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='utf8',
    )
    helm_out, err = p.communicate(smiles)
    if p.returncode != 0:
        raise Exception(f'smiles_to_helm returned {p.returncode}:\n{err}')
    return helm_out


def test_smiles_to_helm():
    """
    This test has two SMILES, one with a name, and only uses standard monomers.
    """
    smiles = "NCC(=O)NCC(=O)O gly 1\nN[C@@H](C)C(=O)O\n"
    helm = smiles_to_helm(smiles)
    expected_helm = 'PEPTIDE1{G.G}$$$$V2.0 gly 1\nPEPTIDE1{A}$$$$V2.0\n'
    assert helm == expected_helm


def test_X1_no_db():
    """
    X1 is a custom monomer from ChEMBL. Here we run without a DB, so we
    get an inline SMILES.
    """
    smiles = 'O=C(C1(N)CCCCC1)O\n'
    helm = smiles_to_helm(smiles)
    assert helm == 'CHEM1{[NC1(C(=O)O)CCCCC1]}$$$$V2.0\n'


def test_X1_json_db(json_db):
    smiles = 'O=C(C1(N)CCCCC1)O\n'
    helm = smiles_to_helm(smiles, db=json_db)
    assert helm == 'PEPTIDE1{[X1]}$$$$V2.0\n'


def test_X1_sqlite_db(sqlite_db):
    smiles = 'O=C(C1(N)CCCCC1)O\n'
    helm = smiles_to_helm(smiles, db=sqlite_db)
    assert helm == 'PEPTIDE1{[X1]}$$$$V2.0\n'


def test_bad_smiles():
    with pytest.raises(Exception, match='Error processing SMILES'):
        smiles_to_helm('x\n')


def test_bad_db():
    with pytest.raises(Exception, match='Error loading database'):
        smiles_to_helm('O\n', db='does_not_exist.json')
