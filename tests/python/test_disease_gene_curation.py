"""Disease-gene curation list is read from the given path and a missing or empty list is an error."""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src/python'))
from parsing.modules.disease_hpo import apply_disease_phenotype_curation, load_cardiac_disease_gene_curation

DISEASE_GENE_LIST = str(Path(__file__).resolve().parents[1] / 'fixtures/mini_disease_gene_list.csv')


def test_listed_genes_are_curated_and_others_are_background():
    out = apply_disease_phenotype_curation(pd.DataFrame({'SYMBOL': ['MYH7', 'A1BG', 'NOTAGENE']}), DISEASE_GENE_LIST)
    assert out['DISEASE_GENE_EVIDENCE'].iloc[0] == 'Definitive / Primary Evidence'
    assert out['DISEASE_GENE_EVIDENCE'].iloc[1] == 'Moderate / Candidate Evidence'
    assert out['DISEASE_GENE_EVIDENCE'].iloc[2] == 'Uncurated / Other Genomic'
    assert out['DISEASE_CURATION_SOURCE'].iloc[2] == 'Genome-Wide Background'


def test_missing_list_fails(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_cardiac_disease_gene_curation(str(tmp_path / 'absent.csv'))
    with pytest.raises(FileNotFoundError):
        apply_disease_phenotype_curation(pd.DataFrame({'SYMBOL': ['MYH7']}), str(tmp_path / 'absent.csv'))


def test_list_without_genes_fails(tmp_path):
    empty = tmp_path / 'empty.csv'
    empty.write_text('GENE,NAME,SOURCE,EVIDENCE\n')
    with pytest.raises(ValueError):
        load_cardiac_disease_gene_curation(str(empty))
