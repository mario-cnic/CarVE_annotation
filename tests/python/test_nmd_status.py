"""Synthetic checks for NMD_STATUS / NMD_RULE / NMD_PTC_BASIS judged at the premature stop codon."""
import gzip
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src/python'))
import build_transcript_exon_table as bte
from parsing.modules.predictors import load_transcript_exons, parse_nmd

# TX1: 4 exons of 100, 200, 150, 120 nt (cDNA ends 100, 300, 450, 570); 5'UTR 50 nt, so CDS c.1 = cDNA 51.
# Penultimate exon ends at cDNA 450: a stop starting at cDNA 401..450 (c.351..400) is in its last 50 nt.
EXONS = {
    'ENST01': ('1', (100, 300, 450, 570), 51, False, False),
    'ENST02': ('3', (900,), 11, False, False),                  # intronless
    'ENST03': ('1', (100, 300, 450, 570), 51, False, True),     # 3' structure incomplete
    'ENST04': ('1', (100, 300, 450, 570), 51, True, False),     # 5' end of the CDS incomplete
}


def row(conseq, hgvsc, hgvsp, cds, feature='ENST01', version='1', utr=50):
    cds_first = int(str(cds).split('-')[0])
    cdna = f"{cds_first + utr}-{cds_first + utr + 1}" if '-' in str(cds) else str(cds_first + utr)
    return {'Consequence': conseq, 'Feature': feature, 'HGVSc': f'{feature}.{version}:{hgvsc}',
            'HGVSp': f'ENSP01.1:{hgvsp}' if hgvsp else None, 'CDS_position': str(cds), 'cDNA_position': cdna}


def judge(rows):
    out = parse_nmd(pd.DataFrame(rows), EXONS)
    return list(zip(out.NMD_STATUS, out.NMD_RULE, out.NMD_PTC_BASIS))


def test_stop_gained_rules_at_the_stop_codon():
    assert judge([
        row('stop_gained', 'c.298C>T', 'p.Arg100Ter', 298),     # cDNA 348: exon 3, far from the last junction
        row('stop_gained', 'c.349C>T', 'p.Gln117Ter', 349),     # cDNA 399: 52 nt before the junction
        row('stop_gained', 'c.352C>T', 'p.Gln118Ter', 352),     # cDNA 402: inside the last 50 nt
        row('stop_gained', 'c.403C>T', 'p.Gln135Ter', 403),     # cDNA 453: last exon
        row('stop_gained', 'c.58C>T', 'p.Gln20Ter', 58),        # first 100 coding nt
    ]) == [
        ('nmd_predicted', 'none', 'ptc'),
        ('nmd_predicted', 'none', 'ptc'),
        ('escape_predicted', 'penultimate_last_50nt', 'ptc'),
        ('escape_predicted', 'last_exon', 'ptc'),
        ('nmd_predicted', 'start_proximal', 'ptc'),
    ]


def test_frameshift_judged_at_new_stop_not_variant():
    # deletion in exon 2 (cDNA 250); the new stop codon is aa 67 + 70 - 1 = 136, c.406 in the variant transcript,
    # c.407 in the reference (one base deleted upstream) = cDNA 457, in the last exon
    assert judge([row('frameshift_variant', 'c.200del', 'p.Lys67ArgfsTer70', 200)]) == [
        ('escape_predicted', 'last_exon', 'ptc')]
    # same deletion, stop soon after: aa 67 + 5 - 1 = 71 -> c.212 - (-1) = c.213 = cDNA 263, exon 2
    assert judge([row('frameshift_variant', 'c.200del', 'p.Lys67ArgfsTer5', 200)]) == [
        ('nmd_predicted', 'none', 'ptc')]
    # insertion of one base: stop at aa 114 + 5 - 1 = 118, c.352 in the variant transcript, c.351 = cDNA 401
    assert judge([row('frameshift_variant', 'c.340_341insA', 'p.Leu114ProfsTer5', '340-341')]) == [
        ('escape_predicted', 'penultimate_last_50nt', 'ptc')]


def test_variant_start_proximal_but_stop_downstream_is_nmd():
    # the VEP plugin flags any variant in c.1..101; the new stop (aa 30 + 60 - 1 = 89, c.265) is far from it
    assert judge([row('frameshift_variant', 'c.88del', 'p.Ala30GlnfsTer60', 88)]) == [('nmd_predicted', 'none', 'ptc')]


def test_no_stop_in_new_frame():
    assert judge([
        row('frameshift_variant', 'c.420del', 'p.Ala140GlnfsTer?', 420),   # cDNA 470: last exon
        row('frameshift_variant', 'c.200del', 'p.Lys67ArgfsTer?', 200),    # exon 2
    ]) == [('escape_predicted', 'last_exon', 'variant'), ('not_evaluated', 'no_stop_in_new_frame', '')]


def test_other_ptc_forms_and_unknowns():
    assert judge([
        row('stop_gained', 'c.28_33delinsTGGTAA', 'p.Gly10_Leu11delinsTrpTer', 28),     # stop at aa 11, c.31
        row('stop_gained&frameshift_variant', 'c.298del', 'p.Arg100Ter', 298),
        row('frameshift_variant&start_lost', 'c.1del', 'p.Met1?', 1),
        row('frameshift_variant', 'c.200+1_200+2del', 'p.Lys67fs', 200),
    ]) == [
        ('nmd_predicted', 'start_proximal', 'ptc'),
        ('nmd_predicted', 'none', 'ptc'),
        ('not_evaluated', 'ptc_position_unknown', ''),
        ('not_evaluated', 'ptc_position_unknown', ''),
    ]


def test_transcript_model_checks():
    bad_utr = row('stop_gained', 'c.298C>T', 'p.Arg100Ter', 298, utr=49)
    no_hgvs = row('stop_gained', 'c.298C>T', 'p.Arg100Ter', 298)
    no_hgvs['HGVSc'] = None
    assert judge([
        row('stop_gained', 'c.298C>T', 'p.Arg100Ter', 298, feature='ENST99'),
        row('stop_gained', 'c.298C>T', 'p.Arg100Ter', 298, version='2'),
        bad_utr,
        no_hgvs,
        row('stop_gained', 'c.298C>T', 'p.Arg100Ter', 298, feature='ENST03'),
        row('stop_gained', 'c.58C>T', 'p.Gln20Ter', 58, feature='ENST04'),
        row('stop_gained', 'c.298C>T', 'p.Arg100Ter', 298, feature='ENST02', version='3', utr=10),
        {'Consequence': 'stop_gained', 'Feature': 'ENST01', 'HGVSc': 'ENST01.1:c.298C>T', 'HGVSp': 'ENSP01.1:p.Arg100Ter',
         'CDS_position': '?-12', 'cDNA_position': '?-62'},
    ]) == [
        ('not_evaluated', 'transcript_not_in_model', ''),
        ('not_evaluated', 'transcript_version_mismatch', ''),
        ('not_evaluated', 'model_mismatch', ''),
        ('not_evaluated', 'hgvs_missing', ''),
        ('not_evaluated', 'incomplete_3prime', ''),
        ('nmd_predicted', 'none', 'ptc'),          # start-proximal not judged when the CDS start is unknown
        ('escape_predicted', 'intronless', 'ptc'),
        ('not_evaluated', 'position_unknown', ''),
    ]


def test_scope_splice_and_missing_table():
    rows = [row('splice_donor_variant', 'c.200+1G>A', None, 200), row('missense_variant', 'c.10A>G', 'p.Lys4Glu', 10),
            row('stop_gained', 'c.298C>T', 'p.Arg100Ter', 298)]
    df = pd.DataFrame(rows, index=[7, 3, 5])
    out = parse_nmd(df.copy(), EXONS)
    assert out.NMD_STATUS.tolist() == ['not_evaluated', 'not_applicable', 'nmd_predicted']
    assert out.NMD_RULE.tolist() == ['splice_outcome_unknown', '', 'none']
    assert out.index.equals(df.index)
    pd.testing.assert_frame_equal(out[df.columns], df)
    out = parse_nmd(df.copy(), None)
    assert out.NMD_STATUS.tolist() == ['not_annotated', 'not_applicable', 'not_annotated']
    assert out.NMD_RULE.tolist() == ['no_transcript_model', '', 'no_transcript_model']


def test_build_exon_table_minus_strand(tmp_path):
    attrs = 'gene_id "G"; transcript_id "ENST05.2";'
    lines = [
        ('transcript', 500, 1099, attrs + ' tag "basic"; tag "cds_end_NF";'),
        ('exon', 500, 649, attrs), ('exon', 800, 899, attrs), ('exon', 1000, 1099, attrs),
        ('CDS', 560, 649, attrs), ('CDS', 800, 899, attrs), ('CDS', 1000, 1049, attrs),
        ('stop_codon', 557, 559, attrs),
    ]
    gtf = tmp_path / 'mini.gtf.gz'
    with gzip.open(gtf, 'wt') as fh:
        for feat, s, e, a in lines:
            fh.write(f'chr1\tTEST\t{feat}\t{s}\t{e}\t.\t-\t.\t{a}\n')
    out = tmp_path / 'exons.tsv'
    bte.main(['--gtf', str(gtf), '--out', str(out)])
    t = pd.read_csv(out, sep='\t', dtype=str).iloc[0].to_dict()
    # transcript order 1000-1099, 800-899, 500-649; CDS starts at 1049 (cDNA 51), stop codon ends at 557 (cDNA 293)
    assert t == {'transcript_id': 'ENST05', 'transcript_version': '2', 'chrom': 'chr1', 'strand': '-',
                 'exon_cdna_ends': '100,200,350', 'cds_start_cdna': '51', 'cds_end_cdna': '293',
                 'start_incomplete': '0', 'end_incomplete': '1'}
    assert load_transcript_exons(str(out)) == {'ENST05': ('2', (100, 200, 350), 51, False, True)}
