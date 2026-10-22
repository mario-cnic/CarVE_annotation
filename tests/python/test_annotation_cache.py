"""Predictor score cache (annotation_cache.py): version identity, split, segments, verify, seed.

Synthetic fixtures only; nothing here touches the real cache or cluster data.
"""
import copy
import gzip
import json
import os
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import annotation_cache as ac

GENOMICS_BIN = "/data_lab_PGP/shared/utils/conda_envs/genomics/bin"


def find(name):
    return shutil.which(name) or os.path.join(GENOMICS_BIN, name)


BCFTOOLS, TABIX, BGZIP = find("bcftools"), find("tabix"), find("bgzip")
needs_bcftools = pytest.mark.skipif(
    not all(os.access(t, os.X_OK) for t in (BCFTOOLS, TABIX, BGZIP)), reason="bcftools/tabix/bgzip not available")

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FILTER=<ID=LowQual,Description="low">
##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="--input S1.cram">
##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">
##INFO=<ID=Pangolin,Number=.,Type=String,Description="Pangolin splice score predictions">
##INFO=<ID=PangolinTissue,Number=.,Type=String,Description="tissue">
##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""


def vcf(path, records):
    """records: (pos, ref, alt, info) on chr1 with one sample."""
    text = HEADER + "".join(f"chr1\t{p}\trs{p}\t{r}\t{a}\t50\tPASS\t{i or '.'}\tGT\t0/1\n" for p, r, a, i in records)
    plain = str(path)[:-3]
    with open(plain, "w") as f:
        f.write(text)
    subprocess.run([BCFTOOLS, "view", plain, "-Oz", "-o", str(path)], check=True)
    subprocess.run([TABIX, "-f", "-p", "vcf", str(path)], check=True)
    return path


def body(path):
    return [l.rstrip("\n").split("\t") for l in gzip.open(path, "rt") if not l.startswith("#")]


def header(path):
    return [l.rstrip("\n") for l in gzip.open(path, "rt") if l.startswith("##")]


def segment(tmp_path, records, name="seg"):
    segdir = tmp_path / "segments"
    segdir.mkdir(exist_ok=True)
    src = vcf(tmp_path / f"{name}_src.vcf.gz", records)
    ac.write_segment("pangolin", src, segdir / name, {"version_id": "v", "run_id": name}, BCFTOOLS, TABIX, BGZIP)
    return segdir


# ---------------------------------------------------------------- version identity

NF = """// two processes, as in older module files
process OTHER {
    script:
    \"\"\"
    other --flag
    \"\"\"
}

process PANGOLIN_ANNOTATE {
    label 'process_long'
    script:
    \"\"\"
    pangolin -d ${params.pangolin_distance}
    \"\"\"

    stub:
    \"\"\"
    cp a b
    \"\"\"
}
"""


def test_script_block_takes_the_named_process_only():
    block = ac.script_block(NF, "PANGOLIN_ANNOTATE")
    assert "pangolin -d" in block and "other --flag" not in block and "cp a b" not in block
    with pytest.raises(ac.CacheError):
        ac.script_block(NF, "MISSING")


def git(repo, *args):
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", *args],
                          cwd=repo, capture_output=True, check=True, text=True).stdout.strip()


@pytest.fixture
def repo(tmp_path):
    root = tmp_path / "repo"
    (root / "src/external/Pangolin-main/pangolin/models").mkdir(parents=True)
    (root / "modules/local").mkdir(parents=True)
    (root / "src/external/Pangolin-main/pangolin/pangolin.py").write_text("score()\n")
    (root / "src/external/Pangolin-main/pangolin/models/final.1").write_bytes(b"weights")
    (root / "src/external/Pangolin-main/pangolin/.fuse_hidden01").write_text("junk\n")
    (root / "modules/local/pangolin.nf").write_text(NF)
    git(root, "init", "-q")
    git(root, "add", "-A")
    git(root, "commit", "-qm", "init")
    return root


def manifest(commit):
    return {"launch": {
        "annotation_pipeline_repo": {"commit": commit, "dirty": False, "dirty_files": []},
        "pangolin_distance": "4999",
        "external_provenance": {
            "resources": {"fasta": {"size_bytes": 100, "sha256": None, "sidecars": {".fai": {"sha256": "fai1"}}},
                          "pangolin_db": {"sha256": "db1"}},
            "envs": {"pangolin_python": {"fingerprint_sha256": "env1"}}}}}


def test_version_identity_is_deterministic_and_skips_hidden_files(repo):
    m = manifest(git(repo, "rev-parse", "HEAD"))
    a = ac.version_identity(m, "pangolin", repo)
    assert ac.version_id(a) == ac.version_id(ac.version_identity(copy.deepcopy(m), "pangolin", repo))
    assert not any(".fuse_hidden" in k for k in a["code"])
    assert "src/external/Pangolin-main/pangolin/models/final.1" in a["code"]


@pytest.mark.parametrize("path,value", [
    (("external_provenance", "resources", "fasta", "sidecars", ".fai", "sha256"), "fai2"),
    (("external_provenance", "resources", "pangolin_db", "sha256"), "db2"),
    (("external_provenance", "envs", "pangolin_python", "fingerprint_sha256"), "env2"),
    (("pangolin_distance",), "50"),
])
def test_each_identity_field_changes_the_version(repo, path, value):
    m = manifest(git(repo, "rev-parse", "HEAD"))
    before = ac.version_id(ac.version_identity(m, "pangolin", repo))
    d = m["launch"]
    for k in path[:-1]:
        d = d[k]
    d[path[-1]] = value
    assert ac.version_id(ac.version_identity(m, "pangolin", repo)) != before


def test_code_and_command_changes_change_the_version_but_hidden_files_do_not(repo):
    v0 = ac.version_id(ac.version_identity(manifest(git(repo, "rev-parse", "HEAD")), "pangolin", repo))
    (repo / "src/external/Pangolin-main/pangolin/.fuse_hidden01").write_text("other junk\n")
    git(repo, "commit", "-qam", "hidden")
    assert ac.version_id(ac.version_identity(manifest(git(repo, "rev-parse", "HEAD")), "pangolin", repo)) == v0
    (repo / "modules/local/pangolin.nf").write_text(NF.replace("pangolin -d", "pangolin -m True -d"))
    git(repo, "commit", "-qam", "flag")
    assert ac.version_id(ac.version_identity(manifest(git(repo, "rev-parse", "HEAD")), "pangolin", repo)) != v0


def test_unknown_provenance_or_dirty_code_is_refused(repo):
    commit = git(repo, "rev-parse", "HEAD")
    m = manifest(commit)
    del m["launch"]["external_provenance"]["resources"]["pangolin_db"]
    with pytest.raises(ac.CacheError, match="pangolin_db"):
        ac.version_identity(m, "pangolin", repo)
    m = manifest(commit)
    m["launch"]["external_provenance"]["envs"] = {}
    with pytest.raises(ac.CacheError, match="env fingerprint"):
        ac.version_identity(m, "pangolin", repo)
    m = manifest(commit)
    m["launch"]["annotation_pipeline_repo"]["dirty_files"] = [" M modules/local/pangolin.nf"]
    with pytest.raises(ac.CacheError, match="uncommitted"):
        ac.version_identity(m, "pangolin", repo)


def test_version_dir_is_created_once_and_checked(tmp_path, repo):
    ident = ac.version_identity(manifest(git(repo, "rev-parse", "HEAD")), "pangolin", repo)
    vdir = ac.ensure_version_dir(tmp_path / "cache", ident)
    assert vdir.name == ac.version_id(ident) and (vdir / "segments").is_dir()
    assert ac.ensure_version_dir(tmp_path / "cache", ident) == vdir
    data = json.loads((vdir / "VERSION.json").read_text())
    data["identity"]["distance"] = "1"
    (vdir / "VERSION.json").write_text(json.dumps(data))
    with pytest.raises(ac.CacheError, match="differs"):
        ac.ensure_version_dir(tmp_path / "cache", ident)


# ---------------------------------------------------------------- segments and split

@needs_bcftools
def test_segment_is_sites_only_with_minimal_header(tmp_path):
    segdir = segment(tmp_path, [(100, "A", "G", "DP=5;Pangolin=G1|10:0.5;PangolinTissue=G1|Heart|Liver"),
                                (200, "C", "T", "DP=7")])
    seg = segdir / "seg.vcf.gz"
    rows = body(seg)
    assert [r[:8] for r in rows] == [
        ["chr1", "100", ".", "A", "G", ".", ".", "Pangolin=G1|10:0.5;PangolinTissue=G1|Heart|Liver"],
        ["chr1", "200", ".", "C", "T", ".", ".", "."]]
    assert all(len(r) == 8 for r in rows)
    hdr = header(seg)
    assert not any(l.startswith(("##GATKCommandLine", "##INFO=<ID=DP", "##FORMAT")) for l in hdr)
    assert "##annotation_cache_predictor=pangolin" in hdr
    side = json.loads((segdir / "seg.json").read_text())
    assert side["n_records"] == 2 and side["sha256"] == ac.sha256_file(seg)


@needs_bcftools
def test_split_without_cache_sends_everything_to_misses(tmp_path):
    inp = vcf(tmp_path / "in.vcf.gz", [(100, "A", "G", "DP=5"), (200, "C", "T", "DP=7")])
    res = ac.split("pangolin", inp, tmp_path / "none", tmp_path / "hits.vcf.gz", tmp_path / "misses.vcf.gz",
                   tmp_path / "r.json", BCFTOOLS, TABIX)
    assert (res["hits"], res["misses"]) == (0, 2)
    assert body(tmp_path / "hits.vcf.gz") == [] and len(body(tmp_path / "misses.vcf.gz")) == 2


@needs_bcftools
def test_split_hits_misses_unscored_and_same_position_alleles(tmp_path):
    segdir = segment(tmp_path, [(100, "A", "G", "Pangolin=G1|10:0.5"), (300, "C", "T", "DP=1")])
    inp = vcf(tmp_path / "in.vcf.gz", [
        (100, "A", "G", "DP=5"),   # cached with a score
        (100, "A", "T", "DP=5"),   # same position, other allele: not cached
        (300, "C", "T", "DP=9"),   # cached, computed without a score
        (400, "G", "A", "DP=2")])  # not cached
    res = ac.split("pangolin", inp, segdir, tmp_path / "hits.vcf.gz", tmp_path / "misses.vcf.gz",
                   tmp_path / "r.json", BCFTOOLS, TABIX)
    assert (res["input_records"], res["hits"], res["misses"]) == (4, 2, 2)
    hits = {(r[1], r[4]): r for r in body(tmp_path / "hits.vcf.gz")}
    assert hits[("100", "G")][7] == "DP=5;Pangolin=G1|10:0.5" and hits[("100", "G")][9] == "0/1"
    assert hits[("300", "T")][7] == "DP=9"
    assert {(r[1], r[4]) for r in body(tmp_path / "misses.vcf.gz")} == {("100", "T"), ("400", "A")}
    assert json.loads((tmp_path / "r.json").read_text())["segments_used"] == ["seg.vcf.gz"]


@needs_bcftools
def test_damaged_segment_is_skipped(tmp_path):
    segdir = segment(tmp_path, [(100, "A", "G", "Pangolin=G1|10:0.5")])
    side = segdir / "seg.json"
    meta = json.loads(side.read_text())
    meta["sha256"] = "0" * 64
    side.write_text(json.dumps(meta))
    good, skipped = ac.list_segments(segdir)
    assert good == [] and "sha256" in skipped[0]
    inp = vcf(tmp_path / "in.vcf.gz", [(100, "A", "G", "DP=5")])
    res = ac.split("pangolin", inp, segdir, tmp_path / "h.vcf.gz", tmp_path / "m.vcf.gz",
                   tmp_path / "r.json", BCFTOOLS, TABIX)
    assert (res["hits"], res["misses"]) == (0, 1)


@needs_bcftools
def test_make_segment_keeps_only_computed_misses(tmp_path):
    annotated = vcf(tmp_path / "ann.vcf.gz", [(100, "A", "G", "Pangolin=hit"), (200, "C", "T", "Pangolin=new")])
    misses = vcf(tmp_path / "misses.vcf.gz", [(200, "C", "T", "DP=1")])
    (tmp_path / "out").mkdir()
    info = ac.make_segment("pangolin", annotated, misses, tmp_path / "out" / "run1", {"version_id": "v", "run_id": "run1"},
                           BCFTOOLS, TABIX, BGZIP)
    assert info["n_records"] == 1
    assert [r[1] for r in body(tmp_path / "out" / "run1.vcf.gz")] == ["200"]


@needs_bcftools
def test_two_segments_combine_and_verify_flags_conflicts(tmp_path):
    segment(tmp_path, [(100, "A", "G", "Pangolin=x"), (200, "C", "T", "Pangolin=y")], "a")
    segdir = segment(tmp_path, [(200, "C", "T", "Pangolin=y"), (300, "G", "A", "Pangolin=z")], "b")
    inp = vcf(tmp_path / "in.vcf.gz", [(100, "A", "G", ""), (200, "C", "T", ""), (300, "G", "A", "")])
    res = ac.split("pangolin", inp, segdir, tmp_path / "h.vcf.gz", tmp_path / "m.vcf.gz",
                   tmp_path / "r.json", BCFTOOLS, TABIX)
    assert (res["hits"], res["misses"]) == (3, 0)
    ok = ac.verify(tmp_path, "pangolin", BCFTOOLS)
    assert (ok["segments"], ok["duplicates"], ok["conflicts"]) == (2, 1, 0)
    segment(tmp_path, [(300, "G", "A", "Pangolin=DIFFERENT")], "c")
    bad = ac.verify(tmp_path, "pangolin", BCFTOOLS)
    assert bad["conflicts"] == 1 and bad["examples"][0]["site"] == "chr1:300:G:A"


def test_find_conflicts_groups_by_position_and_allele():
    rows = [["chr1", "1", "A", "G", "x"], ["chr1", "1", "A", "T", "y"], ["chr1", "1", "A", "G", "x"],
            ["chr1", "2", "C", "T", "p"], ["chr1", "2", "C", "T", "q"]]
    dups, conflicts, ex = ac.find_conflicts(rows)
    assert (dups, conflicts) == (2, 1) and ex[0]["site"] == "chr1:2:C:T"


# ---------------------------------------------------------------- seed

@needs_bcftools
def test_seed_writes_new_sites_and_refuses_conflicts(tmp_path, repo):
    commit = git(repo, "rev-parse", "HEAD")
    cache = tmp_path / "cache"

    def run_dir(name, records):
        d = tmp_path / name
        (d / "pipeline_info").mkdir(parents=True)
        (d / "predictors").mkdir()
        (d / "pipeline_info" / "RUN_MANIFEST.json").write_text(json.dumps(manifest(commit)))
        vcf(d / "predictors" / f"{name}.annPangolin.vcf.gz", records)
        return d

    r1 = run_dir("R1", [(100, "A", "G", "Pangolin=x"), (200, "C", "T", "DP=1")])
    first = ac.seed(r1, "pangolin", cache, None, None, BCFTOOLS, TABIX, BGZIP, repo)
    assert first["n_records"] == 2

    r2 = run_dir("R2", [(100, "A", "G", "Pangolin=x"), (300, "G", "A", "Pangolin=z")])
    second = ac.seed(r2, "pangolin", cache, None, None, BCFTOOLS, TABIX, BGZIP, repo)
    assert second["n_records"] == 1 and second["overlap_with_cache"] == {"shared_sites": 1, "conflicts": 0}

    r3 = run_dir("R3", [(300, "G", "A", "Pangolin=OTHER")])
    with pytest.raises(ac.CacheError, match="differ from the cache"):
        ac.seed(r3, "pangolin", cache, None, None, BCFTOOLS, TABIX, BGZIP, repo)
    vdir = cache / "pangolin" / first["version_id"]
    assert len(list((vdir / "segments").glob("*.json"))) == 2

    m = manifest(commit)
    m["launch"]["pangolin_distance"] = "50"
    (r3 / "pipeline_info" / "RUN_MANIFEST.json").write_text(json.dumps(m))
    with pytest.raises(ac.CacheError, match="distance"):
        ac.seed(r3, "pangolin", cache, first["version_id"], None, BCFTOOLS, TABIX, BGZIP, repo)


@needs_bcftools
def test_segment_defines_fields_even_when_the_source_header_lacks_them(tmp_path):
    src = vcf(tmp_path / "src.vcf.gz", [(100, "A", "G", "DP=5")])
    plain = tmp_path / "noinfo.vcf"
    plain.write_text("".join(l for l in gzip.open(src, "rt") if not l.startswith("##INFO=<ID=Pangolin")))
    nohdr = tmp_path / "noinfo.vcf.gz"
    subprocess.run([BCFTOOLS, "view", str(plain), "-Oz", "-o", str(nohdr)], check=True)
    (tmp_path / "segments").mkdir()
    ac.write_segment("pangolin", nohdr, tmp_path / "segments" / "s", {}, BCFTOOLS, TABIX, BGZIP)
    hdr = header(tmp_path / "segments" / "s.vcf.gz")
    assert any(l.startswith("##INFO=<ID=Pangolin,") for l in hdr)
    assert any(l.startswith("##INFO=<ID=PangolinTissue,") for l in hdr)
    inp = vcf(tmp_path / "in.vcf.gz", [(100, "A", "G", "DP=5")])
    res = ac.split("pangolin", inp, tmp_path / "segments", tmp_path / "h.vcf.gz", tmp_path / "m.vcf.gz",
                   tmp_path / "r.json", BCFTOOLS, TABIX)
    assert (res["hits"], res["misses"]) == (1, 0)


@needs_bcftools
def test_segment_list_names_every_file_of_valid_segments_only(tmp_path):
    segment(tmp_path, [(100, "A", "G", "Pangolin=x")], "good")
    segdir = segment(tmp_path, [(200, "C", "T", "Pangolin=y")], "bad")
    os.remove(segdir / "bad.vcf.gz.tbi")
    n, skipped = ac.write_segment_list(tmp_path, tmp_path / "list.txt")
    lines = (tmp_path / "list.txt").read_text().splitlines()
    assert n == 1 and "bad.vcf.gz" in skipped[0]
    assert [os.path.basename(l) for l in lines] == ["good.vcf.gz", "good.vcf.gz.tbi", "good.json"]
    assert all(os.path.isabs(l) and os.path.exists(l) for l in lines)
