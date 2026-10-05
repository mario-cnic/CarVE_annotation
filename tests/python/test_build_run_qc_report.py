import json
import os
import stat
import sys

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import build_run_qc_report as qc

pa = pytest.importorskip("pyarrow")
import pyarrow.parquet as pq

TRACE_HEADER = "task_id\thash\tnative_id\tprocess\ttag\tname\tstatus\texit\tattempt\tsubmit\tstart\tcomplete\tduration\trealtime\tcpus\t%cpu\tmemory\tpeak_rss\tpeak_vmem\thostname\tworkdir\n"


def trace_line(name, status, exit_code, attempt, realtime, memory, peak, host, workdir):
    return "\t".join(["1", "ab/cdef", "1", name.split(" ")[0], "r1", name, status, exit_code, str(attempt), "", "", "", realtime,
                      realtime, "1", "100%", memory, peak, "-", host, workdir]) + "\n"


def fake_bcftools(tmp_path, count=5):
    path = tmp_path / "bcftools"
    path.write_text(f"#!/bin/sh\necho {count}\n")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return str(path)


def make_run(tmp_path, layout="subfolders", genotype_mode="multi_sample", with_trace=True):
    run = tmp_path / "run1"
    d = {k: run / k for k in ("annotation", "predictors", "tables", "genotypes", "reports", "pipeline_info")} \
        if layout == "subfolders" else dict.fromkeys(("annotation", "predictors", "tables", "genotypes", "reports", "pipeline_info"), run)
    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)
    (d["reports"] / "run1.assembly_check.tsv").write_text(
        f"key\tvalue\nstatus\tPASS\nrecords_scanned\t5\ngenotype_mode\t{genotype_mode}\nsample_names\tA,B\n")
    for tag in ("annVEP", "annSPiP"):
        (d["predictors"] / f"run1.{tag}.vcf.gz").write_bytes(b"x")
    (d["annotation"] / "run1.annotated.vcf.gz").write_bytes(b"x")
    pq.write_table(pa.table({"Locus": ["1:1", "1:1", "1:2"], "SPiP_status": ["scored", "not_covered", "scored"],
                             "TRANSCRIPT_PRIORITY_TIER": ["1", "2", "3"]}), d["tables"] / "run1.parsed.clean.pq")
    pq.write_table(pa.table({"x": [1, 2, 3, 4, 5, 6]}), d["tables"] / "run1.parsed.clean.withGT.pq")
    if genotype_mode == "multi_sample":
        pq.write_table(pa.table({"x": [1] * 6}), d["genotypes"] / "run1.genotypes.pq")
    (d["pipeline_info"] / "RUN_MANIFEST.json").write_text(json.dumps({
        "launch": {"generated_at": "2026-01-01T00:00:00+00:00", "cli_flags": {"run_id": "run1"},
                   "annotation_pipeline_repo": {"commit": "abc", "dirty": False}, "spliceai_distance": 4999},
        "completed": {"generated_at": "2026-01-01T03:00:00+00:00", "exit_code": 0}}))
    if with_trace:
        wd = tmp_path / "wd_fail"
        wd.mkdir()
        (wd / ".command.err").write_text("some noise\nMemoryError\n")
        (d["pipeline_info"] / "execution_trace_2026-01-01_00-00-00.txt").write_text(
            TRACE_HEADER
            + trace_line("WF:TAG_TRANSCRIPT_PRIORITY (run1)", "FAILED", "137", 1, "30s", "12 GB", "11 GB", "n1", str(wd))
            + trace_line("WF:TAG_TRANSCRIPT_PRIORITY (run1)", "COMPLETED", "0", 2, "1h 30m", "24 GB", "20 GB", "n2", str(tmp_path))
            + trace_line("WF:VEP_ANNOTATE (run1)", "COMPLETED", "0", 1, "343ms", "8 GB", "-", "n1", str(tmp_path)))
    return run


def test_unit_parsers():
    assert qc.to_ms("343ms") == 343 and qc.to_ms("1d 2h 3m 4.5s") == 93784500 and qc.to_ms("-") is None and qc.to_ms("1500") == 1500
    assert qc.to_bytes("1.5 GB") == int(1.5 * 1024 ** 3) and qc.to_bytes("-") is None and qc.to_bytes("2048") == 2048
    assert qc.process_of("A:B:PROC (tag x)") == "PROC"


def test_full_report_from_subfolder_layout(tmp_path):
    rep = qc.build(make_run(tmp_path), None, fake_bcftools(tmp_path, 5), True)
    assert rep["summary"]["wall_hours"] == 3.0 and rep["summary"]["commit"] == "abc" and rep["summary"]["genotype_mode"] == "multi_sample"
    assert rep["stage_counts"]["input records"] == 5 and rep["stage_counts"]["annVEP records"] == 5
    assert rep["stage_counts"]["annPangolin records"] is None  # file absent: unavailable, not 0
    assert rep["table"]["rows"] == 3 and rep["table"]["distinct_locus"] == 2
    assert rep["table"]["status"]["SPiP_status"] == {"scored": 2, "not_covered": 1}
    assert rep["table"]["tier"] == {"1": 1, "2": 1, "3": 1}
    assert rep["genotype_rows"] == 6 and rep["genotype_note"] is None
    tag = next(p for p in rep["processes"] if p["process"] == "TAG_TRANSCRIPT_PRIORITY")
    assert tag["tasks_ok"] == 1 and tag["failed_attempts"] == 1 and tag["total_h"] == 1.5 and tag["peak_rss_gb_max"] == 20.0
    vep = next(p for p in rep["processes"] if p["process"] == "VEP_ANNOTATE")
    assert vep["peak_rss_gb_max"] is None and vep["nodes"] == "n1"
    assert len(rep["failed_attempts"]) == 1 and rep["failed_attempts"][0]["last_error"] == "MemoryError"


def test_flat_layout_and_missing_sources(tmp_path):
    rep = qc.build(make_run(tmp_path, layout="flat", genotype_mode="no_samples", with_trace=False), None, fake_bcftools(tmp_path), True)
    assert rep["table"]["rows"] == 3 and rep["stage_counts"]["annotated (merged) records"] == 5
    assert rep["genotype_rows"] is None and "not applicable" in rep["genotype_note"]
    assert rep["processes"] == [] and rep["task_source"] is None
    assert any("no execution trace" in n for n in rep["notes"])


def test_no_manifest_is_reported_unavailable(tmp_path):
    run = make_run(tmp_path)
    (run / "pipeline_info" / "RUN_MANIFEST.json").unlink()
    rep = qc.build(run, None, fake_bcftools(tmp_path), False)
    assert rep["summary"]["commit"] is None and rep["summary"]["wall_hours"] is None
    assert rep["table"]["distinct_locus"] is None
    assert any("no RUN_MANIFEST.json" in n for n in rep["notes"])


def test_log_fallback_uses_workdir_times_and_counts_retries(tmp_path):
    def task_dir(name, seconds, exit_file=True, node="nodeX"):
        d = tmp_path / name
        d.mkdir()
        (d / ".command.log").write_text(f"header\n{node}\n")
        (d / ".command.begin").write_text("")
        os.utime(d / ".command.begin", (1000, 1000))
        end = d / (".exitcode" if exit_file else ".command.log")
        if exit_file:
            end.write_text("0")
        os.utime(end, (1000 + seconds, 1000 + seconds))
        return d
    ok, killed = task_dir("ok", 7200), task_dir("killed", 3600, exit_file=False)
    line = "Oct-02 14:28:48.875 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobId: 1; id: 1; name: WF:P (r1); status: COMPLETED; exit: {}; error: -; workDir: {} started: 1; exited: x; ]\n"
    log = tmp_path / "nf.log"
    log.write_text("unrelated line\n" + line.format("-", killed) + line.format("0", ok))
    rows = qc.read_log(log)
    assert [r["attempt"] for r in rows] == [1, 2] and rows[0]["duration_ms"] == 3600000 and rows[0]["host"] == "nodeX"
    table = qc.process_table(rows)[0]
    assert table["tasks_ok"] == 1 and table["failed_attempts"] == 1 and table["median_h"] == 2.0 and table["failed_attempts_h"] == 1.0


def test_reports_are_never_overwritten(tmp_path):
    rep = qc.build(make_run(tmp_path), None, fake_bcftools(tmp_path), False)
    files = qc.write_reports(rep, tmp_path / "out")
    assert len(files) == 3 and all(f.exists() for f in files)
    assert "unavailable" in files[0].read_text()  # HTML renders missing values as such
    with pytest.raises(FileExistsError):
        qc.write_reports(rep, tmp_path / "out")  # same second: same names
