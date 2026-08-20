import sys
import os
import pytest
import pysam
import tempfile
import shutil

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from split_vcf_chunks import split_vcf_into_chunks


def create_mock_vcf(vcf_path: str, num_records: int = 50):
    header = pysam.VariantHeader()
    header.add_line('##fileformat=VCFv4.2')
    header.add_line('##INFO=<ID=TEST_INFO,Number=1,Type=String,Description="Test Info">')
    header.add_line('##contig=<ID=chr1,length=248956422>')
    
    with pysam.VariantFile(vcf_path, "w", header=header) as out_vcf:
        for i in range(1, num_records + 1):
            rec = out_vcf.new_record(
                contig="chr1",
                start=1000 + (i * 10),
                stop=1001 + (i * 10),
                alleles=("A", "G"),
                info={"TEST_INFO": f"val_{i}"}
            )
            out_vcf.write(rec)


def test_split_vcf_into_chunks():
    temp_dir = tempfile.mkdtemp()
    try:
        mock_vcf = os.path.join(temp_dir, "input.vcf")
        create_mock_vcf(mock_vcf, num_records=45)
        
        chunks_dir = os.path.join(temp_dir, "chunks")
        chunk_files = split_vcf_into_chunks(mock_vcf, chunks_dir, chunk_size=10)
        
        assert len(chunk_files) == 5
        assert os.path.exists(chunk_files[0])
        assert os.path.exists(chunk_files[0] + ".tbi")
        
        # Verify chunk record counts
        records_per_chunk = []
        total_records = 0
        for cf in chunk_files:
            with pysam.VariantFile(cf) as vf:
                recs = list(vf)
                records_per_chunk.append(len(recs))
                total_records += len(recs)
                
        assert records_per_chunk == [10, 10, 10, 10, 5]
        assert total_records == 45
    finally:
        shutil.rmtree(temp_dir)


def test_split_vcf_header_fidelity():
    temp_dir = tempfile.mkdtemp()
    try:
        mock_vcf = os.path.join(temp_dir, "input.vcf")
        create_mock_vcf(mock_vcf, num_records=15)
        
        chunks_dir = os.path.join(temp_dir, "chunks")
        chunk_files = split_vcf_into_chunks(mock_vcf, chunks_dir, chunk_size=10)
        
        with pysam.VariantFile(chunk_files[0]) as vf:
            assert "TEST_INFO" in vf.header.info
            assert vf.header.contigs["chr1"].length == 248956422
    finally:
        shutil.rmtree(temp_dir)
