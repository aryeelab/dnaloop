import pytest
from click.testing import CliRunner
from dnaloop import cli
import md5

def file_checksums_equal(file1, file2):
    with open(file1) as f:
        checksum1 = md5.new(f.read()).digest()
    with open(file2) as f:
        checksum2 = md5.new(f.read()).digest()
    return checksum1==checksum2    

#@pytest.fixture
  
def test_preproc_run():
    runner = CliRunner()
    result = runner.invoke(cli.main, ['--out', 'naive_vs_primed', '--bwa-index', 'test_genome.fa', '--keep-temp-files', 'samples.txt'])
    assert not result.exception
    assert result.exit_code == 0

def test_peaks():
    assert file_checksums_equal('correct_output/peaks/anchor_peaks.narrowPeak', 'naive_vs_primed/peaks/anchor_peaks.narrowPeak')
    assert file_checksums_equal('correct_output/peaks/anchor_peaks.merged.bed', 'naive_vs_primed/peaks/anchor_peaks.merged.bed')
    
def test_loop_counts():
    assert file_checksums_equal('correct_output/naive_esc_1.loop_counts.bedpe', 'naive_vs_primed/naive_esc_1.loop_counts.bedpe')

def test_preproc_run_mergegap():
    runner = CliRunner()
    result = runner.invoke(cli.main, ['--out', 'naive_vs_primed', '--bwa-index', 'test_genome.fa', '--merge-gap', '1000', 'samples.txt'])
    assert not result.exception
    assert result.exit_code == 0

def test_parse_yaml_manifest():
    samples = [ {   'name': 'primed_esc',
                    'read1': 'fastq/primed_esc_1.r1.fastq.gz,fastq/primed_esc_2.r1.fastq.gz', 
                    'read2': 'fastq/primed_esc_1.r2.fastq.gz,fastq/primed_esc_2.r2.fastq.gz'

                }, 
                {   'name': 'naive_esc',
                    'read1': 'fastq/naive_esc_1.r1.fastq.gz,fastq/naive_esc_2.r1.fastq.gz', 
                    'read2': 'fastq/naive_esc_1.r2.fastq.gz,fastq/naive_esc_2.r2.fastq.gz'
                }
               ]
    assert cli.parse_manifest('samples.yaml') == samples
