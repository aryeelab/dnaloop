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
