import click
import os
import shutil
import yaml
from pkg_resources import get_distribution
from subprocess import call

def parse_manifest(manifest):
    samples = []
    if manifest.endswith(('.yaml', '.yml')):
        with open(manifest, 'r') as f: 
            m = yaml.load(f)
        sample_names = m['samples'].keys()
        sample_names.sort()
        for sample_name in sample_names:
            runs = m['samples'][sample_name]
            read1 = []
            read2 = []
            for run in runs:
                fastq1, fastq2 = run.split(" ")
                read1.append(fastq1)
                read2.append(fastq2)
            d = {
                'name': sample_name, 
                'read1': ','.join(read1),
                'read2': ','.join(read2)
                }
            samples.append(d)
        return samples
    else:
        with open(manifest) as f:
            lines = f.readlines()
        for line in lines:
            fields = line.strip().split("\t")
            if len(fields)==3:
                d = {   'name': fields[0], 
                        'read1': fields[1],
                        'read2': fields[2]
                        }
                samples.append(d)
            else:
                if line != "\n":
                    click.echo ("Skipping line: " + line)
        return samples

@click.command()
#@click.option('--cluster-command', '-c', default='bsub', help='Cluster submit command')
@click.option('--out', default=".", required=True, help='Output directory name')
@click.option('--bwa-index', required=True, help='BWA index location')
@click.option('--merge-gap', default="1500", help='Max gap size for merging peaks')
@click.option('--keep-temp-files', is_flag=True, help='Keep temporary files?')
@click.argument('manifest')
#def main(manifest, cluster):
def main(manifest, out, bwa_index, merge_gap, keep_temp_files):
    """A preprocessing and QC pipeline for ChIA-PET data."""
    __version__ = get_distribution('dnaloop').version
    click.echo("Starting dnaloop pipeline v%s" % __version__)
    os.mkdir(out)
    os.mkdir(os.path.join(out, 'log'))
    with open(os.path.join(out, 'log', 'VERSION.txt'), 'w') as f: 
        f.write(__version__ + '\n')
    script_dir = os.path.dirname(os.path.realpath(__file__))
    out = os.path.abspath(out)    
    bwa_index = os.path.abspath(bwa_index)    
    click.echo("Output folder: %s" % out) 
    click.echo("BWA index: %s\n" % bwa_index)     
    # Preprocess individual samples
    samples = parse_manifest(manifest)
    i = 0
    for sample in samples:
        i += 1
        click.echo("\nProcessing sample %d of %d: %s" % (i, len(samples), sample['name']))
        click.echo("    Read 1: %s" % sample['read1']) 
        click.echo("    Read 2: %s" % sample['read2'])    
        preproc_fastq = os.path.join(script_dir, 'preprocess_chiapet_fastq.sh')
        cmd = [preproc_fastq, os.path.join(out, 'samples', sample['name']), bwa_index, sample['read1'], sample['read2']]
        click.echo("    Executing: %s" % " ".join(cmd))
        call(cmd)
    # Create the ChIA-PET analysis set
    preproc_set = os.path.join(script_dir, 'preprocess_chiapet_set.sh')
    cmd = [preproc_set, out, merge_gap] + [os.path.join(out, 'samples', x['name']) for x in samples]
    click.echo("Creating ChIA-PET set")
    click.echo("    Executing: %s\n" % " ".join(cmd))
    call(cmd)
    if keep_temp_files:
        click.echo("Temporary files not deleted since --keep-temp-files was specified")
    else:
        click.echo("Deleting temporary files")
        shutil.rmtree(os.path.join(out, 'peaks'))
        shutil.rmtree(os.path.join(out, 'samples'))
    click.echo("Done")

