import click
import os
import shutil
from subprocess import call

def parse_manifest(manifest):
    samples = []
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
@click.option('--keep-temp-files', is_flag=True, help='Keep temporary files?')
@click.argument('manifest')
#def main(manifest, cluster):
def main(manifest, out, bwa_index, keep_temp_files):
    """A preprocessing and QC pipeline for ChIA-PET data."""
    script_dir = os.path.dirname(os.path.realpath(__file__))
    out = os.path.abspath(out)    
    bwa_index = os.path.abspath(bwa_index)    
    click.echo("Output folder: %s" % out) 
    click.echo("BWA index: %s\n" % bwa_index)     
    # Preprocess individual samples
    samples = parse_manifest("samples.txt")
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
    cmd = [preproc_set, out] + [os.path.join(out, 'samples', x['name']) for x in samples]
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

