import click
import os
import shutil
import yaml
import shutil
import random
import string
from pkg_resources import get_distribution
from subprocess import call, check_call
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

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

#Grab Sample Names
def get_subdirectories(dir):
    return [name for name in os.listdir(dir)
            if os.path.isdir(os.path.join(dir, name))]

def qc_report(dir):
    #Read in the summary statistics            
    statslist = []
    sampleNames = get_subdirectories(os.path.join(dir, "samples"))
    for s in sampleNames:
        sfile = open(os.path.join(dir, "samples/" + s + "/read_stats.txt"), 'r')
        stats = []
        header = sfile.readline()
        for l in sfile.readlines():
            stats.append(int(l.split(": ")[1].strip()))	
        statslist.append(stats)
    
    readstats = pd.DataFrame(statslist)
    readstats = readstats.transpose()
    readstats.columns = sampleNames
    
    # Now make the plots
    with PdfPages(os.path.join(dir, 'qc-report.pdf')) as pdf:    
        plt.figure(figsize=(8, 11))
        readstats.ix[0].plot(kind='bar')
        plt.title('Total Reads Per Sample')
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()


@click.command()
#@click.option('--cluster-command', '-c', default='bsub', help='Cluster submit command')
@click.option('--out', default=".", required=True, help='Output directory name')
@click.option('--bwa-index', required=True, help='BWA index location')
@click.option('--peak-pad', default="0", help='Peak padding width (applied on both left and right)')
@click.option('--merge-gap', default="1500", help='Max gap size for merging peaks')
@click.option('--use-lsf', is_flag=True, help='Submit jobs to an LSF cluster?')
@click.option('--bsub-opts', default="", help='LSF bsub options')
@click.option('--keep-temp-files', is_flag=True, help='Keep temporary files?')
@click.option('--no-qc-report', is_flag=True, help='Skip QC report generation? (Requires R)')
@click.argument('manifest')
#def main(manifest, cluster):
def main(manifest, out, bwa_index, peak_pad, merge_gap, use_lsf, bsub_opts, keep_temp_files, no_qc_report):
    """A preprocessing and QC pipeline for ChIA-PET data."""
    __version__ = get_distribution('dnaloop').version
    click.echo("Starting dnaloop pipeline v%s" % __version__)
    if os.path.exists(out):
        shutil.rmtree(out)
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
    if use_lsf:
        job_ids = []
        lsf_id = ''.join(random.SystemRandom().choice(string.ascii_lowercase + string.digits) for _ in range(5))
    for sample in samples:
        i += 1
        click.echo("\nProcessing sample %d of %d: %s" % (i, len(samples), sample['name']))
        click.echo("    Read 1: %s" % sample['read1']) 
        click.echo("    Read 2: %s" % sample['read2'])    
        preproc_fastq = os.path.join(script_dir, 'preprocess_chiapet_fastq.sh')
        cmd = [preproc_fastq, os.path.join(out, 'samples', sample['name']), bwa_index, merge_gap, sample['read1'], sample['read2']]
        if use_lsf:
            job_id = 'dnaloop_sample_%s_%d' % (lsf_id, i)
            cmd = "bsub -J %s %s %s" % (job_id, bsub_opts,  " ".join(cmd))
            click.echo("    Submitting job to LSF: %s" % cmd)    
            check_call(cmd, shell=True)
            job_ids.append(job_id)        
        else:
            click.echo("    Executing: %s" % " ".join(cmd))
            call(cmd)
    if use_lsf:
        conditions = ["ended(%s)" % job_id for job_id in job_ids]
        depend_cond = " && ".join(conditions)
        cmd = 'bsub -K -o /dev/null -q long -w "%s" exit 0' % depend_cond
        print("Waiting for %d dnaloop sample preprocessing jobs to finish" % len(job_ids))
        check_call(cmd, shell=True)        
    # Create the ChIA-PET analysis set
    preproc_set = os.path.join(script_dir, 'preprocess_chiapet_set.sh')
    cmd = [preproc_set, out, peak_pad, merge_gap] + [os.path.join(out, 'samples', x['name']) for x in samples]
    click.echo("Creating ChIA-PET set")
    click.echo("    Executing: %s\n" % " ".join(cmd))
    call(cmd)
    if no_qc_report:
        click.echo("Skipping QC report since --no-qc-report was specified")
    else:
        click.echo("Creating QC report")
        cmd = ['Rscript', os.path.join(script_dir, 'qc-report.R'), out] 
        call(cmd)
    if keep_temp_files:
        click.echo("Temporary files not deleted since --keep-temp-files was specified")
    else:
        click.echo("Deleting temporary files")
        shutil.rmtree(os.path.join(out, 'peaks'))
        shutil.rmtree(os.path.join(out, 'samples'))
    click.echo("Done")


