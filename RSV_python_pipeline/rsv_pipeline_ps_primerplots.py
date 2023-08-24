#!/usr/bin/python3

import pickle
import gzip
from Bio import SeqIO
from pathlib import Path
import os
base_path = Path(os.path.abspath(os.path.curdir))


### set up logging
import os
import sys
sys.path.insert(0, '../AdV_wgs2')


import logbook
from logbook.more import ColorizedStderrHandler
def get_console_logger(level = 10, format_string = None, name = 'main'):
    '''
    Returns an application logbook logger with simple format
    
    Levels
    ------
    10 'DEBUG'
    11 'INFO'
    12 'NOTICE'
    13 'WARNING'
    14 'ERROR'
    15 'CRITICAL'
    '''
    # prevents +/- 1 hour offset issue that logging has?
    logbook.set_datetime_format('local')
    handler = ColorizedStderrHandler(
            bubble = False,
            level = level,
    )
    if not format_string:
        handler.format_string = \
                '{record.time:%Y-%m-%d %H:%M} {record.level_name}: {record.message}'
        # NB more conventional to have channel included (name to logbook.Logger):
        # '{record.time:%Y-%m-%d %H:%M} {record.level_name} {record.channel}: {record.message}'
    else:
        handler.format_string = format_string
    handler.push_application()
    return logbook.Logger(name, level = level)

log = get_console_logger(10)

log.debug('debug')       # grey
log.info('info')         # grey
log.notice('notice')     # yellow
log.warning('warning')   # yellow
log.error('error')       # red
log.critical('critical') # red
######################################################################################
######################################################################################
######################################################################################

### 0.4 select samples and references

import pathlib
from glob import glob
from nisvirology import illumina
import os
import subprocess
import pysam
import pandas as pd
import ast

run_folder = '230426_NS552083_0414_AH37WWAFX5'

base_path = Path(os.path.abspath(os.path.curdir))

base_data_path = Path(
        os.path.expanduser('~') + (
            f'/AdV_wgs2/{run_folder}/{run_folder}'
        )
)

sample_fastqs = {}
for f in sorted(glob(f'{base_data_path}/*fastq.gz')):
    path = pathlib.Path(f)
    try:
        sample_fastqs[path.name.split('.')[0]] += [path]
    except KeyError:
        sample_fastqs[path.name.split('.')[0]] = [path]

sample_fastqs = {
        sample_name:(R1, R2)
        for sample_name,(R1, R2) in sample_fastqs.items()
        #if 'VTM' not in sample_name # or 'NEG' in sample_name
        #if sample_name in sample_selection
}

len(sample_fastqs) == 21

targets = {
    #'Adenovirus F41 MU22': 'MG925782.1',
    'England174820313': 'LR699737',
    'England174880237': 'MZ515970',
    #'England144960079': 'KU950564',
    #'England144960102': 'KM517573',
    #'England145060047': 'JX576736',
    #'England154640011': 'KU950540',
    #'England155000459': 'KU950523',
    
}

use_genomes = targets
use_genomes = {
        #'Adenovirus F41': 'MG925782.1',
        #'Adenovirus F40': 'NC_001454.1',
        'England174820313': 'LR699737',
        'England174880237': 'MZ515970',
        #'England144960079': 'KU950564',
        #'England144960102': 'KM517573',
        #'England145060047': 'JX576736',
        #'England154640011': 'KU950540',
        #'England155000459': 'KU950523',
}

sorted(sample_fastqs)



sample_to_ref_combos = {
    
        
        '01121231_313_1_a-1': {'England174820313'},
        '01121232_313_2_397_4_a-1': {'England174820313'},
        '01121233_313_3_397_3_a-1': {'England174820313'},
        '01121234_313_4_397_2_a-1': {'England174820313'},
        '01121235_397_1_a-1': {'England174820313'},
        '01121236_237_2_b-1': {'England174880237'},
        '01121237_237_3_410_5_b-1': {'England174880237'},
        '01121238_237_4_410_4_b-1': {'England174880237'},
        '01121239_237_5_410_3_b-1': {'England174880237'},
        '01121240_410_2_b-1': {'England174880237'},
       
       
}
#print(sample_fastqs)

#################################################################################################
#################################################################################################
#################### 0.45 get reference sequences

import nisvirology

## Or reload if already downloaded
seq_result = {}
for ref_genome_type,ref_genome_acc in targets.items():
    seq_result[ref_genome_acc] = SeqIO.read(
            f'{base_path}/{ref_genome_acc}.gbk',
            'genbank'
    )

###############################################################################################
############################################################################################


### 0.5 map primers

import os
import toml
from Bio import SeqIO
from Bio.Seq import Seq

import importlib

#importlib.reload(nisvirology)
from nisvirology import align
#importlib.reload(align)
from nisvirology.align import (
        SinglePCRAssay,
        MappedOrientationError,
        DidNotMapError,
)
from collections import defaultdict
from Bio.Seq import Seq

parameter_file = (
        os.path.expanduser('~') + (
        f'/AdV_wgs2/{run_folder}/RSV Orthopneumovirus.toml'
        )
)

log.info(f"Parsing '{parameter_file}'")
try:
    parameters = toml.load(parameter_file)
except toml.TomlDecodeError:
    log.exception(
            f"Duplicate keys were found while parsing {parameter_file}. "
            "This can be caused by higher and lower ranks in the same "
            "lineage both having IDs assigned in the taxonomy section. "
            "Else check the error message from the parser.")
    sys.exit(1)

pcr_assays = {}
for ref_genome,ref_genome_acc in targets.items():
    ref_genome_path = f'{base_path}/{ref_genome_acc}.fna'
    ref_genome_rec = SeqIO.read(ref_genome_path, 'fasta')
    pcr_assays[ref_genome] = {}
    for assay_name,assay_info in parameters['Assays'].items():
        if not (assay_name.startswith('HRSVA1')
               or
               'HRSVA2' in assay_name
           #     or
          #     'HRSVA3' in assay_name
                or
               'HRSVB1' in assay_name
                or
               'HRSVB2' in assay_name
                
        ):
            continue
        log.notice(f"{ref_genome} vs {assay_name}")
        if 'multiplex' in assay_info:
            for r_name,reaction in assay_info['multiplex'].items():
                fwd_name, fwd_seq = reaction['forward']
                rev_name, rev_seq = reaction['reverse']
                if assay_info.get('reverse_strand_same_as_forward'):
                    rev_seq = str(Seq(rev_seq).reverse_complement())
                ## NB this is inconsistent in toml
                if 'description' in assay_info:
                    description = assay_info['description']
                    sop_id = assay_info['sop_id']
                else:
                    description = reaction['description']
                    sop_id = reaction['sop_id']
                try:
                    pcr_assay = SinglePCRAssay(
                            f'{assay_name}, {r_name}',
                            ref_genome_rec,
                            forward_primer = (fwd_name, fwd_seq),
                            reverse_primer = (rev_name, rev_seq),
                            description = description,
                            sop_id = sop_id,
                            min_pID = 0.75,
                            check_slices = False,
                            allow_partial_end_mapping = True,
                    )
                except MappedOrientationError as exc:
                    log.warning(exc)
                    raise
                except DidNotMapError as exc:
                    continue
                if not pcr_assay.fwd_primer.mapping.end < \
                        pcr_assay.rev_primer.mapping.start:
                    ## NB until actual pairings provided, just omit
                    ## apparent mispairings
                    # "Primer pair not provided in conventional way."
                    # should probably raise a MappedOrientationError?
                    continue
                pcr_assays[ref_genome][assay_name,r_name] = pcr_assay

        else:
            fwd_name, fwd_seq = assay_info['forward']
            rev_name, rev_seq = assay_info['reverse']
            if assay_info.get('reverse_strand_same_as_forward'):
                rev_seq = str(Seq(rev_seq).reverse_complement())
            try:
                pcr_assays[ref_genome][assay_name] = SinglePCRAssay(
                        assay_name,
                        ref_genome_rec,
                        forward_primer = (fwd_name, fwd_seq),
                        reverse_primer = (rev_name, rev_seq),
                        description = assay_info['description'],
                        sop_id = assay_info['sop_id'],
                        min_pID = 0.75,
                        check_slices = False,
                        allow_partial_end_mapping = True,
                )
            except MappedOrientationError as exc:
                log.warning(exc)
                raise
            except DidNotMapError as exc:
                continue

#pcr_assay = pcr_assays['Adenovirus F41 Tak']['Aden-groupF-only']
pcr_assay.fwd_primer.mapping
# Out[428]: Mapping(30370, 30390, 1, target_id = MG925782.1, seq_identity = 1.0)
pcr_assay.rev_primer.mapping
# Out[429]: Mapping(30502, 30522, -1, target_id = MG925782.1, seq_identity = 1.0)

pcr_assay.fwd_primer.mapping.extract(pcr_assay.target)
# Out[432]: Seq('CACTTAATGCTGACACGGGC')
pcr_assay.fwd_primer.sequence
# NB no annotations
pcr_assay.target.features

## some OK some not . . . means mode of presentation is OK but some sequences not?
## need at least 'own' reference to work against all
## these need checking else errors in sequences <== TODO

###################################################################################
###################################################################################
###################################################################################

### 0.6 amplicon ranges

## sort into reactions for multiplexed reactions
## no explicit info in toml schema: parse ids for now?
# pcr_assays[ref_genome][assay_name]

amplicon_ranges = {}
for ref_genome,assays in pcr_assays.items():
    amplicon_ranges[ref_genome] = {}
    for assay_name,assay in assays.items():
        # NB multiplexed reactions have multiple regions (amplicons)
        # encoded in a 2-tuple key here - checked.
        if isinstance(assay_name, tuple):
            use_assay_name,region = assay_name
        else:
            use_assay_name,region = assay_name,'region1'
        print(use_assay_name, region)
        if use_assay_name not in amplicon_ranges[ref_genome]:
            amplicon_ranges[ref_genome][use_assay_name] = {}

        amplicon_start = assay.fwd_primer.mapping.start
        amplicon_end = assay.rev_primer.mapping.end
        assert amplicon_start < amplicon_end
        amplicon_id = assay.fwd_primer.name, assay.rev_primer.name
        amplicon_ranges[ref_genome][use_assay_name][region, amplicon_id] = \
                amplicon_start, amplicon_end


##################################################################################
##################################################################################
##################################################################################


        
### 1. map reads


import pathlib
from glob import glob
from nisvirology import illumina
import os
import subprocess
import pysam

picard_path = '/home/tiina/mambaforge/share/picard-2.27.4-0/picard.jar'

# bwa mem -M -R @RG\tID:reads\tSM:reads\tPL:ILLUMINA
# args_to_bwa : [('-L','10,15'),('-H',), ('T',4)]

## TODO test and improve this mapping short read code
## see nf core iVar pipeline for primer and other trimming etc.
## **are these reads trimmed**?

args_to_bwa = [['-T','4']]
print(args_to_bwa)
path_to_tmp = f'{base_path}/tmp'

#Path(
#        os.path.expanduser('~') + (
 #           '/AdV_wgs2/tmp'
  #      )
#)
print(path_to_tmp)
done = set()
#done.add((ref_genome,sample_name))

for ref_genome,ref_genome_acc in targets.items():
    #if ref_genome not in use_genomes:
     #   continue
    rec = seq_result[ref_genome_acc]
    ref_genome_path = f'{base_path}/{ref_genome_acc}.fna'
    ref_gnm_fname = ref_genome.replace(' ','_')
    for sample_name,(R1, R2) in sorted(sample_fastqs.items()):
        print(sample_name)
        if ref_genome not in sample_to_ref_combos[sample_name]:
            continue
            log.info(f"Skipping {ref_genome} vs {sample_name}")
        log.info(f"Selecting {ref_genome} vs {sample_name}")

        if (ref_genome,sample_name) in done:
            log.info(
                    f"Skipping {ref_genome} vs {sample_name} - "
                    "already done"
            )
            continue
        log.info(f"Doing {ref_genome} vs {sample_name}")
        indir = R1.parent
        print(indir)
        # mapping using BWA MEM (VirVarSeq used BWA sampe)
        sample_outdir = f'{indir}/results/{sample_name}'
        #print(sample_outdir)
        ref_sam_path = Path(
                f'{sample_outdir}/'
                f'{sample_name}_v_{ref_gnm_fname}_{ref_genome_acc}.sam'
        )
        ref_bam_path = Path(str(ref_sam_path)[:-3]+"bam")
        bam_path = Path(str(ref_sam_path)[:-3]+"sorted.bam")
        bam_path_dd = Path(str(bam_path)[:-4]+'_dd.bam')
        if bam_path_dd.exists():
            log.info(f"Already done: {bam_path_dd}")
            try:
                ref_sam_path.unlink()
                log.info(f"Deleted {ref_sam_path}")
            except FileNotFoundError:
                log.info(f"Already deleted {ref_sam_path}")
            try:
                ref_bam_path.unlink()
                log.info(f"Deleted {ref_bam_path}")
            except FileNotFoundError:
                log.info(f"Already deleted {ref_bam_path}")
            try:
                bam_path.unlink()
                log.info(f"Deleted {bam_path}")
            except FileNotFoundError:
                log.info(f"Already deleted {bam_path}")
            continue
        os.makedirs(sample_outdir, exist_ok = True)
        if not ref_sam_path.exists():
            # include read group information for e.g., GATK
            # just a place holder "reads" because one group
            # also use -M to mark split hits as secondary for picard compatibility
            # a single arbitrary read group is hard-wired here.
            RGinfo = r"@RG\tID:reads\tSM:reads\tPL:ILLUMINA"
            cmd = ['bwa', 'mem', '-M']
            if args_to_bwa:
                for arg in args_to_bwa:
                    cmd += list(arg)
            cmd += ['-R', RGinfo, ref_genome_path, str(R1), str(R2)]
            print(cmd)
            log.info('Will call "{}"'.format(' '.join(cmd)))
            with open(ref_sam_path, 'w') as sam:
                proc = subprocess.Popen(
                        cmd,
                        stdout = sam,
                        stderr = subprocess.PIPE
                )
                o,bwa_stderr = proc.communicate()
                print(bwa_stderr)
                bwa_stderr = bwa_stderr.decode('utf-8')
                print(bwa_stderr)
            log.info('Processing alignments to reference')

        if not ref_bam_path.exists():
            pysam.view("-b", "-o", str(ref_bam_path), str(ref_sam_path),
                    catch_stdout = False)

        if not bam_path.exists():
            pysam.sort(
                    "-T", path_to_tmp+"/aln.sorted",
                    "-o", str(bam_path),
                    str(ref_bam_path),
            )
            pysam.index(str(bam_path))
        # remove PCR duplicates using Picard
        if not bam_path_dd.exists():
            picard_command = [
                    'MarkDuplicates',
                    'I=', bam_path,
                    'O=', bam_path_dd,
                    'M=', str(bam_path_dd)[:-3]+'log',
            ] #, 'VALIDATION_STRINGENCY=','LENIENT']
            cmd = ['java', '-Xmx4g', '-jar', picard_path] + picard_command
            log.info('Will call "{}"'.format(' '.join(map(str,cmd))))
            subprocess.call(cmd)
            pysam.index(str(bam_path_dd))
        done.add((ref_genome,sample_name))
        try:
            ref_sam_path.unlink()
            log.info(f"Deleted {ref_sam_path}")
        except FileNotFoundError:
            log.info(f"Already deleted {ref_sam_path}")
        try:
            ref_bam_path.unlink()
            log.info(f"Deleted {ref_bam_path}")
        except FileNotFoundError:
            log.info(f"Already deleted {ref_bam_path}")
        try:
            bam_path.unlink()
            log.info(f"Deleted {bam_path}")
        except FileNotFoundError:
            log.info(f"Already deleted {bam_path}")

##############################################################################
##############################################################################
##############################################################################

### 1.2 iVar filter mapped reads for primer sequences

# https://andersen-lab.github.io/ivar/html/manualpage.html

# ivar trim -b test_primers.bed -p test.trimmed -i test.bam -q 15 -m 50 -s 4

'''
in terms of primer trimming I would personally start by using iVar trim on an aligned, sorted BAM file with a command such as: 
"ivar trim -e -i {mapped and sorted bam} -b {bedfile} -m {minimum expected fragment length} -q {minimum illumina quality threshold} -p {output prefix}"


Which all being well should yield a version of the bam where any primers inside your reads should be softclipped
meaning that the positions will still be present in the BAM file but the cigar string for said read will have an "S" for
any primer positions meaning that it will be excluded by any later pileup / variant calling you perform on that BAM file.
'''

samples = {
        p.split('/')[7]:Path(p)
        for p in sorted(glob(f'{base_data_path}/results/*/'))
}

#############################################################################
#############################################################################
#############################################################################


## 1. generate the BED data

# first 100 bases of chromosome 1 are defined as chrom=1, chromStart=0, chromEnd=100
# i.e., last base is number 99 and 0:100 is Python slice-like.

# Required fields

 # chrom - name of the chromosome or scaffold. Any valid seq_region_name can be used, and chromosome names can be given with or without the 'chr' prefix.
 # chromStart - Start position of the feature in standard chromosomal coordinates (i.e. first base is 0).
 # chromEnd - End position of the feature in standard chromosomal coordinates
# Optional fields
# Nine additional fields are optional. Note that columns cannot be empty - lower-numbered fields must always be populated if higher-numbered ones are used.
 # name - Label to be displayed under the feature, if turned on in "Configure this page".
 # score - A score between 0 and 1000. See track lines, below, for ways to configure the display style of scored data.
 # strand - defined as + (forward) or - (reverse).


bed_data = {}
primerpair_data = {}
for ref_genome,assays in pcr_assays.items():
    bed_data[ref_genome] = {}
    primerpair_data[ref_genome] = {}
    for assay_name,assay in assays.items():
        print(ref_genome,assay_name)
        # NB multiplexed reactions have multiple regions (amplicons)
        # encoded in a 2-tuple key here - checked.
        if isinstance(assay_name, tuple):
            use_assay_name,region = assay_name
        else:
            use_assay_name,region = assay_name,'region1'
        if not use_assay_name in bed_data[ref_genome]:
            bed_data[ref_genome][use_assay_name] = []
        if not use_assay_name in primerpair_data[ref_genome]:
            primerpair_data[ref_genome][use_assay_name] = []

        # print(use_assay_name, region)
        # e.g., Puerto  28  52  400_1_out_L 60  +
        # 60 is apparently an arbitrary value for an unused 'score' column
        for oligo in assay.iter_oligos():
            line = "{}\t{}\t{}\t{}\t60\t{}".format(
                    assay.target.id,
                    oligo.mapping.start,
                    oligo.mapping.end,
                    oligo.name,
                    '+' if oligo.mapping.strand == 1 else '-',
            )
            print(line)
            bed_data[ref_genome][use_assay_name] += [line]
        primerpair_data[ref_genome][use_assay_name] += [
                f"{assay.fwd_primer.name}\t{assay.rev_primer.name}"
        ]

#################################################################################
#################################################################################
#################################################################################


## 2. Merge per PCR-based assay BED data according to Illumina sequencing pooling

# 220519_NS552079_0310_AHWWHGAFX2

## NEXT make v2 only BED and process

sorted(bed_data[ref_genome])
assay_names = set()
for ref_genome,assays in bed_data.items():
    use_lines = []
    use_lines_pp = []
    for assay_name, lines in assays.items():
        # if (
                # '2000jh_v1' in assay_name
                # or
                # '1200jh_v2' in assay_name
        # ):
            use_lines += lines
            assay_names.add(assay_name)
            use_lines_pp += primerpair_data[ref_genome][assay_name]
    ref_genome_acc = targets[ref_genome]
    ref_genome_nm = ref_genome.replace(' ','_')
    bed_path = Path(
            f'{ref_genome_nm}_{ref_genome_acc}_RSV.bed'
    )
    pp_path = Path(
            f'{ref_genome_nm}_{ref_genome_acc}_RSV.pp'
    )
    bed_path.open('w').write('\n'.join(use_lines))
    pp_path.open('w').write('\n'.join(use_lines_pp))
    
for sample_name,(R1, R2) in sorted(sample_fastqs.items()):
    indir = R1.parent
    break

print(line)

done = set()
for ref_genome,ref_genome_acc in targets.items():
    if ref_genome not in use_genomes:
        continue
    rec = seq_result[ref_genome_acc]
    ref_genome_path = f'{base_path}/{ref_genome_acc}.fna'
    ref_genome_nm = ref_genome.replace(' ','_')
    bed_path = Path(
            f'{ref_genome_nm}_{ref_genome_acc}_RSV.bed'
    )
    pp_path = Path(
            f'{ref_genome_nm}_{ref_genome_acc}_RSV.pp'
    )
    for sample_name,sample_outdir in samples.items():
        if (ref_genome,sample_name) in done:
            log.info(
                    f"Skipping {ref_genome} vs {sample_name} - "
                    "already done"
            )
            continue
        try:
            if ref_genome not in sample_to_ref_combos[sample_name]:
                continue
                log.info(f"Skipping {ref_genome} vs {sample_name}")
        except NameError:
            pass
        # log.info(f"Selecting {ref_genome} vs {sample_name}")
        # for dev:
        # if '_jhv2' in sample_name:
            # bed_path = Path(
                    # f'{ref_genome_nm}_{ref_genome_acc}_1200jh_v2.bed'
            # )
            # pp_path = Path(
                    # f'{ref_genome_nm}_{ref_genome_acc}_1200jh_v2.pp'
            # )
        # else:
            # bed_path = Path(
                    # f'{ref_genome_nm}_{ref_genome_acc}_2000jh_v1_1200jh_v2.bed'
            # )
            # pp_path = Path(
                    # f'{ref_genome_nm}_{ref_genome_acc}_2000jh_v1_1200jh_v2.pp'
            # )
        bam_path_dd = Path(
                f'{sample_outdir}/'
                f'{sample_name}_v_{ref_genome_nm}_{ref_genome_acc}.sorted_dd.bam'
        )
        if not bam_path_dd.exists():
            log.warning(f"Could not find {bam_path_dd}")
            continue
        log.info(f"Doing {ref_genome} vs {sample_name}")
        #raise Exception()
        expected_out_path = Path(
                f'{sample_outdir}/'
                f'ivar_trim__{sample_name}_v_{ref_genome_nm}_{ref_genome_acc}.sorted_dd.bam'
        )
        try:
            expected_out_path.unlink()
        except FileNotFoundError:
            pass
        if not expected_out_path.exists():
            cmd = [
                'ivar', 'trim',
                '-e', # include reads with no primer sequences
                '-i', str(bam_path_dd),
                '-b', str(bed_path),
                '-p', str(pp_path),
                '-m', '50', # minimum expected fragment length
                '-q', '15', # minimum illumina quality threshold
                '-p', f'{expected_out_path.parent}/{expected_out_path.stem}',
            ]
            log.info(cmd)
            subprocess.call(
                    cmd,
                    #env = {'LD_LIBRARY_PATH':  os.path.expanduser('~')+'/bamtools/lib/'},
            )
        bam_path_ivar_sorted = Path(str(expected_out_path)[:-4]+'_sorted2.bam')
        try:
            bam_path_ivar_sorted.unlink()
        except FileNotFoundError:
            pass
        if not bam_path_ivar_sorted.exists():
            pysam.sort(
                    "-T", path_to_tmp+"/aln.sorted",
                    "-o", str(bam_path_ivar_sorted),
                    str(expected_out_path),
            )
            pysam.index(str(bam_path_ivar_sorted))
        assert expected_out_path.exists()
        assert bam_path_ivar_sorted.exists()
        done.add((ref_genome,sample_name))

## Still no "Amplicons detected" despite primer pairs <--


###########################################################################################
###########################################################################################
###########################################################################################


### 1.5 QuasiBAM for variants

# LD_LIBRARY_PATH=~/bamtools/lib/ ~/git_repos/Quasi_bam/quasi_bam

# -m [30]         #mapq minimum
# -mc [50]        #min mapped read length
# -f [1]          #minimum percentage of variants reported
# -c [20]         #minimum variation for inclusion in consensus sequence - can be a comma separated list
# -q [33]         #asci string cutof for fastq
# -d [100]        #minimum number of reads required to report a consensus position - can be a comma separated list
# -cd [off]       #report all combinations of -c and -d, [on] report all seqs, [dr] report non-gap coverage only
# -i [on]         #insert optimisation
# -p [on]         #remove pcr duplicates
# -e [on]         #perform error checking based on position of minority variants in a read
# -e1 [0.1]       #fraction at each end of a read where errors are expected
# -e2 [0.9]       #fraction of reads containing a variant that are within the read region e1
# -st [off]       #do analaysis on fwd and reverse strands and print seperate text files
# -k [off]        #report fasta header as key value pairs
# -r [off]        #report mapped reference in sequence header
# -n [off]        #convert N's in consensus to gaps '-'
# -mg [off]       #mask gaps in mixed base/gap position where both frequencies are higher than -c
# -dr [off]       #report % of non-gap positions in fasta header at depth -d
# -pm [off]       #if multiple frags in ref - print mapping fragment
# -vb [off]       #Use [on] to generate stdout
# -o1 [stem.txt]  #variant outfile name
# -o2 [stem.fas]  #consensus fasta sequence
# -o3 [stem.err]  #error reporting
# -o6 [stem.stat.csv]     #mapped/filtered reads
# -o4 [stem.fwd.txt]      #forward variant outfile
# -o5 [stem.rev.txt]      #reverse variant outfile


import subprocess, os
from glob import glob
from pathlib import Path

# FOLDER = "220506_NS552079_0303_AHKY3HAFX3"
# BAMs = [Path(p) for p in sorted(glob(os.path.expanduser('~')+f"/Illumina_data/RVU/{FOLDER}/results/*/*sorted_dd.bam"))]
# ref_genome_filename = os.path.expanduser('~')+"/VRD/RVU/Adenovirus/2020-04 Adv WGS for infant hepatitis/MG925782.1.fna"

# ls /home/phe.gov.uk/david.williams/Illumina_data/RVU/220506_NS552079_0303_AHKY3HAFX3/results/1626321_Ad41v2HS_2-1/

# QuasiBAM bug: does not handle dots in path i.e., in folder names prior to file name.
# Silently fails to write to /home/phe ?

#!mkdir QuasiBAM

samples = {
        p.split('/')[7]:Path(p)
        for p in sorted(glob(f'{base_data_path}/results/*/'))
#        if 'WB_41F-' in p or 'NEG' in p
#        if 'Ad41v2SP' in p.split('/')[8] or 'NegMBv2SP' in p
#        if p.split('/')[8] in sample_selection
}
samples

from pathlib import Path

depth_for_cov = 15
min_mapping_q = 30

done = set()
for ref_genome,ref_genome_acc in targets.items():
    if ref_genome not in use_genomes:
        continue
    rec = seq_result[ref_genome_acc]
    ref_genome_path = f'{base_path}/{ref_genome_acc}.fna'
    ref_genome_nm = ref_genome.replace(' ','_')
    for sample_name,sample_outdir in samples.items():
        if (ref_genome,sample_name) in done:
            log.info(
                    f"Skipping {ref_genome} vs {sample_name} - "
                    "already done"
            )
            continue
        try:
            if ref_genome not in sample_to_ref_combos[sample_name]:
                continue
                log.info(f"Skipping {ref_genome} vs {sample_name}")
            log.info(f"Selecting {ref_genome} vs {sample_name}")
        except (KeyError, NameError):
            pass
        # raise ValueError()

        log.info(f"Doing {ref_genome} vs {sample_name}")
        bam_path_dd = Path(
                f'{sample_outdir}/'
                #f'ivar_trim__{sample_name}_v_{ref_genome_nm}_{ref_genome_acc}.sorted_dd.bam'
                f'ivar_trim__{sample_name}_v_{ref_genome_nm}_{ref_genome_acc}.sorted_dd_sorted2.bam'
                #f'{sample_name}_v_{ref_genome_nm}_{ref_genome_acc}.sorted_dd.bam'
        )
        # if bam_path_dd.exists():
            # raise Exception()
        stem = (
                f"{sample_outdir}/"
                f"QuasiBAM__ivar_trim__{sample_name}_v_{ref_genome_nm}_{ref_genome_acc}"
                #f"QuasiBAM__{sample_name}_v_{ref_genome_nm}_{ref_genome_acc}"
        )
        if not Path(f"{stem}.txt").exists():
            cmd = [
                os.path.expanduser('~')+'/mambaforge/pkgs/bamtools-2.5.2-hd03093a_0/Quasi_bam-master/quasi_bam',
                bam_path_dd,
                ref_genome_path,
                '-d', str(depth_for_cov),  # minimum depth to be included in a consensus
                '-c', '20',                # minimum proportion to affect an ambiguity
                '-m', str(min_mapping_q),  # minimum mapping quality for inclusion
                '-o1', f"{stem}.txt",      # variant outfile name
                '-o2', f"{stem}.fas",      # consensus fasta sequence
                '-o3', f"{stem}.err",      # error reporting
                '-o6', f"{stem}.stat.csv", # mapped/filtered reads
                '-o4', f"{stem}.fwd.txt",  # forward variant outfile
                '-o5', f"{stem}.rev.txt",  # reverse variant outfile
            ]
            log.info(cmd)
            subprocess.call(
                    cmd,
                    env = {'LD_LIBRARY_PATH':  os.path.expanduser('~')+'/mambaforge/pkgs/bamtools-2.5.2-hd03093a_0/lib/'},
            )
        done.add((ref_genome,sample_name))

# /home/phe,mapped:1887680,filtered:1887680
# this looks like a bug: "s"

## Generate archives of QuasiBam output from among shild folders
# !find ~/Illumina_data/RVU/220809_NS552079_0330_AHCCHTAFX3/results/ -type f -name "Quas*" -print0 | xargs -0 tar cfv ~/temp.tar
# !tar --strip-components 8 -xvf ~/temp.tar
# !zip ./QB_220809_NS552079_0330_AHCCHTAFX3.zip ./QuasiBAM__*
# !rm ./QuasiBAM__*



####################################################################################################################
####################################################################################################################
####################################################################################################################


### 1.6 parse QuasiBAM into dict

def parse_qbam_variants(filename, depth_resolution = 10):
    '''Parse variant file from QuasiBAM in mode 1 for SNPs: single sequence

    No codon information is returned. This is useful for comparing against
    a sequenced control (known sequence).

    Mode 2 is where coding regions also provided in the referernce FASTA input
    file: "Quasi_bam will assume that fragment 1 and fragment 2 sequences are
    in-frame"
    '''
    header = next(open(filename)).rstrip().split('\t')
    contig_raw = open(filename).read()
    freqs = []
    consensus = []
    depths = []
    for line in contig_raw.split('\n')[1:-1]:
        bits = line.split('\t')
        depths += [int(bits[2])]
        pos0 = int(bits[0])-1
        # raise Exception
        # A, C, G, T, Gap, Ins
        freqs += [{
                char:float(freq)
                for freq,char in zip(bits[3:9],header[3:9])
                if float(freq)
        }]
        consensus += [bits[10]]

    return freqs,consensus,depths[::depth_resolution]

# !less -S QuasiBAM/1626321_Ad41v2HS_2-1.txt

# sample_to_ref_combos = { sample_name : set() for sample_name in samples}
# sample_name = '1635099_WB_41F-2_JQ-1'
# ref_genome = 'Adenovirus F41 Tak'
# sample_to_ref_combos[sample_name] = {ref_genome,}

def parse_qb(targets, samples, stem_prefix):
    freqs = {}
    consensus = {}
    qbam_depths = {}
    all_variants = {}
    for ref_genome,ref_genome_acc in targets.items():
        if ref_genome not in use_genomes:
            continue
        ref_genome_path = f'{base_path}/{ref_genome_acc}.fna'
        ref_sequence = SeqIO.read(ref_genome_path,'fasta')
        ref_genome_nm = ref_genome.replace(' ','_')
        freqs[ref_genome] = {}
        consensus[ref_genome] = {}
        qbam_depths[ref_genome] = {}
        all_variants[ref_genome] = {}
        for sample_name,sample_outdir in samples.items():
            #print(sample_name)
            try:
                if ref_genome not in sample_to_ref_combos[sample_name]:
                    log.info(f"Skipping {ref_genome} vs {sample_name}")
                    continue
            except (KeyError, NameError):
                pass
            log.info(f"Selecting {ref_genome} vs {sample_name}")
            log.info(f"Doing {ref_genome} vs {sample_name}")
            # ref_sam_path = (
                    # f'{sample_outdir}/'
                    # f'{sample_name}_v_{ref_genome_nm}_{ref_genome_acc}.sam'
            # )
            # bam_path = ref_sam_path[:-3]+"sorted.bam"
            # bam_path_dd =  Path(bam_path[:-4]+'_dd.bam')
            stem = (
                    f"{sample_outdir}/"
                    f"{stem_prefix}{sample_name}_v_{ref_genome_nm}_{ref_genome_acc}"
            )
            filename = f"{stem}.txt"
            log.info(f'Parsing {filename} for {ref_genome} vs {sample_name}')
            #if not Path(filename).exists():
            #    continue
                #raise Exception()
            (
                    freqs[ref_genome][sample_name],
                    consensus[ref_genome][sample_name],
                    qbam_depths[ref_genome][sample_name],
            ) = parse_qbam_variants(filename, depth_resolution = 1)

            SNV_positions = []
            for pos, (r, s, f, d) in enumerate(zip(
                    ref_sequence.seq,
                    consensus[ref_genome][sample_name],
                    freqs[ref_genome][sample_name],
                    qbam_depths[ref_genome][sample_name],
            )):
                if r != s and d and f:
                    # Sometimes depth >1000 but no frequency info.
                    # Ignore those for now.
                    print(pos, f, d)
                    SNV_positions += [(pos,r,s,d,f)]
            all_variants[ref_genome][sample_name] = SNV_positions

    return freqs, consensus, qbam_depths, all_variants


#freqs, consensus, qbam_depths, all_variants = parse_qb(targets, samples,
 #       stem_prefix = 'QuasiBAM__')

freqs_ivar, consensus_ivar, qbam_depths_ivar, all_variants_ivar = parse_qb(targets, samples,
        stem_prefix = 'QuasiBAM__ivar_trim__')

#assert len()

##############################################################################################
##############################################################################################
##############################################################################################


### 1.8 collect BAM depths

import pysam
from array import array
#from glob import glob

resolution = 10

# Initial analysis
# bam_depths = pickle.load(gzip.open('adv_wgs_depths.p.gz','rb'))
# Mapping 3 samples to various 6 Adv references
# bam_depths = pickle.load(gzip.open('adv_wgs_depths.p.gz','rb'))
# len(bam_depths) # 20

from nisvirology.illumina import getCoverageDepths

### This already mapped well - why this time so patchy??
# ref_genome = 'Adenovirus F41'
# ref_genome_acc = targets[ref_genome] # 'MG925782.1'
# sample_name = '1626321_Ad41v2HS_2-1'

bam_depths = {}

# TODO implement softclipping aware depth collecting
# bam_depths_ivar = {}

for ref_genome,ref_genome_acc in targets.items():
    if ref_genome not in use_genomes:
        continue
    if not ref_genome in bam_depths:
        bam_depths[ref_genome] = {}
    # if not ref_genome in bam_depths_ivar:
        # bam_depths_ivar[ref_genome] = {}
    ref_genome_nm = ref_genome.replace(' ','_')
    for sample_name,sample_outdir in samples.items():
        try:
            if ref_genome not in sample_to_ref_combos[sample_name]:
               continue
        except KeyError:
            pass
        except NameError:
            pass
        bam_path = Path(
                f'{sample_outdir}/'
                f'{sample_name}_v_{ref_genome_nm}_{ref_genome_acc}'
                '.sorted_dd.bam'
        )
        # bam_path = Path(
                # f'{sample_outdir}/'
                # f'{sample_name}_v_{ref_genome_nm.replace("F41_MU22","F41")}_{ref_genome_acc}'
                # '.sorted_dd.bam'
        # )
        assert bam_path.exists()
        # bam_path_ivar = Path(
                # f'{sample_outdir}/'
                # f'ivar_trim__{sample_name}_v_{ref_genome_nm}_{ref_genome_acc}'
                # '.sorted_dd.bam'
        # )
        # assert bam_path_ivar.exists()
        if not sample_name in bam_depths[ref_genome]:
            log.info(
                    f"getting depths for {ref_genome} vs {sample_name}"
            )
            # Returns tuple of dicts - we only need value of first dict: one
            # replicon in this case so one dict entry. Second dict is proper
            # pairs only.
            read_depths, _ = getCoverageDepths(
                            bam_path,
                            resolution,
                            min_mapping_quality = min_mapping_q,
                            max_depth = 30000,
            )
            bam_depths[ref_genome][sample_name] = list(read_depths.values())[0]
        log.info(f'done {sample_name} for {ref_genome}: {bam_path}')

        # if not sample_name in bam_depths_ivar[ref_genome]:
            # log.info(
                    # f"getting depths for {ref_genome} vs {sample_name} iVar filtered"
            # )
            # read_depths, _ = getCoverageDepths(
                            # bam_path_ivar_sorted,
                            # resolution,
                            # min_mapping_quality = min_mapping_q,
                            # max_depth = 30000,
            # )
            # bam_depths_ivar[ref_genome][sample_name] = list(read_depths.values())[0]
        # log.info(f'done {sample_name} for {ref_genome}: {bam_path} iVar filtered')


len(bam_depths)
len(bam_depths['England174820313'])
pickle.dump(
    bam_depths,
    gzip.open(f'rsv1_wgs_depths__{run_folder}.p.gz','wb'),
)
assert bam_depths == pickle.load(
    gzip.open(f'rsv1_wgs_depths__{run_folder}.p.gz','rb'),
)
bam_depths = pickle.load(
    gzip.open(f'rsv1_wgs_depths__{run_folder}.p.gz','rb'),
)

len(bam_depths['England174880237'])
pickle.dump(
    bam_depths,
    gzip.open(f'rsv2_wgs_depths__{run_folder}.p.gz','wb'),
)
assert bam_depths == pickle.load(
    gzip.open(f'rsv2_wgs_depths__{run_folder}.p.gz','rb'),
)
bam_depths = pickle.load(
    gzip.open(f'rsv2_wgs_depths__{run_folder}.p.gz','rb'),
)



###########################################################################################
###########################################################################################
###########################################################################################


### 1.9 check variant proximity to gaps

from nisvirology.utils import scan_threshold

resolution = 10
depth_for_cov = 15

pad = 100

freqs_padded = {}
consensus_padded = {}
all_variants_padded = {}

for ref_genome, sample_depths in bam_depths.items():
    freqs_padded[ref_genome] = {}
    consensus_padded[ref_genome] = {}
    all_variants_padded[ref_genome] = {}
    for sample_name, read_depths in sample_depths.items():
        log.info(f"{sample_name} vs {ref_genome} depths")
        # if 'H19108076801' in sample_name:
        #     raise ValueError()
        # Calculate coverage at specified depth
        deep_ranges = scan_threshold(
                read_depths,
                depth_for_cov,
                offset = 0,
                test = '>=',
                filter_name = 'threshold',
                resolution = resolution,
                buffer_distance = 0,
        )
        # Shrink the ranges considered deep enough
        deep_ranges_padded = [
                (s+pad,e-pad)
                for (s,e) in deep_ranges
                if e-s > pad*2
        ]

        these_freqs_padded = []
        these_consensus_padded = []
        for pos,(f,c) in enumerate(zip(
                freqs_ivar[ref_genome][sample_name],
                consensus_ivar[ref_genome][sample_name],
        )):
            selected = False
            for s,e in deep_ranges_padded:
                if s <= pos < e:
                    these_freqs_padded += [f]
                    these_consensus_padded += [c]
                    selected = True
                    break
            if not selected:
                these_freqs_padded += [{}]
                these_consensus_padded += ['N']
                log.info(f"{sample_name} filtered at {pos}")

        these_all_variants_padded = []
        for pos,r,sq,d,f in all_variants_ivar[ref_genome][sample_name]:
            selected = False
            for s,e in deep_ranges_padded:
                if s <= pos < e:
                    these_all_variants_padded += [(pos,r,sq,d,f)]
                    selected = True
                    break
            if not selected:
                log.info(
                        f"{sample_name} variant filtered at "
                        f"{pos} {r}, {sq}, {d}, {f}"
                )

        freqs_padded[ref_genome][sample_name] = \
                these_freqs_padded
        consensus_padded[ref_genome][sample_name] = \
                these_consensus_padded
        all_variants_padded[ref_genome][sample_name] = \
                these_all_variants_padded


###############################################################################################
###############################################################################################
###############################################################################################


### 2. plot coverage

import pickle
import gzip
from Bio import SeqIO

# NB other potential references
# https://www.ncbi.nlm.nih.gov/genome/browse/#!/viruses/4693/

# ref_genome_filename = 'MG925782.1.fna'
# ref_sequence = SeqIO.to_dict(SeqIO.parse(ref_genome_filename,'fasta'))
# Annoying identifier.
# use_id = 'gi|1511893241|gb|MG925782.1|'

from nisvirology import illumina
from glob import glob
import os
import pathlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from nisvirology.utils import scan_threshold

depth_for_cov = 15

from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text

def scan_for_pGC(ref_replicon_seq, window_length = 100, step_size = 10):
    '''Return proportion GC / GCAT
    
    Returns
    -------
    pGC_y : list
        proportion
    pGC_x : list
        postion on 
    '''
    pGC_y = []
    pGC_x = []
    for s in range(0, len(ref_replicon_seq)-window_length, step_size):
        window_seq = ref_replicon_seq[s:s+window_length]
        window_pGC = \
                (window_seq.count('G') + \
                window_seq.count('C')) / window_length
        log.debug(f"{s} {s+window_length} {window_pGC}")
        pGC_y += [window_pGC]
        pGC_x += [s+(window_length/2)]
    return pGC_x, pGC_y

pGC_by_ref = {}
for ref_genome,ref_genome_acc in targets.items():
    if ref_genome not in use_genomes:
        continue
    ref_genome_path = f'{base_path}/{ref_genome_acc}.fna'
    ref_sequence = SeqIO.read(ref_genome_path,'fasta')
    pGC_by_ref[ref_genome] = scan_for_pGC(ref_sequence.seq)

## B. plotting without amplicon labels

plot_lim = 40
plot_num = 0
mapping_q = 30
out_path = Path(f'plots_{run_folder}')
if not out_path.exists():
    out_path.mkdir()

sorted(samples)

## sample v primer set would be useful too
sample_to_assay_title = {}
sample_to_assays = {}
for sample_name in samples:
    if '_a' in sample_name:
        sample_to_assay_title[sample_name] = 'RSVA_1'
        sample_to_assays[sample_name] = {
                'HRSVA1_pool1':'ALL REGIONS',
                'HRSVA1_pool2':'ALL REGIONS',
                'HRSVA2_pool1':'ALL REGIONS',
                'HRSVA2_pool2':'ALL REGIONS'
        }
    elif '-a_' in sample_name:
        sample_to_assay_title[sample_name] = 'RSVA_1'
        sample_to_assays[sample_name] = {
                'HRSVA1_pool1':'ALL REGIONS',
                'HRSVA1_pool2':'ALL REGIONS',
                'HRSVA2_pool1':'ALL REGIONS',
                'HRSVA2_pool2':'ALL REGIONS'
        }
    elif '_a-' in sample_name:
        sample_to_assay_title[sample_name] = 'RSVA_1'
        sample_to_assays[sample_name] = {
                'HRSVA1_pool1':'ALL REGIONS',
                'HRSVA1_pool2':'ALL REGIONS',
                'HRSVA2_pool1':'ALL REGIONS',
                'HRSVA2_pool2':'ALL REGIONS'
        }
    elif '_a_' in sample_name:
        sample_to_assay_title[sample_name] = 'RSVA_1'
        sample_to_assays[sample_name] = {
                'HRSVA1_pool1':'ALL REGIONS',
                'HRSVA1_pool2':'ALL REGIONS',
                'HRSVA2_pool1':'ALL REGIONS',
                'HRSVA2_pool2':'ALL REGIONS'
        }
    elif '_b' in sample_name:
        sample_to_assay_title[sample_name] = 'RSVB_1'
        sample_to_assays[sample_name] = {
                'HRSVB1_pool1':'ALL REGIONS',
                'HRSVB1_pool2':'ALL REGIONS',
                'HRSVB2_pool1':'ALL REGIONS',
                'HRSVB2_pool2':'ALL REGIONS'
        }
    elif '-b_' in sample_name:
        sample_to_assay_title[sample_name] = 'RSVB_1'
        sample_to_assays[sample_name] = {
                'HRSVB1_pool1':'ALL REGIONS',
                'HRSVB1_pool2':'ALL REGIONS',
                'HRSVB2_pool1':'ALL REGIONS',
                'HRSVB2_pool2':'ALL REGIONS'
        }
    elif '_b-' in sample_name:
        sample_to_assay_title[sample_name] = 'RSVB_1'
        sample_to_assays[sample_name] = {
                'HRSVB1_pool1':'ALL REGIONS',
                'HRSVB1_pool2':'ALL REGIONS',
                'HRSVB2_pool1':'ALL REGIONS',
                'HRSVB2_pool2':'ALL REGIONS'
        }
    elif '_b_' in sample_name:
        sample_to_assay_title[sample_name] = 'RSVB_1'
        sample_to_assays[sample_name] = {
                'HRSVB1_pool1':'ALL REGIONS',
                'HRSVB1_pool2':'ALL REGIONS',
                'HRSVB2_pool1':'ALL REGIONS',
                'HRSVB2_pool2':'ALL REGIONS'
        }
    
plot_samples = sorted(samples)

#plot_samples = []
#plot_samples = df['Sample:'].tolist()

plot_samples = [
       
        '01121231_313_1_a-1',
        '01121232_313_2_397_4_a-1',
        '01121233_313_3_397_3_a-1',
        '01121234_313_4_397_2_a-1',
        '01121235_397_1_a-1',
        '01121236_237_2_b-1',
        '01121237_237_3_410_5_b-1',
        '01121238_237_4_410_4_b-1',
        '01121239_237_5_410_3_b-1',
        '01121240_410_2_b-1',
      

      
]

#1 A PDF per sample against each ref for easier comparison
for sample_name in plot_samples:
    use_name = '_'.join(sample_name.split('_')[1:])
    file_name = (
            f'rsv_wgs_q{mapping_q}_d{depth_for_cov}_'
            f'{use_name}_vs_{len(use_genomes)}refs.pdf'
    )
    log.info(f"Plotting to '{out_path}/{file_name}'")
    with PdfPages(f'{out_path}/{file_name}') as pdf_pages:
        figure, subplots = plt.subplots(
                nrows=2, ncols=1,
                sharex = 'all', sharey = 'none',
                gridspec_kw = {'height_ratios': [0.85,0.15]},
        )
        # Upper
        main_plot = subplots[0]
        # Lower
        extra_plot = subplots[1]
        # Fix inconsistent layout on first plot (page) by
        # inserting a blank one.
        figure.suptitle('Blank (first axes slightly different size to others)')
        figure.tight_layout()
        pdf_pages.savefig(figure)



y_max = 30000
amp_line_interval = 3 # /100
text_offset = 2.5  #  1.5 or 6
amplicon_min_dist = 0.10 # proportion of x range
rctn_init_pos_y = 0

resolution = 10

plot_text = False

# default for QB is 30

# as used for quasibam -d, which combined with -c determines variants
# called, and used in plotting to visualise and calculate coverage
depth_for_cov = 15

# TODO implement soft-clipped aware getCoverage()
include_ivar_depths = False

# all to one PDF multi samples and refs, usually one sample to one ref
file_name = (
        f'{run_folder}__{len(plot_samples)}samples_vs_{len(use_genomes)}refs_'
        f'mq{mapping_q}_d{depth_for_cov}_'
        #'ivar_'
        'RSV_WGS_amplicons_plot.pdf'
)
log.info(f"Plotting to '{out_path}/{file_name}'")
with PdfPages(f'{out_path}/{file_name}') as pdf_pages:
    figure, subplots = plt.subplots(
            nrows=2, ncols=1,
            sharex = 'all', sharey = 'none',
            gridspec_kw = {'height_ratios': [0.85,0.15]},
    )
    # Upper
    main_plot = subplots[0]
    # Lower
    extra_plot = subplots[1]
    # Fix inconsistent layout on first plot (page) by
    # inserting a blank one.
    figure.suptitle('Blank (ignore - first axes slightly different size to others, a bug)')
    figure.tight_layout()
    pdf_pages.savefig(figure)
    for sample_name in plot_samples:
        use_name = '_'.join(sample_name.split('_')[1:])


        # Plot each sample for which we have mapped read depths
        for ref_genome,ref_acc in targets.items():
            if ref_genome not in use_genomes:
                continue
            try:
                if ref_genome not in sample_to_ref_combos[sample_name]:
                    continue
            except (NameError, KeyError):
                pass
            pGC_x, pGC_y = pGC_by_ref[ref_genome]
            ref_sequence = SeqIO.read(f"{ref_acc}.gbk", 'genbank')

            # Plot the depth data as a line
            read_depths = bam_depths[ref_genome][sample_name]
            #depths = info[0][accession]
            positions = np.arange(0,len(read_depths)) * resolution
            main_plot.plot(positions, read_depths, lw = 0.5, color = '#2200ee')

            if include_ivar_depths:
                read_depths_ivar = bam_depths_ivar[ref_genome][sample_name]
                positions = np.arange(0,len(read_depths)) * resolution
                main_plot.plot(positions, read_depths, lw = 0.5, color = '#ee22ee')


            # Calculate coverage at specified depth
            deep_ranges = scan_threshold(read_depths, depth_for_cov,
                                    offset = 0,
                                    test = '>=',
                                    filter_name = 'threshold',
                                    resolution = 10,
                                    buffer_distance = 0,
            )
            total_deep = sum([(e-s) for s,e in deep_ranges])
            pc_covered = total_deep / len(ref_sequence.seq)

            # Plot regions with sufficient coverage as translucent
            # boxes over the depth
            for s,e in deep_ranges:
                if 'water' in sample_name:
                    main_plot.add_patch(
                    Rectangle(
                            (s,1), e-s, depth_for_cov,
                            alpha = 0.5,
                            color = 'chartreuse',
                            lw = 0,
                    )
                )
                elif 'WATER' in sample_name:
                    main_plot.add_patch(
                    Rectangle(
                            (s,1), e-s, depth_for_cov,
                            alpha = 0.5,
                            color = 'chartreuse',
                            lw = 0,
                    )
                )
                elif 'wtr' in sample_name:
                    main_plot.add_patch(
                    Rectangle(
                            (s,1), e-s, depth_for_cov,
                            alpha = 0.5,
                            color = 'chartreuse',
                            lw = 0,
                    )
                )
                elif 'VTM' in sample_name:
                    main_plot.add_patch(
                    Rectangle(
                            (s,1), e-s, depth_for_cov,
                            alpha = 0.5,
                            color = 'chartreuse',
                            lw = 0,
                    )
                )
                elif 'NEG' in sample_name:
                    main_plot.add_patch(
                    Rectangle(
                            (s,1), e-s, depth_for_cov,
                            alpha = 0.5,
                            color = 'chartreuse',
                            lw = 0,
                    )
                )
                elif 'neg' in sample_name:
                    main_plot.add_patch(
                    Rectangle(
                            (s,1), e-s, depth_for_cov,
                            alpha = 0.5,
                            color = 'chartreuse',
                            lw = 0,
                    )
                )
                elif 'POS' in sample_name:
                    main_plot.add_patch(
                    Rectangle(
                            (s,1), e-s, depth_for_cov,
                            alpha = 0.5,
                            color = 'orange',
                            lw = 0,
                    )
                )
                elif 'pos' in sample_name:
                    main_plot.add_patch(
                    Rectangle(
                            (s,1), e-s, depth_for_cov,
                            alpha = 0.5,
                            color = 'orange',
                            lw = 0,
                    )
                )
                elif 'RSVAB_B' in sample_name:
                    main_plot.add_patch(
                    Rectangle(
                            (s,1), e-s, depth_for_cov,
                            alpha = 0.5,
                            color = 'orange',
                            lw = 0,
                    )
                )
                else:
                    main_plot.add_patch(
                    Rectangle(
                            (s,1), e-s, depth_for_cov,
                            alpha = 0.5,
                            color = 'mediumslateblue',
                            lw = 0,
                    )
                )

            # Could this change after tight_layout() below?
            # Establish x-axis range manually and apply and refer to later on
            # (because .get_xlim() was having weird side effects on the range)
            x_max = resolution*len(positions)
            main_plot.set_xlim(0, x_max)
            main_plot.tick_params(
                axis = 'both',        # Changes apply to the x-axis.
                which = 'major',      # Major ticks are affected.
                labelsize = 'small',  # Smaller than default to be less
                                      # intrusive on a detailed plot.
            )
            extra_plot.tick_params(
                axis = 'both',        # Changes apply to the x-axis.
                which = 'major',      # Major ticks are affected.
                labelsize = 'small',  # Smaller than default to be less
                                      # intrusive on a detailed plot.
            )

            # Plot ORFs (not in use - too much on one axis)
            # main_plot.vlines(max(positions),0,5000,linestyle='dotted')
            # for i,(ORF,y) in enumerate(sorted(segment2coding[x].items())):
                # for h,(z,(s,e)) in enumerate(y.items()):
                    # main_plot.add_line(
                        # Line2D( [s,e], [i*m,i*m], linestyle='solid', color = 'red')
                    # )
                # main_plot.set_title('main')
            # plot amplicons TODO colour code by reaction and amplicon

            # Instantiate a second axes that shares the same x-axis for adding
            # "artists" (annotations) onto a simple linear 0-100 y scale.
            ax2 = main_plot.twinx()
            ax2.set_ylim(0,100)
            # No ticks on second axis (on the right)
            ax2.tick_params(
                axis = 'y',          # changes apply to the x-axis
                which = 'both',      # both major and minor ticks are affected
                right = False,       # ticks along the bottom edge are off
                labelright = False,  # ticks along the bottom edge are off
            )

            # Plot the amplicons in each reaction (multiplexed have more than one)

            #for n,(reaction,amplicons) in enumerate(use_amplicons.items()):

            ### HERE: allow for specified amplicons: sample_to_assays
            ### newer toml layout is more specific per amplicon within a reaction

            n = 0
            for assay_name,amplicons in amplicon_ranges[ref_genome].items():
                #continue
                if not assay_name in sample_to_assays[sample_name]:
                    # No data for this assay with this sample: skip
                    continue
                rctn_base_pos_y = (
                        rctn_init_pos_y + n *
                        ((amp_line_interval + text_offset) * 2 * 1.0)
                )
                if n:
                    # Horizontal break line between mutliplexed reactions
                    # (except first iteration when n == 0)
                    ax2.hlines(rctn_base_pos_y, 0, x_max, linestyles='dotted',
                            lw=0.5, color = 'black')


                #log.debug(sample_name,reaction,primer_set_name)
                log.debug(
                        f"Plotting amplicons for {reaction} at "
                        f"y={rctn_base_pos_y}"
                )
                # Offset adjacent amplicons if they are too close
                last_e = 0
                amplicon_y = amp_line_interval
                amplicon_y_raised = False
                last_amplicon_y_raised = False
                min_dist_positions = amplicon_min_dist * ax2.get_xlim()[1]
                for (region, amp_id), (amp_start, amp_end) in amplicons.items():
                    if not (
                        region in sample_to_assays[sample_name][assay_name]
                        or
                        sample_to_assays[sample_name][assay_name] == 'ALL REGIONS'
                    ):
                        log.info(
                                f"For {assay_name} with {sample_name}, "
                                f"omitting region {region} from plot"
                        )
                        continue

                    log.info(
                            f"For {assay_name} with {sample_name}, "
                            f"plotting region {region}"
                    )

                    # If distance to previous plotted amplicon is shorter
                    # than amplicon_min_dist, plot slightly higher.
                    if last_e and amp_start < last_e + min_dist_positions:
                        if not last_amplicon_y_raised:
                            amplicon_y += amp_line_interval
                            amplicon_y_raised = True
                            last_amplicon_y_raised = True
                        else:
                            # Last one is too close but was already raised, so
                            # plot this one on lower level
                            last_amplicon_y_raised = False

                    # ax2.add_line(
                        # Line2D(
                                # [amp_start, amp_end],
                                # [rctn_base_pos_y+amplicon_y, rctn_base_pos_y+amplicon_y],
                                # linestyle='solid', color = 'red',
                                # lw = 3,
                                # solid_capstyle = 'butt',
                        # )
                    # )

                    ax2.add_patch(
                        Rectangle(
                                #(s,1), e-s, depth_for_cov,
                                (amp_start, rctn_base_pos_y+amplicon_y), amp_end-amp_start, 1.5,
                                alpha = 0.6,
                                facecolor = 'purple',
                                edgecolor = 'black',
                                linestyle = '-',
                                lw = 0,
                        )
                    )

                    if plot_text:
                        f_name, r_name = amp_id
                        ax2.add_artist(
                            Text(
                                x = (amp_start+amp_end)/2,
                                y = rctn_base_pos_y+amplicon_y+text_offset,
                                text = (
                                        f"{region}\n{f_name}+{r_name}\n"
                                        f"{amp_start}-{amp_end}"
                                ),
                                fontsize = 3,
                                color=None,
                                verticalalignment='baseline', horizontalalignment='center',
                                multialignment=None,
                                fontproperties=None,
                                rotation=None, linespacing=None, rotation_mode=None,
                                usetex=None, wrap=False,
                                transform_rotates_text=False, parse_math=True,
                            )
                        )
                    if amplicon_y_raised:
                        amplicon_y -= amp_line_interval
                        amplicon_y_raised = False
                    last_e = amp_end
                n += 1

            # if 'H19108076801' in sample_name:
                # raise ValueError()
            # Plot variants as vertical lines as indicator of quality for
            # sequencing of known sequence and mapping back to that
            # reference sequence.
            try:
                SNV_positions = all_variants[ref_genome][sample_name]

                fixed_positions = [
                        pos for (pos,r,s,d,f) in SNV_positions
                        if max(f.values()) >= 99
                ]
                allele_positions = [
                        pos for (pos,r,s,d,f) in SNV_positions
                        if max(f.values()) < 99
                ]
            except NameError:
                SNV_positions = []
                fixed_positions = []
                allele_positions = []

            try:
                SNV_positions_ivar = all_variants_ivar[ref_genome][sample_name]
                ivar_title_suffix = ' iVar trim'
                fixed_positions_ivar = [
                        pos for (pos,r,s,d,f) in SNV_positions_ivar
                        if max(f.values()) >= 99
                ]
                allele_positions_ivar = [
                        pos for (pos,r,s,d,f) in SNV_positions_ivar
                        if max(f.values()) < 99
                ]
                fixed_positions_filtered = sorted(
                        set(fixed_positions) - set(fixed_positions_ivar)
                )
                allele_positions_filtered = sorted(
                        set(allele_positions) - set(allele_positions_ivar)
                )
            except NameError:
                SNV_positions_ivar = []
                ivar_title_suffix = ''
                fixed_positions_ivar = []
                allele_positions_ivar = []
                fixed_positions_filtered = []
                allele_positions_filtered = []

            try:
                SNV_positions_padded = all_variants_padded[ref_genome][sample_name]
                fixed_positions_padded = [
                        pos for (pos,r,s,d,f) in SNV_positions_padded
                        if max(f.values()) >= 99
                ]
                allele_positions_padded = [
                        pos for (pos,r,s,d,f) in SNV_positions_padded
                        if max(f.values()) < 99
                ]
                fixed_positions_filtered2 = sorted(
                        set(fixed_positions_ivar) - set(fixed_positions_padded))
                allele_positions_filtered2 = sorted(
                        set(allele_positions_ivar) - set(allele_positions_padded)
                )
            except NameError:
                SNV_positions_padded = []
                fixed_positions_padded = []
                allele_positions_padded = []
                fixed_positions_filtered2 = []
                allele_positions_filtered2 = []


            if fixed_positions_ivar:
                main_plot.vlines(
                        fixed_positions_ivar,
                        0, y_max,
                        linestyles='-',
                        lw = 0.5,
                        color = 'red',
                )
            if allele_positions_ivar:
                main_plot.vlines(
                        allele_positions_ivar,
                        0, y_max,
                        linestyles='dotted',
                        lw = 0.5,
                        color = 'red',
                )
            if fixed_positions_filtered:
                main_plot.vlines(
                        fixed_positions_filtered,
                        0, y_max,
                        linestyles='-',
                        lw = 0.5,
                        color = 'green',
                )
            if allele_positions_filtered:
                main_plot.vlines(
                        allele_positions_filtered,
                        0, y_max,
                        linestyles='dotted',
                        lw = 0.5,
                        color = 'green',
                )
            if fixed_positions_filtered2:
                main_plot.vlines(
                        fixed_positions_filtered2,
                        0, y_max,
                        linestyles='-',
                        lw = 0.5,
                        color = 'dimgrey',
                )
            if allele_positions_filtered2:
                main_plot.vlines(
                        allele_positions_filtered2,
                        0, y_max,
                        linestyles='dotted',
                        lw = 0.5,
                        color = 'dimgrey',
                )

            # raise ValueError()

            # Label axes etc
            main_plot.set_ylabel('Read depth')
            # PCR amplification characteristics usually mean log scale is useful
            # for depths.
            main_plot.set_yscale('log')
            main_plot.set_ylim(1, y_max)
            x_min, x_max = main_plot.get_xlim()
            main_plot.set_xlim(-500, x_max+500)

            # Plot %GC content
            extra_plot.set_yscale('linear')
            extra_plot.plot(pGC_x, pGC_y, lw = 0.5, color = '#292')

            extra_plot.set_ylabel('GC/ATGC')
            extra_plot.set_xlabel('position (bp)')

            assay_title = sample_to_assay_title[sample_name] + ivar_title_suffix
            figure.suptitle(
                    f"{sample_name} + {assay_title}\n"
                    f"vs {ref_genome} ({ref_acc}).\n"
                    f"{pc_covered:.1%} at depth {depth_for_cov}",
                    fontsize = 'medium',
            )

            # Some other padding or margin keeping them apart
            plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                    wspace=None, hspace=0)

            figure.tight_layout()
            # savefig() writes the page.
            pdf_pages.savefig(figure)
            # After writing the page, clear the axes.
            main_plot.clear()
            ax2.clear()
            extra_plot.clear()
            plot_num += 1
            ### Is plt.close(figure) needed? clearing axes was implying
            ### closing didn't do anything.
            # Closing the figure writes a page (but without .clear() plotting is
            # retained in the axes)
            # plt.close(figure)

