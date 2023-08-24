#!/usr/bin/python3

'''
rsv_pipeline_ivar.py

DW 20220620
- initial script
TT 20220920
- adjusted for RSV
- input: .bam file from bwa aligner, RSV primer .toml file
- output: .bam file with filtered mapped reads

'''

#DEPENDENCIES:

from Bio import SeqIO
from pathlib import Path
import os
import sys
import pathlib
from glob import glob
import subprocess
import pysam
import toml
from Bio.Seq import Seq
import importlib
from nisvirology.align import SinglePCRAssay
from nisvirology import align
from nisvirology.align import (
        SinglePCRAssay,
        MappedOrientationError,
        DidNotMapError,
)


#Hard coded bits - locations etc, needs to be adjusted to users system:

picard_path = '/home/tiina/mambaforge/share/picard-2.27.4-0/picard.jar'
base_path = Path(os.path.abspath(os.path.curdir))
sys.path.insert(0, '../AdV_wgs2')
run_folder = '220816_VL00115_120_AAAWW3YM5'
parameter_file = (
        os.path.expanduser('~') +
        "/AdV_wgs2/220816_VL00115_120_AAAWW3YM5_NF/RSV Orthopneumovirus.toml"
)
base_data_path = Path(
        os.path.expanduser('~') + (
            f'/AdV_wgs2/220816_VL00115_120_AAAWW3YM5_NF/{run_folder}'
        )
)

path_to_tmp = f'{base_path}/tmp'

#mapping reference:
targets = {
    'RSVB MZ515970': 'MZ515970'
}









#Code start:

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
}


seq_result = {}
for ref_genome_type,ref_genome_acc in targets.items():
    seq_result[ref_genome_acc] = SeqIO.read(
            f'{base_path}/{ref_genome_acc}.gbk',
            'genbank'
    )
try:
    parameters = toml.load(parameter_file)
except toml.TomlDecodeError:
    sys.exit(1)

for ref_genome,ref_genome_acc in targets.items():
    rec = seq_result[ref_genome_acc]
  

samples = {
        p.split('/')[7]:Path(p)
        for p in sorted(glob(f'{base_data_path}/results/*/'))
}


pcr_assays = {}
for ref_genome,ref_genome_acc in targets.items():
    ref_genome_path = f'{base_path}/{ref_genome_acc}.fna'
    ref_genome_rec = SeqIO.read(ref_genome_path, 'fasta')
    pcr_assays[ref_genome] = {}
    for assay_name,assay_info in parameters['Assays'].items():
        if not (assay_name.startswith('HRSVA1')
               or
               'HRSVA2' in assay_name
                or
               'HRSVA3' in assay_name
                or
               'HRSVB1' in assay_name
                or
               'HRSVB2' in assay_name
        ):
            continue
        if 'multiplex' in assay_info:
            for r_name,reaction in assay_info['multiplex'].items():
                fwd_name, fwd_seq = reaction['forward']
                rev_name, rev_seq = reaction['reverse']
                if assay_info.get('reverse_strand_same_as_forward'):
                    rev_seq = str(Seq(rev_seq).reverse_complement())
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
                    raise
                except DidNotMapError as exc:
                    continue
                if not pcr_assay.fwd_primer.mapping.end < \
                        pcr_assay.rev_primer.mapping.start:
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
                raise
            except DidNotMapError as exc:
                continue

bed_data = {}
primerpair_data = {}
primer_set = []
for ref_genome,assays in pcr_assays.items():
    bed_data[ref_genome] = {}
    primerpair_data[ref_genome] = {}
    for assay_name,assay in assays.items():
        if isinstance(assay_name, tuple):
            use_assay_name,region = assay_name
        else:
            use_assay_name,region = assay_name,'region1'
        if not use_assay_name in bed_data[ref_genome]:
            bed_data[ref_genome][use_assay_name] = []
        if not use_assay_name in primerpair_data[ref_genome]:
            primerpair_data[ref_genome][use_assay_name] = []
        primer_set.append(assay_name[0])
        for oligo in assay.iter_oligos():
            line = "{}\t{}\t{}\t{}\t60\t{}".format(
                    assay.target.id,
                    oligo.mapping.start,
                    oligo.mapping.end,
                    oligo.name,
                    '+' if oligo.mapping.strand == 1 else '-',
            )
            bed_data[ref_genome][use_assay_name] += [line]
        primerpair_data[ref_genome][use_assay_name] += [
                f"{assay.fwd_primer.name}\t{assay.rev_primer.name}"
        ]
primer_set = str(set(primer_set))
primer_set = ''.join(n for n in primer_set if n.isalnum())
sorted(bed_data[ref_genome])
assay_names = set()
for ref_genome,assays in bed_data.items():
    use_lines = []
    use_lines_pp = []
    for assay_name, lines in assays.items():
            use_lines += lines
            assay_names.add(assay_name)
            use_lines_pp += primerpair_data[ref_genome][assay_name]
    ref_genome_acc = targets[ref_genome]
    ref_genome_nm = ref_genome.replace(' ','_')
    bed_path = Path(
            f'{ref_genome_nm}_{ref_genome_acc}_{primer_set}.bed'
    )
    pp_path = Path(
            f'{ref_genome_nm}_{ref_genome_acc}_{primer_set}.pp'
    )
    bed_path.open('w').write('\n'.join(use_lines))
    pp_path.open('w').write('\n'.join(use_lines_pp))
    
for sample_name,(R1, R2) in sorted(sample_fastqs.items()):
    indir = R1.parent
    break

depth_for_cov = 15
min_mapping_q = 30

for ref_genome,ref_genome_acc in targets.items():
    rec = seq_result[ref_genome_acc]
    ref_genome_path = f'{base_path}/{ref_genome_acc}.fna'
    ref_genome_nm = ref_genome.replace(' ','_')
    for sample_name,sample_outdir in samples.items():
        bam_path_dd = Path(
                f'{sample_outdir}/'
                f'{sample_name}_v_{ref_genome_nm}_{ref_genome_acc}.sorted_dd.bam'
        )
        if not bam_path_dd.exists():
            continue
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
            subprocess.call(
                    cmd
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
   
