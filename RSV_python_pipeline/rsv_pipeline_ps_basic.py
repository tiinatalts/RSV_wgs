#!/usr/bin/python3

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
from nisvirology.align import SinglePCRAssay

#Hard coded bits - locations etc, needs to be adjusted to users system:
#adjust line 317 and 332 within the code
#also 137 to 145 according to primer set used for the samples in analysis 

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
targets = {
    'RSVB MZ515970': 'MZ515970'
}
args_to_bwa = [['-T','4']]
path_to_tmp = '/home/tiina/AdV_wgs2/220816_VL00115_120_AAAWW3YM5_RB2/tmp'

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

len(sample_fastqs) == 21

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
    ref_genome_path = f'{base_path}/{ref_genome_acc}.fna'
    ref_gnm_fname = ref_genome.replace(' ','_')
    for sample_name,(R1, R2) in sorted(sample_fastqs.items()):
        indir = R1.parent
        # mapping using BWA MEM (VirVarSeq used BWA sampe)
        sample_outdir = f'{indir}/results/{sample_name}'
        ref_sam_path = Path(
                f'{sample_outdir}/'
                f'{sample_name}_v_{ref_gnm_fname}_{ref_genome_acc}.sam'
        )
        ref_bam_path = Path(str(ref_sam_path)[:-3]+"bam")
        bam_path = Path(str(ref_sam_path)[:-3]+"sorted.bam")
        bam_path_dd = Path(str(bam_path)[:-4]+'_dd.bam')
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
            with open(ref_sam_path, 'w') as sam:
                proc = subprocess.Popen(
                        cmd,
                        stdout = sam,
                        stderr = subprocess.PIPE
                )
                o,bwa_stderr = proc.communicate()
                bwa_stderr = bwa_stderr.decode('utf-8')
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
            subprocess.call(cmd)
            pysam.index(str(bam_path_dd))

samples = {
        p.split('/')[7]:Path(p)
        for p in sorted(glob(f'{base_data_path}/results/*/'))
}
print(samples)

pcr_assays = {}
for ref_genome,ref_genome_acc in targets.items():
    ref_genome_path = f'{base_path}/{ref_genome_acc}.fna'
    ref_genome_rec = SeqIO.read(ref_genome_path, 'fasta')
    pcr_assays[ref_genome] = {}
    for assay_name,assay_info in parameters['Assays'].items():
        if not (assay_name.startswith('HRSVB2')
               #or
               #'HRSVA2' in assay_name
                #or
               #'HRSVA3' in assay_name
                #or
               #'HRSVB1' in assay_name
                #or
               #'HRSVA1' in assay_name
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
        bam_path_dd = Path(
                f'{sample_outdir}/'
                f'ivar_trim__{sample_name}_v_{ref_genome_nm}_{ref_genome_acc}.sorted_dd_sorted2.bam'
        )
        stem = (
                f"{sample_outdir}/"
                f"QuasiBAM__ivar_trim__{sample_name}_v_{ref_genome_nm}_{ref_genome_acc}"
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
            subprocess.call(
                    cmd,
                    env = {'LD_LIBRARY_PATH':  os.path.expanduser('~')+'/mambaforge/pkgs/bamtools-2.5.2-hd03093a_0/lib/'}
            )
