#!/bin/env python3
import os
import os.path
import copy
from Bio import SeqIO, AlignIO
import numpy as np

from defaults import reference_sequences

#def parse_args():
 #   import argparse
  #  parser = argparse.ArgumentParser(description='Annotate sequences using a genbank reference')
   # parser.add_argument('--reference', help='Genbank accession of reference sequence (will be fetched from genbank)')
  #  parser.add_argument('--virus', help='Virus name, default reference for the virus will be fetched from genbank')
   # parser.add_argument('--sequences', required=True, help='Fasta file of sequences to annotate')
  #  parser.add_argument('--output-format', choices=['genbank', 'gff3', 'tbl'], default='genbank', type=str, help='Output format')
   # parser.add_argument('--output-dir', required=True, type=str, help='Output directory')
  #  return parser.parse_args()

def parse_args(sequences, output_format):
    
    sequences = 'C:\\Users\\Tiina.Talts\\Desktop\\EFAR_0037\\8 - GISAID_uploads\\20230325_England_RSV_GISAID.fas'
    path = os.path.dirname(sequences)
    output_format = output_format
    output_dir = os.path.join(path, r'/Annotated_outputs')
    return sequences, output_format, output_dir

def get_reference_sequence(accession):
    from Bio import Entrez
    Entrez.email = "ninni004@hotmail.com"
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    print(f"Fetching reference sequence {accession} from genbank")
    return SeqIO.read(handle, "genbank")

def get_coordinate_map(reference, record):
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as tmp_dir:
        SeqIO.write([reference, record], f"{tmp_dir}/sequences.fasta", "fasta")

        os.system(f"mafft --quiet {tmp_dir}/sequences.fasta > {tmp_dir}/aligned.fasta")
        ref_aligned, qry_aligned = AlignIO.read(f"{tmp_dir}/aligned.fasta", "fasta")

    aln_to_ref = np.cumsum(np.array(ref_aligned.seq) !='-')
    aln_to_qry = np.cumsum(np.array(qry_aligned.seq) !='-')

    ref_to_qry = aln_to_qry[np.array(ref_aligned.seq) !='-'] - 1
    qry_to_ref = aln_to_ref[np.array(qry_aligned.seq) !='-'] - 1

    positions_mapping_to_end_of_query = np.where(ref_to_qry == len(record.seq)-1)[0]
    if len(positions_mapping_to_end_of_query) > 1:
        ref_to_qry[positions_mapping_to_end_of_query[1]:] = len(record.seq)

    positions_mapping_to_end_of_ref = np.where(qry_to_ref == len(reference.seq)-1)[0]
    if len(positions_mapping_to_end_of_ref) > 1:
        qry_to_ref[positions_mapping_to_end_of_ref[1]:] = len(reference.seq)

    return np.concatenate([ref_to_qry, [min(ref_to_qry[-1]+1, len(record.seq))]]), np.concatenate([qry_to_ref, [min(qry_to_ref[-1]+1, len(reference.seq))]])


def lift_location(location, ref_to_qry, query):
    # check if the feature is entirely outside of the query sequence
    if ref_to_qry[location.start]==ref_to_qry[location.end]:
        return None

    from Bio import SeqFeature
    return SeqFeature.FeatureLocation(
            max(0,int(ref_to_qry[location.start])),
            min(int(ref_to_qry[location.end-1]+1), len(query.seq)))


def annotate_sequence(reference, record):
    from Bio import SeqFeature

    ref_to_qry, qry_to_ref = get_coordinate_map(reference, record)

    new_features = []
    # annotate features
    for feature in reference.features:
        new_feat = copy.deepcopy(feature)
        if new_feat.type == 'source':
            new_feat.location = SeqFeature.FeatureLocation(0, len(record.seq))
            new_feat.qualifiers = {}
            for v in ['organism', 'mol_type']:
                if v in feature.qualifiers:
                    new_feat.qualifiers[v] = feature.qualifiers[v]
        else:
            # for compound locations, loop over parts
            if isinstance(feature.location, SeqFeature.CompoundLocation):
                parts = []
                for loc in feature.location.parts:
                    new_loc = lift_location(loc, ref_to_qry, record)
                    if new_loc:
                        parts.append(new_loc)
                if len(parts) > 0:
                    new_feat.location = SeqFeature.CompoundLocation(parts)
                elif len(parts) == 1:
                    new_feat.location = parts[0]
                else: # continue on empty feature
                    continue
            else:
                new_loc = lift_location(feature.location, ref_to_qry, record)
                if new_loc:
                    new_feat.location = new_loc
                else: # continue on empty feature
                    continue

            # exclude short features less that 10 amino acids long that often correspond to alignment problems
            if new_feat.type in  ['CDS', 'gene']:
                if (new_feat.location.end - new_feat.location.start < 30) and (feature.location.end - feature.location.start > 30):
                    print("Skipping feature", new_feat.qualifiers['gene'], "because it is too short, old location", feature.location, "new location", new_feat.location)
                    continue

            if new_feat.type == 'CDS':
                # check of CDS starts before the beginning of the sequence
                if ref_to_qry[feature.location.start]==-1:
                    frame=(feature.location.start-qry_to_ref[new_feat.location.start]) % 3
                    new_feat.qualifiers['codon_start'][0] = f"{int(feature.qualifiers['codon_start'][0]) + frame}"

                new_feat.qualifiers['translation'] = new_feat.translate(record.seq, cds=False)

            if 'locus_tag' in new_feat.qualifiers:
                new_feat.qualifiers.pop('locus_tag')

        new_features.append(new_feat)

    record.features = new_features
    record.annotations = {k:v for k,v in reference.annotations.items() if k in ['organism', 'taxonomy', 'molecule_type']}

    return record

def tbl_entry(start, end, key, qualifiers):
    return f"{start+1}\t{end}\t{key}\t\t\n" + "\n".join([f"\t\t\t{q[0]}\t{q[1]}" for q in qualifiers]) + "\n"

def feature_to_tbl(feat):
    if feat.type == 'source':
        return tbl_entry(feat.location.start, feat.location.end, 'source',
                         [('mol_type', feat.qualifiers['mol_type'][0])])
    elif feat.type == 'CDS':
        return tbl_entry(feat.location.start, feat.location.end, 'CDS',
                        [('product', feat.qualifiers['product'][0]),
                         ('protein_id', feat.qualifiers['protein_id'][0])])
    elif feat.type == 'gene':
        return tbl_entry(feat.location.start, feat.location.end, 'gene',
                        [('gene', feat.qualifiers['gene'][0])])
    else:
        return ''

def write_tbl(record, filename):
    with open(filename, 'w') as f:
        f.write(f">Feature {record.id}\n")
        for feat in record.features:
            f.write(feature_to_tbl(feat))

if __name__=="__main__":
    
    args = parse_args()
    if not os.path.isdir(args.output_dir):
        print("Output directory does not exist!")
        exit()

    if args.virus is None and args.reference is None:
        print("Please specify a reference sequence accession or a virus name.")
        exit()

    if args.virus is not None and args.virus not in reference_sequences:
        print("Virus name not found in reference sequences.")
        exit()

    #ref_id = reference_sequences[args.virus] if args.virus else args.reference

    #reference_sequence = get_reference_sequence(ref_id)

    annotated_sequences = []

    for record in SeqIO.parse(args.sequences, "fasta"):
        ref_id0 = record.id.split('/')
        if ref_id0[1] == 'A':
          ref_id = "rsv-a"
        else:
          ref_id = "rsv-b"
        reference_sequence = get_reference_sequence(ref_id)
        annotated_sequences.append(annotate_sequence(reference_sequence, record))
        print(f"Annotated {record.id}")
        file_base_name = record.id.replace('/', '_')
        if args.output_format == 'genbank':
            SeqIO.write(record, f"{args.output_dir}/{file_base_name}.gb", "genbank")
        elif args.output_format == 'gff3':
            from BCBio.GFF import GFF3Writer
            gffwriter = GFF3Writer()
            with open(f"{args.output_dir}/{file_base_name}.gff3", "w") as handle:
                gffwriter.write([record], handle)
        elif args.output_format == 'tbl':
            write_tbl(record, f"{args.output_dir}/{file_base_name}.tbl")

