# This script modifies VR sequence genbank file in two main ways. First it
# uses the pacbio barcoding primers binding sites to define the sequence used
# (the whole plasmid is not used just the amplicon) and secondly it
# converts the amplicon sequence into a fasta file with an rlooper safe
# header. This script is not part of the workflow and I just ran it once
# on the VR genbank files to create the amplicon files.

from Bio import SeqIO
from Bio.SeqUtils import nt_search
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from rotate_genbank import rotate


# Primer binding sites, just hardcoded them in for simplicty they should
# not change anyway. These are not the complete primers, just the regions
# which actually have homology to the plasmid (no barcodes)
FWD_PRIMER = 'AAGCACTAAATCGGAACCCTAAAG'
REV_PRIMER= 'GGCGAACTACTTACTCTAGCTTCC'

SEQ_PATH = 'resources/VRs'
OUTPUT_PATH = 'resources/sequences'



def read_gb_files(gb_dir):
    records = []
    for each_record in Path(gb_dir).iterdir():
        if each_record.suffix == '.gb':
            records.append(
                SeqIO.read(each_record, 'genbank')
            )
    return records


def locate_primer_binding_sites(sequence):
    # find forward primer
    fwd_site = str(sequence).find(FWD_PRIMER)
    # find reverse primer
    rev_site = str(sequence).find(str(Seq(REV_PRIMER).reverse_complement()))

    # return the indexes where the sites were found
    return fwd_site, rev_site


def truncate_gb_to_amplicon(gb_record, fwd_site, rev_site):


    def make_rlooper_safe_header():
        # get just the name of the VR
        name = gb_record.name.split('_')[-1].split('.')[0]
        template = f"{name} range={name}FIXED:1-{rev_site-fwd_site} 5'pad=0 3'pad=0 strand=+ repeatMasking=none"
        return template


    amplicon_seq = gb_record.seq[fwd_site:rev_site]
    return SeqRecord(
        seq=amplicon_seq,
        id=make_rlooper_safe_header(),
        description=''
    )


def convert_gb_records_to_amplicons(gb_records):

    amplicons = []

    for each_record in gb_records:
        sites = locate_primer_binding_sites(each_record.seq)
        print(sites)
        amplicon = truncate_gb_to_amplicon(each_record, *sites)
        amplicons.append(amplicon)
    
    return amplicons



def write_amplicon_seqs(output_path, amplicon_records):
    for each_record in amplicon_records:
        filename = f"{each_record.id.split(' ')[0]}.amplicon.fa"
        filepath = Path(output_path).joinpath(filename)
        SeqIO.write(
            each_record, str(filepath), 'fasta'
        )
    
    return 0


def main():


    gb_records = read_gb_files(SEQ_PATH)
    amplicon_records = convert_gb_records_to_amplicons(gb_records)
    write_amplicon_seqs(OUTPUT_PATH, amplicon_records)



if __name__ == '__main__':
    main()
