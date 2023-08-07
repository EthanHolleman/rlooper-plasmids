# One aspect of R-looper behavior I would like to test is how providing more or less
# of a region of interest changes the predictions about that region. Towards this,
# this script takes in a fasta file, a center point and a step distance and outputs 
# a number of fasta files that include slices of the provided sequence expanding outwards
# from the center point at the rate of the step parameter. This script will also
# ensure that the sliced fasta files have R-looper ready fasta headers (assuming that the
# input fasta file has a valid R-looper header). 


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re


def find_replace_capture_group(r, sub, text):
    match = re.search(r, text)
    if not match:
        raise Exception('No regex match found you suck lol loser!!!!')
    start, end = match.span(1)  # assume only 1 capture group
    return text[:start] + sub + text[end:]


def slice(seq, center, step, min_size=300):
    r = 'range=.+:\\d-(\\d+)'
    slices = []
    outside_range, num_steps = False, 0


    while outside_range == False:
        dist = step*num_steps
        left_cord, right_cord = center - dist, center + dist
        seq_str = str(seq.seq)

        if left_cord < 0:  # reached end of seq on left side
            left_cord = 0 
        if right_cord >= len(seq_str):  # reached on right
            right_cord = len(seq_str) - 1
        if left_cord == 0 and right_cord == len(seq_str)-1:  # at the end of both sides
            slice = seq_str
            outside_range = True
        
        slice = seq_str[left_cord:right_cord]

        print(left_cord, right_cord)
        
        slice_id = find_replace_capture_group(
            r, str(len(slice)), seq.description
        )

        slices.append(
            SeqRecord(str(slice), slice_id)
        )
        num_steps += 1
    
    return slices



def write_slices():
    # write slices to invidiviual fasta files as this is format
    # rlooper is expecting
    pass

    

def main():

    fa = SeqIO.read(
        '/home/ethollem/workflows/rlooper-runs/resources/sequences/vasquezLabSeqs/pfc8_genbank_full_coding_gene_start.fa',
        'fasta')
    slices = slice(
            fa, int(len(fa.seq) / 2),   # for now always use midpoint of the provided sequence
           snakemake.params['step'],
    )
    SeqIO.write(snakemake.output[0], slices, 'fasta')


if __name__ == '__main__':
    main()