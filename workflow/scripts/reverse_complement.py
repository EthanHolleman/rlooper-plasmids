# Input fasta file with Rlooper formatted header and output a fasta
# with the same header but with the sequence reverse complemented.

from Bio import SeqIO


def main():
    print(snakemake.input['fasta'])
    record = SeqIO.read(snakemake.input['fasta'], 'fasta')
    rc = record.seq.reverse_complement()
    record.seq = rc
    SeqIO.write([record], snakemake.output[0], 'fasta')


if __name__ == '__main__':
    main()