

rule reverse_complement_seqs:
    conda:
        '../envs/py.yml'
    input:
        fasta=lambda w: seq_paths[w.plasmid]
    output:
        'output/reverseComplements/{plasmid}.rc.fa'
    script:'../scripts/reverse_complement.py'
