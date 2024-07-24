
rule run_rlooper:
    input:
        rlooper_exe = 'submodules/rlooper/bin/rlooper',
        fasta=lambda w: seq_paths[w.plasmid]
    output:
        directory(f'output/rlooper_runs/fwd/{{plasmid}}/{paramspace.wildcard_pattern}')
    params:
        sigma = lambda w: paramspace.instance(w)['sigma'],
        N = lambda w: paramspace.instance(w)['N'],
        a = lambda w: paramspace.instance(w)['a'],
        plasmid = lambda w: w.plasmid

    shell:'''
    mkdir {output}
    {input.rlooper_exe} {input.fasta} {output}/{params.plasmid} --sigma {params.sigma} --N {params.N} --a {params.a} --localaverageenergy
    '''

rule run_rlooper_rc:
    input:
        rlooper_exe = 'submodules/rlooper/bin/rlooper',
        fasta='output/reverseComplements/{plasmid}.rc.fa'
    output:
        directory(f'output/rlooper_runs/rc/{{plasmid}}/{paramspace.wildcard_pattern}')
    params:
        sigma = lambda w: paramspace.instance(w)['sigma'],
        N = lambda w: paramspace.instance(w)['N'],
        a = lambda w: paramspace.instance(w)['a'],
        plasmid = lambda w: w.plasmid

    shell:'''
    mkdir {output}
    {input.rlooper_exe} {input.fasta} {output}/{params.plasmid} --sigma {params.sigma} --N {params.N} --a {params.a} --localaverageenergy
    '''
    

rule wig_to_long:
    conda:
        '../envs/py.yml'
    input:
        f'output/rlooper_runs/{{orrientation}}/{{plasmid}}/{paramspace.wildcard_pattern}', 
    output:
        f'output/rlooper_probs/{{orrientation}}/{{plasmid}}/{{plasmid}}.{paramspace.wildcard_pattern}.tsv'
    params:
        plasmid = lambda wildcards: wildcards.plasmid,
        orrientation = lambda wildcards: wildcards.orrientation,
        sigma = lambda w: paramspace.instance(w)['sigma'],
        N = lambda w: paramspace.instance(w)['N'],
        a = lambda w: paramspace.instance(w)['a'],
    script:'../scripts/wig_to_long.py'


rule agg_prob_files:
    conda:
        '../envs/py.yml'
    input:
        expand(
            'output/rlooper_probs/{orrientation}/{plasmid}/{plasmid}.{params}.tsv',
            params=paramspace.instance_patterns,
            plasmid=sequence_names,  orrientation=['fwd']
        ),
        expand(
            'output/rlooper_probs/{orrientation}/{plasmid}/{plasmid}.{params}.tsv',
            params=paramspace.instance_patterns,
            plasmid=sequence_names,  orrientation=['rc']
        )
    output:
        'output/agg_probs/aggregated.rlooper.probs.all.runs.tsv'
    script:'../scripts/concat_tsv.py'


rule plot_probs:
    conda:
        '../envs/R.yml'
    input:
        'output/aggProbs/{plasmid}/{plasmid}.agg.probs.a.params.tsv'
    output:
        'output/plotProbs/{plasmid}/{plasmid}.agg.probs.a.params.png'
    script:'../scripts/plot_a_probs.R'

    

