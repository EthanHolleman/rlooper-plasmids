

rule calc_max_probs:
    conda:
        '../envs/py.yml'
    resources:
        mem_mb=20000
    input:
        'output/agg_probs/aggregated.rlooper.probs.all.runs.tsv'
    output:
        'output/agg_probs/aggregated.rlooper.max.probs.tsv'
    script:
        '../scripts/calc_max_prob.py'


rule plot_sigma_vals:
    conda:
        '../envs/R.yml'
    resources:
        mem_mb=20000
    input:
        'output/agg_probs/aggregated.rlooper.max.probs.tsv'
    output:
        'output/plots/sigma.vs.rloop_prob.pdf'
    script:
        '../scripts/plot_max_prob.R'
