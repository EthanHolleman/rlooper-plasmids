
rule complie_rlooper:
    output:
        'submodules/rlooper/bin/rlooper'
    shell:'''
    cd submodules/rlooper
    make all
    '''




