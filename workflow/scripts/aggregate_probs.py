import pandas as pd
from pathlib import Path


def make_wig_path(rlooper_dir, plasmid_name):
    return Path(rlooper_dir).joinpath(plasmid_name + '_bpprob.wig')


def read_wig(filepath):
    # pass
    df = pd.read_csv(str(filepath), skiprows=4, sep='\t')
    df.columns = ['prob']

    # pull out a value (last dir name)
    a = int(filepath.parent.parent.name.split('~')[-1])
    df['a'] = a
    df['position'] = list(range(0, len(df)))

    return df


def read_all_wigs(rlooper_dirs, plasmid_name):

    dfs = []
    for each_dir in rlooper_dirs:
        wig_path = make_wig_path(each_dir, plasmid_name)
        df = read_wig(wig_path)
        df['plasmid'] = plasmid_name
        dfs.append(df)
    
    return pd.concat(dfs)

def main():


    rlooper_dirs = snakemake.input
    plasmid_name = snakemake.params['plasmid']
    df = read_all_wigs(rlooper_dirs, plasmid_name)
    df['orrientation'] = snakemake.params['orrientation']
    df.to_csv(snakemake.output[0], sep='\t', index=False)


if __name__ == '__main__':
    main()

