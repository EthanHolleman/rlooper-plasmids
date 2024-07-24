from pathlib import Path
import pandas as pd


def make_wig_path(rlooper_dir, plasmid_name):
    return Path(rlooper_dir).joinpath(plasmid_name + "_bpprob.wig")


def read_wig(filepath):
    df = pd.read_csv(str(filepath), skiprows=4, sep="\t")
    df.columns = ['rloop_prob']
    df["position"] = list(range(0, len(df)))

    return df

def make_energy_path(rlooper_dir, plasmid_name):
    return Path(rlooper_dir).joinpath(plasmid_name + '_avgG.wig')


def add_parameter_columns(wig_df, **kwargs):
    for each_param in kwargs:
        wig_df[each_param] = kwargs[each_param]

    return wig_df


def main():
    rlooper_dir = snakemake.input[0]
    
    wig_path = make_wig_path(rlooper_dir, snakemake.params['plasmid'])
    energy_path = make_energy_path(rlooper_dir, snakemake.params['plasmid'])
    
    wig_df = read_wig(wig_path)
    energy_df =read_wig(energy_path)
    
    wig_df['avgG'] = energy_df['rloop_prob']  # called rloop prob but actually energy
    
    wig_df_params = add_parameter_columns(wig_df, **dict(snakemake.params))

    wig_df_params.to_csv(snakemake.output[0], sep="\t", index=False)


if __name__ == "__main__":
    main()
