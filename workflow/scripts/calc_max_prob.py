import pandas as pd


def add_integer_name(row):
    return int(row.plasmid.split('-')[-1])


def main():
    
    df_path = snakemake.input[0]
    agg_df = pd.read_csv(df_path, sep='\t')
    agg_df['plot_name'] = agg_df.apply(lambda row: add_integer_name(row), axis=1)
    groups = agg_df.groupby(
        ['plasmid', 'orrientation', 'sigma', 'N', 'a', 'plot_name']
        ).max().reset_index()

    groups.to_csv(snakemake.output[0], sep='\t', index=None)


if __name__ == '__main__':
    main()