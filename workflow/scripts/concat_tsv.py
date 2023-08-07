import pandas as pd


def main():

    tsvs = snakemake.input
    dfs = [
        pd.read_csv(tsv, sep='\t') for tsv in tsvs
    ]
    concat_df = pd.concat(dfs)

    concat_df.to_csv(snakemake.output[0], sep='\t', index=False)


if __name__ == '__main__':
    main()