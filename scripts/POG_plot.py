#!/usr/bin/env python
"""Generate HTML plot from table of contig clusters to POG hits as generated by
POG_table.py"""""
import argparse
import pandas as pd
import seaborn as sns

from bokeh.charts import CategoricalHeatMap


def plot_pogs(clu_pog_table, output_html, title):
    clu_pog = pd.read_csv(clu_pog_table, sep="\t")

    # get only the pog counts of the table
    pogs = clu_pog.ix[:, range(3, len(clu_pog.columns))]
    pogs.set_index(pogs.index.astype(str), inplace=True)
    pogs.columns = [str(int(c[3:])) for c in pogs.columns]

    # generate enough colors using seaborn to represent all the different counts of the POGs
    # TODO: should maybe use a more common library
    chpal = sns.cubehelix_palette(pogs.max().max() + 1, start=2, rot=0, dark=0, light=.95)
    palette = ['#%02x%02x%02x' % tuple([round(255 * rgbi) for rgbi in rgb]) for rgb in chpal]

    hm = CategoricalHeatMap(pogs.transpose(), palette=palette,
            xlabel="Cluster", ylabel="POG", title=title, filename=output_html)
    hm.width(30 * len(pogs)).height(20 * len(pogs.columns)).show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-c', '--clu_pog_table', required=True,
            help=('Cluster to POG table as generated by POG_table.py'))
    parser.add_argument('-o', '--output_html', required=True,
            help=('HTML plot output file.'))
    parser.add_argument('--title',
            help=('default is POG hits per cluster'), default='POG hits per cluster')
    args = parser.parse_args()
    plot_pogs(args.clu_pog_table, args.output_html, args.title)