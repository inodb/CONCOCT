#!/usr/bin/env python
"""Generate table of contigs clusters to POG hits given blastp output file and
target to POG mapping. Resulting table looks similar to COG_table.py e.g.:

Cluster Contigs Num_contigs POG1234 POG0005 ETC..

Only clusters with hits are listed. Only POG columns that have a hit are
included."""""
import argparse
import pandas as pd
import sys


def get_id(sseqid):
    """Get only the gid of a subject sseqid in case it has one"""
    if "gi" in sseqid:
        return sseqid.split("|")[1]
    else:
        return sseqid


def create_pog_row(protein, pro_pogdf):
    """Create a Series with nr of columns equal to the amount of POGs with one
    row representing the POGs the given protein belongs to. If it belongs to
    the POG the value is True, otherwise False."""
    s = pd.Series([False] * len(set(pro_pogdf.pog)), index=sorted(list(set(pro_pogdf.pog))))
    s[set(pro_pogdf.pog[pro_pogdf.protein == protein])] = True
    return s


def count_pogs_per_cluster(blastoutfile, cluster_file, protein_pog_file,
        separator):
    """Count POGs per cluster and print result"""
    # parse blast output
    blastdf = pd.read_csv(blastoutfile, sep="\t",
        names=["qseqid", "sseqid" , "evalue", "pident", "score", "qstart",
            "qend", "sstart", "send", "length", "slen"],
        converters={"sseqid": get_id, })

    # A certain fraction of the pog should be covered to avoid the same pog
    # to be counted twice in the case when a pog is split across two or more contigs.
    #aln_reqs = ((blastdf.send - blastdf.sstart + 1) / blastdf.slen * 100.0 >=
    #        scovs_threshold) & (blastdf.pident >= pident_threshold)

    # parse protein id to pog mapping
    pro_pogdf = pd.read_csv(protein_pog_file, sep="\t",
            names=["protein", "pog"])

    # assign each protein in the blast output to one or more POGs
    blastdf.pogs = blastdf.sseqid.apply(create_pog_row, args=(pro_pogdf,))
    # get the contig name of the protein
    blastdf.contig = blastdf.qseqid.apply(lambda x: x.split(separator)[0])

    # Load clustering
    clusterdf = pd.read_csv(cluster_file, sep=",", names=["contig", "cluster"])

    # get cluster to pog mappings
    clusterdf.pogs = clusterdf.groupby('cluster').apply(
        lambda x: blastdf.pogs[blastdf.contig.isin(x.contig)].sum())
    # select only cluster rows that have a pog hit
    clu_pog_sub = clusterdf.pogs[clusterdf.pogs.sum(axis=1) > 0]
    # select only pog columns that have a hit
    clu_pog_sub2 = clu_pog_sub.ix[:, [i for i in
        range(len(clu_pog_sub.columns)) if sum(clu_pog_sub.ix[:, [i]].sum()) > 0]]

    # Keep format consistent with COG_table.py
    # add all contig names of cluster
    clu_pog_sub2.insert(0, "Contigs", clu_pog_sub2.apply(lambda x:
        "|".join(clusterdf.contig[clusterdf.cluster == x.name]), axis=1))
    # add number of contigs
    clu_pog_sub2.insert(1, "Num_contigs", clu_pog_sub2.apply(lambda x:
        len(clusterdf.contig[clusterdf.cluster == x.name]), axis=1))
    clu_pog_sub2.index.name = "Cluster"

    sys.stdout.write(clu_pog_sub2.to_csv(sep="\t"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-b', '--blastoutfile', required=True,
            help=('Output of blastp run, assumed to be in tabular format whith '
                'columns: qseqid sseqid evalue pident score qstart qend sstart send length slen. '
                'The contigs ids are assumed to be recoverable by removing the last underscore '
                'and the characters following it from the qseqid column.'))
    parser.add_argument('-c', '--cluster_file', required=True,
            help=('Clusters file from concoct execution.'))
    parser.add_argument('--protein_pog_file',
            help=('Supply a protein to pog mapping file in a tsv format'))
    parser.add_argument('--separator', default="_",
            help=('Character that is used to separate the contig id from the '
                    'protein identifier. Everything before the last occurence '
                    'of this character will be used as the contig id. Default '
                    'value is "_"'))
    #parser.add_argument('-s', '--scovs-threshold', type=float, default=50.0,
    #        help='Threshold covered in percent, default=50.0')
    #parser.add_argument('-p', '--pident-threshold', type=float, default=0.0,
    #        help='Threshold identity in percent, default=0.0')
    args = parser.parse_args()
    count_pogs_per_cluster(args.blastoutfile, args.cluster_file,
            args.protein_pog_file, args.separator)
