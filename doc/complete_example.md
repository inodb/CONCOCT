Complete Example V0.3
=====================
This documentation page aims to be a complete example walk through for the usage of the CONCOCT package version 0.3.
It assumes you have successfully gone through the installation description found in the README. 

Required software
----------------------
To run the entire example you need to install all dependencies as stated in the [README dependencies](../README.md#dependencies). This includes all the optional dependencies. You can also look at [doc/Dockerfile](Dockerfile) to help you install these packages on your server.

Another way to get everything set up is to use our full Docker image (binnisb/concoct_0.3.0) as suggested in the [README docker](../README.md#using-docker).

It is not required to run all steps. The output files for each step are in the test data repository. At the end of this example the results should be the same as the results in the corresponding test data repository: https://github.com/BinPro/CONCOCT-test-data/releases. The version numbers listed above are the ones used to generate the results in that repository. Using newer versions will probably not be a problem, but your results may be different in that case.

Downloading test data
-----------------------
First download the test data repository of CONCOCT corresponding to the version of CONCOCT, that you have installed. The test data repository can be downloaded [here](https://github.com/BinPro/CONCOCT-test-data/releases). Then extract it in a suitable location.

If you are running the current unstable master branch of concoct, you need to clone the latest version of the test-data-repository as well.

Setting up the test environment
-------------------------------
###Using Docker###
On your HOST machine create the following folder structure below (Data, Data/CONCOCT-complete-example, Data/CONCOCT-test-data):

    mkdir -p /HOST/path/to/Data
    mkdir /HOST/path/to/Data/CONCOCT-complete-example
    
Move the test data that was downloaded and extracted (CONCOCT-test-data) to the Data folder
    
    # Move the test data you extracted in the download part into the Data folder
    mv /HOST/extracted/test/data/CONCOCT-test-data /HOST/path/to/Data/CONCOCT-test-data

Now you want to execute the following command to log into our Docker image and to map the ```/HOST/path/to/Data``` to your image and the Data folder will be accessable in /opt/Data:

    sudo docker run -v /HOST/path/to/Data:/opt/Data/ -i -t binnisb/concoct_0.3.0 bash

This will download the 2G image to your machine and then leaves you in a BASH shell. In the Docker imgage, the following environmental variables have been set. So if you have your folders set up differently in the steps above you need to alter these variables accordingly:

    CONCOCT=/opt/CONCOCT-0.3.0
    CONCOCT_TEST=/opt/Data/CONCOCT-test-data
    CONCOCT_EXAMPLE=/opt/Data/CONCOCT-complete-example

###Your own setup###
After obtaining the test data, create a folder where you want all the output from this example to go:

    mkdir CONCOCT-complete-example
    cd CONCOCT-complete-example

Set three variables with full paths. One pointing to the root directory of the ```CONCOCT``` software, one pointing to the test data repository, named ```CONCOCT_TEST``` and one to the directory we just created. If you now have these in the folder ```/home/username/src/```, for instance, then use:

    CONCOCT=/home/username/src/CONCOCT
    CONCOCT_TEST=/home/username/src/CONCOCT-test-data
    CONCOCT_EXAMPLE=/home/username/CONCOCT-complete-example

You can see the full path of a directory you are located in by running the command ```pwd```.

Assembling Metagenomic Reads
----------------------------
The first step in the analysis is to assemble all reads into contigs, here we use the software [Velvet](http://www.ebi.ac.uk/~zerbino/velvet/) for this. This step can be computationaly intensive but for this small data set comprising a synthetic community of four species and 16 samples (100,000 reads per sample) it can be performed in a few minutes. If you do not wish to execute this step, the resulting contigs are already in the test data repository, and you can copy them from there insted. The commands for running Velvet are:

    cd $CONCOCT_EXAMPLE
    cat $CONCOCT_TEST/reads/Sample*_R1.fa > All_R1.fa
    cat $CONCOCT_TEST/reads/Sample*_R2.fa > All_R2.fa
    velveth velveth_k71 71 -fasta -shortPaired -separate All_R1.fa All_R2.fa
    velvetg velveth_k71 -ins_length 400 -exp_cov auto -cov_cutoff auto	

After the assembly is finished create a directory with the resulting contigs and copy the result of Velvet there (this output is also in ```$CONCOCT_TEST/contigs```):

    mkdir contigs
    cp velveth_k71/contigs.fa contigs/velvet_71.fa
    rm All_R1.fa
    rm All_R2.fa

Map the Reads onto the Contigs
------------------------------
After assembly we map the reads of each sample back to the assembly using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and remove PCR duplicates with [MarkDuplicates](http://picard.sourceforge.net/command-line-overview.shtml#MarkDuplicates). The coverage histogram for each bam file is computed with [BEDTools](https://github.com/arq5x/bedtools2) genomeCoverageBed. The script that calls these programs is provided with CONCOCT. 

If you are not using the Docker image, then one does need to set an environment variable with the full path to the MarkDuplicates jar file. ```$MRKDUP``` which should point to the MarkDuplicates jar file e.g.

    #NOTE not necessary if using the Docker image
    export MRKDUP=/home/username/src/picard-tools-1.77/MarkDuplicates.jar

It is typically located within your picard-tools installation.

The following command is to be executed in the ```$CONCOCT_EXAMPLE``` dir you created in the previous part. First create the index on the assembly for bowtie2:

    cd $CONCOCT_EXAMPLE
    bowtie2-build contigs/velvet_71.fa contigs/velvet_71.fa
    
Then run this for loop, which for each sample creates a folder and runs ```map-bowtie2-markduplicates.sh```:

    for f in $CONCOCT_TEST/reads/*_R1.fa; do
        mkdir -p map/$(basename $f);
        cd map/$(basename $f);
        bash $CONCOCT/scripts/map-bowtie2-markduplicates.sh -ct 1 -p '-f' $f $(echo $f | sed s/R1/R2/) pair $CONCOCT_EXAMPLE/contigs/velvet_71.fa asm bowtie2;
        cd ../..;
    done

The parameters used for `map-bowtie2-markduplicates.sh` are:

* `-c` option to compute coverage histogram with genomeCoverageBed
* `-t` option is number of threads
* `-p` option is the extra parameters given to bowtie2. In this case `-f`.

The five arguments are:
* pair1, the fasta/fastq file with the #1 mates
* pair2, the fasta/fastq file with the #2 mates
* pair_name, a name for the pair used to prefix output files
* assembly, a fasta file of the assembly to map the pairs to
* assembly_name, a name for the assembly, used to postfix outputfiles
* outputfolder, the output files will end up in this folder

Generate coverage table
------------------------
Use the bam files of each sample to create a table with the coverage of each contig per sample.

    cd $CONCOCT_EXAMPLE/map
    python $CONCOCT/scripts/gen_input_table.py --isbedfiles \
        --samplenames <(for s in Sample*; do echo $s | cut -d'_' -f1; done) \
        ../contigs/velvet_71.fa */bowtie2/asm_pair-smds.coverage \
    > concoct_inputtable.tsv
    mkdir $CONCOCT_EXAMPLE/concoct-input
    mv concoct_inputtable.tsv $CONCOCT_EXAMPLE/concoct-input/

Generate linkage table
------------------------
The same bam files can be used to give linkage per sample between contigs:

    cd $CONCOCT_EXAMPLE/map
    python $CONCOCT/scripts/bam_to_linkage.py -m 8 \
        --regionlength 500 --fullsearch \
        --samplenames <(for s in Sample*; do echo $s | cut -d'_' -f1; done) \
        ../contigs/velvet_71.fa Sample*/bowtie2/asm_pair-smds.bam \
    > concoct_linkage.tsv
    mv concoct_linkage.tsv $CONCOCT_EXAMPLE/concoct-input/
    

Run concoct
-----------

To see possible parameter settings with a description run

    $CONCOCT/bin/concoct --help

We will only run concoct for some standard settings here. First we need to parse the input table to just contain the mean coverage for each contig in each sample:

    cd $CONCOCT_EXAMPLE
    cut -f1,3-26 concoct-input/concoct_inputtable.tsv > concoct-input/concoct_inputtableR.tsv

Then run concoct with 40 as the maximum number of cluster `-c 40`, that we guess is appropriate for this data set:

    cd $CONCOCT_EXAMPLE
    concoct -c 40 --coverage_file concoct-input/concoct_inputtableR.tsv --composition_file contigs/velvet_71.fa -b concoct-output/

When concoct has finished the message "CONCOCT Finished, the log shows how it went." is piped to stdout. The program generates a number of files in the output directory that can be set with the `-b` parameter and will be the present working directory by default. 

Evaluate output
---------------

This will require that you have Rscript with the R packages [gplots](http://cran.r-project.org/web/packages/gplots/index.html), [reshape](http://cran.r-project.org/web/packages/reshape/index.html), [ggplot2](http://cran.r-project.org/web/packages/ggplot2/index.html), [ellipse](http://cran.r-project.org/web/packages/ellipse/index.html), [getopt](http://cran.r-project.org/web/packages/getopt/index.html) and [grid](http://cran.r-project.org/web/packages/grid/index.html) installed.

First we can visualise the clusters in the first two PCA dimensions:

    cd $CONCOCT_EXAMPLE
    mkdir evaluation-output
    Rscript $CONCOCT/scripts/ClusterPlot.R -c concoct-output/clustering_gt1000.csv -p concoct-output/PCA_transformed_data_gt1000.csv -m concoct-output/pca_means_gt1000.csv -r concoct-output/pca_variances_gt1000_dim -l -o evaluation-output/ClusterPlot.pdf

<https://github.com/BinPro/CONCOCT-test-data/tree/master/evaluation-output/ClusterPlot.pdf>

We can also compare the clustering to species labels. For this test data set we know these labels, they are given in the file ```clustering_gt1000_s.csv```. For real data labels may be obtained through taxonomic classification, e.g. using:

<https://github.com/umerijaz/TAXAassign>

In either case we provide a script Validate.pl for computing basic metrics on the cluster quality:

    cd $CONCOCT_EXAMPLE
    cp $CONCOCT_TEST/evaluation-output/clustering_gt1000_s.csv evaluation-output/
    $CONCOCT/scripts/Validate.pl --cfile=concoct-output/clustering_gt1000.csv --sfile=evaluation-output/clustering_gt1000_s.csv --ofile=evaluation-output/clustering_gt1000_conf.csv

This script requires the clustering output by concoct ```concoct-output/clustering_gt1000.csv``` these have a simple format of a comma separated file listing each contig id followed by the cluster index and the species labels that have the same format but with a text label rather than a cluster index. The script should output:

    N	M	S	K	Rec.	Prec.	NMI	Rand	AdjRand
    88	88	4	4	0.920455	0.988636	0.757418	0.871473	0.695022

This gives the no. of contigs N clustered, the number with labels M, the number of unique labels S, the number of clusters K, the recall, the precision, the normalised mutual information (NMI), the Rand index, and the adjusted Rand index. It also generates a file called a `confusion matrix` with the frequencies of each species in each cluster. We provide a further script for visualising this as a heatmap:

    $CONCOCT/scripts/ConfPlot.R  -c evaluation-output/clustering_gt1000_conf.csv -o  evaluation-output/clustering_gt1000_conf.pdf

This generates a file with normalised frequencies of contigs from each cluster across species:

<https://github.com/BinPro/CONCOCT-test-data/tree/master/evaluation-output/clustering_gt1000_conf.pdf>

Validation using single-copy core genes
---------------------------------------

We can also evaluate the clustering based on single-copy core genes. You first need to find genes on the contigs and functionally annotate these. Here we used prokka (http://www.vicbioinformatics.com/software.prokka.shtml) for gene prediction and annotation, but you can also use for example prodigal. The corresponding protein sequences are here:

    $CONCOCT_TEST/annotations/proteins/velvet_71.faa

and GFF3 file:

    $CONCOCT_TEST/annotations/proteins/velvet_71.gff
    
And we used RPS-Blast to COG annotate the protein sequences using (PROKKA_RPSBLAST.sh). With the following command on eight cores:

    $CONCOCT/scripts/PROKKA_RPSBLAST.sh -f annotations/proteins/velvet_71.faa -p -c 8 -r 1

To run this yourself the file ```velvet_71.faa``` will have to be copied into the test directory i.e.

    mkdir $CONCOCT_EXAMPLE/annotations
    mkdir $CONCOCT_EXAMPLE/annotations/proteins
    mkdir $CONCOCT_EXAMPLE/annotations/cog-annotations
    cp $CONCOCT_TEST/annotations/proteins/* $CONCOCT_EXAMPLE/annotations/proteins/
The blast output has been placed in:

    $CONCOCT_TEST/annotations/cog-annotations/velvet_71.out
    
Finally, we filtered for COGs representing a majority of the subject to ensure fragmented genes are not over-counted and generated a table of counts of single-copy core genes in each cluster generated by CONCOCT. Remember to use a real email adress, this is supplied since information is fetched from ncbi using their service eutils, and the email is required to let them know who you are.

    cd $CONCOCT_EXAMPLE
    $CONCOCT/scripts/COG_table.py -g annotations/proteins/velvet_71.gff -b annotations/cog-annotations/velvet_71.out -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c concoct-output/clustering_gt1000.csv -e mail@example.com > evaluation-output/clustering_gt1000_scg.tab

The script requires the clustering output by concoct ```concoct-output/clustering_gt1000.csv``` and a file listing a set of SCGs (e.g. a set of COG ids) to use ```scgs/scg_cogs_min0.97_max1.03_unique_genera.txt```.

Since these protein sequences are generated by Prokka, the names of the contig ids need to be recovered from the gff file. If prodigal would have been used, the contig ids would instead have been recovered from the protein ids using a separator character, in which case only the string before (the last instance of) the separator will be used as contig id in the annotation file. In the case of prodigal the separator that should be used is _ and this is the default value, but other characters can be given through the '--separator' argument.

The output file is a tab-separated file with basic information about the clusters (cluster id, ids of contigs in cluster and number of contigs in cluster) in the first three columns, and counts of the different SCGs in the following columns.

This can also be visualised graphically using the R script:

    cd $CONCOCT_EXAMPLE
    $CONCOCT/scripts/COGPlot.R -s evaluation-output/clustering_gt1000_scg.tab -o evaluation-output/clustering_gt1000_scg.pdf

The plot is downloadable here:

<https://github.com/BinPro/CONCOCT-test-data/tree/master/evaluation-output/clustering_gt1000_scg.pdf>

Incorporating linkage information
---------------------------------

To perform a hierarchical clustering of the clusters based on linkage we simply run:

    $CONCOCT/scripts/ClusterLinkNOverlap.pl --cfile=concoct-output/clustering_gt1000.csv --lfile=concoct-input/concoct_linkage.tsv --covfile=concoct-input/concoct_inputtableR.tsv --ofile=concoct-output/clustering_gt1000_l.csv

The output indicates that the clusters have been reduced from four to three. The new clustering is given by ```concoct-output/clustering_gt1000_l.csv```. This is a significant improvement in recall:

    $CONCOCT/scripts/Validate.pl --cfile=concoct-output/clustering_gt1000_l.csv --sfile=evaluation-output/clustering_gt1000_s.csv --ofile=evaluation-output/clustering_gt1000_conf.csv

    N	M	S	K	Rec.	Prec.	NMI	Rand	AdjRand
    88	88	4	3	1.000000	0.988636	0.976458	0.999478	0.998515
