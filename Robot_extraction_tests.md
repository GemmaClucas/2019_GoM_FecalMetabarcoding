Robot extraction tests, June 2022
================
Gemma Clucas
6/21/2022

## 1. Import the data into Qiime2

Load qiime environment and cd to correct directory.

    cd /Users/gemmaclucas/GitHub/Fecal_metabarcoding/2019_GoM_FecalMetabarcoding/Robot_extraction_test
    conda activate qiime2-2021.4

The raw reads are saved on my solid state hard drive, Data_SS1, in
`/Volumes/Data_SS1/MiFish/Robot_tests/`. There is the robot extraction
plate, and then the first half of plate 12 and the second half of plate
14 in that folder.

As I import them I am just going to name them plate 1, 2, 3 to make it
easier to loop through them in subsequent steps.

    qiime tools import\
      --type 'SampleData[PairedEndSequencesWithQuality]'\
      --input-path /Volumes/Data_SS1/MiFish/Robot_tests/June2022_robotplate/reads \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt\
      --output-path demux_plate1.qza
     
    qiime tools import\
      --type 'SampleData[PairedEndSequencesWithQuality]'\
      --input-path /Volumes/Data_SS1/MiFish/Robot_tests/Plate12_forcomparison/ \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt\
      --output-path demux_plate2.qza
      
    qiime tools import\
      --type 'SampleData[PairedEndSequencesWithQuality]'\
      --input-path /Volumes/Data_SS1/MiFish/Robot_tests/Plate14_forcomparison/ \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt\
      --output-path demux_plate3.qza

## 2. Trim primers with cutadapt plugin

F primer: GTCGGTAAAACTCGTGCCAGC (21 bp)  
R primer: CATAGTGGGGTATCTAATCCCAGTTTG (27 bp)

### Trim 3’ ends first

At the 3’ end of the read, the primer will have been read through after
reading the MiFish sequence. I need to be looking for the reverse
complement of the reverse primer in R1 (`—p-adapter-f`) and the reverse
complement of the forward primer in R2 (`—p-adapter-r`)

F primer reverse complement: GCTGGCACGAGTTTTACCGAC  
R primer reverse complement: CAAACTGGGATTAGATACCCCACTATG

    for K in {1..3}; do
      qiime cutadapt trim-paired \
        --i-demultiplexed-sequences demux_Plate$K.qza \
        --p-adapter-f CAAACTGGGATTAGATACCCCACTATG \
        --p-adapter-r GCTGGCACGAGTTTTACCGAC \
        --o-trimmed-sequences trimd_Plate$K.qza \
        --verbose > cutadapt_out_Plate$K.txt
    done

To see how much data passed the filter for each sample:

    grep "Total written (filtered):" cutadapt_out_Plate$K.txt 

About 77% is passing the filters for the robot plate and the second half
of plate 14. Plate 12 is much more variable.

### Trim 5’ ends of reads

All R1 should begin with the forward primer: GTCGGTAAAACTCGTGCCAGC (21
bases).  
All R2 should begin with the reverse primer: CATAGTGGGGTATCTAATCCCAGTTTG
(27 bases).

Trim these with the following commands:

    for K in {1..3}; do
      qiime cutadapt trim-paired \
        --i-demultiplexed-sequences trimd_Plate$K.qza \
        --p-front-f GTCGGTAAAACTCGTGCCAGC \
        --p-front-r CATAGTGGGGTATCTAATCCCAGTTTG \
        --o-trimmed-sequences trimd2_Plate$K.qza \
        --verbose > cutadapt_out2_Plate$K.txt
    done

About 87% of the data is passing the filter very consistently this time.

    grep "Total written (filtered):" cutadapt_out2_Plate$K.txt

## 3. Denoise with dada2

I am going to use the same settings that I used for the 2017 and 2018
tern fecal samples here, except I need to add the `--p-min-overlap`
parameter, otherwise I seem to be getting a load of rubbish reads which
are 250bp long and start with long strings of Cs. This is only available
in qiime2-2021.4 and later. I initially tried specifying an overlap of
30, but that didn’t seem like enough as I was still getting junk
sequences, but I think 50 is working well now.

Note, this step is pretty slow to run, a whole plate takes about 40 mins
(but that changes depending on sequencing depth).

    for K in {1..3}; do
      qiime dada2 denoise-paired \
        --i-demultiplexed-seqs trimd2_Plate$K.qza \
        --p-trunc-len-f 133 \
        --p-trunc-len-r 138 \
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
        --p-min-overlap 50 \
        --p-n-threads 12 \
        --o-representative-sequences rep-seqs_Plate$K \
        --o-table table_Plate$K \
        --o-denoising-stats denoise_Plate$K
    done

Create visualizations for the denoising stats.

    for K in {1..3}; do  
      qiime metadata tabulate\
        --m-input-file denoise_Plate$K.qza\
        --o-visualization denoise_Plate$K.qzv
    done

This looks really good. Most samples in the robot plate have about 90%
of the sequences passing through the denoising steps, which is similar
to the amount that passed for plate 12 and 14 when doing them by hand.
Most of the blanks had all data filtered out here, which is great.

## 4. Merge across plates

I need to merge both the feature tables, which contain the counts of
each feature, and the rep-seqs, which contain the actual sequence for
each feature.

    qiime feature-table merge \
      --i-tables table_Plate1.qza \
      --i-tables table_Plate2.qza \
      --i-tables table_Plate3.qza \
      --o-merged-table table_merged.qza
      
    qiime feature-table summarize \
        --i-table table_merged.qza \
        --m-sample-metadata-file metadata.txt \
        --o-visualization table_merged
      
    qiime feature-table merge-seqs \
      --i-data rep-seqs_Plate1.qza \
      --i-data rep-seqs_Plate2.qza \
      --i-data rep-seqs_Plate3.qza \
      --o-merged-data rep-seqs_merged.qza
      
    qiime feature-table tabulate-seqs \
      --i-data rep-seqs_merged.qza \
      --o-visualization rep-seqs_merged

## 5. Assign taxonomy

I’m going to use the database I created for Will’s ATPU project, which
contains all fish 12S sequences on GenBank in 2021, and Devin’s blast
method. The blast method uses an older version of qiime.

    conda activate qiime2-2019.4

    ./mktaxa.py ncbi-refseqs-withHuman.qza \
      ncbi-taxonomy-withHuman.qza \
      rep-seqs_merged.qza
      
    qiime metadata tabulate \
      --m-input-file superblast_taxonomy.qza \
      --o-visualization superblast_taxonomy

## 6. Make barplots

    qiime taxa barplot \
      --i-table table_merged.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization barplot_before_filtering.qzv

## 7. Rarefy to a sampling depth of 400

This is what I found was deep enough in the 2017 and 2018 tern samples,
and so I’ll use it here too.

    qiime feature-table rarefy \
      --i-table table_merged.qza \
      --p-sampling-depth 400 \
      --o-rarefied-table table_merged_rarefied400 
      
    qiime feature-table rarefy \
      --i-table table_merged.qza \
      --p-sampling-depth 1000 \
      --o-rarefied-table table_merged_rarefied1000 

Recreate the barplots for the rarefied table.

    qiime taxa barplot\
      --i-table table_merged_rarefied400.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization baprlot_rarefied400
      
    qiime taxa barplot\
      --i-table table_merged_rarefied1000.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization baprlot_rarefied1000

Test for composition and diversity changes uisng qiime.
