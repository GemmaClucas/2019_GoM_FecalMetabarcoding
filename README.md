# 2019 GoM Fecal Metabarcoding
Started 21st September 2021

Fecal metabarcoding analysis for samples collected in the Gulf of Maine, mostly 2019 samples.

Load qiime environment and cd to correct directory.
```
cd /Users/gemmaclucas/GitHub/Fecal_metabarcoding/2019_GoM_FecalMetabarcoding
conda activate qiime2-2021.4
```

## Import the data into Qiime2

It's saved on my back-up hardrive. This is now full, so the puffin MiFish test reads are on my solid state hardrive "Data_SS1".
```
qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path /Volumes/DATA_BACKUP/2021_metabarcoding/Plate7/reads/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt\
  --output-path MiFish/demux_Plate7.qza
  
qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path /Volumes/DATA_BACKUP/2021_metabarcoding/Plate8/reads/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt\
  --output-path MiFish/demux_Plate8.qza
  
qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path /Volumes/DATA_BACKUP/2021_metabarcoding/Plate9/reads/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt\
  --output-path MiFish/demux_Plate9.qza
  
qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path /Volumes/DATA_BACKUP/2021_metabarcoding/Plate10/reads/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt\
  --output-path MiFish/demux_Plate10.qza
  
qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path /Volumes/Data_SS1/Puffin_MiFish_tests/reads/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt\
  --output-path MiFish/demux_Puffin_tests.qza  
```

Move to the diectory where all the next steps will be. Summarise read quality and number of reads for each plate.
```
cd MiFish/

for K in {7..10}; do
  qiime demux summarize \
    --i-data demux_Plate$K.qza \
    --o-visualization demux_Plate$K.qzv
done

qiime demux summarize \
  --i-data demux_Puffin_tests.qza \
  --o-visualization demux_Puffin_tests.qzv
```

## 2. Trim primers with cutadapt plugin

F primer: GTCGGTAAAACTCGTGCCAGC (21 bp)  
R primer: CATAGTGGGGTATCTAATCCCAGTTTG (27 bp)

### Trim 3' ends first
At the 3' end of the read, the primer will have been read through after reading the MiFish sequence. I need to be looking for the reverse complement of the reverse primer in R1 (```—p-adapter-f```) and the reverse complement of the forward primer in R2 (```—p-adapter-r```)

F primer reverse complement: GCTGGCACGAGTTTTACCGAC    
R primer reverse complement: CAAACTGGGATTAGATACCCCACTATG

```
for K in {7..10}; do
  qiime cutadapt trim-paired \
    --i-demultiplexed-sequences demux_Plate$K.qza \
    --p-adapter-f CAAACTGGGATTAGATACCCCACTATG \
    --p-adapter-r GCTGGCACGAGTTTTACCGAC \
    --o-trimmed-sequences trimd_Plate$K.qza \
    --verbose > cutadapt_out_Plate$K.txt
done

for K in {7..10}; do
  qiime demux summarize \
    --i-data trimd_Plate$K.qza \
    --o-visualization trimd_Plate$K.qzv
done    
```

And the same for the puffins
```
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux_Puffin_tests.qza \
  --p-adapter-f CAAACTGGGATTAGATACCCCACTATG \
  --p-adapter-r GCTGGCACGAGTTTTACCGAC \
  --o-trimmed-sequences trimd_Puffin_tests.qza \
  --verbose > cutadapt_out_Puffin_tests.txt

qiime demux summarize \
  --i-data trimd_Puffin_tests.qza \
  --o-visualization trimd_Puffin_tests.qzv
```

To see how much data passed the filter for each sample:

```
grep "Total written (filtered):" cutadapt_out_Plate$K.txt 
```

Looks like not very much passed the filters in the puffin samples. 50 - 75% seems to be about the average for the tern samples; some drop as low as 20% though.


### Trim 5' ends of reads

All R1 should begin with the forward primer: GTCGGTAAAACTCGTGCCAGC (21 bases).  
All R2 should begin with the reverse primer: CATAGTGGGGTATCTAATCCCAGTTTG (27 bases).

Trim these with the following commands:

```
for K in {7..10}; do
  qiime cutadapt trim-paired \
    --i-demultiplexed-sequences trimd_Plate$K.qza \
    --p-front-f GTCGGTAAAACTCGTGCCAGC \
    --p-front-r CATAGTGGGGTATCTAATCCCAGTTTG \
    --o-trimmed-sequences trimd2_Plate$K.qza \
    --verbose > cutadapt_out2_Plate$K.txt
done

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences trimd_Puffin_tests.qza \
  --p-front-f GTCGGTAAAACTCGTGCCAGC \
  --p-front-r CATAGTGGGGTATCTAATCCCAGTTTG \
  --o-trimmed-sequences trimd2_Puffin_tests.qza \
  --verbose > cutadapt_out2_Puffin_tests.txt
```

About 90% of the data seems to pass this filter for both puffins and terns. 

Not going to make the qzv files for now.

## 3. Denoise with dada2

I am going to use the same settings that I used for the 2017 and 2018 tern fecal samples here, except I need to add the ```--p-min-overlap``` parameter, otherwise I seem to be getting a load of rubbish reads which are 250bp long and start with long strings of Cs. This is only available in qiime2-2021.4 and later. I initally tried specifying an overlap of 30, but that didn't seem like enough as I was still getting junk sequences, but I think 50 is working well now.

Note, this step is pretty slow to run, a whole plate takes about 40 mins (but that changes depending on sequencing depth).

```
for K in {7..10}; do
  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs trimd2_Plate$K.qza \
    --p-trunc-len-f 133 \
    --p-trunc-len-r 138 \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-min-overlap 50 \
    --p-n-threads 16 \
    --o-representative-sequences rep-seqs_Plate$K \
    --o-table table_Plate$K \
    --o-denoising-stats denoise_Plate$K
done

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimd2_Puffin_tests.qza \
  --p-trunc-len-f 133 \
  --p-trunc-len-r 138 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-min-overlap 50 \
  --p-n-threads 16 \
  --o-representative-sequences rep-seqs_Puffin_tests \
  --o-table table_Puffin_tests \
  --o-denoising-stats denoise_Puffin_tests 

```

Create visualizations for the denoising stats.
```
qiime metadata tabulate\
  --m-input-file denoise_Puffin_tests.qza\
  --o-visualization denoise_Puffin_tests.qzv

for K in {7..10}; do  
  qiime metadata tabulate\
    --m-input-file denoise_Plate$K.qza\
    --o-visualization denoise_Plate$K.qzv
done
```
This looks good. It seems like 80% of sequences tend to get through all the filters/denoising from the actual samples, with <1% getting through in the blanks and negatives.

To view the rep-seqs (only running this for one plate, will view all after merging plates)
```
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_Plate7.qza \
  --o-visualization rep-seqs_Plate7
```
This looks good. Using 50 as a minimum overlap seems to get rid of junk sequences.


## 4. Merge across plates

I need to merge both the feature tables, which contain the counts of each feature, and the rep-seqs, which contain the actual sequence for each feature. The overlap method isn't important here as no samples were sequenced on multiple plates. I'm going to include the puffin samples, so that I can do the taxonomy assigment all together later.

```
qiime feature-table merge \
  --i-tables table_Plate7.qza \
  --i-tables table_Plate8.qza \
  --i-tables table_Plate9.qza \
  --i-tables table_Plate10.qza \
  --i-tables table_Puffin_tests.qza \
  --o-merged-table table_merged.qza
  
qiime feature-table merge-seqs \
  --i-data rep-seqs_Plate7.qza \
  --i-data rep-seqs_Plate8.qza \
  --i-data rep-seqs_Plate9.qza \
  --i-data rep-seqs_Plate10.qza \
  --i-data rep-seqs_Puffin_tests.qza \
  --o-merged-data rep-seqs_merged.qza
  
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_merged.qza \
  --o-visualization rep-seqs_merged
```

## 5. Assign taxonomy

I think for my WSC3 presentation, I will use my old database and Devin's blast method, for comparability between 2017/18 and 2019 samples. Before publishing, I should probably reclassify everything with an updated database.

```
conda activate qiime2-2019.4

./mktaxa.py 12SnMito.qza full_taxonomy_strings.qza rep-seqs_merged.qza
```

## 6. Make barplots

```
qiime metadata tabulate \
  --m-input-file superblast_taxonomy.qza \
  --o-visualization superblast_taxonomy
  
qiime taxa barplot \
  --i-table table_merged.qza \
  --i-taxonomy superblast_taxonomy.qza \
  --m-metadata-file metadata_7-10-Puffin.txt \
  --o-visualization superblast-barplots.qzv
```

Looking at the Puffin samples, and comparing the ones which had weak bands where I sent both diluted and raw PCR product, the number of reads is slightly higher from the raw PCR product, but it's only a small difference. The composition of the samples seems almost identical in terms of the relative read abundance of each species, so that's very reassuring.

