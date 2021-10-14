# 2019 GoM Fecal Metabarcoding
Started 21st September 2021

Fecal metabarcoding analysis for samples collected in the Gulf of Maine, mostly 2019 samples.

Load qiime environment and cd to correct directory.
```
cd /Users/gemmaclucas/GitHub/Fecal_metabarcoding/2019_GoM_FecalMetabarcoding
conda activate qiime2-2021.4
```

## Import the data into Qiime2

It's saved on my back-up hardrive. This is now full, so the puffin MiFish test reads are on my solid state hardrive "Data_SS1". I've also added the data from plates 11 - 14 to Data_SS1. 
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


for K in {11..14}; do
  qiime tools import\
    --type 'SampleData[PairedEndSequencesWithQuality]'\
    --input-path /Users/gemmaclucas/Desktop/Fecal_metabarcoding/Plate$K/reads/ \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt\
    --output-path MiFish/demux_Plate$K.qza
done 

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
  
for K in {11..14}; do
  qiime demux summarize \
    --i-data demux_Plate$K.qza \
    --o-visualization demux_Plate$K.qzv
done

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

for K in {11..14}; do
  qiime cutadapt trim-paired \
    --i-demultiplexed-sequences demux_Plate$K.qza \
    --p-adapter-f CAAACTGGGATTAGATACCCCACTATG \
    --p-adapter-r GCTGGCACGAGTTTTACCGAC \
    --o-trimmed-sequences trimd_Plate$K.qza \
    --verbose > cutadapt_out_Plate$K.txt
done


for K in {11..14}; do
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

Looks like not very much passed the filters in the puffin samples. 50 - 75% seems to be about the average for the tern samples; some drop as low as 20% though. Plate 11 looks like it didn't do too well, Plate 13 also looks a bit iffy.


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
  
for K in {11..14}; do
  qiime cutadapt trim-paired \
    --i-demultiplexed-sequences trimd_Plate$K.qza \
    --p-front-f GTCGGTAAAACTCGTGCCAGC \
    --p-front-r CATAGTGGGGTATCTAATCCCAGTTTG \
    --o-trimmed-sequences trimd2_Plate$K.qza \
    --verbose > cutadapt_out2_Plate$K.txt
done
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

for K in {11..14}; do
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

for K in {11..14}; do  
  qiime metadata tabulate\
    --m-input-file denoise_Plate$K.qza\
    --o-visualization denoise_Plate$K.qzv
done
```
This looks good. It seems like 80% of sequences tend to get through all the filters/denoising from the actual samples, with <1% getting through from the blanks and negatives.

To view the rep-seqs (only running this for one plate, will view all after merging plates)
```
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_Plate7.qza \
  --o-visualization rep-seqs_Plate7
```
This looks good. Using 50 as a minimum overlap seems to get rid of junk sequences.


## 4. Merge across plates

I need to merge both the feature tables, which contain the counts of each feature, and the rep-seqs, which contain the actual sequence for each feature. The overlap method isn't important here as no samples were sequenced on multiple plates. 

I'm going to keep the puffin samples separate and I've moved puffin stuff into the puffin folder.

```
qiime feature-table merge \
  --i-tables table_Plate7.qza \
  --i-tables table_Plate8.qza \
  --i-tables table_Plate9.qza \
  --i-tables table_Plate10.qza \
  --i-tables table_Plate11.qza \
  --i-tables table_Plate12.qza \
  --i-tables table_Plate13.qza \
  --i-tables table_Plate14.qza \
  --o-merged-table table_terns7-14_merged.qza
  
qiime feature-table summarize \
    --i-table table_terns7-14_merged.qza \
    --m-sample-metadata-file metadata_terns_7-14.txt \
    --o-visualization table_terns7-14_merged
  
qiime feature-table merge-seqs \
  --i-data rep-seqs_Plate7.qza \
  --i-data rep-seqs_Plate8.qza \
  --i-data rep-seqs_Plate9.qza \
  --i-data rep-seqs_Plate10.qza \
  --i-data rep-seqs_Plate11.qza \
  --i-data rep-seqs_Plate12.qza \
  --i-data rep-seqs_Plate13.qza \
  --i-data rep-seqs_Plate14.qza \
  --o-merged-data rep-seqs_terns7-14_merged.qza
  
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_terns7-14_merged.qza \
  --o-visualization rep-seqs_terns7-14_merged
```

## 5. Assign taxonomy

I think for my WSC3 presentation, I will use my old database and Devin's blast method, for comparability between 2017/18 and 2019 samples. Before publishing, I should probably reclassify everything with an updated database.

```
conda activate qiime2-2019.4

./mktaxa.py 12SnMito.qza full_taxonomy_strings.qza rep-seqs_terns7-14_merged.qza
```

For the puffins, I'll do this in their own folder so that it doesn't get confused with the terns.

```
cd Puffin_tests
./mktaxa.py 12SnMito.qza full_taxonomy_strings.qza rep-seqs_Puffin_tests.qza
```

## 6. Make barplots
Note that my metadata table is labelled plates 7-10 but it actually contains all the plates.

For the terns:
```
cd /Users/gemmaclucas/GitHub/Fecal_metabarcoding/2019_GoM_FecalMetabarcoding/MiFish/

qiime metadata tabulate \
  --m-input-file superblast_taxonomy.qza \
  --o-visualization superblast_taxonomy
  
qiime taxa barplot \
  --i-table table_terns7-14_merged.qza \
  --i-taxonomy superblast_taxonomy.qza \
  --m-metadata-file metadata_terns_7-14.txt \
  --o-visualization superblast-barplots.qzv
```

For the puffins:
```
qiime metadata tabulate \
  --m-input-file Puffin_tests/superblast_taxonomy.qza \
  --o-visualization Puffin_tests/superblast_taxonomy
  
qiime taxa barplot \
  --i-table Puffin_tests/table_Puffin_tests.qza \
  --i-taxonomy Puffin_tests/superblast_taxonomy.qza \
  --m-metadata-file metadata_7-10-Puffin.txt \
  --o-visualization Puffin_tests/superblast-barplots.qzv
```

Looking at the Puffin samples, and comparing the ones which had weak bands where I sent both diluted and raw PCR product, the number of reads is slightly higher from the raw PCR product, but it's only a small difference. The composition of the samples seems almost identical in terms of the relative read abundance of each species, so that's very reassuring.




## 7. Remove non-food reads

I need to filter out any sequences from the bird, mammals, and unnassigned sequences before rarefying.

```
qiime taxa filter-table \
  --i-table table_terns7-14_merged.qza \
  --i-taxonomy superblast_taxonomy.qza \
  --p-exclude Unassigned,Archelosauria,Mammalia \
  --o-filtered-table table_terns7-14_merged_noBirdsMammalsUnassigned.qza
  
qiime feature-table summarize \
    --i-table table_terns7-14_merged_noBirdsMammalsUnassigned.qza \
    --m-sample-metadata-file metadata_terns_7-14.txt \
    --o-visualization table_terns7-14_merged_noBirdsMammalsUnassigned
```


## 8. Rarefy to a sampling depth of 400

This is what I found was deep enough in the 2017 and 2018 tern samples, and so to make this comparable across years, I am going to use the same rarefaction depth. I will lose 7 samples with fewer than 400 food reads.

```
qiime feature-table rarefy \
  --i-table table_terns7-14_merged_noBirdsMammalsUnassigned.qza \
  --p-sampling-depth 400 \
  --o-rarefied-table table_terns7-14_merged_noBirdsMammalsUnassigned_rarefied400 
```

Recreate the barplots for the rarefied table.
```
qiime taxa barplot\
  --i-table table_terns7-14_merged_noBirdsMammalsUnassigned_rarefied400.qza \
  --i-taxonomy superblast_taxonomy.qza \
  --m-metadata-file metadata_terns_7-14.txt \
  --o-visualization terns7-14_merged_noBirdsMammalsUnassigned_rarefied400_barplots
```

STOPPING HERE FOR NOW.

Just out of interest, I am going to make rarefaction curves to see what they look like for these samples. But first I need to collapse the taxa, so that I am rarefying based on taxonomy and not ASVs (which I don't care as much about for these purposes).

```
qiime taxa collapse \
  --i-table table_terns7-14_merged_noBirdsMammalsUnassigned.qza \
  --i-taxonomy superblast_taxonomy.qza \
  --p-level 19 \
  --o-collapsed-table table_terns7-14_merged_noBirdsMammalsUnassigned_collapsed.qza

qiime diversity alpha-rarefaction \
  --i-table table_terns7-14_merged_noBirdsMammalsUnassigned_collapsed.qza \
  --m-metadata-file metadata_terns_7-14.txt \
  --p-min-depth 100 \
  --p-max-depth 50000 \
  --o-visualization alpha-rarefaction-100-50000
  
```

## 9. What to do with my blanks?

Four of my extraction blanks have fish DNA, none of the no template controls did. So I probably just messed these up during my extractions and introduced a tiny bit of cross-contamination.

The blanks with fish DNA were:

* BLANK_10-4 - herring and 2 pollock reads
* BLANK_8-2 - only herring
* BLANK_9-4 - looks like DNA from well G6 got into this (which was in G7), herring, pollock, river herring
* BLANK_9-7 - only herring

Given that I did 32 negatives and only four had cross-contamination, I think I am comfortable with this level of cross-contamination/cross-talk. Herring is by far the dominant prey type in these plates, so it is not surprising that cross-contamination came from herring.

People on the qiime forum do not suggest removing a fixed number of reads from all samples in the feature table, since the negative control is going to be very different to a real sample that gets contaminated, so it is kind of like removing an arbitrary number of reads by doing it that way. So I think I will just be cool with it and leave it for now, but keep an eye out for any other ways to deal with it.

