# 2019 GoM Fecal Metabarcoding
Started 21st September 2021

Fecal metabarcoding analysis for samples collected in the Gulf of Maine, mostly 2019 samples.

Load qiime environment and cd to correct directory.
```
cd /Users/gemmaclucas/GitHub/Fecal_metabarcoding/2019_GoM_FecalMetabarcoding
conda activate qiime2-2021.2
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

# NOT RUN YET:
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

# NOT RUN YET:
qiime demux summarize \
  --i-data trimd_Puffin_tests.qza \
  --o-visualization trimd_Puffin_tests.qzv
```

To see how much data passed the filter for each sample

```
grep "Total written (filtered):" cutadapt_out_Plate$K.txt 
```

Looks like not very much passed the filters in the puffin samples - compare to terns.


