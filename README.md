# 2019 GoM Fecal Metabarcoding
Started 21st September 2021

Fecal metabarcoding analysis for samples collected in the Gulf of Maine, mostly 2019 samples.

Load qiime environment and cd to correct directory.
```
cd /Users/gemmaclucas/GitHub/Fecal_metabarcoding/2019_GoM_FecalMetabarcoding
conda activate qiime2-2021.2
```

## Import the data into Qiime2

It's saved on my back-up hardrive. 
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
  --input-path /Volumes/DATA_BACKUP/2021_metabarcoding/Puffin_MiFish_tests/reads/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt\
  --output-path MiFish/demux_Puffin_tests.qza  
```