From raw Chip-Seq data to results:
<img width="583" height="484" alt="image" src="https://github.com/user-attachments/assets/053625b7-a11c-4965-8b65-481a4c7e9a59" />

### 1) QC
Aligning with bowtie (example for 1 file):
```bash
/rds/bear-apps/2019b/EL8-cas/software/Bowtie2/2.3.5.1-GCC-8.3.0/bin/bowtie2-align-s --wrapper basic-0 -p 8 --very-sensitive-local -x /rds/projects/c/cazierj-ccbservice/Genomes/Mouse/mm10/bowtie2_indices/mm10 -S HB_Tal1_aligned.sam -U ./HB_Tal1.fastq.gz"
```
-x flag specifies index file (used for aligning). You can create the index with bowtie2. You need to download the whole genome of interest, and see bowtie2 documentation.  

Sorting with samtools:
```bash
samtools sort HB_Tal1_aligned.sam
```

Looking at header of bam alignment file:
```bash
samtools head HB_Lmo2_aligned_sorted.bam
samtools head HB_Tal1_aligned_sorted.bam
samtools head HB_input_aligned_sorted.bam
```
This contains info on reference genome, what aligner was used and versions

Duplicates removing:
```bash
java -jar ../../software/picard.jar MarkDuplicates INPUT=HB_Lmo2_aligned_sorted.bam OUTPUT=HB_Lmo2_aligned_sorted_noDups.bam METRICS_FILE=HB_Lmo2_aligned_sorted_noDups.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true
java -jar ../../software/picard.jar MarkDuplicates INPUT=HB_Tal1_aligned_sorted.bam OUTPUT=HB_Tal1_aligned_sorted_noDups.bam METRICS_FILE=HB_Tal1_aligned_sorted_noDups.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true
java -jar ../../software/picard.jar MarkDuplicates INPUT=HB_input_aligned_sorted.bam OUTPUT=HB_input_aligned_sorted_noDups.bam METRICS_FILE=HB_input_aligned_sorted_noDups.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true
```

Text files created have info on duplicates removed, number of paired-end and single-end reads.

Explanation of output example:
```output
LIBRARY = Unknown Library
The read group (@RG) in your BAM has no LB: tag, so Picard can’t name the library. (You can add it with picard AddOrReplaceReadGroups.)

UNPAIRED_READS_EXAMINED = 20,725,662
Number of mapped, primary single-end reads that Picard actually examined for duplicates. (Pairs are zero, so everything is treated as unpaired.)

READ_PAIRS_EXAMINED = 0
No paired-end fragments were examined—this is a single-end BAM or lacks proper pairing flags.

SECONDARY_OR_SUPPLEMENTARY_RDS = 0
Count of secondary (0x100) or supplementary (0x800) alignments Picard saw but excluded from duplicate marking. Hero it’s zero.

UNMAPPED_READS = 25,917,921
Primary reads with the unmapped flag (0x4). Reported but not considered for duplicate detection. (Example - can be mapping rate ≈ 44%.)

UNPAIRED_READ_DUPLICATES = 7,952,974
Among the unpaired reads examined, this many were flagged as duplicates (same start position/strand after optical duplicate filtering).

READ_PAIR_DUPLICATES = 0
Duplicate pairs (for PE data). Zero because you have no pairs.

READ_PAIR_OPTICAL_DUPLICATES = 0
Optical duplicates among pairs (detected by proximity on flowcell). Zero since there are no pairs.

PERCENT_DUPLICATION = 0.383726
Fraction of examined reads marked duplicate. For single-end, roughly
UNPAIRED_READ_DUPLICATES / UNPAIRED_READS_EXAMINED ≈ 7,952,974 / 20,725,662 ≈ 0.3837.

ESTIMATED_LIBRARY_SIZE = (blank)
Picard couldn’t estimate library complexity (often happens with single-end data, high duplication, or insufficient unique reads).
```

### 2) Peak calling:
```bash
#breaking down all option, input and output
#macs2 callpeak 
#  -t HB_Tal1_aligned_sorted_noDups.bam     # treatment: TAL1 ChIP BAM
#  -c HB_input_aligned_sorted_noDups.bam    # control: input DNA BAM (background)
#  -f BAM                                   # file format is BAM (single-end). For PE use -f BAMPE
#  -g mm                                    # effective genome size: mouse (mm10 ~1.87e9)
#  -n tal1                                  # output name prefix “tal1”
#  -q 0.05                                  # FDR threshold 5% for calling peaks
#  --keep-dup auto                          # allow identical start sites up to a data-driven cap
#  -B                                       # write bedGraph signal tracks (pileup & lambda)
#  --trackline                              # add UCSC track lines to bedGraphs

#actual commands for 2 files:
macs2 callpeak -t HB_Tal1_aligned_sorted_noDups.bam -c HB_input_aligned_sorted_noDups.bam -f BAM -g mm -n tal1 -q 0.05 --keep-dup auto -B --trackline
macs2 callpeak -t HB_Lmo2_aligned_sorted_noDups.bam -c HB_input_aligned_sorted_noDups.bam -f BAM -g mm -n lmo2 -q 0.05 --keep-dup auto -B --trackline
```
  
**How MACs work**:  
- **A per-base coverage track** is a genome track where every single nucleotide position has a number: how many sequencing fragments overlap that base.  
    - What it represents: at position i on chr1, the value is the count of reads/fragments covering that base. Do this for base 1, base 2, base 3… across the genome → you get a “track” of coverage.  
    - Why it’s useful: shows exactly where signal is high/low, allows you to see sharp peaks and summits (e.g., TF binding), and is the rawest view of enrichment before binning/smoothing.  
- **Input DNA** (a.k.a. background control): genomic DNA from the same cells, processed the same way but without the IP (no antibody pulldown).  
- **treatment pileup**: It’s a per-base coverage track of your ChIP (treatment) library after MACS2 has reconstructed fragments from reads.   
- At each genomic base, the value is the number of inferred DNA fragments overlapping that base (sometimes a float due to internal scaling).
MACS2 writes this to *_treat_pileup.bdg when you use -B.  
- The ENCODE **blacklist** contains genomic regions which have anomalous, unstructured or high signals in sequencing data in all sequencing experiments, irrespective of cell type. These regions include, for example, repetitive regions, and have to be removed from the data.  
- The MACS peak calling algorithm can call **narrow or broad peaks**. Narrow peaks are suitable for finding short genomic regions such as transcription factor binding sites. Broad peaks are often comparable to a gene’s length and can be used to identify where RNA Pol II of the transcription machinery is on the DNA (if you use antibodies for RNA pol ii)  
    - Narrow mode (default) finds sharp local maxima and keeps peaks small and distinct. Good for motif finding and precise summits.  
    - Broad mode (--broad, --broad-cutoff) merges adjacent enriched bins into longer regions, so you don’t split one gene-body signal into dozens of tiny peaks.  

**Key files**:
- *_peaks.narrowPeak = the regions (only peaks)
    - A list of start–end intervals where MACS2 says “there’s a real peak here.”
    - Each row = one peak with stats (strength, p/q-values).
    - Use when you need the peak regions themselves.
- *_summits.bed = the exact tip (same number of rows as narrowPeak)
    - One point inside each peak: the highest point (the “summit”).
    - Great for motif finding or pinpointing where binding is strongest.
- *_treat_pileup.bdg = the signal track (all nucleotides)
- A continuous curve: how much ChIP signal at each base across the genome.
    - Looks like a waveform in IGV/UCSC; peaks rise above background.

**Steps**:  
1) Filtering out duplicate reads.  
2) Building a background model for the peaks, against which the actual coverages can be compared.  
3) Peak detection against the computed background.  
4) Multiple testing correction using the Benjamini-Hochberg (aka FDR) correction.  

**Options**:
-t: the treatment file name.  
-c: the control (ChIP-seq input) file name.  
-g: the genome used, to set the effective genome size  
--keep-dup: it sets if we should keep all duplicates (all) or remove all (1, default), or something in between.  
-q: Q-value cutoff for calling statistically significant peaks during the FDR correction.  
-B: If given, store the pile-up (i.e. the coverage track) in a bedgraph format.  


A number of new files have been created for the treat pile-up or read coverage (ending in _treat_pileup.bdg), the background model (ending with _control_lambda.bdg and _model.r) and finally the called peaks (ending with  _peaks.xls, _summits.bed and _peaks.narrowPeak). Inspect the peak files to identify peaks.

Removing blacklisted regions:
```bash
bedtools intersect -a lmo2_summits.bed -b annotation_files/mm10_simpleRepeat.bed -v | bedtools intersect -a - -b annotation_files/mm10-blacklist.bed -v > lmo2_noBL_summits.bed
bedtools intersect -a tal1_summits.bed -b annotation_files/mm10_simpleRepeat.bed -v | bedtools intersect -a - -b annotation_files/mm10-blacklist.bed -v > tal1_noBL_summits.bed
```

Get peaks from summit files:
```bash
cut -f4 lmo2_noBL_summits.bed | grep -Fwf - lmo2_peaks.narrowPeak | grep -v chrUn | grep -v chrM | grep -v random > lmo2_filtered_peaks.bed
cut -f4 tal1_noBL_summits.bed | grep -Fwf - tal1_peaks.narrowPeak | grep -v chrUn | grep -v chrM | grep -v random > tal1_filtered_peaks.bed
```

Visualizing in IGV (doesn't work, says image out of bounds for the first 2 tracks):
```python
# Load the igv_notebook package
import igv_notebook

# Initialise a session and load the data
igv_notebook.init()

igv_browser= igv_notebook.Browser(
    {
        "genome": "mm10",
        "tracks": [{
            "name": "Tal1",
            "path": "./results/tal1_treat_pileup.bdg.tdf",
            "format": "tdf"
        },
        {
            "name": "Lmo2",
            "path": "./results/lmo2_treat_pileup.bdg.tdf",
            "format": "tdf"
        }]
    }
)
```

Visualizing filtered peaks (doesn't work, says image out of bounds for the first 2 tracks):
```python
# Load the igv_notebook package
import igv_notebook

# Initialise a session and load the data
igv_notebook.init()

igv_browser= igv_notebook.Browser(
    {
        "genome": "mm10",
        "tracks": [{
            "name": "Tal1",
            "path": "./results/tal1_treat_pileup.bdg.tdf",
            "format": "tdf"
        },
        {
            "name": "Lmo2",
            "path": "./results/lmo2_treat_pileup.bdg.tdf",
            "format": "tdf"
        },
        {
            "name": "Tal1 peaks",
            "path": "tal1_filtered_peaks.bed",
            "format": "bed"
        },
        {
            "name": "Lmo2 peaks",
            "path": "lmo2_filtered_peaks.bed",
            "format": "bed"
        }]
    }
)
```


### 3) Comparative analysis:
Finding which peaks are specific for each cell type (-v option excludes, and -u option combines, -a and -b specify the input files):
```bash
bedtools intersect -a lmo2_filtered_peaks.bed -b tal1_filtered_peaks.bed -v > lmo2_specific.bed
bedtools intersect -a lmo2_filtered_peaks.bed -b tal1_filtered_peaks.bed -u > shared_peaks.bed
bedtools intersect -a tal1_filtered_peaks.bed -b lmo2_filtered_peaks.bed -v > tal1_specific.bed
```
Shared ChIP-seq peaks of files indicate that they bind at the same site. However, without further wet lab experiments, we cannot tell if the TFs compete with each other at these shared binding sites, or if they act collaboratively and one is required for the recruitment of the other TF.  


Annotating separate and common peaks with HOMER:
```bash
annotatePeaks.pl lmo2_specific.bed mm10 -annStats lmo2_peaks_annstats.tsv > lmo2_peaks_annotation.tsv 2> lmo2_peaks_homer.out
annotatePeaks.pl tal1_specific.bed mm10 -annStats tal1_peaks_annstats.tsv > tal1_peaks_annotation.tsv 2> tal1_peaks_homer.out
annotatePeaks.pl shared_peaks.bed mm10 -annStats shared_peaks_annstats.tsv > shared_peaks_annotation.tsv 2> shared_peaks_homer.out
```

Input:   
- bed file with peaks  
- genome version (example - mm10; hg18) - Must be installed in HOMER  
  
Output:  
- lmo2_peaks_annotation.tsv – with annotation (main!)  
- lmo2_peaks_annstats.tsv (from -annStats) – with annotation stats (a summary table of counts/percentages by annotation category and distance bins (e.g., how many peaks in promoters, introns, intergenic, within ±1 kb of TSS, etc.).) In the annStat files, the enrichment of features in the peaks, log2(obs/exp), is reported. If this number is positive, the observed number is higher than expected. The significance of the enrichment is reported as raw log p-values. The sign of the raw log p-value should be negative. As there are 12 annotation categories tested, we should apply a multiple testing correction.

- lmo2_peaks_homer.out - log: progress, genome loading, numbers of peaks processed, warnings (e.g., unmapped chromosomes), parameter echo, any errors.  


GO and KEGG enrichment saved in R file in this directory
