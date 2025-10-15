From raw Chip-Seq data to results:

Aligning with bowtie (example for 1 file):
```bash
/rds/bear-apps/2019b/EL8-cas/software/Bowtie2/2.3.5.1-GCC-8.3.0/bin/bowtie2-align-s --wrapper basic-0 -p 8 --very-sensitive-local -x /rds/projects/c/cazierj-ccbservice/Genomes/Mouse/mm10/bowtie2_indices/mm10 -S HB_Tal1_aligned.sam -U ./HB_Tal1.fastq.gz"
```

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
Text files created have info on duplicates removed, number of paired-end and single-end reads
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

Peak calling:
```bash
macs2 callpeak -t HB_Lmo2_aligned_sorted_noDups.bam -c HB_input_aligned_sorted_noDups.bam -f BAM -g mm -n lmo2 -q 0.05 --keep-dup auto -B --trackline<img width="477" height="80" alt="image" src="https://github.com/user-attachments/assets/f6740387-4c56-4817-a700-406d2174e4be" />
```

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
