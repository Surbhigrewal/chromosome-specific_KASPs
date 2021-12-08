# chromosome-specific_KASPs

## subline


###index genome for mapping

```
bwa index iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta
```

###map paired-end reads to genome and filter

```
bwa mem -t 8 -M iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta muticum_reads_1.fq muticum_reads_2.fq | awk '$1 ~ /^@/ || $2 == 65 || $2 == 129 || $2 == 67 || $2 == 131 || $2 == 113 || $2 == 177 || $2 == 81 || $2 == 161 || $2 == 163 || $2 == 83 || $2 == 97 || $2 == 145 || $2 == 99 || $2 == 147 || $2 == 137 || $2 == 73 {print $0}' | > muticum.sam
```

###keep uniquely mapping reads and sort bam file

```
samtools view -u -q 10 muticum.sam | samtools sort -o muticum_filt_srt.bam -
```

###index bam file

```
samtools index -c muticum_filt_srt.bam
```

###make dictionary of genome

```
java -jar -Xmx10g picard.jar CreateSequenceDictionary R= iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta O= iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.dict
```

###remove PCR duplicates from bam

```
java -jar -Xmx10g picard.jar MarkDuplicates I= muticum_filt_srt.bam O= muticum_remove_dups.bam M=duplication.txt REMOVE_DUPLICATES=true AS=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 VALIDATION_STRINGENCY=LENIENT
```

###index bam file

```
samtools index -c muticum_remove_dups.bam
```

###call variants using bcftools - each chromosome runs on different thread

```
python3 bcftools_call_parallel.py muticum_remove_dups.bam muticum
```

###concatanate variant calling from each chromosome

```
bcftools concat -o muticum.raw.vcf muticum_variants*
```


###filter variants

```
java -jar -Xmx10g GenomeAnalysisTK.jar -R iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta -T VariantFiltration --variant muticum.raw.vcf -o muticum_vfalleleINDELcalls.vcf --clusterSize 3 --clusterWindowSize 10 --filterExpression "DP < 5 || QUAL < 30.0" --filterName "DodgySNPs" --genotypeFilterExpression "isHet == 1" --genotypeFilterName Heterozygote --genotypeFilterExpression "isHomVar == 1" --genotypeFilterName Homozygote --genotypeFilterExpression "isHomRef == 1" --genotypeFilterName HomozygoteRef
```

###remove INDELs

```
grep -v "INDEL" muticum_vfalleleINDELcalls.vcf | grep -v "^#" > muticum_vfallelecalls_filt_flanking_snps.vcf
```

###remove SNPs with SNP within 50bp

```
python3 filter_SNPs_with_nearby_SNPs.py muticum_vfallelecalls.vcf 50 muticum_vfallelecalls_filt_flanking_snps.vcf
```

###keep homozygous SNPs that passed the GATK filters

```
grep "PASS" muticum_vfallelecalls_filt_flanking_snps.vcf | grep "Homo" > muticum_vfallelecalls_filt_flanking_snps_PASS_homo.vcf
```

###extract the sequence 50bp either side of SNP

```
python3 extract_100bp_KASP_sequence_for_blast.py iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta muticum_vfallelecalls_filt_flanking_snps_PASS_homo.vcf > muticum_vfallelecalls_filt_flanking_snps_PASS_homo.fasta
```

###blast sequences to genome

```
blastn -db iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta -query muticum_vfallelecalls_filt_flanking_snps_PASS_homo.fasta -outfmt 6 -num_threads 16 > muticum_vfallelecalls_filt_flanking_snps_PASS_homo_blast.tsv
```

###discard sequences that hit to more than one place in the genome (any non self hits)

```
cat muticum_vfallelecalls_filt_flanking_snps_PASS_homo_blast.tsv | python3 filter_any_blast_hit.py muticum_vfallelecalls_filt_flanking_snps_PASS_homo_single_copy_id.txt
```

###produce final fasta of sequences

```
grep -A1 -f muticum_vfallelecalls_filt_flanking_snps_PASS_homo_single_copy_id.txt muticum_vfallelecalls_filt_flanking_snps_PASS_homo.fasta | grep -v "^-" > Am_muticum_snp_regions_single_copy.fasta
```

###produce final vcf

```
sed 's/-/\.\*/' muticum_vfallelecalls_filt_flanking_snps_PASS_homo_single_copy_id.txt | grep f - muticum_vfallelecalls_filt_flanking_snps_PASS_homo.fasta > Am_muticum_homo_snps_single_copy.vcf
```

###produce SNP file with ref and alt alleles in square brackets

```
python3 produce_final_KASP_file.py muticum_vfallelecalls_filt_flanking_snps_PASS_homo.fasta muticum_vfallelecalls_filt_flanking_snps_PASS_homo.fasta muticum_vfallelecalls_filt_flanking_snps_PASS_homo_single_copy_id.txt
```
