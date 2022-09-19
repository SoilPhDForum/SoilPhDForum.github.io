# Metagenomics pipeline

**_Kankan Zhao (zhaokk@zju.edu.cn)_**

Up to now, there is no state-of-the-art pipeline of metagenomics analysis (especially for environmental samples) because of the huge number of _**microbial dark matters**_ and a wide variety of emerging bioinformatic tools. Hence, this pipeline is a assemble-based method of paired-end metagenomic data. It assumes you have a basic knowledge of metagenomics and Linux commands. Also, it is important to read the references listed to have a better understanding of each parameters (and it is better to read more articles cited the tools you used), because sometimes there is a lack of consensus among different articles. I hope it works for you and please feel free to contact me with any questions or tell me the mistakes of this pipeline.

## 1. Quality control

### _1.1 Evaluate data quality_

```
export PATH=$PATH:/public/home/bma/application/FastQC
fastqc -o 01.rawdata/ -t 28 01.rawdata/S1_raw_R1.fq.gz
```

> [FastQC Manual](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) and [Tutorial](https://rtsf.natsci.msu.edu/sites/_rtsf/assets/File/FastQC_TutorialAndFAQ_080717.pdf).

### _1.2 Quality control_

If the results do not reach the standard, you should do it repeatedly.

```
java -jar <path to trimmomatic.jar> PE [-threads <threads>] | <input 1> <input 2>] |<paired output 1> <unpaired output 1> <paired output 2> <unpaired output 2> | ILLUMINACLIP:<path to adapter-PE.fa>:2:30:10 | LEADING:3 TRAILING:3 | SLIDINGWINDOW:<windowSize>:<requiredQuality> | MINLEN:50

# ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
# SLIDINGWINDOW: Perform a sliding window trimming approach. It starts scanning at the 5‟ end and clips the read once the average quality within the window falls below a threshold.
# LEADING: Cut bases off the start of a read, if below a threshold quality.
# TRAILING: Cut bases off the end of a read, if below a threshold quality.
# CROP: Cut the read to a specified length by removing bases from the end.
# HEADCROP: Cut the specified number of bases from the start of the read.
# MINLEN: Drop the read if it is below a specified length.
```

Example:

```
java -jar /public/home/bma/application/Trimmomatic-0-2.39/trimmomatic-0.39.jar PE -threads 28 01.rawdata/S1_raw_R1.fq.gz 01.rawdata/S1_raw_R2.fq.gz 02.paired.fq.gz/S1_paired_R1.clean.fastq.gz  02.unpaired.fq.gz/S1_unpaired_R1.clean.fastq.gz  02.paired.fq.gz/S1_paired_R2.clean.fastq.gz  02.unpaired.fq.gz/S1_unpaired_R2.clean.fastq.gz ILLUMINACLIP:/public/home/bma/application/Trimmomatic-0-2.39/adapter-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
```

> Bolger, A. M., Lohse, M., & Usadel, B. [Trimmomatic: A flexible trimmer for Illumina Sequence Data](https://academic.oup.com/bioinformatics/article/30/15/2114/2390096). *Bioinformatics*, 2014, 30(15): 2114-2120.

### _1.3 Evaluate data quality again_

```
export PATH=$PATH:/public/home/bma/application/FastQC
fastqc -o 02.paired.fq.gz/ -t 28 02.paired.fq.gz/S1_paired_R1.clean.fastq.gz
```

### _1.4 Remove host reads (Optional)_

```
export PATH=$PATH:/public/home/bma/application/bowtie2-2.3.2
bowtie2-build 03.dehost/01.ref/host_ref_genome.fasta 03.dehost/01.ref/index
bowtie2 -x 03.dehost/01.ref/index -1 02.paired.fq.gz/S1_paired_R1.clean.fastq.gz -2 02.paired.fq.gz/S1_paired_R2.clean.fastq.gz -S 03.dehost/02.sam/S1.sam --threads 28 --un-conc-gz 03.dehost/03.cleandata/S1.dehost.fq.gz

# BWA is also a great tool.
# S1.dehost.fq.gz would be devided into S1.dehost.fq.1.gz and S1.dehost.fq.2.gz automatically. So you could rename them as S1.dehost_R1.fq.gz and S1.dehost_R2.fq.gz for further analysis.

#export PATH=$PATH:/share/home/bmalab/samtools-1.11/
#export PATH=$PATH:/share/home/bmalab/bwa-master/
#bwa index haishen.fasta
#bwa mem -t 36 00.host/haishen.fasta 01.rawdata/SRR8524814_paired_R1.clean.fastq.gz 01.rawdata/SRR8524814_paired_R2.clean.fastq.gz|samtools sort -O bam -@ 36 -o - > 02.cleandata/SRR8524814.bam
#samtools view -@ 36 -b -f 12 -F 256 02.cleandata/SRR8524814.bam > 02.cleandata/SRR8524814_unmap.bam
#samtools sort -n 02.cleandata/SRR8524814_unmap.bam -o 02.cleandata/SRR8524814_unmap.sorted.bam
#bedtools bamtofastq -i SRR8524814_unmap.sorted.bam -fq SRR8524814_R1.fq -fq2 SRR8524814_R2.fq
```

> Langmead B, Salzberg S. [Fast gapped-read alignment with Bowtie 2](http://www.nature.com/nmeth/journal/v9/n4/full/nmeth.1923.html). *Nature Methods*. 2012, 9:357-359.

> Li H. and Durbin R. [Fast and accurate short read alignment with Burrows-Wheeler Transform](https://academic.oup.com/bioinformatics/article/25/14/1754/225615). *Bioinformatics*, (2009), 25:1754-60.

### _1.5 Repair paired reads (Optional)_

After host reads removal, the reads would be unpaired again. You could remove host reads **before** quality control **or** repair them with ``repair.sh`` in ``BBMap`` or other tools like ``Samtools``.

Before repairation, decompression is needed.

```
export PATH=$PATH:/public/home/bma/application/pigz-2.4/
unpigz -k -p 28 *.gz

# -k: keep the package files.
# -p: thread number.
```

```
export PATH=$PATH:/public/home/bma/application/bbmap/
repair.sh in=S1_paired_R1.clean.fq in2=S1_paired_R1.clean.fq out1=S1_R1.fq out2=S1_R2.fq
```

Then, it is better to compress the fastq files again.

```
export PATH=$PATH:/public/home/bma/application/pigz-2.4/
pigz -p 28 *.fq
```

> Bushnell B, Rood J, Singer E. [BBMerge–Accurate paired shotgun read merging via overlap](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0185056). *PloS one*, 2017, 12(10): e0185056.

> Li H, Handsaker B, Wysoker A, et al. [The sequence alignment/map format and SAMtools](https://academic.oup.com/bioinformatics/article/25/16/2078/204688?login=true). *Bioinformatics*, 2009, 25(16): 2078-2079.

---

## 2. Assemble and dereplication

### _2.1 Assemble_

If you have enough computational resources and time, you could choose [metaSPAdes](https://genome.cshlp.org/content/27/5/824), believed to be the best metagenome assembler. Besides, to improve reads retrieval rate, it is better to assemble contigs by single sample, divided groups and all samples, respectively. [Mash](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x) could be used to separated samples into groups based on microbial community similarities like _[Jin et al. (2018)](https://www.nature.com/articles/s41467-018-07343-2#Sec7)_.

```
export PATH=$PATH:/public/home/bma/application/MEGAHIT-1.2.9-Linux-x86_64-static/bin/
megahit -t 28 -m 0.95 --min-contig-len 1500 --k-min 21 --k-step 10 -1 S1_R1_fq.gz -2 S1_R2_fq.gz -o 04.assemble/S1/
```

To avoid unnecessary problems, you could remove the characters after space of contain name.

```
./seqkit replace -p " .+" -i assemble_contigs.fa > contigs.fa
```

> Li D, Liu C M, Luo R, et al. [MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph](https://academic.oup.com/bioinformatics/article/31/10/1674/177884). _Bioinformatics_, 2015, 31(10): 1674-1676.

> Xu J, Zhang Y, Zhang P, et al. [The structure and function of the global citrus rhizosphere microbiome](https://www.nature.com/articles/s41467-018-07343-2#Sec7). _Nature communications_, 2018, 9(1): 1-10.

### _2.2 Dereplication (Optional)_

If you need the information of the relationship between contigs and genes, you should skip this step.

---

## 3. Abundance calculation of contigs and taxanomy annotation

### _3.1 Abundance calculation of contigs_

I apply a stratagy based on depth-of-coverage by dividing the summed depth per base by the length of the respective sequence, which is same as *[Woodcroft et al.(2018)](https://www.nature.com/articles/s41586-018-0338-1#Sec8)*, _[Salazar et al. (2019)](https://www.sciencedirect.com/science/article/pii/S009286741931164X#sec3)_ and _[Herold et al. (2020)](https://www.nature.com/articles/s41467-020-19006-2#Sec8)_. Other methods like [Kraken2](https://link.springer.com/article/10.1186/s13059-019-1891-0), [Kaiju](https://www.nature.com/articles/ncomms11257) and [MetaPhlAn2](https://www.nature.com/articles/nmeth.3589?report=reader) are also feasible. the comparision of their performances could be seen in _[Sun et al. (2021)](https://www.nature.com/articles/s41592-021-01141-3)_. However, they usually perform not well in complex environment samples with a large amount of microbial dark matters.

```
export PATH=$PATH:/share/home/bmalab/bwa-master/
bwa index contigs.fa
```

```
export PATH=$PATH:/share/home/bmalab/samtools-1.11/
bwa mem -t 28 contigs.fa S1_R1.fq S1_R2.fq|samtools sort -O bam -@ 28 -o - > S1_contig.bam
```

BamM is no longer being maintained. Instead try CoverM which is easier to both install and use, and is faster. Parameters is same as _[Woodcroft et al. (2018)](https://www.nature.com/articles/s41586-018-0338-1#Sec8)_, which is much stricter than _[Salazar et al. (2019)](https://www.sciencedirect.com/science/article/pii/S009286741931164X#sec3)_ and _[Herold et al. (2020)](https://www.nature.com/articles/s41467-020-19006-2#Sec8)_. The relative abundance of each contig in each sample could be calculated as its coverage divided by the total coverage of all contigs. Other [coverage measure methods](https://github.com/wwood/CoverM/#calculation-methods) could also be applied.

```
export PATH=$PATH:/public/home/bma/application/coverm-x86_64-unknown-linux-musl-0.2.0-alpha7/
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b S1_contig.bam -o S1_contig_filter.bam -t 28
coverm contig --trim-max 90 --trim-min 10 --bam-files S1_contig_filter.bam  > S1_contig_coverage.csv -t 28
```

> Sun Z, Huang S, Zhang M, et al. [Challenges in benchmarking metagenomic profilers](https://www.nature.com/articles/s41592-021-01141-3). _Nature methods_, 2021, 18(6): 618-626.

> Woodcroft B J, Singleton C M, Boyd J A, et al. [Genome-centric view of carbon processing in thawing permafrost](https://www.nature.com/articles/s41586-018-0338-1#Sec8). _Nature_, 2018, 560(7716): 49-54.

> Salazar G, Paoli L, Alberti A, et al. [Gene expression changes and community turnover differentially shape the global ocean metatranscriptome](https://www.sciencedirect.com/science/article/pii/S009286741931164X#sec3). *Cell*, 2019, 179(5): 1068-1083. e21.

> Herold M, Arbas S M, Narayanasamy S, et al. [Integration of time-series meta-omics data reveals how microbial ecosystems respond to disturbance](https://www.nature.com/articles/s41467-020-19006-2#Sec8). *Nature communications*, 2020, 11(1): 1-14.

### _3.2 Taxanomy annotation_

---

## 4. ORFs prediction and dereplication

### _4.1 ORFs (open reading frames) prediction_

Here, you could choose the raw contigs or non-redundant contigs to predict ORFs. If you need the information of the source contigs of ORFs, you should choose the raw contigs. To improve the efficiency, you could split the contig file into several parts and combine the results after prediction.

```
export PATH=$PATH:/public/home/bma/application/Prodigal/
prodigal -i contigs.fa -a amino.faa -d nucl.fnn -o genes.txt -p meta

# faa file contains amino acids
# fnn file contains nucleotides
```

To avoid unnecessary problems, you could remove the characters after space of gene name.

```
./seqkit replace -p " .+" -i nucl.fnn > nucl_rename.fnn
./seqkit replace -p " .+" -i amino.faa > amino_rename.faa
```

> Hyatt D, Chen G L, LoCascio P F, et al. [Prodigal: prokaryotic gene recognition and translation initiation site identification](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119). *BMC bioinformatics*, 2010, 11(1): 1-11.

### _4.2 Dereplication (Optional)_

If you use the non-redundant contigs before and need the information of the source contigs of ORFs, you could skip this step.

CD-HIT is also a great tool, but much slower than MMseqs2. In MMseqs2, ``--cluster-mode 1 or 2`` is similar to greedy incremental clustering strategy to cluster sequences in CD-HIT. Here, I choose parameters ``--min-seq-id 0.95 -c 0.9 --cluster-mode 2 --cov-mode 1``, which is same as *[Salazar et al. (2019)](https://www.sciencedirect.com/science/article/pii/S009286741931164X#sec3)*.

```
export PATH=$PATH:/share/home/bmalab/mmseqs/bin/
mmseqs easy-linclust amino_rename.faa clusterRes tmp --min-seq-id 0.95 -c 0.9 --threads 36 --cluster-mode 2 --cov-mode 1

# clusterRes_rep_seq.fasta is representative amino acids
```

> [MMseqs2 User Guide](https://mmseqs.com/latest/userguide.pdf)

> Steinegger M, Söding J. [MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets](https://www.nature.com/articles/nbt.3988). _Nature biotechnology_, 2017, 35(11): 1026-1028.

> Salazar G, Paoli L, Alberti A, et al. [Gene expression changes and community turnover differentially shape the global ocean metatranscriptome](https://www.sciencedirect.com/science/article/pii/S009286741931164X#sec3). *Cell*, 2019, 179(5): 1068-1083. e21.

---

## 5. Abundance calculation of ORFs and gene annotation

If you already have target gene databases, you could mapping the reads to databases directly by [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)-like tools (such as [DIAMOND](https://www.nature.com/articles/s41592-021-01101-x), [MMseqs2](https://www.nature.com/articles/nbt.3988) and [USEARCH11](https://www.drive5.com/usearch/)) or [HMMER](http://hmmer.org/download.html). Here is a [comparasion](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07132-6) of BLAST-like tools. Briefly, MMseqs2 had the lowest error rates among the programs tested and DIAMOND (ultra-sensitive mode) offered the best balance between speed and quality. Then you could filter the mapping reads by some parameters like e-value, identity and so on. The remaining reads could be used to calculate gene abundance by [RPKM, TPM or other methods](https://github.com/wwood/CoverM/#calculation-methods).

### _5.1 Abundance calculation of ORFs_

I extract representative genes' name and nucleotide sequences, because I apply a stratagy based on depth-of-coverage by dividing the summed depth per base by the length of the respective sequence, which is same as _[Salazar et al. (2019)](https://www.sciencedirect.com/science/article/pii/S009286741931164X#sec3)_ and _[Herold et al. (2020)](https://www.nature.com/articles/s41467-020-19006-2#Sec8)_.

To generate gene index, extract representative nucleotide sequence first.

```
./seqkit seq -n clusterRes_rep_seq.fasta > gene_rep.txt
sed -i 's/ //g' gene_rep.txt
./seqkit grep -f gene_rep.txt nucl_rename.fnn > gene_rep.fa
```

```
export PATH=$PATH:/share/home/bmalab/bwa-master/
bwa index gene_rep.fa
```

```
export PATH=$PATH:/share/home/bmalab/samtools-1.11/
bwa mem -t 28 gene_rep.fa S1_R1.fq S1_R2.fq|samtools sort -O bam -@ 28 -o - > S1_gene.bam
```

Parameters is same as _[Woodcroft et al. (2018)](https://www.nature.com/articles/s41586-018-0338-1#Sec8)_, which is much stricter than _[Salazar et al. (2019)](https://www.sciencedirect.com/science/article/pii/S009286741931164X#sec3)_ and _[Herold et al. (2020)](https://www.nature.com/articles/s41467-020-19006-2#Sec8)_. The relative abundance of each gene in each sample could be calculated as its coverage divided by the total coverage of all genes. Other [coverage measure methods](https://github.com/wwood/CoverM/#calculation-methods) could also be applied.

```
export PATH=$PATH:/public/home/bma/application/coverm-x86_64-unknown-linux-musl-0.2.0-alpha7/
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b S1_gene.bam -o S1_gene_filter.bam -t 28

coverm contig --trim-max 90 --trim-min 10 --bam-files S1_gene_filter.bam  > S1_gene_coverage.csv
```

> Woodcroft B J, Singleton C M, Boyd J A, et al. [Genome-centric view of carbon processing in thawing permafrost](https://www.nature.com/articles/s41586-018-0338-1#Sec8). _Nature_, 2018, 560(7716): 49-54.

> Salazar G, Paoli L, Alberti A, et al. [Gene expression changes and community turnover differentially shape the global ocean metatranscriptome](https://www.sciencedirect.com/science/article/pii/S009286741931164X#sec3). *Cell*, 2019, 179(5): 1068-1083. e21.

> Herold M, Arbas S M, Narayanasamy S, et al. [Integration of time-series meta-omics data reveals how microbial ecosystems respond to disturbance](https://www.nature.com/articles/s41467-020-19006-2#Sec8). *Nature communications*, 2020, 11(1): 1-14.

> Hernández-Salmerón J E, Moreno-Hagelsieb G. [Progress in quickly finding orthologs as reciprocal best hits: comparing blast, last, diamond and MMseqs2](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07132-6). _BMC genomics_, 2020, 21(1): 1-9.

### _5.2 Gene annotation_

---

## 6. Binning

---

## 7. Virus prediction
