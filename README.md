### Identification of sequence blocks

The following pipeline was used to identify B and sex chromosomes sequences on the genome of the species *Holochilus sciureus* (2n = 56, NF = 56), a Neotropical rodent of Oryzomyini tribe. This same pipeline was already used to identify sequence blocks of B chromosome on the cichlid fish *Astatotilapia latifasciata* ([Valente et al. 2014](https://pubmed.ncbi.nlm.nih.gov/24770715/)), and on the *Astyanax* fish ([Ahmad et al. 2020](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07072-1)).

- Softwares used:

[Python 2.7](https://www.python.org/download/releases/2.7/).

[bedtools-v2.29.2](https://bedtools.readthedocs.io/en/latest/).

[Sushi: An R/Bioconductor package for visualizing genomic data](https://www.bioconductor.org/packages/release/bioc/vignettes/Sushi/inst/doc/Sushi.pdf).

Pipeline developed by [ivanrwolf](https://github.com/ivanrwolf/CovDetect/blob/master/LICENSE). Valente GT, Conte MA, Fantinatti BEA, Cabral-de-Mello DC, Carvalho RF, Vicari MR, Kocher TD, Martins C (2014) Origin and evolution of B chromosomes in the cichlid fish *Astatotilapia latifasciata* based on integrated genomic analyses. Mol Biol Evol, 31(8):2061-2072.

**Input data:**

- The *.sorted.bam* file resultant of the [align of filtered libraries against the assembled genome](https://github.com/MoreiraCN/Genomic_alignment).

### Step 1 > Alignment:

- [Align the filtered libraries against the assembled genome](https://github.com/MoreiraCN/Genomic_alignment).

### Step 2 > Coverage calculation per base:

- /[bedtools-v2.29.2](https://bedtools.readthedocs.io/en/latest/) genomecov -ibam alignment.sorted.bam -d > per_base_coverage.bed

### Step 3 > Merge *.bed* files into a single file:

- awk 'NR==FNR{a[$1,$2]=$3;next} ($1,$2) in a{print $0, a[$1,$2]}' per_base_coverage1.bed per_base_coverage2.bed > merged_bed.txt

### Step 4 > Discard genomic sites with less than 15X coverage:

- awk '$3 > 15 && $4 > 15' merged_bed.txt > merged_bed_greater15x.txt

### Step 5 > Calculate the sequence blocks coverage ratio:

- awk -v OFS='\t' '{$5 = sprintf("%.3f", $4 / $3)}1' merged_bed_greater15x.txt > merged_bed_ratio.txt

### Step 6 > Extract genomic regions with sequence blocks having at least 2x greater coverage:

- awk '$5 > 2' merged_bed_ratio.txt > merged_bed_ratio2xcoverage.txt

### Step 7 > Format the output file for identification of sequence blocks:

- awk -v OFS='\t' '{print $1, $2, $3, $4}' merged_bed_ratio2xcoverage.txt > merged_bed_inputfile.txt

### Step 8 > Identification of sequence blocks:

- python2.7 [CovDetect.py](https://github.com/ivanrwolf/CovDetect/blob/master/CovDetect.py) -bp 100 -stdv 2 merged_bed_inputfile.txt

### Step 9 > Select sequence blocks larger than 200bp:

- awk '$4 > 200' merged_bed_inputfile_STDV2_BP100.blocks.txt > blocks_larger_than_200bp.txt

The file *blocks_larger_than_200bp.txt* contain a list of scaffolds with sequence blocks larger than 200bp. In order to vizualize the sequence blocks, the following steps were performed:

### Step 10 > Obtaining of *.bg* files:

- /[bedtools-v2.29.2](https://bedtools.readthedocs.io/en/latest/) genomecov -ibam [alignment.sorted.bam](https://github.com/MoreiraCN/Genomic_alignment) -g [assembly.fa](https://github.com/MoreiraCN/Assembling_Illumina_sequences) -bg > sample_name.bg

### Step 11 > extrair os sccafolds de interesse o arquivo .bg:

### Step 10 > View the graphs with the Sushi library of Rstudio:

- Command line used:

