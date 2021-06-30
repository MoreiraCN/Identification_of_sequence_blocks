### Identification of sequence blocks

The following pipeline was used to identify B and sex chromosomes sequences on the genome of the species *Holochilus sciureus* (2n = 56, NF = 56), a Neotropical rodent of Oryzomyini tribe. This same pipeline was already used to identify sequence blocks of B chromosome on the cichlid fish *Astatotilapia latifasciata* ([Valente et al. 2014](https://pubmed.ncbi.nlm.nih.gov/24770715/)), and on the *Astyanax* fish ([Ahmad et al. 2020](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07072-1)).

- Softwares used:

[Python 2.7](https://www.python.org/download/releases/2.7/).

[bedtools-v2.29.2](https://bedtools.readthedocs.io/en/latest/).

[Sushi: An R/Bioconductor package for visualizing genomic data](https://www.bioconductor.org/packages/release/bioc/vignettes/Sushi/inst/doc/Sushi.pdf).

Pipeline developed by [Ivan R Wolf](https://github.com/ivanrwolf/CovDetect/blob/master/LICENSE). Valente GT, Conte MA, Fantinatti BEA, Cabral-de-Mello DC, Carvalho RF, Vicari MR, Kocher TD, Martins C (2014) Origin and evolution of B chromosomes in the cichlid fish *Astatotilapia latifasciata* based on integrated genomic analyses. Mol Biol Evol, 31(8):2061-2072.

**Step 1 > Alignment (filtered libraries 0B, 1B, and probes against assembled genome 0B):**

[Align the filtered libraries against the assembled genome](https://github.com/MoreiraCN/Genomic_alignment).

**Step 2 > Coverage calculation per base:**

/[bedtools-v2.29.2](https://bedtools.readthedocs.io/en/latest/) genomecov -ibam alignment_0B_0B.sorted.bam -d > per_base_coverage_0B_0B.bed

/[bedtools-v2.29.2](https://bedtools.readthedocs.io/en/latest/) genomecov -ibam alignment_0B_1B.sorted.bam -d > per_base_coverage_0B_1B.bed

/[bedtools-v2.29.2](https://bedtools.readthedocs.io/en/latest/) genomecov -ibam alignment_0B_probe.sorted.bam -d > per_base_coverage_0B_probe.bed

**Step 3 > Merge two *.bed* files, row by row, into a single file:**

awk 'NR==FNR{a[$1,$2]=$3;next} ($1,$2) in a{print $0, a[$1,$2]}' per_base_coverage_0B_1B.bed per_base_coverage_0B_0B.bed > merged_1B_0B.txt

awk 'NR==FNR{a[$1,$2]=$3;next} ($1,$2) in a{print $0, a[$1,$2]}' per_base_coverage_0B_probe.bed per_base_coverage_0B_0B.bed > merged_probe_0B.txt

**Step 4 > Discard genomic sites with less than 15x coverage ($3 and $4 are the columns which represent the per base coverage of 0B and 1B or  0B and probe, respectively):**

awk '$3 > 15 && $4 > 15' merged_1B_0B.txt > merged_1B_0B_greater15x.txt

awk '$3 > 15 && $4 > 15' merged_probe_0B.txt > merged_probe_0B_greater15x.txt

**Step 5 > Calculate 1B/0B and probe/0B ratio:**

'awk -v OFS='\t' '{$5 = sprintf("%.3f", $4 / $3)}1' merged_1B_0B_greater15x.txt > merged_1B_0B_greater15x_ratio.txt'

awk -v OFS='\t' '{$5 = sprintf("%.3f", $4 / $3)}1' merged_probe_0B_greater15x.txt > merged_probe_0B_greater15x_ratio.txt







**Step 6 > Extract genomic regions with sequence blocks having at least 2x greater coverage:**
Calculate the 2B/0B ratio and extract the genomic regions with 2B genome having atleast 2x times greater coverage than 0B this means the extra copies present on the B chromosome.

awk '$5 > 2' merged_bed_ratio.txt > merged_bed_ratio2xcoverage.txt

awk '$5 > 2' merged_bed_ratio.txt > merged_bed_ratio2xcoverage.txt




### Step 7 > Format the output file for identification of sequence blocks:

- awk -v OFS='\t' '{print $1, $2, $3, $4}' merged_bed_ratio2xcoverage.txt > merged_bed_inputfile.txt





### Step 8 > Identification of sequence blocks:

- python2.7 [CovDetect.py](https://github.com/ivanrwolf/CovDetect/blob/master/CovDetect.py) -bp 100 -stdv 2 merged_bed_inputfile.txt





### Step 9 > Select sequence blocks larger than 200bp:

- awk '$4 > 200' merged_bed_inputfile_STDV2_BP100.blocks.txt > blocks_larger_than_200bp.txt

The file *blocks_larger_than_200bp.txt* contain a list of scaffolds with sequence blocks larger than 200bp. In order to vizualize the sequence blocks, the following steps were performed:




### Step 10 > Obtaining the bed graph file:

- /[bedtools-v2.29.2](https://bedtools.readthedocs.io/en/latest/) genomecov -ibam [alignment.sorted.bam](https://github.com/MoreiraCN/Genomic_alignment) -g [assembly.fa](https://github.com/MoreiraCN/Assembling_Illumina_sequences) -bg > sample_name.bg




### Step 11 > Extract interesting scaffolds from bed graph file:

grep 'scaffold_number' sample_name.bg > samplename_scaffoldnumber.bg



### Step 10 > View the graphs with the Sushi library of Rstudio:

- Command line used:

#call sushi
library(Sushi)

####open_sequence_blocks
referencesample_scaffoldnumber = read.table(file="/referencesample_scaffoldnumber.bg")
targetsample_scaffoldnumber = read.table(file="/targetsample_scaffoldnumber.bg")
#plot_interesting_region
chrom = "scaffoldnumber"
chromstart = number_of_base_pair_to_start
chromend = number_of_base_pair_to_end
#plot_graph
plotBedgraph(targetsample_scaffoldnumber,chrom,chromstart,chromend,transparency=.70,color=SushiColors(2)(2)[2])
plotBedgraph(referencesample_scaffoldnumber,chrom,chromstart,chromend,transparency=.70,color=SushiColors(2)(2)[1],overlay=TRUE,rescaleoverlay=FALSE)
labelgenome(chrom,chromstart,chromend,n=5,scale="Mb")
mtext("Read Depth",side=5,line=4,cex=1,font=5)> axis(side=2,las=1,tcl=.2)

