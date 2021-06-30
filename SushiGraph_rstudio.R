#call_sushi
library(Sushi)

####open_sequence_blocks
0B_0B_scaffoldnumber = read.table(file="/path_to_bg_file/0B_0B_scaffoldnumber.bg")
1B_0B_scaffoldnumber = read.table(file="/path_to_bg_file/1B_0B_scaffoldnumber.bg")
#plot_interesting_region
chrom = "scaffoldnumber"
chromstart = number_of_base_pair_to_start
chromend = number_of_base_pair_to_end
#plot_graph
plotBedgraph(1B_0B_scaffoldnumber,chrom,chromstart,chromend,transparency=.70,color=SushiColors(2)(2)[2])
plotBedgraph(0B_0B_scaffoldnumber,chrom,chromstart,chromend,transparency=.70,color=SushiColors(2)(2)[1],overlay=TRUE,rescaleoverlay=FALSE)
labelgenome(chrom,chromstart,chromend,n=5,scale="Mb")
mtext("Read Depth",side=5,line=4,cex=1,font=5)> axis(side=2,las=1,tcl=.2)


#call_sushi
library(Sushi)

####open_sequence_blocks
0B_0B_scaffoldnumber = read.table(file="/path_to_bg_file/0B_0B_scaffoldnumber.bg")
probe_0B_scaffoldnumber = read.table(file="/path_to_bg_file/probe_0B_scaffoldnumber.bg")
#plot_interesting_region
chrom = "scaffoldnumber"
chromstart = number_of_base_pair_to_start
chromend = number_of_base_pair_to_end
#plot_graph
plotBedgraph(probe_0B_scaffoldnumber,chrom,chromstart,chromend,transparency=.70,color=SushiColors(2)(2)[2])
plotBedgraph(0B_0B_scaffoldnumber,chrom,chromstart,chromend,transparency=.70,color=SushiColors(2)(2)[1],overlay=TRUE,rescaleoverlay=FALSE)
labelgenome(chrom,chromstart,chromend,n=5,scale="Mb")
mtext("Read Depth",side=5,line=4,cex=1,font=5)> axis(side=2,las=1,tcl=.2)
