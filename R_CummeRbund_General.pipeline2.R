#source("https://bioconductor.org/biocLite.R")
#biocLite("cummeRbund")

library(cummeRbund)
library(gridExtra)
library(pheatmap)

args = commandArgs(trailingOnly=TRUE)


flags = args[which(substr(args,1,1)=="-")]

File_Check.bool = F
GTF_Check.bool = F
Gene_List.bool = F
Output_Check.bool = F
Genome_Check.bool = F
Rebuild.bool = T
Sample_Pairs.bool = F

if(any(grepl("-i",tolower(flags))))
{
  File_Check.bool = T
  filename = args[which(tolower(args)=="-i")+1]
}
if(any(grepl("-g",tolower(flags))))
{
  GTF_Check.bool = T
  GTF_filename = args[which(tolower(args)=="-g")+1]
}
if(any(grepl("-l",tolower(flags))))
{
  Gene_List.bool = T
  Gene_List = args[which(tolower(args)=="-l")+1]
}
if(any(grepl("-p",tolower(flags))))
{
  Sample_Pairs.bool = T
  Sample_Pairs_List = args[which(tolower(args)=="-p")+1]
  if(substr(as.character(Sample_Pairs_List),1,1)=="-") Sample_Pairs_List = ""
}
if(any(grepl("-o",tolower(flags))))
{
  Output_Check.bool = T
  out_folder = args[which(tolower(args)=="-o")+1]
}
if(any(grepl("-r",tolower(flags))))
{
  Genome_Check.bool = T
  Set_Genome = tolower(args[which(tolower(args)=="-r")+1])
}
if(any(grepl("-nb",tolower(flags))))
{
  Rebuild.bool = F
}


if(File_Check.bool == T)
{
  if(substr(filename,nchar(filename),nchar(filename)) == '/') filename = substr(filename,1,(nchar(filename)-1))
}

if((File_Check.bool == T) & (GTF_Check.bool == T) & (GTF_Check.bool == T))
{
  if(Rebuild.bool==T)  cuff_data<-readCufflinks(filename,gtfFile = GTF_filename ,genome = Set_Genome,rebuild=T)
  if(Rebuild.bool==F)  cuff_data<-readCufflinks(filename,gtfFile = GTF_filename ,genome = Set_Genome)
}

if(File_Check.bool == T & (GTF_Check.bool == F | GTF_Check.bool == F))
{
  if(Rebuild.bool==F) cuff_data <- readCufflinks(filename)
  if(Rebuild.bool==T) cuff_data <- readCufflinks(filename,rebuild=T) 
}


if(File_Check.bool == T)
{

if(Output_Check.bool == F) out_folder = paste("filename_",output,sep="")

sysCommand = paste("mkdir",out_folder)
system(sysCommand)

pdf(paste(out_folder,"/Summary.pdf",sep=""))

print("Performing General Comparison")

Dendro = csDendro(genes(cuff_data)) 
plot(Dendro)
csScat =  csScatterMatrix(genes(cuff_data))
plot(csScat)
csVol =csVolcanoMatrix(genes(cuff_data))
plot(csVol)
dispP = dispersionPlot(genes(cuff_data)) + ggtitle("Dispersion Plot (Quality Check)")
plot(dispP)
csDH = csDistHeat(genes(cuff_data)) + ggtitle("Distance heat map between samples")
plot(csDH)
PCAP = PCAplot(genes(cuff_data),"PC1","PC2") + ggtitle("PCA Plot with all genes")
plot(PCAP)
MDSP = MDSplot(genes(cuff_data)) + ggtitle("MDS Plot")
plot(MDSP)

SMGenes = sigMatrix(cuff_data,level='genes',alpha=0.05)
plot(SMGenes)
SMIsoforms = sigMatrix(cuff_data,level='isoforms',alpha=0.05)
plot(SMIsoforms)

Sample_names = samples(cuff_data)[,2]


print("Determining Expressed Genes")

gene_diff_data <- diffData(genes(cuff_data))
sig_gene_data <- subset(gene_diff_data,(significant=='yes'))
gene_diff_data.FC2.5rpkm = subset(gene_diff_data,((value_1 > 5 | value_2 > 5) & abs(log2_fold_change) > 1)) 

#gene_diff_data.FC2.5rpkm[,1]
myGenes.FC2.5rpkm = getGenes(cuff_data,gene_diff_data.FC2.5rpkm[,1])
myGenes.FC2.5rpkm.annotations = annotation(myGenes.FC2.5rpkm)

myGenes.FC2.5rpkm.names = myGenes.FC2.5rpkm.annotations$gene_short_name
myGenes.FC2.5rpkm.names = subset(myGenes.FC2.5rpkm.names,!is.na(myGenes.FC2.5rpkm.names))

myGenes.FC2.5rpkm.names = paste(myGenes.FC2.5rpkm.names,collapse = " ")
myGenes.FC2.5rpkm.names = gsub(","," ",myGenes.FC2.5rpkm.names)
myGenes.FC2.5rpkm.names = unlist(strsplit(myGenes.FC2.5rpkm.names,split=" "))
myGenes.FC2.5rpkm.names = unique(myGenes.FC2.5rpkm.names)


Temp.AllGenes = data.frame(unique(gene_diff_data$gene_id))
for(i in 1:(length(Sample_names)-1))
{
  Temp.RPKM = subset(gene_diff_data,(sample_1 == Sample_names[i]))[c(1,2,5)]
  Temp.AllGenes = data.frame(Temp.AllGenes,unique(Temp.RPKM)[,3])
}
Temp.RPKM = subset(gene_diff_data,(sample_2 == Sample_names[(i+1)]))[c(1,3,6)]
Temp.AllGenes = data.frame(Temp.AllGenes,unique(Temp.RPKM)[,3])
colnames(Temp.AllGenes) = c("Gene_ID",Sample_names)
annotated.genes = annotation(genes(cuff_data))
annotated.genes.no.dup = annotated.genes[!duplicated(annotated.genes[1]),]
gene_symbols = annotated.genes.no.dup$gene_short_name
Data.AllGenes = cbind(gene_symbols,Temp.AllGenes)

print("check1")
Data.ExpressedGenes = Data.AllGenes[apply(Data.AllGenes[, c(-1,-2)], MARGIN = 1, function(x) any(x > 5)), ]
print("check2")

print("Writing Files")
write.table(Data.AllGenes, file = paste(out_folder,'/Data.All.Genes.txt',sep=""),sep = "\t")
write.table(Data.ExpressedGenes, file = paste(out_folder,'/Data.Expressed.Genes.txt',sep=""),sep = "\t")
write.csv(Data.AllGenes, file = paste(out_folder,'/Data.All.Genes.csv',sep=""))
write.csv(Data.ExpressedGenes, file = paste(out_folder,'/Data.Expressed.Genes.csv',sep=""))
write(myGenes.FC2.5rpkm.names, file = paste(out_folder,'/All.Genes.Names.FC2.5rpkm.txt',sep=""))

dev.off()


####################-Pair_wise Comparison-##############
print("Performing Pairwise Comparison")


if(Sample_Pairs.bool == T)
{
  Compare_input = as.character(Sample_Pairs_List)
  if(nchar(Sample_Pairs_List)<3)
  {
    print("Which Samples Would you like to Compare? (example: 2,3;5,1;2,3)")
    for(i in 1:length(Sample_names))print(paste(i,": ",Sample_names[i],sep=""))

    Compare_input <- readline("Enter samples You would like to compare using a ';' between each pair and a ',' between each sample. Example: 1,2;1,3;1,2")
  }

  Compare_input.pairs = unlist(strsplit(Compare_input, ";"))
  Compare_input.pairs.set.total =c()  
  for(i in 1:length(Compare_input.pairs))
  {
    Compare_input.pairs.set = unlist(strsplit(Compare_input.pairs[i], ","))
    if(Compare_input.pairs.set > 2)
    {
       print("Paired List Input Error")
       stop("Program Executed. Please Check paired Lists")
    } 
    Compare_input.pairs.set.total = cbind(Compare_input.pairs.set.total,Compare_input.pairs.set)
  }
  paired.combs = Compare_input.pairs.set.total


}
if(Sample_Pairs.bool == F) paired.combs = combn(length(Sample_names),2)





for(mm in 1:length(paired.combs[1,]))
{
  ## Do a pairwise comparison between experimental conditions ##


  Sample1 = Sample_names[paired.combs[1,mm]]
  Sample2 = Sample_names[paired.combs[2,mm]]



#pdf(paste(out_folder,"/Pairwise.",Sample1,".vs.",Sample2,".pdf",sep=""))

#gene_diff_data.FC2.5rpkm

 #csVolcano(genes(cuff_data),Sample1,Sample2)
 #csScatter(genes(cuff_data),Sample1,Sample2)


 #Paired.gene_diff_data.FC2.5rpkm = subset(gene_diff_data.FC2.5rpkm,((sample_1 == Sample1 & sample_2 == Sample2)|(sample_1 == Sample2 & sample_2 == Sample1)))
 Paired.gene_diff_data.FC2.5rpkm = subset(gene_diff_data.FC2.5rpkm,(((sample_1 == Sample1 & sample_2 == Sample2)|(sample_1 == Sample2 & sample_2 == Sample1))&abs(log2_fold_change)>1 & p_value < 0.05))

 Paired.genes.diff = getGenes(cuff_data,Paired.gene_diff_data.FC2.5rpkm$gene_id)
 Paired.gene_diff_data.FC2.5rpkm.annotation = annotation(Paired.genes.diff)
 Paired.gene.FC2.5rpkm.names = Paired.gene_diff_data.FC2.5rpkm.annotation$gene_short_name
 Paired.gene.FC2.5rpkm.names = subset(Paired.gene.FC2.5rpkm.names,!is.na(Paired.gene.FC2.5rpkm.names))
 Paired.gene.FC2.5rpkm.names = unique(unlist(strsplit(gsub(","," ",paste(Paired.gene.FC2.5rpkm.names,collapse = " ")),split=" ")))

 

 if(length(Paired.gene.FC2.5rpkm.names)>0)
 {

  pdf(paste(out_folder,"/Pairwise.",Sample1,".vs.",Sample2,".pdf",sep=""))

  unique.paired.Genes = getGenes(cuff_data,Paired.gene.FC2.5rpkm.names)
  unique.paired.Genes.diff = diffData(unique.paired.Genes)
  unique.paired.Genes.diff = subset(unique.paired.Genes.diff,(((sample_1 == Sample1 & sample_2 == Sample2)|(sample_1 == Sample2 & sample_2 == Sample1))&abs(log2_fold_change)>1 & p_value < 0.05))


 if(nrow(unique.paired.Genes.diff) > 0)
 {
  unique.paired.Genes.diff.GENES = getGenes(cuff_data,unique.paired.Genes.diff$gene_id)
  unique.paired.Genes.diff.GENES.ANN = annotation(unique.paired.Genes.diff.GENES)
 
  

  unique.paired.Genes.diff.GENES.ANN = subset(unique.paired.Genes.diff.GENES.ANN,(!is.na(gene_short_name)))

 if(nrow(unique.paired.Genes.diff.GENES.ANN) > 0)
 {

  Diff.Genes.FC2.RPKM5.pval0.05 = cbind(unique.paired.Genes.diff.GENES.ANN$gene_short_name,unique.paired.Genes.diff[match(unique.paired.Genes.diff.GENES.ANN$gene_id,unique.paired.Genes.diff$gene_id),])


  colnames(Diff.Genes.FC2.RPKM5.pval0.05)[1] = "gene_symbol"
  #annotated.genes = annotation(unique.paired.Genes)
  #annotated.genes.no.dup = annotated.genes[!duplicated(annotated.genes[1]),]

  #gene_symbol = annotated.genes.no.dup[which(unique.paired.Genes.diff$gene_id %in% annotated.genes.no.dup$gene_id),]$gene_short_name
  #gene_symbol = annotated.genes.no.dup$gene_short_name

  #Diff.Genes.FC2.RPKM5.pval0.05 = cbind(gene_symbol,unique.paired.Genes.diff)

  NEW.pseudo.log2FC = log2((Diff.Genes.FC2.RPKM5.pval0.05$value_2+.1)/(Diff.Genes.FC2.RPKM5.pval0.05$value_1+.1))
  NEW.abs.diff = abs((Diff.Genes.FC2.RPKM5.pval0.05$value_2)-(Diff.Genes.FC2.RPKM5.pval0.05$value_1))
  row.names(Diff.Genes.FC2.RPKM5.pval0.05) = paste(Diff.Genes.FC2.RPKM5.pval0.05$gene_symbol,"|",Diff.Genes.FC2.RPKM5.pval0.05$gene_id,sep="")
  Diff.Genes.FC2.RPKM5.pval0.05 = cbind(Diff.Genes.FC2.RPKM5.pval0.05,NEW.pseudo.log2FC,NEW.abs.diff)
  Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC = Diff.Genes.FC2.RPKM5.pval0.05[order(abs(Diff.Genes.FC2.RPKM5.pval0.05$NEW.pseudo.log2FC),decreasing = T),]
  Diff.Genes.FC2.RPKM5.pval0.05.sorted.absDiff = Diff.Genes.FC2.RPKM5.pval0.05[order(abs(Diff.Genes.FC2.RPKM5.pval0.05$NEW.abs.diff),decreasing = T),]

  colnames(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC)[6] = as.character(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC[1,3])
  colnames(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC)[7] = as.character(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC[1,4])
  

  pheatmap(Diff.Genes.FC2.RPKM5.pval0.05.sorted.absDiff[,c(6:7)],cluster_rows=F,cluster_cols=F,annotation_names_row=T,main = paste(Sample1,"vs",Sample2,"(Sorted by FC)\n",Sample2,"Up-regulated genes: ",
             length(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change[Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change>0]),
             "\n",Sample2,"Down-regulated genes:",length(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change[Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change<0])))

   pheatmap(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC[,c(6:7)],cluster_rows=F,cluster_cols=F,annotation_names_row=T,main = paste(Sample1,"vs",Sample2,"(Sorted by FPKM difference)\n",Sample2,"Up-regulated genes: ",
            length(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change[Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change>0]),
            "\n",Sample2,"Down-regulated genes:",length(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change[Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change<0])))


  if(length(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$gene_symbol)>100)
  {
       pheatmap(Diff.Genes.FC2.RPKM5.pval0.05.sorted.absDiff[c(1:100),c(6:7)],cluster_rows=F,cluster_cols=F,annotation_names_row=T,main = paste(Sample1,"vs",Sample2,"(Sorted by FC)\nTop 100\n",Sample2,"Up-regulated genes: ",
             length(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change[Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change>0]),
             "\n",Sample2,"Down-regulated genes:",length(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change[Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change<0])))

       pheatmap(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC[c(1:100),c(6:7)],cluster_rows=F,cluster_cols=F,annotation_names_row=T,main = paste(Sample1,"vs",Sample2,"(Sorted by FPKM difference)\nTop 100\n",Sample2,"Up-regulated genes: ",
            length(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change[Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change>0]),
            "\n",Sample2,"Down-regulated genes:",length(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change[Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change<0])))

        pheatmap(Diff.Genes.FC2.RPKM5.pval0.05.sorted.absDiff[c(1:50),c(6:7)],cluster_rows=F,cluster_cols=F,annotation_names_row=T,main = paste(Sample1,"vs",Sample2,"(Sorted by FC)\nTop 50\n",Sample2,"Up-regulated genes: ",
             length(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change[Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change>0]),
             "\n",Sample2,"Down-regulated genes:",length(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change[Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change<0])))

       pheatmap(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC[c(1:50),c(6:7)],cluster_rows=F,cluster_cols=F,annotation_names_row=T,main = paste(Sample1,"vs",Sample2,"(Sorted by FPKM difference)\nTop 50\n",Sample2,"Up-regulated genes: ",
            length(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change[Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change>0]),
            "\n",Sample2,"Down-regulated genes:",length(Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change[Diff.Genes.FC2.RPKM5.pval0.05.sorted.FC$log2_fold_change<0])))



  }




  write.table(Diff.Genes.FC2.RPKM5.pval0.05,file= paste(out_folder,"/",Sample1,".vs.",Sample2,".Pval0.05.rpkm5.FC2.txt",sep=""),sep = "\t",quote=F )
  write.csv(Diff.Genes.FC2.RPKM5.pval0.05, file= paste(out_folder,"/",Sample1,".vs.",Sample2,".Pval0.05.rpkm5.FC2.csv",sep=""))

  if(nrow(unique.paired.Genes.diff.GENES.ANN)>=2)
  {
   print(paste(Sample1,Sample2))
   tryCatch({
    csHM = csHeatmap(unique.paired.Genes,cluster='both') + ggtitle(paste("Heatmap:",Sample1,"vs",Sample2))
    plot(csHM)
    csVP = csVolcano(unique.paired.Genes,Sample1,Sample2)
    plot(csVP)
    csScat_smooth = csScatter(unique.paired.Genes,Sample1,Sample2,smooth=T) + ggtitle(paste("Scatterplot:",Sample1,"vs",Sample2))
    plot(csScat_smooth)
    csScat_norm =  csScatter(unique.paired.Genes,Sample1,Sample2)+ ggtitle(paste("Scatterplot:",Sample1,"vs",Sample2))
    plot(csScat_norm)
          }) 
 }
     
  dev.off()

 }#f(nrow(unique.paired.Genes.diff.GENES.ANN) > 0)
 }#if(nrow(unique.paired.Genes.diff) > 0)
 }#if(length(Paired.gene.FC2.5rpkm.names)>0)

}#for(mm in 1:length(paired.combs[1,]))

########################### Genes of Interest ###########################


if(Gene_List.bool == T)
{

pdf(paste(out_folder,"/GeneList.analysis.pdf",sep=""))


 Gene_List2 = unlist(strsplit(Gene_List,split=",")) 
 Gene_List2 = toupper(Gene_List2)

 myGeneList = getGenes(cuff_data,Gene_List2)


 csHeatmap(myGeneList,cluster='both')
 csScatter(myGeneList,samples(cuff_data)[1,2],samples(cuff_data)[2,2])


 PCAplot(myGeneList,"PC1","PC2")

 if(length(Gene_List2)>5)
 {
   ic = csCluster(myGeneList,k=4)
   csClusterPlot(ic)
 }

 expressionBarplot(myGeneList,showErrorbars = F)
 ID_Table=cbind(annotation(myGeneList)$gene_id,annotation(myGeneList)$gene_short_name)
 plot1 = expressionPlot(myGeneList,showErrorbars = F)
 plot2 = tableGrob(ID_Table)
 grid.arrange(plot1, plot2, nrow=2,padding=0.1)



 for(i in 1:length(Gene_List2))
 {
  mygene <- getGene(cuff_data,Gene_List2[i])
 
  expBPT = expressionBarplot(mygene,showErrorbars = F)
  expBPTiso =  expressionBarplot(isoforms(mygene),showErrorbars = F)
  expPT = expressionPlot(mygene,showErrorbars = F)

  plot(expBPT)
  plot(expBPTiso)
  plot(expPT)

  
  mySimilar<-findSimilar(cuff_data,Gene_List2[i],n=10)
  ID_Similar_Genes =cbind(annotation(mySimilar)$gene_id,annotation(mySimilar)$gene_short_name)
  plot1 =expressionPlot(mySimilar,logMode=T,showErrorbars=F) + ggtitle(paste("Genes wiith similar expression pattern to",Gene_List2[i]))
  plot2 = tableGrob(ID_Similar_Genes)
  grid.arrange(plot1, plot2, ncol=2,padding=0.1,widths = c(10,5))
  
 
  if(GTF_Check.bool == T  & Genome_Check.bool == T)
  {
    trackList<-list()
    myStart<-min(features(mygene)$start)
    myEnd<-max(features(mygene)$end)
    myChr<-unique(features(mygene)$seqnames)
    genome<-Set_Genome
    ideoTrack <- IdeogramTrack(genome = genome, chromosome = myChr)
    trackList<-c(trackList,ideoTrack)
    axtrack<-GenomeAxisTrack()
    trackList<-c(trackList,axtrack)
    genetrack<-makeGeneRegionTrack(mygene)
    genetrack
  
    trackList<-c(trackList,genetrack)
    biomTrack<-BiomartGeneRegionTrack(genome=genome,chromosome=as.character(myChr),start=myStart,end=myEnd,name="ENSEMBL",showId=T)
    trackList<-c(trackList,biomTrack)
    conservation <- UcscTrack(genome = genome, chromosome = myChr,
                            track = "RefSeq Genes", table = "refGene",
                            from = myStart-2000, to = myEnd+2000, trackType = "DataTrack",
                            start = "start", end = "end", data = "score",
                            type = "hist", window = "auto", col.histogram = "darkblue",
                            fill.histogram = "darkblue", ylim = c(-3.7, 4),
                            name = "Conservation")
    trackList<-c(trackList,conservation)
    plotTracks(trackList,from=myStart-2000,to=myEnd+2000)
    plotTracks(trackList,from=myStart-20000,to=myEnd+20000)

  }#if(GTF_Check.bool == T  & Genome_Check.bool = T)
 }#for(i in 1:length(Gene_List2))

 dev.off()
}#if(Gene_List.bool == T) 
 
  

}#if(File_Check.bool == T)



