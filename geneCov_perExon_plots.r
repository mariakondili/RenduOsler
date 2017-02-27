#!usr/bin/R
#setwd("~/Documents/myR_projects/RenduOsler_git/")

##"" Receives genes of interest as input (found mutated in RenduOsler) and 
##   finds exons regions in SQL dB, compiles data of coverage measured from
##   GATK depth-of-coverage and shows a histogram with coverage of gene per exon ""

source("http://bioconductor.org/biocLite.R")
library(biomaRt)
library(DBI)
library(RMySQL)


CovData <- "/media/maritilia/DATA/#2015/ICAN_RenduOsler/myData/DepthOfCov/"
Exome_files  <- list.files(path=CovData, pattern="^Sorted_MDup_Exome_39_", all.files=TRUE)
Genome_files <- list.files(path=CovData, pattern="^Sorted_MDup_Genome_39_",all.files=TRUE)

WE_covtab <-lapply(paste(CovData,Exome_files,sep=""), read.table, as.is=TRUE,header=TRUE)
WG_covtab <-lapply(paste(CovData,Genome_files,sep=""),read.table, as.is=TRUE,header=TRUE)
#Each elem of list WE,WG refers to one Patient

geneID <-  c("ALK","ENG","SMAD4","ACVRL1")
dB <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") 
filters <- listFilters(dB)
attrib  <- listAttributes(dB)

for ( g in geneID) { 
  # Search the description (given by : attrib[5,1] )
  des <- as.character(getBM(filters="external_gene_name", attributes=attrib[5,1],values=g, mart=dB)) 
  des <- strsplit(des,'\\[\\Source')[[1]][1]
  
  cat("Gene Description :",des,"\n\n")
  
  ## mySQL Connection with 1000Genomes-> Extract start-end position from description
  mydb<- dbConnect(MySQL(), user='anonymous',password='',dbname='homo_sapiens_core_73_37',
                   host='mysql-db.1000genomes.org',port=4272 )
  
  ##> Search and get output directly ( sendQuery needs 'Fetch' to get output)
  command<- paste('select description, seq_region_start,seq_region_end from gene where description like "%',des,'%"', sep="")
  query  <- suppressWarnings (dbGetQuery(mydb,command)) ##** Attention: Warnings
  cat( "~~> The query in the 1000Genomes Database gives the following results :\n")
  print(query)
  
  chosen_query <- ifelse(nrow(query)==1, 1, 1)  #Always first query corresponds to gene
  
  gene.start <- query[chosen_query,2]
  gene.end   <- query[chosen_query,3]

  ##>Extract exons regions(same for both tables WG,WE) 
  exPairs  <- gsub( "(^\\d|\\d\\d):","",WE_covtab[[1]][,1] ) # ex: "114503825-114503944"
  
  exStarts <- as.integer( unlist(lapply(strsplit (exPairs,"-"),"[",1)) ) 
  exEnds   <- as.integer( unlist(lapply(strsplit (exPairs,"-"),"[",2)) )
  s <- which(exStarts >= as.integer(gene.start))
  e <- which(exEnds   <= as.integer(gene.end  )) 
  
  exonsOfGene <- intersect(s,e) # Exons(pos in exStarts) inside gene Region
  
  ##>> Coverage of gene by Exons (each bar = an exon)
  
  if (length(exonsOfGene) >20) {
    #chosenExons <- readline("The Gene contains >20 exons.Please choose which exons you wish to visualise.
    #                         Write a-b for an interval or a-a for the exon a): ")
    chosenExons <- "1-20"
    chosenExons <- as.integer(unlist(strsplit(chosenExons,"-")))
    chosenExons <- seq(chosenExons[1], chosenExons[2])
    
    .create.plots(p, g,WG_covtab, WE_covtab, chosenExons,exonsOfGene)
    
  }else {
   .create.plots(p,g, WG_covtab, WE_covtab, seq(1,length(exonsOfGene)),exonsOfGene)
  }
  
  .density.dsb(p,WE_covtab, WG_covtab)

} #end.for ea gene
  
                          ##   F U N C T I O N S   ## 
.create.plots <- function(p, WG_tab, WE_tab,chosenExons) {
  Plot_dir <- "/media/maritilia/DATA/#2015/ICAN_RenduOsler/myData/Plots/"
  png(filename=paste(Plot_dir,"GeneCoverage/",g,"_WG-WE_Subject_",p,".png",sep=""))
  par(mfrow=c(1,2))
  .create.exon.barplot(WG_tab, whatData="WG",chosenExons,g,exonsOfGene,"mediumblue")
  .create.exon.barplot(WE_tab, whatData="WE",chosenExons,g,exonsOfGene,"olivedrab")
  mtext(paste("Gene Coverage of Subject ",p,sep=""),side = 3, at =-2, line = 0.8, cex=1.5)
  dev.off()
}

.create.exon.barplot <- function(seqData, whatData="",chosenExons, g, exonsOfGene, color) {
    mc <-  sum(seqData[exonsOfGene,3]) / length(exonsOfGene)
    mc <-  format(round(mc,2),nsmall=2)
    barplot (seqData[exonsOfGene[chosenExons],3],col=color, ylab="Coverage",
             names.arg=paste("Exons of Gene",g,"in", whatData, sep=" "),
             sub=paste("Mean Cov=",mc, sep="") ) 
}

.common.exons <- function(exom_dat, genom_dat) {
  # find how many exons are common in 1st dataset from 2nd.Return "which" of 2nd found in 1st.
  cov_ex0_E <- which(exom_dat[,3] == 0) # Coverage==0
  exons_E   <- exom_dat[cov_ex0_E,1]
  ex.of.genom <- which(genom_dat[,1] %in% exons_E )
  return(ex.of.genom)
}

.density.dsb <- function(p, seqData1, seqData2) {
  
  ##>> DENSITY DISTRIBUTION of WG_Cov for the exons where WE_Cov=0 (For each Patient) 
  
  ex.of.dat2 <- .common.exons(seqData1, seqData2)
  ex.of.dat1 <- .common.exons(seqData2, seqData1)
  
  x11() 
  plot(density(seqData2[ex.of.dat2,3]),col="darkblue",xlim=range(0,150),xlab="",
       main=paste("Distribution of Exons Coverage (Subject No",p,")",sep=""))
  lines(density(seqData1[ex.of.dat1,3]),col="red") # Overlapped curves, BUT Exome curve is very thin around 0 !!
  legend("topright", c("Cov on WGS data","Cov on WES data"), col=c("darkblue","red"),pch=19)
}



