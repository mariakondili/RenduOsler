#__________________ Nb ofBases covered by certain Fold Coverage(by chr)____________________# 

### Another Interesting Output of GATK is a file which gives Nb.bases covered by specific Cut-off Coverage value.
### sample_nterval_statistics is a table which has for many different Depth values :0-500 the Nb of bases that are covered.
### These files are created by chromosome per Subject because origin bam file and exons were cut by chromosome.
### So, here we are going to search the Depth of Coverage for reads that are on the chromosomes where the causal genes are.
### From script: Calc_Gene_Cover_AllPlots_V1.R

Plot_dir <- "/media/maritilia/DATA/#2015/ICAN_RenduOsler/myData/Plots/"
data_dir <- "/media/maritilia/DATA/#2015/ICAN_RenduOsler/myData/"

Pat<-seq(1,5) 
D=50 # a by default value for plotting (can be given by user, too) 
Chrom_Depth_WGS <- list(matrix(0,nrow=4, ncol=5))
Chrom_Depth_WES <- list(matrix(0,nrow=4, ncol=5))

chrom <- c(7,9,12,18)

for (p in Pat ){
  Chrom_Depth_WGS[[p]] <- .calc.chrom.depth(chrom,D,p, Chrom_Depth, data_dir, data_type="WGS")
  Chrom_Depth_WES[[p]] <- .calc.chrom.depth(chrom,D,p, Chrom_Depth, data_dir, data_type="WES")
  
}

.calc.chrom.depth <-function(chrom,D,p,Chrom_Depth, data_dir, data_type="WGS") {
  for (k in 1:length(chrom)) {
      interval_stat <- read.table(paste(data_dir,"/DepthOfCov/",data_type,"_Sample_interval_statistics/39_",p,
                                  "_MarkDup_CovbyEx_",chrom[k],".sample_interval_statistics",sep=""), 
                                  as.is=TRUE)
    
      bases_cov <- interval_stat_g[1,D+1]
      Chrom_Depth[k,p]<- bases_cov
  }
  rownames(Chrom_Depth)= c("chr7", "chr9", "chr12","chr18")
  
  png(paste(Plot_dir,"/BaseCov_D>50_perCHR.png", sep=""))
    par(mfrow=c(1,2),oma=c(1,0,2,0) ,cex=1 ) #mar=c(4,2,2,1),
    barplot(t(Chrom_Depth), beside=TRUE, 
            names.arg=rownames(Chrom_Depth), 
            main=data_type, ylab="per Base Coverage")   
    title("Nucleotide Coverage for chromosomes of interest",outer=TRUE)
  dev.off()
}

