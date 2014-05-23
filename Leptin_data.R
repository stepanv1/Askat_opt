#!/usr/bin/Rscript
#script to load snp and methilation data in chromosome1 into iBMQ-required 
#data objects
#cuts LINKAGE-style data (sex, affect, ets)

#gen_replace(x, rep) - a function to replace element x in with rep in gp$geno
#x- genotype to replace, rep- replacement of that genotype
#geno - matrix of genotypes
# function ids needed to recode function codeGeno to work
#on our data
##TODO: remove access to the global object "gp"  
gen_replace<-function(x, rep){      
  idx<-which(grepl(x, gp$geno), arr.ind=TRUE)   
  #print(idx)
  gp$geno[idx] <<- rep                     
}



######################################################
#1. Create kinship matrix from pedigree information, using 
Pheno <- as.data.frame(read.csv("IDs.Leptin.Phenotypes.csv"))[1:3960,]#[1:1000,]
Pheno$fam<-as.numeric(Pheno$fam)
Pheno$actual_zygosity<-as.numeric(Pheno$actual_zygosity)
#Pheno<-Pheno[-unique(Pheno$fam),]

nSample<-nrow(Pheno)

#if zigocity is 3 I put 0.25
kinPair<-function(i,j){
  if (i==j) return(0.5)
  if (Pheno[i,"fam"]==Pheno[j,"fam"]) {
    if (Pheno[i,"actual_zygosity"]==1) {return(0.5)}
    else {return(0.25)}
  }  
  return(0)
}

kinEst<-outer(1:nSample, 1:nSample, FUN = Vectorize(kinPair))
#indexes of the individuals missing phenotypes and missing zygosity
IDX_miss<-which((is.na(Pheno$leptin) | (Pheno$actual_zygosity == 3)), arr.ind=TRUE)

kinEst<-kinEst[-IDX_miss, -IDX_miss]




###################################################################
#preparing the snp data

haps<-read.table("TwinsUK.chr1.posLEPR.phasing.impute2_haps")
#remove five 5 st column with id and position information for creation of MAF-based genotype data
h<-haps[, -c(1:5)]
hapMAF<-matrix(NA, nrow(h), ncol(h)/2)
#hapMAF<-hapMAF[1:100, ]
for (i in 1:nrow(h)){
  for (j in 1:(ncol(h)/2)){
    hapMAF[i,j]<-abs(h[i,2*j-1]+h[i,2*j]-2)
  }
}
hapMAF<-cbind(haps[, c(1:5)], hapMAF )
#save(hapMAF, file="hapMAR.RData")

library(itertools)
library(doParallel)

#parallel version of the same code
hapMAFp<-matrix(NA, nrow(h), ncol(h)/2)
#register cluster
cl <- makeCluster(4)
registerDoParallel(cl, cores= 4)
#separate data in chunks
chunks <- getDoParWorkers()
hapMAFp<-foreach (these = isplitIndices(nrow(h),
                                        chunks=chunks),
                  .combine=rbind) %dopar% {
                    IDX<-these
                    a <- matrix(0, length(IDX), (ncol(h)/2))
                    for (i in 1:length(IDX)){  
                      for (j in 1:(ncol(h)/2)){
                        #a[i,j] = i*j
                        a[i,j]<-abs(h[IDX[i],2*j-1]+h[IDX[i],2*j]-2)
                      }
                    }  
                    return(a)
                  }
stopCluster(cl)
hapMAFp<-cbind(haps[, c(1:5)], hapMAFp )

#creating data file consisiting of phenotypes and genotypes
hapMAF_transposed<-t((hapMAF[, -c(1:5)])[, 1:3960])
PED<-cbind(c(1:nrow(hapMAF_transposed)), as.numeric(Pheno$leptin), hapMAF_transposed)
save()
PED<-as.matrix(PED[-IDX_miss, ])

#remove monozygous alleles or those with having just one minor allele in single sample 
low_MAF_IDX<-apply(PED, 2, function(col) !(all(col ==0 ) || (sum(col)==1) ) )
PED<-PED[, low_MAF_IDX]


write.table(file="PED_leptin.dat", x=PED, sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
save(kinEst, file="kinEst.RData")


#msnp<-snps[1:22, 7:26] #snippet fo a data set for testing purposes
#coding msnp AA, GT, T.. into 1, 2, 3 
library("plyr")
library("synbreed")
library("gtools")
# separate first and second allele
omsnp=(msnp)[seq(1, length(msnp), by = 2)]
emsnp=(msnp)[seq(2, length(msnp), by = 2)]
# remove numbers of alleles (e.g.: "1", "2")
names(omsnp)<-gsub("\\.[0-9]", "", names(omsnp))
names(emsnp)<-gsub("\\.[0-9]", "", names(emsnp))
# join homologous alleles back together now under the same SNP id
# form genotypeslooking like AT, GC etc. 
temp=cbind(omsnp, emsnp)
len=dim(temp)[2]
len
genotypes = data.frame(t(apply(temp, 1, function(x) paste(x[1:(len/2)], x[(len/2+1):len], sep="/"))))
names(genotypes)<-names(omsnp)
# adding rownames for synbreed package
row.names(genotypes) <- snps[,2] 
#checking if there is more than 4  values per column SNP
gt_stats<-lapply(1:(dim(genotypes)[2]), function(x) (table(genotypes[, x])))
gt_tbls<-lapply(1:length(gt_stats),  function(x) {!("0/0" %in% unlist(dimnames(gt_stats[[x]]))) & (length(gt_stats[[x]])>=4)})
if (length(which(unlist(gt_tbls)))>0) 
  cat("Warning: possibly more than 3 genotypes per SNP")
# create object of class ’gpData’ from synbreed
gp <- create.gpData(geno=genotypes)

#linkage disequilibrium calculation
library(genetics)

listofargs=list(list("0/0", NA))
gen_replace("0/0", NA)

#Here is the encoded genotype data
#check results of this with plink!
gp.coded <- codeGeno(gp, label.heter = "alleleCoding")

#correcting from 0,1,2 code to 1,2,3 required by iBMQ
gp.coded$geno <- gp.coded$geno+1
#gp$geno
#gp.coded$geno
summary(gp.coded)
#add the id column 
genom<-cbind(strtoi(row.names(gp.coded$geno)), gp.coded$geno)
colnames(genom)[1]<-"clientid"
genom<-as.data.frame(genom)

###########################################################
#preparing the methylation data
meth<-d.all[,1:96]

###########################################################
#match clientid's to have equal number of individuals in genotype and probes sets
#clean genom: if there are more tham 20% of empty snps -remove row
l<-floor(dim(genom)[2]/2*0.2)
genom<-as.data.frame(genom[rowSums(is.na(genom))<=l,])


idx <- intersect(genom[,1], meth[,1])
st= length(idx)
meth<-meth[meth$clientid %in% idx[1:st], ]#[, 1:10]
genom<-genom[genom$clientid %in% idx[1:st], ]#[, 1:10]
row.names(genom)<-NULL
############################################################
#applying iBMQ functions
library(iBMQ)
#goofy imputation for testing purposes
#g<-apply(genom, c(1,2), function(x){ if (is.na(x))  1 else x})
# Imputation is done for now with fastPHASE algorithm

# This is to put imputed genotypes back into genom variable.
# imputed_genoms_genotypes.out contains two haplotypes per individual,
# coded as 0 and 1. 
# Each genotype is transformed into the minor allele frequency notation:
# E.g., if if we have 0/0 genotype we map it into 0+0 + 1 (+1 is because
# of iBMQ's "123" coding)
tl <- scan("imputed_genoms_genotypes.out", what="", sep="\n")
#remove header?
tl<-tl[-c(1:19)]
#remove line separators
tl<-tl[tl!='# ']
#remove tail
tl<-tl[-length(tl)]

len<-length(tl)/2
for (x in 1:len){
  h1<-as.numeric(unlist(strsplit(tl[2*x-1]," ")))
  h2<-as.numeric(unlist(strsplit(tl[2*x]," ")))
  genom[x,-1] <- h1+h2+1
}
#impute the methilation data
library(mice) 
imp<-mice(meth[,-1])
temp<-complete(imp)
meth[,-1]<-temp

#meth<-apply(meth, c(1,2), function(x){ if (is.na(x))  0.5 else x})
#scale methilation levels with inverse logit function
meth[,-1]<-apply(meth[,-1], c(1,2), function(x) (logit(x, min=-0.0001, max=1.0001)))
rownames(genom)<-genom[,1]
genom<-genom[,-1]
#minor allel frequency threshold, discard those with MAF less tham 0.05
genom<-genom[sapply(genom, function(x) sum(x-1)/length(x)/2 )>0.05]
genom<-genom[sapply(genom, function(x) length(unique(x))>1 )]
hist(sapply(genom, function(x) var(x)))
rownames(meth)<-meth[,1]
meth<-meth[,-1]

#TEMP!!! to check contimuous covariate idea
#genom<-cbind(genom[-length(genom)], (rnorm(1175,2,1)))

genomQ <-new ('SnpSet',  call = t(genom))
methQ <- new("ExpressionSet", exprs = t(meth))
save(genomQ, methQ, file = "eqtlMcmc.RData")
#apply those to check empty rows/columns
#all(apply(meth[!complete.cases(meth),], 2, function(x) all(is.na(x)))==FALSE)
#all(apply(meth[!complete.cases(meth),], 1, function(x) all(is.na(x)))==FALSE)

PPA <- eqtlMcmc(genomQ, methQ, n.iter=200, burn.in=100, n.sweep=20, mc.cores=4, RIS=FALSE)
source("./eqtl.mcmc.R")
#PPA <- eqtl.mcmc(genomQ, methQ, n.iter=200, burn.in=100, n.sweep=20, nproc=4, RIS=FALSE)
###############################################################
#adding mapping data
#gene.data<-annot.cgs[,c(2,3,8,9)]
#gene.data[,4]<-gsub(';(.*)', "", gene.data[,4])
require(data.table)
gene.data<-annot.cgs[,c(2,3,8)]
#no gene names are needed for now
#gene.data[,4]<-gsub(';(.*)', "", gene.data[,4])
#experimenting with data.table format - powerful and fast sql-like data access
gene.data<-data.table(gene.data)
setkey(gene.data, TargetID)
#extract coordinates of the probes in meth 
z<-data.table(colnames(meth[ ,-1]))
cc<-gene.data[z, nomatch=0]

genemap<-data.frame(cc$TargetID, cc$CHR, cc$MAPINFO, cc$MAPINFO + 1)
snpmap2<-data.frame(snpmap$marker, snpmap$chr, snpmap$position)


#############################################################################

cutoff <- calculateThreshold(PPA, 0.1)
eqtl <- eqtlFinder(PPA, cutoff)
eqtltype <- eqtlClassifier(eqtl, snpmap2, genemap, 1000000)
library(ggplot2)
ggplot(eqtltype, aes(y=GeneStart, x=MarkerPosition)) +
  geom_point(aes(y=GeneStart, x=MarkerPosition,color = PPA))+
  facet_wrap(~PPA) + scale_colour_gradientn(colours=rainbow(3))+
  facet_grid(GeneChrm~MarkerChrm)+theme_bw(base_size = 12, base_family = "")  
#theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())+scale_x_reverse()

#taking subset of the data for one probe
#18432992
eqtltype2<-eqtltype[eqtltype$GeneStart>186653723-5000000,]
ggplot(eqtltype2, aes(y=GeneStart, x=MarkerPosition)) +
  geom_point(aes(y=GeneStart, x=MarkerPosition, color = PPA))+
  facet_wrap(~PPA) + scale_colour_gradientn(colours=rainbow(3))+
  facet_grid(GeneChrm~MarkerChrm)+theme_bw(base_size = 12, base_family = "")  
#theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())+scale_x_reverse()

#useful snippet for quering the result
eqtlDT<-data.table(eqtltype, key=c("SNP", "Gene"))
tt<-eqtlDT[PPA>0.2][, count:= .N, by=SNP][order(-count)][,list(SNP, Gene, count, PPA, MarkerPosition, GeneStart)]

#visualisation of data
#boxplot(meth[c("cg15898116","cg25147026")], main = "Methylation on two top probes",   ylab="Methylation logit")
library(LDheatmap)

cor_gen <- cor(genom)^2
LDheatmap(cor_gen, distances="genetic", add.map=F)

cor_meth <- cor(logit(meth))
LDheatmap(abs(cor_meth), distances="genetic", add.map=F, title="Correlation between probes")
