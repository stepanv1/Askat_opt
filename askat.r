#!/usr/bin/Rscript
#-------------------------------------------------------------------------------
# ASKAT Package
#
# By: Karim Oualkacha
# Optimized: Stepan Grinek      																	
#-------------------------------------------------------------------------------

library(GenABEL)
library(nFactors)
library(CompQuadForm)
library(Matrix)
#
#-------------------------------------------------------------------------------
# STEPs1-2 Checking for missing data and running fastlmm under the null model
#
# Ped    : Is the Pedigree data file. It has subject IDs as first column (IDs 
#          should be differents for all subjects), phenotype as second column 
#          and region-based SNPs that will be analized together  
#
# kin1   : The kinship matrix. If the kinship matrix noted is calculated from 
#          GenABEL we should convert it to the kinship matrix using this 
#          command: 
#              kin1 = 2*diagReplace(kin1, upper=TRUE) ############## Notice that now the kinship matrix is 2 TIMES the GenABEL kinship matrix 
#
# Missing: A logical parameter, if TRUE, it means that there is missing values 
#          in the Pedigree data set: ASKAT function moves out subjects with 
#          missing data
#
#-------------------------------------------------------------------------------
STEPs12 <- function(Ped, kin1, Missing, tmpDir="." ) {
	#####  STEP 1: Check for Missing Data and construction of Pedigree without Missing Data   #####
	##### Also creation of the file 'file.name.data' file, which will be needed by FaST-LMM #####

	file.name.data = tmpFile("kin.FaST.data.miss.txt");		# Data file for Fast-LMM
	eigenDir <- tmpFile("ASKAT_FaSTLMM");						# Directory for eigenvalue output (Fast-LMM)

	# Create file.name.data 
	if( Missing ) {
		# Create file.name.data when having missing values
		data.Without.NA = pheno.geno.kin.Without.NA(Ped, kin1, file.name.data)
		Ped = data.Without.NA$pheno.Geno.No.NA
	} else {
		# Create file.name.data (no missing values)
		kin.FaST = kin1
		n.col = Ped[,1]
		kin.FaST = cbind(n.col,kin.FaST)
		kin.FaST = rbind(c("var",n.col),kin.FaST)

		if( debug )	{ cat('Writing kinship matrix to file: ', file.name.data , '\n' ); }
		write(t(kin.FaST[ ,1:(dim(kin.FaST)[1])]), file = file.name.data, ncolumns = dim(kin.FaST)[1], sep= "\t")
	}

#-------------------------------------------------------------------------------
	##### STEP 2: Under the NULL we call FaST-LMM to estimate the VCs ####
	res.VC.FaST = VC.FaST.LMM(Ped, eigenDir, file.name.data)
	Estim.Sigma.RG = res.VC.FaST$Polygenic
	Estim.Sigma.e = res.VC.FaST$Env
	S = res.VC.FaST$S
	U = res.VC.FaST$U

	return( list(Estim.Sigma.RG = Estim.Sigma.RG, Estim.Sigma.e = Estim.Sigma.e, U = U, S = S) );
}
#-------------------------------------------------------------------------------
# ASKAT Main function
#
# Ped    : Is the Pedigree data file. It has subject IDs as first column (IDs 
#          should be differents for all subjects), phenotype as second column 
#          and region-based SNPs that will be analized together.  
#
# Estim.Sigma.RG  : The polygenic variance component that is calculated from 
#                   the null model using fastlmm. This parameter is estimated 
#                   in STEPs1-2.
#
# Estim.Sigma.e   : The environment variance component that is calculated from 
#                   the null model using fastlmm. This parameter is estimated 
#                   in STEPs1-2.
#
# U      : Orthogonale matrix obtained from the spectrale decomposition of the 
#          kinship matrix:  Kinship = U*S*Ut. This matrix is calculated also in STEPs1-2.
#          This matrix is obtained from fastlmm once we adjust the null model in STEPs1-2.
#
# S      : Diagonal matrix obtained from the spectrale decomposition of the 
#          kinship matrix:  Kinship = U*S*Ut. This matrix is calculated also in STEPs1-2.
#          This matrix is obtained from fastlmm once we adjust the null model in STEPs1-2.
#
#-------------------------------------------------------------------------------
ASKAT <- function(Ped, Estim.Sigma.RG, Estim.Sigma.e, S, U, pACC=1e-9) {
	#####
	Y.trait = Ped[,2]
	X = as.matrix(Ped[,3:dim(Ped)[2]])

##### STEP 1: Calculation of weights matrix W and the matrix K =GWG #####
freq.MAF = apply(X, 2, mean)/2
    Geno.no.sparse = which(!freq.MAF==0)
    X = X[,Geno.no.sparse]
    freq.MAF = apply(X, 2, mean)/2


#start of the trick to reduce the dimensionality of the matrix
if( length(freq.MAF) == 1){
  w <- dbeta(freq.MAF, 1, 25)
  K.sqrt <- w * t(X)
} else {
  w <- vector(length = length(freq.MAF))
  for (i in 1:length(freq.MAF)){
    w[i] <- dbeta(freq.MAF[i], 1, 25)
  }
  w <- diag(w)
  K.sqrt <- w %*% t(X)
}

##### STEP 2: ASKAT score test statistic calculations #####


Gamma = Estim.Sigma.RG / Estim.Sigma.e
UT<-t(U)
D.0 <- (Gamma * S)  + diag(1, dim(X)[1], dim(X)[1])
inv.sqrt.D.0 <- diag(1/sqrt(diag(D.0)))
inv.D.0 = diag(1/diag(D.0))
un.n <- c(rep(1,dim(U)[1]))
Z <- 1/(( t(un.n) %*% (U %*% (inv.D.0 %*% (UT %*% un.n))))[1,1]) 
X.tilde <- inv.sqrt.D.0 %*% (UT %*% un.n)
Y.tilde <- inv.sqrt.D.0 %*% (UT %*% Y.trait)
s2 = Estim.Sigma.e
P.0.tilde = (diag(1, dim(U)[1], dim(U)[2])) - (X.tilde %*% Z) %*% ((t(X.tilde)))

RM <- ((K.sqrt %*% U) %*% inv.sqrt.D.0) %*%  P.0.tilde# 44 s
W <- RM %*% t(RM)

Q <- (t(Y.tilde) %*% t(RM) %*% ((RM %*% Y.tilde)))/(2 * s2)




out = Get_PValue.Modif(W/2, Q, pACC)
#print("new"); print(out$p.value)
pvalue.davies <- out$p.value
lambda <- out$lambda
is_converge <- out$is_converge

	return( list(pvalue.ASKAT = pvalue.davies, Q.ASKAT = Q, Polygenic.VC = Estim.Sigma.RG, Env.VC = Estim.Sigma.e, lambda = lambda, is_converge = is_converge) );
}

#-------------------------------------------------------------------------------
# Write a matrix for FastLmm
#-------------------------------------------------------------------------------
Geno.FaST.LMM <- function(Geno, tpedFile){
	Geno.FaST = vector(length = 2*length(Geno))
	for (i in 1:length(Geno)){
		if (Geno[i] == 0){
			Geno.FaST[((2*(i-1))+1)] = 4
			Geno.FaST[(2*i)] = 4
		}

		if (Geno[i] == 1){
			Geno.FaST[((2*(i-1))+1)] = 1
			Geno.FaST[(2*i)] = 4
		}

		if (Geno[i] == 2){
			Geno.FaST[((2*(i-1))+1)] = 1
			Geno.FaST[(2*i)] = 1
		}
	}
	geno.test = c(1, "snp0", 0, 0, 1, Geno.FaST)
	geno.test = rbind(geno.test)

	if( debug )	{ cat('Writing tped matrix to file: ', tpedFile , '\n' ); }
	write.matrix(geno.test, file = tpedFile, sep="")
}

#-------------------------------------------------------------------------------
# Get lambda UPD: heuristics: "IDX2 <- which(lambda1 > mean(lambda1[IDX1])/1e+05)"
# is removed to address accuracy issues
#-------------------------------------------------------------------------------
Get_Lambda <- function (K) {
  out.s <- eigen(K, symmetric = TRUE, only.values = TRUE)
  lambda1 <- out.s$values
  #filter out eigenvalues with numbers higher than rank of matrix (they are actually zeroes),
  #since the number of non-zero eigenvalues must be less than or equal to the rank
  lambda1<-lambda1[1:rankMatrix(K)[1]]
  #cat("rank:",rankMatrix(K),"\n")
  #print(lambda1)
  IDX1 <- which(lambda1 >= 0)
  #IDX2 <- which(lambda1 > mean(lambda1[IDX1])/1e+05)
  if (length(IDX1) == 0) { fatalError("No Eigenvalue is bigger than 0!!") }
  lambda <- lambda1[IDX1]
  #print(lambda)
  return( lambda );
}

#-------------------------------------------------------------------------------
# Get p-value
# Note: This function was taken from SKAT package, modified as needed
#-------------------------------------------------------------------------------
Get_PValue.Modif <- function(K, Q, pACC){
  #print(K)
  lambda <- Get_Lambda(K)
  p.val<-0
  p.val.liu<-0
  is_converge <- 0
  
  out <- davies(Q, lambda, acc = pACC)
  p.val <- out$Qq
  p.val.liu <- liu(Q, lambda)
  
  if (out$ifault==0){
    is_converge <- 1} 
  else{
    is_converge <- 0}
  
  if ((length(lambda) == 1) | (is_converge == 0)) {
    p.val <- p.val.liu
  } 
  
  if (p.val > 1 || p.val < 0) {
    is_converge <- 0
  }
  
#-------------------------------------------------------------------------------
# Fatal error
#-------------------------------------------------------------------------------
fatalError <- function(errStr)	{ 
  errStr <- paste('FATAL ERROR:',  errStr, "\n");
  stop( errStr ); 
}

	return(list(p.value = p.val, p.val.liu = p.val.liu, is_converge = is_converge, lambda = lambda))
}

#-------------------------------------------------------------------------------
# Phenotype for kinship matrix
#-------------------------------------------------------------------------------
pheno.geno.kin.Without.NA <- function(Ped.All, kin2, file.name.data){
	pheno.geno0 = cbind(Ped.All, c(1:dim(Ped.All)[1]))
	cont1 = 0
	for (i in 2:dim(Ped.All)[2]){
		pheno.geno1 = pheno.geno0
		pheno.geno0 = pheno.geno0[which(!is.na(pheno.geno0[,(i-cont1)])),(1:dim(pheno.geno0)[2])]
		if (dim(pheno.geno0)[1] < 100){
			pheno.geno0 = pheno.geno1[-(i-cont1)]
			cont1 = cont1 + 1
			}
		}
	pheno.geno = pheno.geno0[,1:(dim(Ped.All)[2] - cont1)]

	kin.FaST = 2 * kin2

	ord.kin2.No.Miss = as.integer(pheno.geno0[,dim(pheno.geno0)[2]])
	kin2.data.miss0 = kin2[ord.kin2.No.Miss,ord.kin2.No.Miss]

	n.col = Ped.All[ord.kin2.No.Miss,1]
	kin.FaST.data.miss = kin.FaST[ord.kin2.No.Miss,ord.kin2.No.Miss]
	kin.FaST.data.miss1 = cbind(n.col,kin.FaST.data.miss)
	kin.FaST.data.miss1 = rbind(c("var",n.col),kin.FaST.data.miss1)

	if( debug )	{ cat('Writing kinship matrix to file: ', file.name.data , '\n' ); }
	write(t(kin.FaST.data.miss1[ ,1:(dim(kin.FaST.data.miss1)[1])]), file = file.name.data, ncolumns = dim(kin.FaST.data.miss1)[1], sep= "\t")

	return( list(kin2.data.miss = kin2.data.miss0, pheno.Geno.No.NA = pheno.geno, ord.kin2.No.Miss = ord.kin2.No.Miss, W = kin.FaST) );
}

#-------------------------------------------------------------------------------
# Invoke FastLmm and get VC
#
# Note: In order to invoke Fast-LMM we must create some temporal files
#       After invoking the program, we have to read and parse the result files
#-------------------------------------------------------------------------------
VC.FaST.LMM <- function(Ped, eigenDir, kinshipFileName ){
	# TMP File names to be used
	phenoFileName <- tmpFile("pheno.txt");
	genoName <- tmpFile("geno_test");
	genoTfamFileName <- paste(genoName,".tfam", sep="");
	genoTpedFileName <- paste(genoName,".tped", sep="");
	genoOutFileName <- paste(genoName,".out.txt", sep="");
	fastlmmOutFileName <- tmpFile("OUTFaST-LMM.txt");

	#---
	# Create TMP files used by the program
	#---

	# Create TPED file
	Y.trait = Ped[,2]
	IID = Ped[,1]
	Geno = Ped[,3]
	FID = c(rep(1,length(IID)))
	Geno.FaST.LMM(Geno, genoTpedFileName)

	# Create 'pheno' file
	pheno.test = cbind(FID, IID, Y.trait)

	if( debug )	{ cat('Writing phenotype to file: ', phenoFileName , '\n' ); }
	write.table(pheno.test, file = phenoFileName, sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)

	# Create TFAM file
	geno.tfam = cbind(pheno.test[,1], pheno.test[,2], c(rep(0,length(Y.trait))), c(rep(0,length(Y.trait))), c(rep(0,length(Y.trait))), Y.trait)

	if( debug )	{ cat('Writing tfam to file: ', genoTfamFileName , '\n' ); }
	write.table(geno.tfam, file = genoTfamFileName, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

	#---
	# Invoke Fast-LMM
	#---

	# Create Fast-LMM command line and execute it
	fastlmmcCmd <- paste( path.FastLmm		# Full path to fastlmm binary
		, "-tfile", genoName				# basename for PLINK's transposed .tfam and .tped files
		, "-sim" , kinshipFileName			# file containing the genetic similarity matrix
		, "-eigenOut", eigenDir				# save the spectral decomposition object to the directoryname
		, "-pheno", phenoFileName			# name of phenotype file
		, "-mpheno 1"						# index for phenotype in -pheno file to process, starting at 1 for the first phenotype column
		);
	if( debug )	{ cat('Execute system command: ', fastlmmcCmd , '\n' ); }
	system(fastlmmcCmd)

	#---
	# Read Fast-LMM results (from TMP files)
	#---

	# Read results from 'genoOutFileName'
	if( debug )	{ cat('Reading VC.FaST.LMM table from: ', genoOutFileName , '\n' ); }
	VC.FaST.LMM = read.table(genoOutFileName, header=TRUE)
	Estim.Sigma.RG = VC.FaST.LMM[1,19]
	Estim.Sigma.e = VC.FaST.LMM[1,20]
	pvalue.FaST.LMM = VC.FaST.LMM[1,9]

	# Read 'S' matrix
	p = dim(Ped)[1]
	sFileName <- paste(eigenDir,"S.bin", sep="/")
	if( debug )	{ cat('Reading read.SVD.bin S-matrix from: ', sFileName , '\n' ); }
	read.SVD.bin = file(sFileName, "rb")
	S = readBin(read.SVD.bin, "numeric", n=p, endian="little")
	close(read.SVD.bin)
	S = diag(sort(S, decreasing = T))

	# Read 'U' matrix
	uFileName <- paste(eigenDir,"U.bin", sep="/")
	if( debug )	{ cat('Reading upper read.SVD.bin U-matrix from: ', uFileName , '\n' ); }
	read.SVD.bin = file(uFileName, "rb")
	U = readBin(read.SVD.bin, "numeric", n=p*p, endian="little")
	close(read.SVD.bin)

	U = matrix(U,p,p, byrow=F)
	U = U[ ,ncol(U):1]

	# Remove TMP files and dirs
	unlink( c(fastlmmOutFileName, genoOutFileName, genoTpedFileName , genoTfamFileName, phenoFileName, kinshipFileName) );
	unlink( eigenDir, recursive = TRUE)
	unlink( kinshipFileName );				# WARNING: A function should NOT have secondary effect (e.g. deleting a file created somewhere else)

	# Create results list
	return( list(Polygenic = Estim.Sigma.RG, Env = Estim.Sigma.e, pvalue.FaST.LMM = pvalue.FaST.LMM, S = S, U = U) );
}

#-------------------------------------------------------------------------------
# Create a name for a temporal file
#-------------------------------------------------------------------------------
tmpFile <- function(name)	{ return( paste(tmpDir, name, sep="/") ); }

#-------------------------------------------------------------------------------
# Main program
#-------------------------------------------------------------------------------

# Example of command line arguments
# You can just uncomment one of these lines for debugging.
#
#---
# Default values
#---

#path.FastLmm <- "fastlmmc"
path.FastLmm<-"./fastlmmc"
debug <- FALSE
tmpDir <- '.'
pACC<-1e-06

#load("./kin2.Rdata"); dataFile = "PedB.dat"
#load("./kin1.Rdata"); dataFile = "Ped_EX_ASKAT_NOMissdata.dat"
load("./kinEst.RData"); dataFile = "PED_leptin.dat"; kin1<-kinEst

Ped  <- read.csv(dataFile, sep="", header=FALSE );
Ped <- Ped



#------------             In thid STEP we Calculate:            ------------#
#------------ the variance components under the null model   ---------------#
#------- and we get the spectral decomposition of the kinship matrix -------#
#------- It is The part of the program that should be done only once -------#
Missing = FALSE # this is means that there is no missing data neither in the trait nor in the genotypes

ptm <- proc.time()
results.STEPs12 = STEPs12(Ped, kin1, Missing, tmpDir="." )
print("Time for steps12")
print(proc.time() - ptm)

Estim.Sigma.RG = results.STEPs12$Estim.Sigma.RG
Estim.Sigma.e = results.STEPs12$Estim.Sigma.e

S = results.STEPs12$S
U = results.STEPs12$U

#--------- This is the main ASKAT function for a window of 5 SNPs--------#
#------------- This is what you should loop for your windows-------------#

#ptm <- proc.time()
#Rprof("a1.out")
#res.ASKAT = ASKAT(Ped, Estim.Sigma.RG, Estim.Sigma.e, S, U)
#Rprof(NULL, interval=1/1000)
#res.ASKAT
#print("Time for ASKAT")
#print(proc.time() - ptm)
#library(profr)
#plot(parse_rprof("a1.out", interval=0.001))
#k1<-parse_rprof("a1.out", interval=0.001)
results <- list()
#for (i in 1:753){
for (i in 1:10){
  res.ASKAT <- ASKAT(Ped[, c(1,2, c((4*i-3):(4*i)))], Estim.Sigma.RG, Estim.Sigma.e, S, U, pACC)#UPD add freq.MAF argument
  results <-append(results, res.ASKAT)
  print(i)
  print(unlist(res.ASKAT))
}




