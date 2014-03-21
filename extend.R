#extend Ped_EX_ASKAT.. file
Ped<-rbind(Ped, Ped, Ped, Ped)
Ped[,1]<-c(1:2400)
addData<-data.frame(replicate(20,sample(c(0,1,2), 2400, rep=TRUE)))
PedB<-cbind(Ped,addData)
write.table(file="PedB.dat", x=PedB, sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)

dataFile = "PedB.dat"
Ped  = read.csv(dataFile, sep="", header=FALSE );

#extend kinship data simply by duplicating initial matrix
load("./kin1.Rdata")
k<-cbind(kin1,kin1,kin1,kin1)
kin1<-rbind(k,k,k,k)
save(kin1, file = "kin2.Rdata")
