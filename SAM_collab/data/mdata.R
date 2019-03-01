
setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\data_EDA")


catch = read.table(file='catch.dat',header=FALSE) 

catch<-catch[,-c(1,15,16)]

setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\SAM\\data")
write.table(catch,"cn.txt",sep="\t",row.names=FALSE)


setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\data_EDA")

c_w = read.table(file='comm_wt.txt',header=TRUE) 
c_w<-c_w[,-c(1)]

setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\SAM\\data")
write.table(c_w,"cw.txt",sep="\t",row.names=FALSE, col.names=FALSE)

setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\data_EDA")

mat = read.table(file='mat.txt',header=TRUE) 
mat<-mat[-c(1,2,3,4),]
mat<-mat[-c(57,58,59,60,61),]
mat<-mat[,-c(1)]

setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\SAM\\data")
write.table(mat,"mo.txt",sep="\t",row.names=FALSE, col.names=FALSE)


setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\data_EDA")

rv = read.table(file='RV_OFF.dat',header=TRUE) 
rv[is.na(rv)] <- -1

rv<-rv[-c(10),]
rv<-rv[,-c(1)]

setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\SAM\\data")
write.table(rv,"survey.txt",sep="\t",row.names=FALSE, col.names=FALSE)


setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\data_EDA")

rv = read.table(file='DFO_stock_wt.txt',header=TRUE) 

rv<-rv[,-c(1)]

setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\SAM\\data")
write.table(rv,"sw.txt",sep="\t",row.names=FALSE, col.names=FALSE)
