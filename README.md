# Reading an input cancer data file
data_file=read.csv("C:/Users/91973/Desktop/GSE149650_Read_counts.csv",sep=",",header=T,row.names = 1)

# Create a count per matrix
#Counts per million mapped reads (also known as CPM) is a basic gene expression unit that normalizes only for sequencing depth .

cpmatrix=data_file
for(i in 1:ncol(data_file)){
  cpmatrix[,i]=(data_file[,i]/sum(data_file[,i]))*1000000
}

# To Calculate a log of cpm
# Log pf CPM removes genes that are lowly expressed
logcpm=log2(cpmatrix+1)
saveRDS(logcpm,file="logCPM.rds")
summary(logcpm)

# Calculate a z score 
# Z score indicates how many standard deviations a value is above or below the mean.

library(matrixStats)

z_score = (logcpm - rowMeans(logcpm))/rowSds(as.matrix(logcpm))[row(logcpm)]

# Calculate variance using log 

variance = apply(logcpm, 1, var)

variance = sort(variance,decreasing = T)

top50 = variance[1:50]

pmat = z_score[names(top50),]

# Create a heatmap

library(ComplexHeatmap)

Heatmap(pmat)
![Heatmap](https://user-images.githubusercontent.com/110582335/197956414-265ff60c-1ce3-48b4-92bd-34dead63b887.jpeg)
