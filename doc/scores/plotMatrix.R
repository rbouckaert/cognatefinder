


args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
	fileName = "asjp.mat"
}else{
	fileName = args[1]
}



mat.df = read.table(fileName, sep = "\t", header=F)
symbols = mat.df$V1
mat = as.matrix(mat.df[,-1])
rownames(mat) = symbols
colnames(mat) = symbols




png(gsub("[.]mat", ".png", fileName), width=1000, height=1000, res=100)
heatmap(mat, main=fileName)
dev.off()












