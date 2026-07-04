install.packages('e1071')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("preprocessCore")


inputFile="merge.normalize.txt"      #表达数据文件
setwd("C:\\Users\\fan'yun'bin\\Desktop\\网络毒理学数据\\机器学习\\免疫浸润分析")      #设置工作目录
source("smrML20.CIBERSORT.R")       #引用包

#免疫细胞浸润分析
outTab=CIBERSORT("ref.txt", inputFile, perm=1000)

#对免疫浸润结果过滤，并且保存免疫细胞浸润结果
outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


