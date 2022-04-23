# DEGman

DEGman is a R package for detecting differentially-expressed genes (DEGs) from two cell groups in single cell RNA-Seq data

Version: 0.1.1

Depends: R(>3.6)

## Install the DEGman package
```
install.packages("devtools")
devtools::install_github("shaoqiangzhang/DEGman")
```

# Usage example

**Step 1**: load the package, read file, and preprocess data

```
library(DEGman)
expressdata=read.csv("data/example.csv",header = T,row.names = 1)  ##read csv file
data=datapreprocess(expressdata) ## preprocess data (filter genes, do log-transformation) 
```

**Instead of step1**, we also recommand you use the function in **"Seurat"** to preprocess data because Seurat is a standard pipline for analysing single-cell data 
```
library(DEGman)
library(Seurat)
expressdata=read.csv("data/example.csv",header = T,row.names = 1)  ##read csv file
expressdata=CreateSeuratObject(counts = expressdata, min.cells = 3, min.genes = 200)
expressdata<- NormalizeData(object =expressdata, normalization.method = "LogNormalize", scale.factor = 100000)
data=expressdata@assays$RNA@data
```

**Step 2**: run DEGman

```
result=DEGman(data,116,64) # 116 and 64 are the cell numbers of two groups, respectively. 
```

**Step 3**: save the result into a csv file

```
write.csv(result,"result.csv");
```
