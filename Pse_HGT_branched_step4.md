<!-- TOC -->

- [1.数据信息](#1数据信息)
  - [1.1选择基因表达量数据](#11选择基因表达量数据)
  - [1.2名词解释](#12名词解释)
- [2.数据下载和获取（手动下载15个数据集)](#2数据下载和获取手动下载15个数据集)
  - [2.1处理原始数据，比如GSE8408_RAW.tar，并用RMA标准化或mas5标准化（举例)](#21处理原始数据比如gse8408_rawtar并用rma标准化或mas5标准化举例)
  - [2.3批量化处理cel文件](#23批量化处理cel文件)
- [3.利用WGCNA包进行分析](#3利用wgcna包进行分析)
  - [3.1载入表达量数据](#31载入表达量数据)
  - [3.2样本层级聚类和检测离群值](#32样本层级聚类和检测离群值)
  - [3.3通过函数pickSoftThreshold执行网络拓扑分析，并帮助用户选择合适的软阈值。](#33通过函数picksoftthreshold执行网络拓扑分析并帮助用户选择合适的软阈值)
  - [3.4power值分布图，并根据分布图挑选合适的软阈值](#34power值分布图并根据分布图挑选合适的软阈值)
  - [3.5一步法网络构建和模块识别](#35一步法网络构建和模块识别)
  - [3.6层级聚类树展示各个模块](#36层级聚类树展示各个模块)
  - [3.7绘制模块之间相关性状图](#37绘制模块之间相关性状图)
  - [3.8可视化基因网络](#38可视化基因网络)
  - [3.9导出网络用于Cytoscape](#39导出网络用于cytoscape)
- [4.直接使用R脚本运行命令（批量化)](#4直接使用r脚本运行命令批量化)
  - [4.1拆分条件后两两组合（不能用！！！！！！）](#41拆分条件后两两组合不能用)
  - [4.2拆分条件后三三组合](#42拆分条件后三三组合)
- [5.统计基因频率](#5统计基因频率)
  - [5.1拆分条件后两两组合(条件拆分后两两组合，共有1600个组合，共1524个有结果文件) （不能用！！！！！！）](#51拆分条件后两两组合条件拆分后两两组合共有1600个组合共1524个有结果文件-不能用)
  - [5.2拆分条件后三三组合(条件拆分后三三组合，共有64000个组合，共56663个有结果文件)](#52拆分条件后三三组合条件拆分后三三组合共有64000个组合共56663个有结果文件)
- [6.KEGG和GO分析(不能用!!!!!!)](#6kegg和go分析不能用)
  - [6.1使用david对gene_id进行转换(方法一) (不能用!!!!!!)](#61使用david对gene_id进行转换方法一-不能用)
  - [6.2下载注释数据库 (不能用!!!!!!)](#62下载注释数据库-不能用)
    - [6.2.1GO注释库 (不能用!!!!!!)](#621go注释库-不能用)
    - [6.2.2KEGG注释库 (不能用!!!!!!)](#622kegg注释库-不能用)
  - [6.3构建非模式物种的OrgDb  (不能用!!!!!!)](#63构建非模式物种的orgdb--不能用)
  - [6.4进行GO和KEGG富集分析  (不能用!!!!!!)](#64进行go和kegg富集分析--不能用)
  - [6.5GO富集分析结果可视化  (不能用!!!!!!)](#65go富集分析结果可视化--不能用)
- [7.使用mclust进行聚类分析](#7使用mclust进行聚类分析)
  - [7.1mclust使用举例](#71mclust使用举例)
  - [7.2三三组合的结果文件划分bins为100然后进行mclust分析](#72三三组合的结果文件划分bins为100然后进行mclust分析)
  - [7.3 mclust计算均值和平方差](#73-mclust计算均值和平方差)
  - [7.4ggplot正态分布直方密度图](#74ggplot正态分布直方密度图)
- [8.eggnog在线注释](#8eggnog在线注释)
  - [8.1绘制三三组合的两个sigma和三个sigma的braB和braZ的韦恩图](#81绘制三三组合的两个sigma和三个sigma的brab和braz的韦恩图)
  - [8.2下载注释文件](#82下载注释文件)
  - [8.3改变注释文件的格式](#83改变注释文件的格式)
  - [8.4GO富集分析](#84go富集分析)
  - [8.5KEGG富集分析](#85kegg富集分析)
- [9.绘制三三组合的braZ和braB的基因表达谱](#9绘制三三组合的braz和brab的基因表达谱)
  - [9.1绘制braB的基因表达谱](#91绘制brab的基因表达谱)
  - [9.2绘制braZ的基因表达谱](#92绘制braz的基因表达谱)
- [10.绘制40个条件的braZ和braB的基因表达谱](#10绘制40个条件的braz和brab的基因表达谱)
  - [10.1绘制braB的基因表达谱](#101绘制brab的基因表达谱)
  - [10.2绘制braZ的基因表达谱](#102绘制braz的基因表达谱)
  - [10.3以14,13,13分为低中高组合](#103以141313分为低中高组合)
  - [10.4plotr绘制40个条件的braZ和braB的密度分布图](#104plotr绘制40个条件的braz和brab的密度分布图)
- [11.绘制40个条件中braB和braZ的基因表达量热图](#11绘制40个条件中brab和braz的基因表达量热图)
  - [11.1绘制40个条件中braB和braZ的平均基因表达量热图，不聚类](#111绘制40个条件中brab和braz的平均基因表达量热图不聚类)
  - [11.2绘制40个条件中braB和braZ和相关基因的平均基因表达量热图，对条件进行聚类，基因不聚类,对基因标准化](#112绘制40个条件中brab和braz和相关基因的平均基因表达量热图对条件进行聚类基因不聚类对基因标准化)
  - [11.3使两个热图的条件顺序一致](#113使两个热图的条件顺序一致)
  - [11.3将热图的样本条件改成样本号码](#113将热图的样本条件改成样本号码)

<!-- /TOC -->

# 1.数据信息
## 1.1选择基因表达量数据
1.在PAO1数据库中输入PA号，根据PA号显示对应的braB/Z序列注释及下方的"profiles from geo expression database at NCBI)
2.选择Expression---选择GEO Dataset Details---选择GDS如GDS2870（非GPL,是因为GPL是原始数据，GDS是整理后的数据,质量较好)(Value type:transformed count)---选择GSE如GSE8408---随机点击biosample(查看Data table header descriptions,其处理方式为GC-RMA calculated signal intensity)---返回GSE界面，下载原始数据GSE8408_RAW.tar，这是因为处理后的23个数据集处理方式各不相同，比如有RMA，有log2 qcRNA，有mas5,有GC-RMA算法等等，所以需要重新标准化
3.查看23个数据集的测序平台Platforms  Affymetrix Pseudomonas aeruginosa Array,其expression的处理方式为RMA MAS5等 其value分别为count和transformed count
4.因此需要下载原始数据(.soft或者series.maritx都是经过处理的)_RAW.tar文件，每个tar压缩包里是单个样本的原始数据CEL文件，其芯片平台都是Affymetrix
5.根据统计,共有16个数据集拥有原始文件

| GDS号 | GSE号|value type| 测序平台 | 样本处理方式 | 有无原始文件|max expression|
| :------ | :------ | :------ | :------ | :------ | :------ |  :------ | 
|GDS2870 | GSE8408 | transformed count |Affymetrix Pseudomonas aeruginosa Array |GC-RMA calculated signal intensity|有| 3|  
|GDS4457 | GSE27674| transformed count |Affymetrix Pseudomonas aeruginosa Array | 	log2 qcRNA| 有| 6.4| 
|GDS1469| GSE3090| count| Affymetrix Pseudomonas aeruginosa Array| GCOS-calculated signal intensity| 无| 120| 
|GDS3572| GSE13252| count|Affymetrix Pseudomonas aeruginosa Array |GeneSpring normalized data| 有| 120| 
|GDS3562 | GSE12738 | count|Affymetrix Pseudomonas aeruginosa Array |MAS5 signal|有|600| 
|GDS3251| GSE10030 | count|Affymetrix Pseudomonas aeruginosa Array |value|有|300| 
|GDS2869 | GSE7704| count|Affymetrix Pseudomonas aeruginosa Array |Signal Intensity|有|1000| 
|GDS2893| GSE6741 | count|Affymetrix Pseudomonas aeruginosa Array |Signal Intensity|有|600| 
|GDS2317|GSE5443|count|Affymetrix Pseudomonas aeruginosa Array |Signal value|无|400|
|GDS4244|GSE34141|transformed count|Affymetrix Pseudomonas aeruginosa Array |RMA|有|3| 
|GDS4250|GSE18594|transformed count|Affymetrix Pseudomonas aeruginosa Array |RMA signal|有|9.3| 
|GDS3174|GSE9621| count|Affymetrix Pseudomonas aeruginosa Array |dChip signal|有|750| 
|GDS2377|GSE4152 |count|Affymetrix Pseudomonas aeruginosa Array |MAS5-calculated Signal intensity|无|3000|
|GDS2764|GSE3836| count|Affymetrix Pseudomonas aeruginosa Array |MAS 5.0-processed and scaled signal intensities|无|250|
|GDS1910|GSE2885|count|	Affymetrix Pseudomonas aeruginosa Array |MAS5-calculated Signal intensity|无|300|
|GDS2502|GSE4614| count|Affymetrix Pseudomonas aeruginosa Array |MAS5-calculated Signal intensity|无|6000|
|GDS4254|GSE35248| count|Affymetrix Pseudomonas aeruginosa Array |Signal|有|780| 
|GDS2111|GSE2430|count|Affymetrix Pseudomonas aeruginosa Array |MAS5.0|有|500|
|GDS4479|GSE23007|transformed count|Affymetrix Pseudomonas aeruginosa Array |log2 normalized signal intensity|有|8| 
|GDS4480|GSE25945|transformed count|Affymetrix Pseudomonas aeruginosa Array |log2 normalized signal intensity|有|9| 
|GDS4249|GSE21966|transformed count|Affymetrix Pseudomonas aeruginosa Array |log2 RMA signal intensity|有|9| 
|GDS4959|GSE48982(8个样本)|transformed count|Affymetrix Pseudomonas aeruginosa Array |log2 RMA|有|6.5| 
|GDS4958|GSE48982(同上4个样本)(该数据集不需要)|transformed count|Affymetrix Pseudomonas aeruginosa Array |log2 RMA|有|6.5|

## 1.2名词解释
* RMA标准化和mas5标准化
1.RMA 是一种用于从 Affymetrix 数据创建表达式矩阵的算法。原始强度值经过背景校正、log2 转换，然后分位数归一化。接下来，将线性模型拟合到标准化数据以获得每个阵列上每个探针组的表达测量
2.芯片表达定量过程分为background correct、quantile normalize、summarize 三个主要步骤，可以通过R包 affy 或 oligo 中的一体化算法 RMA（Robust Multiarray Average）直接完成上述四步的数据处理
3.Affy芯片数据的预处理一般有三个步骤：背景处理（background adjustment），归一化处理（normalization，或称为“标准化处理”），汇总（summarization）。最后一步获取表达水平数据。需要说明的是，每个步骤都有很多不同的处理方法（算法），选择不同的处理方法对最终结果有非常大的影响
4.affy包的rma算法返回值经过log2变换，mas5算法在最后一步才进行normalize，而且输出的数据并没有经过log2变换，而大多数的后续分析都要求数据进行log2变换
5.Affymetrix表达谱芯片数据读取的方法分3种：
1)使用affy包读取。（HGU95/HGU133芯片）
2)使用oligo包读取。（Whole Transcriptome 芯片/ NimbleGen 芯片/ SNP芯片等）
3)使用simpleaffy包读取。（HGU95/HGU133芯片）

* 判断GEO芯片数据表达矩阵是否需要log2转换
1.肉眼识别：如果表达量的数值在50以内，通常是经过log2转化后的。如果数字在几百几千，则是未经转化的。因为2的几十次方已经非常巨大，如果2的几百次方，则不符合实际情况。
2.GSE数据下载界面中的SOFT文件和Series Matrix File(s)文件中均有描述该系列的数据是如何进行标准化处理的，常见的标准化处理方法有3种：RMA算法、GC-RMA算法、MAS5算法，其中前两中算法的返回值已经经过log2转换，可直接进行差异表达分析，第三种算法返回值未经过log2转换，需要自行进行log2转换


* 去除批次效应(batch effects)
1.GEO数据挖掘分析中经常会遇到会遇到样品数据不足，需要联合分析多个芯片数据进行分析，那么能将这些数据直接混合分析吗？如果贸然混合，会有什么问题？这个问题就是batch effect。
2.不同平台的数据，同一平台的不同时期的数据，同一个样品不同试剂的数据，以及同一个样品不同时间的数据等等都会产生一种batch effect 。
3.如何检测是否存在这种效应呢.最简单的就是记录实验中时间这个变量，然后对差异表达的基因进行聚类，看是否都和时间相关，如果相关就证明存在batch effect。
4.标准化只能减弱，不能消除批次效应.矫正批次效应：sva中的combat包和limma的removeBatchEffect 函数

* WGCNA分析对样本有什么要求？
推荐以下的样品数：
1）不含生物学重复的独立样本组：样本数≥8
2）包含生物学重复的样本组：样本数≥15
主要考虑的因素：
1）样本必须包含丰富的变化信息，才能区分为多个有
意义的生物学模块（需要多个独立处理组）。
2）必须保证有多个样本，才能保证相关系数计算的准
确性。

* WGCNA分析中需要找出模块中也就共表达网络中的hub基因：
这里可以根据：KME (eigengene connectivity)值来筛选hub基因：Hub genes are those that show most connections in the network as indicated by their high KME (eigengene connectivity) value
kme 值也叫Module Membership,顾名思义就是模块关系，是基因与模块之间的关系。如何计算基因与模块之间的关系，就是计算基因的表达量与模块特征值之间的相关性。其实Kme值 只是衡量基因与模块之间的相关性。绝对值越大，越接近1代表这个基因和这个模块越相关。所以会被认为是这个模块内的。但是这个参数也只能衡量基因与模块之间的关系。和 hubgene 的概念并没有什么直接的交集

* GO和KEGG
1. 基因本体论是对基因在不同维度和不同层次上的描述
对基因的描述一般从三个层面进行：
Cellular component，CC 细胞成分
Biological process， BP 生物学过程
Molecular function，MF 分子功能 
模式生物 ---> 有标准的注释数据库；
非模式生物 ---> 自己搜注释数据库(怎们搜后面具体介绍)，搜不到就用blast的办法解决
2. KEGG enrichment analysis？
把生物体中所有的pathway都要进行富集分析
DO enrichment analysis？
看目标基因是否在某个疾病或某一类疾病当中富集


# 2.数据下载和获取（手动下载15个数据集)
## 2.1处理原始数据，比如GSE8408_RAW.tar，并用RMA标准化或mas5标准化（举例)
```r
BiocManager::install("affy")
BiocManager::install("paeg1a.db")
BiocManager::install("paeg1acdf")
library(affy)
library(paeg1a.db)
library(paeg1acdf)
# 读取cel文件
setwd("D:/WGCNA/rawdata")
rawData <- ReadAffy(celfile.path = './GSE48982_RAW')
# 使用RMA方法进行预处理：依次进行背景处理，log2转化，分位数标化和探针表达量计算
eset.rma <- rma(rawData)
# 提取探针水平表达矩阵
ex.mat.rma <- exprs(eset.rma)
# Write RMA-normalized, mapped data to file
write.table(ex.mat.rma, file = "GSE48982.rma.txt", quote = FALSE, sep = "\t")

# 使用MAS5方法进行预处理
eset.mas5 <- mas5(rawData)
# 提取探针水平表达矩阵
ex.mat.mas5 <- exprs(eset.mas5)
#取对数
ex.mas5<-log(ex.mat.mas5,2)
# Write RMA-normalized, mapped data to file
write.table(ex.mas5, file = "GSE48982.mas5.txt", quote = FALSE, sep = "\t")
```

## 2.3批量化处理cel文件
```bash
cd /mnt/d/WGCNA/rawdata
for filename in *.tar
do
base=$(basename $filename .tar)
echo $base
Rscript rawdata_normalize.R $base
done
```

# 3.利用WGCNA包进行分析
## 3.1载入表达量数据
```R
#安装R包
install.packages("BiocManager")
BiocManager::install("WGCNA")
BiocManager::install("reshape2")
BiocManager::install("DESeq2")
#stringr包是专门用于字符处理的R包
BiocManager::install("stringr")
library(WGCNA)
library(reshape2)
library(stringr)
setwd("D:/braB")
exprMat <- "GDS3562.expression.tsv"
dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=T,
                     quote="", comment="", check.names=F)
#查看数据的前五列，前五行
head(dataExpr)[,1:5]
#查看数据名称
names(dataExpr) 
#查看数据维度 5891   10
dim(dataExpr)

#得到的datExpr一定是基因名在列，样本名在行
#去除第一列即标识符的这一列，并且将矩阵数据转置，并且转换为数据框
dataExpr <- as.data.frame(t(dataExpr))
#查看数据的前五列，前五行
head(dataExpr)[,1:5]
#查看数据维度 9  5891
dim(dataExpr)

``` 
## 3.2样本层级聚类和检测离群值
```r
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
```
## 3.3通过函数pickSoftThreshold执行网络拓扑分析，并帮助用户选择合适的软阈值。
* blockwiseModules()函数构建网络。官方推荐 "signed" 或 "signed hybrid" 。需要注意的是，参数networkType 有"unsigned"、"signed"等选择，"signed"与"unsigned"方法的不同会导致阈值和网络的不同。不过即使该步与之后TOMsimilarityFromExpr()函数均选择"signed"，最终输出用于Cytoscape可视化的网络时仍显示"undirected"。
* 相关性计算 官方推荐 biweight mid-correlation & bicor corType: pearson or bicor 这里选择的是bicor双重相关
* 软阈值的筛选原则是使构建的网络更符合无标度网络特征。
* 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
```r
type = "signed"
corType = "bicor"
#设置一个power区间，一般为c(c(1:10), seq(from = 12, to=20, by=2))。 利用pickSoftThreshold()函数对datExpr在该power区间内筛选出合适的阈值
powers = c(c(1:10), seq(from = 12, to=30, by=2))
#计算有尺度分布拓扑矩阵
# 获得各个阈值下的 R方 和平均连接度
sft = pickSoftThreshold(dataExpr, powerVector=powers,
                        networkType=type, verbose=5)
#一般选择函数估测的阈值(powerEstimate)作为最优值
#挑选软阈值时存在一个很大的问题！！！！比如样本GSE10030.GSE18594.soft，由于其R2<-0.8，会将其误判为R2>+0.8,从而并不能将power定为NA
#此外当R2为负值时，连通度非常大，根本不约等于100左右
data1<-data.frame(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2])
if(min(abs(data1[,2]-0.9)>0.1)){
power=NA}
if(min(abs(data1[,2]-0.9))<=0.1){power=data1[which(abs(data1[,2]-0.9)==min(abs(data1[,2]-0.9))),1]}
#列为基因，列为样本
nGenes = ncol(dataExpr)
#经验power (无满足条件的power时选用)
nSamples = nrow(dataExpr)
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
          ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
          ifelse(type == "unsigned", 6, 12))
          )
          )
}
```

## 3.4power值分布图，并根据分布图挑选合适的软阈值
```r
#sft$fitIndices[,1]为Power这一列，sign函数如果为正数，返回1，如果为负数，返回-1，如果为0，返回0
#如果想知道WGCNA建议的阈值，可以输入sft，列表的第一项值就是估计的最优值
#绘制power图
par(mfrow = c(1,2))
cex1=0.85
pdf(file=paste0(exprMat, ".power1.pdf"),width =12,height = 9)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
dev.off()

pdf(file=paste0(exprMat, ".power2.pdf"),width =12,height = 9)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
dev.off()
```

## 3.5一步法网络构建和模块识别
* 该步骤把输入的表达矩阵的几千个基因组归类成了几十个模块。其中module grey包含的是未被分配的基因，后续研究中可以不关注该module。大体思路：计算基因间的邻接性，根据邻接性计算基因间的相似性，然后推出基因间的相异性系数，并根据此得到基因间的系统聚类树。然后按照混合动态剪切树的标准。
* power :上一步计算的软阈值
* maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)
* corType: pearson or bicor
* numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
* saveTOMs：最耗费时间的计算，存储起来，供后续使用
* mergeCutHeight: 合并模块的阈值，越大模块越少
* deepSplit 参数调整划分模块的敏感度，值越大，越敏感，得到的模块就越多，默认是2
* minModuleSize 参数设置最小模块的基因数，值越小，小的模块就会被保留下来
```R
#一步法网络构建
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.3,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, corType = corType,
                       loadTOMs = TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)
#根据模块中基因数目的多少，降序排列，依次编号为 1-最大模块数。0 (grey)表示未分入任何模块的基因。
table(net$colors)
```
## 3.6层级聚类树展示各个模块
```r
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
```
## 3.7绘制模块之间相关性状图
```R
# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs
### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)

```

## 3.8可视化基因网络
```r
# 如果采用分步计算，或设置的blocksize>=总基因数，直接load计算好的TOM结果
# 否则需要再计算一遍，比较耗费时间
# TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
plotTOM = dissTOM^7
diag(plotTOM) = NA

# 这一部分特别耗时，行列同时做层级聚类
TOMplot(plotTOM, net$dendrograms, moduleColors, 
        main = "Network heatmap plot, all genes")
```

## 3.9导出网络用于Cytoscape
```r
probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)
# threshold 默认为0.5, 可以根据自己的需要调整，也可以都导出后在
# cytoscape中再调整
cyt = exportNetworkToCytoscape(TOM,
             edgeFile = paste(exprMat, ".edges.txt", sep=""),
             nodeFile = paste(exprMat, ".nodes.txt", sep=""),
             weighted = TRUE, threshold = 0,
             nodeNames = probes, nodeAttr = moduleColors)
```

# 4.直接使用R脚本运行命令（批量化)
## 4.1拆分条件后两两组合（不能用！！！！！！）
```bash
cd /mnt/d/WGCNA/rma_WGCNA/
#最后共生成16个RMA和16个MAS5文件
#将16个数据集拆分为40个条件
#对40个条件进行两两组合

#先对文件行列转置，方便后面提取指定行
for filename in *.txt
do
base=$(basename $filename .txt)
echo $base
perl ranks.pl $base.txt  >condition_combine/data/$base.tsv
done

#根据condition文件提取指定行
for filename in  sample_condition2/*.tsv
do
base=$(basename $filename .tsv)
gene=$(echo "$base" | cut -d "_" -f 1)
cat condition_combine/data/$gene.rma.tsv  | grep -f sample_condition2/$base.tsv >condition_combine/$base.condition.tsv
done

#再进行转置，方便两两组合
for filename in condition_combine/*.tsv
do
base=$(basename $filename .tsv)
echo $base
perl ranks.pl condition_combine/$base.tsv  >condition_combine/$base.txt
done

#两两组合
for i in {1..40}
do
name1=$(ls *.condition.txt | cut -d . -f 1 | head -n $i | tail -n 1)
file1=$i.$name1
for s in {1..40}
do
name2=$(ls *.condition.txt | cut -d . -f 1 | head -n $s | tail -n 1)
file2=$s.$name2
echo $file1.$file2 >>combine.tsv
done
done
#去除重复组合的
#cat combine.tsv | tr "." "\t" |tsv-filter --ff-lt 1:3 | cut -f 2,4 >combine_uniq.tsv
#去除重复组合的
cat combine.tsv | tr "." "\t" | cut -f 2,4 >combine_uniq.tsv
#拼接文件
for i in {1..1600}
do
name1=$(cat combine_uniq.tsv| head -n $i | tail -n 1 | cut -f 1)
name2=$(cat combine_uniq.tsv| head -n $i | tail -n 1 | cut -f 2)
echo $name1.$name2
cat $name1.condition.txt | tsv-join --H -f $name2.condition.txt --key-fields ID --append-fields 'GSM*' >combine/$name1.$name2.soft
done

#运行命令
for filename in *.soft
do
base=$(basename $filename .soft)
echo $base
Rscript WGCNA.R $base.soft
done

```
## 4.2拆分条件后三三组合
```bash
cd /mnt/d/WGCNA/rma_WGCNA/
#最后共生成16个RMA和16个MAS5文件
#手动将16个数据集拆分为40个条件，放在在目录sample_conditon2下
#对40个条件进行三三组合

#先对文件行列转置，方便后面提取指定行
for filename in *.txt
do
base=$(basename $filename .txt)
echo $base
perl ranks.pl $base.txt  >condition_three_combine/data/$base.tsv
done

#根据condition文件提取指定行
for filename in  sample_condition2/*.tsv
do
base=$(basename $filename .tsv)
gene=$(echo "$base" | cut -d "_" -f 1)
cat condition_three_combine/data/$gene.rma.tsv  | grep -f sample_condition2/$base.tsv >condition_three_combine/$base.condition.tsv
done

#再进行转置，方便三三组合
for filename in condition_three_combine/*.tsv
do
base=$(basename $filename .tsv)
echo $base
perl ranks.pl condition_three_combine/$base.tsv  >condition_three_combine/$base.txt
done

#三三组合
for i in {1..40}
do
name1=$(ls *.condition.txt | cut -d . -f 1 | head -n $i | tail -n 1)
file1=$i.$name1
for s in {1..40}
do
name2=$(ls *.condition.txt | cut -d . -f 1 | head -n $s | tail -n 1)
file2=$s.$name2
for s in {1..40}
do
name3=$(ls *.condition.txt | cut -d . -f 1 | head -n $s | tail -n 1)
file3=$s.$name3
echo $file1.$file2.$file3 >>combine.tsv
done
done
done

#.替换为换行符
cat combine.tsv | tr "." "\t" | cut -f 2,4,6 >combine_uniq.tsv

#拼接文件
for i in {1..64000}
do
name1=$(cat combine_uniq.tsv | head -n $i | tail -n 1 | cut -f 1)
name2=$(cat combine_uniq.tsv| head -n $i | tail -n 1 | cut -f 2)
name3=$(cat combine_uniq.tsv| head -n $i | tail -n 1 | cut -f 3)
echo $name1.$name2.$name3
cat $name1.condition.txt | tsv-join --H -f $name2.condition.txt --key-fields ID --append-fields 'GSM*' |  tsv-join --H -f $name3.condition.txt --key-fields ID --append-fields 'GSM*'> combine/$name1.$name2.$name3.soft
done

#运行命令
for filename in *.soft
do
base=$(basename $filename .soft)
echo $base
Rscript WGCNA.R $base.soft
done

#批量提交命令 
#把上述批量化命令写入bash脚本，然后用超算运行
#超算安装R包 R CMD INSTALL WGCNA_1.69-81.tar.gz
```

# 5.统计基因频率
## 5.1拆分条件后两两组合(条件拆分后两两组合，共有1600个组合，共1524个有结果文件) （不能用！！！！！！）
```bash
cd /mnt/d/WGCNA/rma_WGCNA/condition_two_combine/result
for filename in whole/*.soft.nodes.txt
do
base=$(basename $filename .soft.nodes.txt)
echo $base
braB=$(cat whole/$base.soft.nodes.txt| grep "PA1590"| cut -f 3)
braZ=$(cat whole/$base.soft.nodes.txt | grep "PA1971"| cut -f 3)
cat whole/$base.soft.nodes.txt | grep "$braB" >> select/$base.select.txt
cat whole/$base.soft.nodes.txt | grep "$braZ" >> select/$base.select.txt
done

#查看颜色种类
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat select/$base.select.txt  | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR" 
done | tsv-summarize -g 2 --count
2       791
6       4
1       643
3       78
4       7
5       1


#统计颜色种类为1的(共643种)
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat select/$base.select.txt  | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR"  | tsv-filter --eq 2:1 
done | cut -f 1 >color_one.tsv

for PROJECT in $(cat color_one.tsv);do
    cat select/$PROJECT.select.txt | sort | uniq > tem &&
    mv tem select/$PROJECT.select.txt
    done

#查看颜色种类为2的(共791种)
for filename in *.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat $base.select.txt | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR"  | tsv-filter --eq 2:2
done | wc -l

#查看颜色为3（共78种）
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat select/$base.select.txt | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR"  | tsv-filter --eq 2:3
done | cut -f 1 >color_three.tsv

#提取含有braB/braZ的颜色
for PROJECT in $(cat color_three.tsv);do
braZ=$( cat select/$PROJECT.select.txt | grep  "braZ" | head -n 1 | cut -f 3 )
braB=$( cat select/$PROJECT.select.txt | grep  "braB" | head -n 1 | cut -f 3 )
echo -e "$braB\n$braZ" | perl -alne '$name=(split/\s/,$_)[0];$name=~s/\r?\n//g;$name=~s/\s+//g;print"\/"."$name"."\/";' >other/color_three_$PROJECT.tsv
done  

#过滤
for f in $(cat color_three.tsv)
do
echo  "$f"
cat select/$f.select.txt |perl -alne '($gene,$num,$name)=(split/\t/,$_)[0,1,2];$name=~s/\r?\n//g;$name=~s/\s+//g;print"$gene"."\t"."$num"."\t"."\/"."$name"."\/";' | grep -f other/color_three_$f.tsv | sed 's/\///g' > tem &&
mv tem select/$f.select.txt
done


#查看颜色为4（共7种）
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat select/$base.select.txt | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR"  | tsv-filter --eq 2:4
done | cut -f 1 >color_four.tsv

#提取含有braB/braZ的颜色
for PROJECT in $(cat color_four.tsv);do
braZ=$( cat select/$PROJECT.select.txt | grep  "braZ" | head -n 1 | cut -f 3 )
braB=$( cat select/$PROJECT.select.txt | grep  "braB" | head -n 1 | cut -f 3 )
echo -e "$braB\n$braZ" | perl -alne '$name=(split/\s/,$_)[0];$name=~s/\r?\n//g;$name=~s/\s//g;print"\/"."$name"."\/";' >other/color_four_$PROJECT.tsv
done  

#过滤
for f in $(cat color_four.tsv)
do
echo -e  "$f"
cat select/$f.select.txt |perl -alne '($gene,$num,$name)=(split/\t/,$_)[0,1,2];$name=~s/\r?\n//g;$name=~s/\s+//g;print"$gene"."\t"."$num"."\t"."\/"."$name"."\/";' | grep -f other/color_four_$f.tsv | sed 's/\///g'  > tem &&
mv tem select/$f.select.txt
done


#查看颜色为5（共1种）
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat select/$base.select.txt | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR"  | tsv-filter --eq 2:5
done | cut -f 1 >color_five.tsv

#提取含有braB/braZ的颜色
for PROJECT in $(cat color_five.tsv);do
braZ=$( cat select/$PROJECT.select.txt | grep  "braZ" | head -n 1 | cut -f 3 )
braB=$( cat select/$PROJECT.select.txt | grep  "braB" | head -n 1 | cut -f 3 )
echo -e "$braB\n$braZ" | perl -alne '$name=(split/\s/,$_)[0];$name=~s/\r?\n//g;$name=~s/\s//g;print"\/"."$name"."\/";'   >other/color_five_$PROJECT.tsv
done  

#过滤
for f in $(cat color_five.tsv)
do
echo -e  "$f"
cat select/$f.select.txt |perl -alne '($gene,$num,$name)=(split/\t/,$_)[0,1,2];$name=~s/\r?\n//g;$name=~s/\s+//g;print"$gene"."\t"."$num"."\t"."\/"."$name"."\/";' | grep -f other/color_five_$f.tsv | sed 's/\///g' > tem &&
mv tem select/$f.select.txt
done


#查看颜色为6（共4种）
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat select/$base.select.txt | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR"  | tsv-filter --eq 2:6
done | cut -f 1 >color_six.tsv

#提取含有braB/braZ的颜色
for PROJECT in $(cat color_six.tsv);do
braZ=$( cat select/$PROJECT.select.txt | grep  "braZ" | head -n 1 | cut -f 3 )
braB=$( cat select/$PROJECT.select.txt | grep  "braB" | head -n 1 | cut -f 3 )
echo -e "$braB\n$braZ" | perl -alne '$name=(split/\s/,$_)[0];$name=~s/\r?\n//g;$name=~s/\s//g;print"\/"."$name"."\/";'   >other/color_six_$PROJECT.tsv
done  

#过滤
for f in $(cat color_six.tsv)
do
echo -e  "$f"
cat select/$f.select.txt |perl -alne '($gene,$num,$name)=(split/\t/,$_)[0,1,2];$name=~s/\r?\n//g;$name=~s/\s+//g;print"$gene"."\t"."$num"."\t"."\/"."$name"."\/";'  | grep -f other/color_six_$f.tsv | sed 's/\///g'  > tem &&
mv tem select/$f.select.txt
done

#重新统计结果
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat select/$base.select.txt  | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR" 
done | tsv-summarize -g 2 --count
2       881
1       643

#合并结果
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
echo $base
braB=$(cat select/$base.select.txt| grep "PA1590"| cut -f 3)
braZ=$(cat select/$base.select.txt | grep "PA1971"| cut -f 3)
cat select/$base.select.txt| grep "$braB" >> PA1590.select.txt
cat select/$base.select.txt | grep "$braZ" >> PA1971.select.txt
done

#统计不同基因出现的次数
cat PA1590.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | grep "PA1590" # PA1590_braB_at  1386
cat PA1590.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | tsv-filter --ge 2:650 | wc -l  #234
cat PA1971.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | grep "PA1971" # PA1590_braB_at  1386
cat PA1971.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | tsv-filter --ge 2:725 | wc -l  #242
```

## 5.2拆分条件后三三组合(条件拆分后三三组合，共有64000个组合，共56663个有结果文件)
```bash
cd /mnt/d/WGCNA/rma_WGCNA/condition_three_combine/result
for filename in whole/*.soft.nodes.txt
do
base=$(basename $filename .soft.nodes.txt)
echo $base
braB=$(cat whole/$base.soft.nodes.txt| grep "PA1590"| cut -f 3)
braZ=$(cat whole/$base.soft.nodes.txt | grep "PA1971"| cut -f 3)
cat whole/$base.soft.nodes.txt | grep "$braB" >> select/$base.select.txt
cat whole/$base.soft.nodes.txt | grep "$braZ" >> select/$base.select.txt
done

#查看颜色种类
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat select/$base.select.txt  | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR" 
done | tsv-summarize -g 2 --count >color_class.tsv
2	31965
4	332
3	1769
1	22498
5	78
6	18
8	3

#统计颜色种类为1的(共22498种)
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat select/$base.select.txt  | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR"  | tsv-filter --eq 2:1 
done | cut -f 1 >color_one.tsv

for PROJECT in $(cat color_one.tsv);do
    cat select/$PROJECT.select.txt | sort | uniq > tem &&
    mv tem select/$PROJECT.select.txt
    done

#查看颜色种类为2的(共31965种)
for filename in *.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat $base.select.txt | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR"  | tsv-filter --eq 2:2
done | wc -l

#查看颜色为3（共1769种）
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat select/$base.select.txt | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR"  | tsv-filter --eq 2:3
done | cut -f 1 >color_three.tsv

#提取含有braB/braZ的颜色
for PROJECT in $(cat color_three.tsv);do
braZ=$( cat select/$PROJECT.select.txt | grep  "braZ" | head -n 1 | cut -f 3 )
braB=$( cat select/$PROJECT.select.txt | grep  "braB" | head -n 1 | cut -f 3 )
echo -e "$braB\n$braZ" | perl -alne '$name=(split/\s/,$_)[0];$name=~s/\r?\n//g;$name=~s/\s+//g;print"\/"."$name"."\/";' >other/color_three_$PROJECT.tsv
done  

#过滤
for f in $(cat color_three.tsv)
do
echo  "$f"
cat select/$f.select.txt |perl -alne '($gene,$num,$name)=(split/\t/,$_)[0,1,2];$name=~s/\r?\n//g;$name=~s/\s+//g;print"$gene"."\t"."$num"."\t"."\/"."$name"."\/";' | grep -f other/color_three_$f.tsv | sed 's/\///g' > tem &&
mv tem select/$f.select.txt
done

#查看颜色为4（共332种）
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat select/$base.select.txt | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR"  | tsv-filter --eq 2:4
done | cut -f 1 >color_four.tsv

#提取含有braB/braZ的颜色
for PROJECT in $(cat color_four.tsv);do
braZ=$( cat select/$PROJECT.select.txt | grep  "braZ" | head -n 1 | cut -f 3 )
braB=$( cat select/$PROJECT.select.txt | grep  "braB" | head -n 1 | cut -f 3 )
echo -e "$braB\n$braZ" | perl -alne '$name=(split/\s/,$_)[0];$name=~s/\r?\n//g;$name=~s/\s//g;print"\/"."$name"."\/";' >other/color_four_$PROJECT.tsv
done  

#过滤
for f in $(cat color_four.tsv)
do
echo -e  "$f"
cat select/$f.select.txt |perl -alne '($gene,$num,$name)=(split/\t/,$_)[0,1,2];$name=~s/\r?\n//g;$name=~s/\s+//g;print"$gene"."\t"."$num"."\t"."\/"."$name"."\/";' | grep -f other/color_four_$f.tsv | sed 's/\///g'  > tem &&
mv tem select/$f.select.txt
done


#查看颜色为5（共78种）
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat select/$base.select.txt | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR"  | tsv-filter --eq 2:5
done | cut -f 1 >color_five.tsv

#提取含有braB/braZ的颜色
for PROJECT in $(cat color_five.tsv);do
braZ=$( cat select/$PROJECT.select.txt | grep  "braZ" | head -n 1 | cut -f 3 )
braB=$( cat select/$PROJECT.select.txt | grep  "braB" | head -n 1 | cut -f 3 )
echo -e "$braB\n$braZ" | perl -alne '$name=(split/\s/,$_)[0];$name=~s/\r?\n//g;$name=~s/\s//g;print"\/"."$name"."\/";'   >other/color_five_$PROJECT.tsv
done  

#过滤
for f in $(cat color_five.tsv)
do
echo -e  "$f"
cat select/$f.select.txt |perl -alne '($gene,$num,$name)=(split/\t/,$_)[0,1,2];$name=~s/\r?\n//g;$name=~s/\s+//g;print"$gene"."\t"."$num"."\t"."\/"."$name"."\/";' | grep -f other/color_five_$f.tsv | sed 's/\///g' > tem &&
mv tem select/$f.select.txt
done


#查看颜色为6（共18种）
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat select/$base.select.txt | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR"  | tsv-filter --eq 2:6
done | cut -f 1 >color_six.tsv

#提取含有braB/braZ的颜色
for PROJECT in $(cat color_six.tsv);do
braZ=$( cat select/$PROJECT.select.txt | grep  "braZ" | head -n 1 | cut -f 3 )
braB=$( cat select/$PROJECT.select.txt | grep  "braB" | head -n 1 | cut -f 3 )
echo -e "$braB\n$braZ" | perl -alne '$name=(split/\s/,$_)[0];$name=~s/\r?\n//g;$name=~s/\s//g;print"\/"."$name"."\/";'   >other/color_six_$PROJECT.tsv
done  

#过滤
for f in $(cat color_six.tsv)
do
echo -e  "$f"
cat select/$f.select.txt |perl -alne '($gene,$num,$name)=(split/\t/,$_)[0,1,2];$name=~s/\r?\n//g;$name=~s/\s+//g;print"$gene"."\t"."$num"."\t"."\/"."$name"."\/";'  | grep -f other/color_six_$f.tsv | sed 's/\///g'  > tem &&
mv tem select/$f.select.txt
done


#查看颜色为8(共3种)
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat select/$base.select.txt | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR"  | tsv-filter --eq 2:8
done | cut -f 1 >color_eight.tsv

#提取含有braB/braZ的颜色
for PROJECT in $(cat color_eight.tsv);do
braZ=$( cat select/$PROJECT.select.txt | grep  "braZ" | head -n 1 | cut -f 3 )
braB=$( cat select/$PROJECT.select.txt | grep  "braB" | head -n 1 | cut -f 3 )
echo -e "$braB\n$braZ" | perl -alne '$name=(split/\s/,$_)[0];$name=~s/\r?\n//g;$name=~s/\s//g;print"\/"."$name"."\/";'  >other/color_eight_$PROJECT.tsv
done  

#过滤
for f in $(cat color_eight.tsv)
do
echo -e  "$f"
cat select/$f.select.txt |perl -alne '($gene,$num,$name)=(split/\t/,$_)[0,1,2];$name=~s/\r?\n//g;$name=~s/\s+//g;print"$gene"."\t"."$num"."\t"."\/"."$name"."\/";' | grep -f other/color_eight_$f.tsv | sed 's/\///g'  > tem &&
mv tem select/$f.select.txt
done

###最后再统计结果发现
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
COLOUR=$(cat select/$base.select.txt  | cut -f 3| sort -n | uniq | wc -l)
echo -e "$base\t$COLOUR" 
done | tsv-summarize -g 2 --count
2	34147
1 22516


#合并结果
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
echo $base
braB=$(cat select/$base.select.txt| grep "PA1590"| cut -f 3)
braZ=$(cat select/$base.select.txt | grep "PA1971"| cut -f 3)
cat select/$base.select.txt| grep "$braB" >> PA1590.select.txt
cat select/$base.select.txt | grep "$braZ" >> PA1971.select.txt
done
```
* 不能用
```bash
#划分区间统计基因出现的次数
cat PA1590.select.txt | grep -v "grey" | tsv-summarize  -g 1 --count >braB_whole_number.tsv
cat braB_whole_number.tsv | tsv-summarize -g 2 --count | tsv-sort -k1,1n >braB_per_gene_num.tsv
for i in $(seq 12500 100 50100)
do
j=$[$i+100]
num=$(cat braB_per_gene_num.tsv | tsv-filter --lt 1:$j| tsv-filter --ge 1:$i| tsv-summarize --sum 2)
echo -e "$num\t$i\t$j" 
done >braB_bins_100.tsv

cat PA1971.select.txt | grep -v "grey" | tsv-summarize  -g 1 --count >braZ_whole_number.tsv
cat braZ_whole_number.tsv | tsv-summarize -g 2 --count | tsv-sort -k1,1n >braZ_per_gene_num.tsv
for i in $(seq 12500 100 50100)
do
j=$[$i+100]
num=$(cat braZ_per_gene_num.tsv | tsv-filter --lt 1:$j| tsv-filter --ge 1:$i| tsv-summarize --sum 2)
echo -e "$num\t$i\t$j" 
done >braZ_bins_100.tsv


#每10个数相加统计基因出现的次数
for f in $(seq 10 10 3920)
do
min=$(cat braB_per_gene_num.tsv | cut -f 1 | head -n $f | tail -n 10| head -n 1)
max=$(cat braB_per_gene_num.tsv | cut -f 1 | head -n $f | tail -n 10| tail -n 1)
num=$(cat braB_per_gene_num.tsv | cut -f 2 | head -n $f | tail -n 10| tsv-summarize --sum 1)
echo -e "$min\t$max\t$num"
done >braB_per_10.tsv

for f in $(seq 10 10 3820)
do
min=$(cat braZ_per_gene_num.tsv | cut -f 1 | head -n $f | tail -n 10| head -n 1)
max=$(cat braZ_per_gene_num.tsv | cut -f 1 | head -n $f | tail -n 10| tail -n 1)
num=$(cat braZ_per_gene_num.tsv | cut -f 2 | head -n $f | tail -n 10| tsv-summarize --sum 1)
echo -e "$min\t$max\t$num"
done >braZ_per_10.tsv

#统计不同基因出现的次数
cat PA1590.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | grep "PA1590"  #PA1590_braB_at  47285
cat PA1590.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | tsv-filter --ge 2:23500 | wc -l  #241
cat PA1971.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | grep "PA1971"  #PA1971_braZ_at  50052
cat PA1971.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | tsv-filter --ge 2:23500 | wc -l  #252
```

# 6.KEGG和GO分析(不能用!!!!!!)
## 6.1使用david对gene_id进行转换(方法一) (不能用!!!!!!)
```bash
#GO和KEGG前需要对gene_id转换成pfam(locus tag)
打开david---shortcut to david tools---gene id conversion
upload---step1(选择文件并上传)---step2(locus tag)---step3(gene list)---step4(submit list)
option1 convert the gene list to(Pfam/locus tag)---for species(Pseudomonas aeruginosa PAO1)---submit to conversion tool---点击结果界面的(convert all) 
#下载braB结果文件
# 提交了171个gene id，共有163个被识别
cd /mnt/d/WGCNA/KEGG_GO/database/
wget https://david.ncifcrf.gov/data/download/conv_0C66F0A958F11651041251276.txt -O braB_Pfam_id.tsv
#结果文件第一列是原ID,第二列是转换后的ID，第三列是物种，第四列是基因名和基因描述
#下载braB结果文件
# 提交了177个gene id,共有174个被识别
wget https://david.ncifcrf.gov/data/download/conv_FBD558AE3C581650714720192.txt -O braZ_Pfam_id.tsv
#结果文件第一列是原ID,第二列是转换后的ID，第三列是物种，第四列是基因名和基因描述
#匹配最后一个括号(.+)(\()
cat braB_Pfam_id.tsv | tsv-select -f 1,4 |grep -v "From" | perl -alne '($name,$gene)=(split/\s+/,$_,2)[0,1];$gene=~/(.+)(\()(.*?)\)$/;$id=$3;print"$name\t$id";' | sed '1iPA\tgene'  >braB_gene_id.tsv
cat braZ_Pfam_id.tsv | tsv-select -f 1,4 |grep -v "From" | perl -alne '($name,$gene)=(split/\s+/,$_,2)[0,1];$gene=~/(.+)(\()(.*?)\)$/;$id=$3;print"$name\t$id";'  | sed '1iPA\tgene' >braZ_gene_id.tsv
```

## 6.2下载注释数据库 (不能用!!!!!!)
### 6.2.1GO注释库 (不能用!!!!!!)
* GOC（Gene Ontology Consortium）提供了41种不同模型生物的GAF格式的注释信息网址[GO的gaf文件](http://current.geneontology.org/products/pages/downloads.html)
* [GO的注释文件](http://geneontology.org/docs/download-ontology)
* GO官网提供了几种常见物种对应的GO注释信息，文件格式为GAF
* EBI对uniprot数据库中的蛋白进行了GO注释分析，这个项目名为gene ontology annotation, 简称GOA, 在FTP也提供了物种对应的注释
* 在NCBI检索基因时，在结果页面会看到该基因对应的很多注释信息，其中就包括了GO注释。gene2go就是基因对应的GO注释文件，这个文件包含了所有物种的GO信息，可以根据物种对应的tax id提取指定物种。NCBI中用Entrez Id 标识每个基因，通过另外的几个文件，可以得到Entrez ID, Ensemble Id, Gene Symbol对应的GO信息
* 从bioconductor获取对饮的注释包(但是只有模式生物的) 通过keys和select函数可以获得基因对应的GO注释信息
* GO 注释主要有两种方法：序列相似性比对（BLAST）和结构域相似性比对（InterProScan）
* gaf文件内容:第一列是DB基因标识的来源数据库,第二列是DB object ID是数据库对应的唯一标识符,第三列是db object symbol是对应的基因名，第四列是qualifier，第五列是GO ID,第六列是DB：Reference，注释的证据来源，一般为文献参考,第十一列是基因symbol ID，第十二是DB object type蛋白产物,第十三列是物种的taxonomic标识符，使用数字编码来代表某个物种。注释参考数据库是PseudoCAP
* 对于非模式生物或者无参考基因组的项目，经常需要进行基因的功能注释，而GO注释是基因功能注释的重要部分。有很多软件能够获得GO注释的信息，例如interproscan、eggnog-mapper和blas2go等。
* 利用interproscan获得非模式生物的GO注释信息，一般需要将注释信息整理成gene2go，go2gene以及wego的格式，有时也要提取GO Level2的GOID进行分析
```bash
cd WGCNA/KEGG_GO/database
#铜绿的gaf文件
 wget http://current.geneontology.org/annotations/pseudocap.gaf.gz
 gunzip pseudocap.gaf.gz
 cat pseudocap.gaf | grep -v "!" | cut -f 2,3,5 | sed '1idb_id\tsymbol\tGO_id' >GO_id_PA_id_Pse.tsv
 
 rm -rf pseudocap.gaf
 #所有的的go注释文件
wget http://current.geneontology.org/ontology/go-basic.obo
python GO_id_descri_abstract.py -i go-basic.obo  -o GO_id_descri.tsv
rm -rf go-basic.obo 
#GO富集拆分成BP，MF，CC三个数据框

#合并PA号，GO号和GO descri等
cat  GO_id_PA_id_Pse.tsv | sed 's/GO_id/GO/g' | tsv-join -H -d GO -f GO_id_descri.tsv -k GO --append-fields Description,level -z | tsv-select -f 1,3,4,5 >GO_PA_descri.tsv

#下载go2gene.tsv对应的gene名(1508共有1497被注释出来了)
wget https://david.ncifcrf.gov/data/download/conv_7134C2C613F91650772068066.txt -O GO_PA_gene.tsv
cat GO_PA_gene.tsv | tsv-select -f 1,4 |grep -v "From" | perl -alne '($name,$gene)=(split/\s+/,$_,2)[0,1];$gene=~/(.+)(\()(.*?)\)$/;$id=$3;print"$name\t$id";'  >GO_PA_sed.tsv
#替换注释信息
cat GO_PA_sed.tsv | 
    perl -nla -e '
        print q{s/^} . quotemeta($F[0]) . q{/} . quotemeta($F[1]) . q{/g;};
    ' \
    > GO_sed.script
cat GO_PA_descri.tsv | sed -f GO_sed.script >GO_PA_replace.tsv

tsv-select -f  2,1,4 GO_PA_replace.tsv >go2gene.tsv
cut -f 2,3,4 GO_PA_replace.tsv >go2name.tsv
```
### 6.2.2KEGG注释库 (不能用!!!!!!)
 * KEGG API是一个连接KEGG各类数据库的应用程序，主要以URL形式进行访问[kegg](https://www.kegg.jp/kegg/rest/keggapi.html)
 * 查看kegg中的铜绿假单胞菌的信息[铜绿信息](https://www.kegg.jp/kegg-bin/show_organism?org=pae)
 * [铜绿标识符](https://www.kegg.jp/entry/gn:T00035)
 * [铜绿假单胞菌的数据库](https://www.genome.jp/dbget-bin/get_linkdb?genome+T00035)
 * [KEGG的ID转换](https://www.genome.jp/kegg/mapper/convert_id.html)

```bash
cd WGCNA/KEGG_GO/database
#直接把连接网址放在浏览器，下载kegg注释的json文件
https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir=
https://www.genome.jp/kegg-bin/download_htext?htext=br08610&format=htext&filedir=
#获得KO与pathway的关系,这会生成KEGG_pathway_ko.txt文件，随后对行去重,最后得到KEGG_pathway_ko_uniq.txt文件，这个文件包含了KO和KEGG pathway的对应关系信息，也包含了pathway的级别分类(KEGG pathway分为3级)，如下所示：
python KO_id_descri_abstract.py
python KO_id_descri_uniq.py
#下载kegg注释信息
wget https://www.genome.jp/ftp/db/kofam/ko_list.gz -O ./Pse_KO_annotation.gz
gunzip Pse_KO_annotation.gz 
mv Pse_KO_annotation KO_id_descri.tsv
#下载PA号对应的KEGG信息(5697)
wget https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+gn:T00035 -O KO_1.tsv
wget https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+-p+2+gn:T00035 -O KO_2.tsv
wget https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+-p+3+gn:T00035 -O KO_3.tsv
wget https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+-p+4+gn:T00035 -O KO_4.tsv
wget https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+-p+5+gn:T00035 -O KO_5.tsv
wget https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+-p+6+gn:T00035 -O KO_6.tsv
#提取想要的注释信息
cat KO_1.tsv | grep "K" | grep "<a" | perl -alne '($pa,$ko,$des)=(split/\s+/,$_,4)[1,2,3];$pa=~m/\>pae\:(.*?)\</;$name=$1;print"$name\t$ko\t$des";' >KO_1_anno.tsv
cat KO_2.tsv | grep "K" | grep "<a" | perl -alne '($pa,$ko,$des)=(split/\s+/,$_,4)[1,2,3];$pa=~m/\>pae\:(.*?)\</;$name=$1;print"$name\t$ko\t$des";' >KO_2_anno.tsv
cat KO_3.tsv | grep "K" | grep "<a" | perl -alne '($pa,$ko,$des)=(split/\s+/,$_,4)[1,2,3];$pa=~m/\>pae\:(.*?)\</;$name=$1;print"$name\t$ko\t$des";' >KO_3_anno.tsv
cat KO_4.tsv | grep "K" | grep "<a" | perl -alne '($pa,$ko,$des)=(split/\s+/,$_,4)[1,2,3];$pa=~m/\>pae\:(.*?)\</;$name=$1;print"$name\t$ko\t$des";' >KO_4_anno.tsv
cat KO_5.tsv | grep "K" | grep "<a" | perl -alne '($pa,$ko,$des)=(split/\s+/,$_,4)[1,2,3];$pa=~m/\>pae\:(.*?)\</;$name=$1;print"$name\t$ko\t$des";' >KO_5_anno.tsv
cat KO_6.tsv | grep "K" | grep "<a" | perl -alne '($pa,$ko,$des)=(split/\s+/,$_,4)[1,2,3];$pa=~m/\>pae\:(.*?)\</;$name=$1;print"$name\t$ko\t$des";' >KO_6_anno.tsv
#合并并删除
cat KO_1_anno.tsv KO_2_anno.tsv KO_3_anno.tsv KO_4_anno.tsv KO_5_anno.tsv KO_6_anno.tsv >KO_id_PA_id_Pse.tsv
rm -rf KO_1.tsv KO_2.tsv KO_3.tsv KO_4.tsv KO_5.tsv  KO_6.tsv KO_1_anno.tsv KO_2_anno.tsv KO_3_anno.tsv KO_4_anno.tsv KO_5_anno.tsv KO_6_anno.tsv

#不需要合并PA号，KO号和KO descri,该文件KO_id_PA_id_Pse.tsv中含有
sed -i '1igene\tKO\tdescri' KO_id_PA_id_Pse.tsv
#下载ko2gene.tsv对应的gene名(5697共有5697被注释出来)
wget https://david.ncifcrf.gov/data/download/conv_7134C2C613F91650773243799.txt -O KO_PA_gene.tsv
cat KO_PA_gene.tsv | tsv-select -f 1,4 |grep -v "From" | perl -alne '($name,$gene)=(split/\s+/,$_,2)[0,1];$gene=~/(.+)(\()(.*?)\)$/;$id=$3;print"$name\t$id";'  >KO_PA_sed.tsv
#替换注释信息
cat KO_PA_sed.tsv |
    perl -nla -e '
        print q{s/^} . quotemeta($F[0]) . q{/} . quotemeta($F[1]) . q{/g;};
    ' \
    > KO_sed.script
cat KO_id_PA_id_Pse.tsv | sed -f KO_sed.script >KO_PA_replace.tsv

tsv-select -f  2,1 KO_PA_replace.tsv >ko2gene.tsv
cut -f 2,3 KO_PA_replace.tsv | cut -d "|" -f 1 | cut  -d "[" -f 1 >ko2name.tsv

 #下载kegg物种列表
wget -c "http://rest.kegg.jp/list/organism" -O species.txt
grep "Pseudomonas aeruginosa" species.txt >Pse_KEGG_info.tsv
#查看物种对应的Ontology (39)
wget https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+gn:T00035 -O ./Pse_KO.tsv

#查看物种对应的pathway(126)
wget https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+gn:T00035 -O ./Pse_pathway.tsv
#查看物种对应的gene号(5697)
wget https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+gn:T00035 -O ./Pse_gene_id_kegg.tsv
```

## 6.3构建非模式物种的OrgDb  (不能用!!!!!!)
* 要进行GO或者KEGG富集分析，就需要知道每个基因对应什么样的GO/KEGG分类，OrgDb就是存储不同数据库基因ID之间对应关系，以及基因与GO等注释的对应关系的 R 软件包.如果自己研究的物种不在[OrgDb库](http://bioconductor.org/packages/release)之列，很大可能就需要自己构建OrgDb，然后用clusterProfiler分析
* 非模式生物要想找到自己的注释包，又分成两类;一类是在AnnotationHub（https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html）中存在的.另一类是在AnnotationHub也不存在相应物种，就需要用AnnotationForge（ https://bioconductor.org/packages/release/bioc/html/AnnotationForge.html ）来自己构建
* 如果说一个data.frame中的元素是factor，你想转化成numeric，你会怎么做？比如d[1,1]是factor。正确答案是 先as.character(x) 再as.numeric(x)。如果直接as.numeric,就不是以前的数字了，原来as.data.frame()有一个参数stringsAsFactors。如果stringAsFactor=F，就不会把字符转换为factor 这样以来，原来看起来是数字变成了character，原来是character的还是character

* GO可以通过读取外部的GO注释文件进行分析.格式:GeneID	GO	GO_Description
* clusterProfiler包只针对含有OrgDb对象，如果是公共数据库中有该物种注释信息，只是未制作成org.db数据库（标准注释库），则可以不需要从头注释，只需手动制作org.db数据库类型，完成后直接使用即可
* 同样地，对于pathway数据库中没有的物种，也支持读取基因的pathway注释文件，然后进行分析，注释文件的格式如下格式:GeneID	Pathway	Path_Description
```r
BiocManager::install("AnnotationHub") 
BiocManager::install("biomaRt")
library(AnnotationHub) 
library(biomaRt)
options(stringsAsFactors = F)
#建立AnnotationHub对象
hub <- AnnotationHub() 
#查看AnonotationHub里面物种
unique(hub$species) 
#看AnonotationHub里是否包含想要的物种
hub$species[which(hub$species=="Pseudomonas aeruginosa")] 
#查看该物种信息
query(hub, "Pseudomonas aeruginosa")  
 #OrgDb属于rdataclass中，因此查看下该物种有没有OrgDb
hub[hub$species=="Pseudomonas aeruginosa" & hub$rdataclass == "OrgDb"]
#AH87086是假单胞菌对应的编号
Pse.OrgDb <- hub[["AH87086"]]
#查看注释信息
Pse.OrgDb
```

## 6.4进行GO和KEGG富集分析  (不能用!!!!!!)
```r
# 安装包
BiocManager::install("clusterProfiler")  #用来做富集分析
BiocManager::install("topGO")  #画GO图用
BiocManager::install("Rgraphviz")#展示结构化信息
BiocManager::install("pathview") #看KEGG pathway
BiocManager::install("ggnewscale") 
# clusterProfiler 是Y叔写的一个功能强大的R包，可以用来做各种富集分析，如GO、KEGG、DO
library(clusterProfiler)
library(topGO)
library(Rgraphviz) 
library(pathview) 
library(enrichplot)
library(ggnewscale)
setwd("D:/WGCNA/KEGG_GO/database")
output.gene_id <- read.table("braB_gene_id.tsv", sep='\t', row.names=1, header=T,
                     quote="", comment="", check.names=F)
gene_list = as.character(output.gene_id$gene)
#通过导入外部注释文件富集分析: 读取GO注释文件
ko2name <- read.delim('ko2name.tsv', stringsAsFactors=FALSE)
ko2gene <- read.delim('ko2gene.tsv', stringsAsFactors=FALSE)
go2name <- read.delim('go2name.tsv', stringsAsFactors=FALSE)
go2gene <-read.delim('go2gene.tsv', stringsAsFactors=FALSE)
#拆分成BP，MF，CC三个数据框
go2gene = split(go2gene , with(go2gene,level))
#以MF为例
enricher_GO_MF <- enricher(gene_list,TERM2GENE=go2gene$molecular_function,TERM2NAME=go2name,pvalueCutoff =1,qvalueCutoff=1,pAdjustMethod = "none",minGSSize=2,maxGSSize=400)
#以BP为例
enricher_GO_BP <- enricher(gene_list,TERM2GENE=go2gene$biological_process,TERM2NAME=go2name,pvalueCutoff =1,qvalueCutoff=1,pAdjustMethod = "none",minGSSize=2,maxGSSize=400)
#以CC为例
enricher_GO_CC <- enricher(gene_list,TERM2GENE=go2gene$cellular_component,TERM2NAME=go2name,pvalueCutoff =1,qvalueCutoff=1,pAdjustMethod = "none",minGSSize=2,maxGSSize=400)
#KEGG富集
enricher_KO <- enricher(gene_list,TERM2GENE = ko2gene,TERM2NAME = ko2name,pvalueCutoff =1,qvalueCutoff=1,pAdjustMethod = "none",minGSSize=2,maxGSSize=400)

#gene差异基因对应的向量;
#keyType指定基因ID的类型，默认为ENTREZID, 可参考keytypes(org.Hs.eg.db)类型 ；
#OrgDb指定该物种对应的org包的名字；
#ont代表GO的3大类别，BP, CC, MF，也可是全部ALL；
#pAdjustMethod指定多重假设检验矫正的方法，有“ holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种；
#cufoff指定对应的阈值；
#readable=TRUE代表将基因ID转换为gene symbol。
```


## 6.5GO富集分析结果可视化  (不能用!!!!!!)
```r
#barplot #默认展示显著富集的top10个，即p.adjust最小的10个，聚类条形图
GO_MF_1<-barplot(enricher_GO_MF, showCategory = 10) 
pdf(file = "GO_MF_1.pdf",width =12,height = 9)
print(GO_MF_1)
dev.off()
#dotplot 按富集的数从大到小的聚类气泡图
GO_MF_2<-dotplot(enricher_GO_MF, showCategory = 10)
pdf(file = "GO_MF_2.pdf",width =12,height = 9)
print(GO_MF_2)
dev.off()
#GO term与差异基因关系网络图
GO_MF_3<-cnetplot(enricher_GO_MF, showCategory = 10,categorySize="pvalue")
pdf(file = "GO_MF_3.pdf",width =12,height = 9)
print(GO_MF_3)
dev.off()
#GO terms关系网络图（通过差异基因关联）
enricher_GO_MF_1<-pairwise_termsim(enricher_GO_MF)
GO_MF_4<-emapplot(enricher_GO_MF_1,showCategory = 10)
pdf(file = "GO_MF_4.pdf",width =12,height = 9)
print(GO_MF_4)
dev.off()
#热图
GO_MF_5<-enrichplot::heatplot(enricher_GO_MF,showCategory =100)
pdf(file = "GO_MF_5.pdf",width =12,height = 9)
print(GO_MF_5)
dev.off()

#barplot #默认展示显著富集的top10个，即p.adjust最小的10个，聚类条形图
GO_BP_1<-barplot(enricher_GO_BP, showCategory = 10) 
pdf(file = "GO_BP_1.pdf",width =12,height = 9)
print(GO_BP_1)
dev.off()
#dotplot 按富集的数从大到小的聚类气泡图
GO_BP_2<-dotplot(enricher_GO_BP, showCategory = 10)
pdf(file = "GO_BP_2.pdf",width =12,height = 9)
print(GO_BP_2)
dev.off()
#GO term与差异基因关系网络图
GO_BP_3<-cnetplot(enricher_GO_BP, showCategory = 10,categorySize="pvalue")
pdf(file = "GO_BP_3.pdf",width =12,height = 9)
print(GO_BP_3)
dev.off()
#GO terms关系网络图（通过差异基因关联）
enricher_GO_BP_1<-pairwise_termsim(enricher_GO_BP)
GO_BP_4<-emapplot(enricher_GO_BP_1,showCategory = 10)
pdf(file = "GO_BP_4.pdf",width =12,height = 9)
print(GO_BP_4)
dev.off()
#热图
GO_BP_5<-enrichplot::heatplot(enricher_GO_BP,showCategory =100)
pdf(file = "GO_BP_5.pdf",width =12,height = 9)
print(GO_BP_5)
dev.off()


GO_CC_1<-barplot(enricher_GO_CC, showCategory = 10) 
pdf(file = "GO_CC_1.pdf",width =12,height = 9)
print(GO_CC_1)
dev.off()
#dotplot 按富集的数从大到小的聚类气泡图
GO_CC_2<-dotplot(enricher_GO_CC, showCategory = 10)
pdf(file = "GO_CC_2.pdf",width =12,height = 9)
print(GO_CC_2)
dev.off()
#GO term与差异基因关系网络图
GO_CC_3<-cnetplot(enricher_GO_CC, showCategory = 10,categorySize="pvalue")
pdf(file = "GO_CC_3.pdf",width =12,height = 9)
print(GO_CC_3)
dev.off()
#GO terms关系网络图（通过差异基因关联）
enricher_GO_CC_1<-pairwise_termsim(enricher_GO_CC)
GO_CC_4<-emapplot(enricher_GO_CC_1,showCategory = 10)
pdf(file = "GO_CC_4.pdf",width =12,height = 9)
print(GO_CC_4)
dev.off()
#热图
GO_CC_5<-enrichplot::heatplot(enricher_GO_CC,showCategory =100)
pdf(file = "GO_CC_5.pdf",width =12,height = 9)
print(GO_CC_5)
dev.off()

#Pathway富集分析结果可视化
#barplot
KEGG_1<-barplot(enricher_KO, showCategory = 10)
pdf(file = "KEGG_1.pdf",width =12,height = 9)
print(KEGG_1)
dev.off()
#dotplot
KEGG_2<-dotplot(enricher_KO, showCategory = 10)
pdf(file = "KEGG_2.pdf",width =12,height = 9)
print(KEGG_2)
dev.off()
#pathway与差异基因关系网络图
KEGG_3<-cnetplot(enricher_KO, showCategory = 10,categorySize="pvalue")
pdf(file = "KEGG_3.pdf",width =12,height = 9)
print(KEGG_3)
dev.off()
#pathway关系网络图（通过差异基因关联）
enricher_KO_1<-pairwise_termsim(enricher_KO)
KEGG_4<-emapplot(enricher_KO_1,showCategory = 10)
pdf(file = "KEGG_4.pdf",width =12,height = 9)
print(KEGG_4)
dev.off()
#热图
KEGG_5<-enrichplot::heatplot(enricher_KO,showCategory =100)
pdf(file = "KEGG_5.pdf",width =13,height =12)
print(KEGG_5)
dev.off()
```

# 7.使用mclust进行聚类分析
## 7.1mclust使用举例
* mclust(Model-based clustering) 能够基于高斯有限混合模型进行聚类，分类以及密度估计(density estimation)。对于具有各种协方差结构的高斯混合模型，它提供了根据EM算法的参数预测函数。它也提供了根据模型进行模拟的函数。还提供了一类函数，整合了基于模型的层次聚类，混合估计的EM算法，用于聚类、密度估计和判别分析中综合性策略的贝叶斯信息判别标准。最后还有一类函数能够对聚类，分类和密度估计结果中的拟合模型进行可视化展示
*  Mclust假设观测数据是一个或多个混合高斯分布的抽样结果，Mclust就需要根据现有数据去推断最优可能的模型参数，以及是由 q几组分布抽样而成。mclust一共提供了14种模型（见下表），可以用mclustModelNames可以查看mclust提供的所有模型
*  继续回到之前的问题，Mclust如何确定模型和确定分组数目。之前我们调用Mclust时，除了必须设置的输入参数，没有修改其他参数。其实Mclust可以设置的参数不少，和问题直接相关的是如下两个参数  G: 分组数，默认情况下是1:9  modelNames: 待拟合的模型，默认使用所有14种  。也就是，Mclust默认得到14种模型1到9组的分析结果，然后根据一定的标准选择最终的模型和分组数
*  Mclust提供了两种方法用于评估不同模型在不同分组下的可能性。BIC( Bayesian Information Criterion ): 贝叶斯信息判别标准和ICL( integrated complete-data likelihood ): 综合完全数据可能性
* Mclust默认用的就是BIC，因此我们可以用plot.Mclust绘制其中BIC变化曲线
```r
# 安装
install.packages("mclust")
# 加载
library(mclust)
#加载数据集
install.packages("gclus")
data("wine", package = "gclus")
dim(wine)
# 第一列和聚类无关
#使用Mclust做聚类分析. Mclust主要功能就是分析当前的提供的数据是由什么统计模型
X <- data.matrix(wine[,-1])
mod <- Mclust(X)
#直接在交互行输入mod会得到如下信息
'Mclust' model object: (VVE,3) 
Available components: 
 [1] "call"           "data"           "modelName"     
 [4] "n"              "d"              "G"             
 [7] "BIC"            "bic"            "loglik"        
[10] "df"             "hypvol"         "parameters"    
[13] "z"              "classification" "uncertainty" 
#第一行告诉我们mclust以WE模型将数据分为3类。第3行开始，它告诉我们'Mclust'的输出结果中包含了如下内容，我们可以通过$来提取。举个例子，我们提取Mclust的聚类结果和已知结果进行比较
table(wine$Class, mod$classification)
# 如下是输出信息
     1  2  3
  1 59  0  0
  2  0 69  2
  3  0  0 48
# adjustedRandIndex:评估聚类效果
adjustedRandIndex(wine$Class, mod$classification)
#从结果中，我们发现仅有2例没有正确聚类，说明Mclust的效果很好。但是随之而来的问题是，Mclust如何挑选模型以及它为什么认为聚成3类比较合适呢？我们可以根据什么信息进行模型选择呢？ -----------  BIC变化曲线
#Mclust默认用的就是BIC，因此我们可以用plot.Mclust绘制其中BIC变化曲线
plot.Mclust(mod, what = "BIC", ylim = range(mod$BIC[,-(1:2)], na.rm = TRUE), legendArgs = list(x = "bottomleft", cex =0.7))
#Mclucst会选择其中BIC最大的模型和分组作为最终的结果
#此外我们可以用MclustBIC和MclustICL分别进行计算
par(mfrow=c(1,2))
BIC <- mclustBIC(X)
ICL <- mclustICL(X)
#从中选择最佳的模型分组和模型作为输入
mod2 <- Mclust(X, G = 3, modelNames = "VVE", x=BIC)
```
## 7.2三三组合的结果文件划分bins为100然后进行mclust分析
```r
cut -f 1,3 braB_bins_100.tsv | perl -alne '($num,$bin)=(split/\s+/,$_)[0,1];@data=($bin)x$num;print join("\n",@data);' | sed '/^$/d'  >braB_bins_mclust.tsv
#braB_bins_100  377行
#braZ_bins_100  377行

###braB选择G和模型这两个参数
library(mclust)
ratio<-read.table("braZ_density_ggplot.tsv",header=T)
attach(ratio)
#mclust分析
BIC <- mclustBIC(ratio$ratio)
ICL <- mclustICL(ratio$ratio)
write.table(BIC, "BIC.tmp.tsv", row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(ICL, "ICL.tmp.tsv", row.names = FALSE, col.names = TRUE, sep = "\t")
#将BIC和ICL合并
cat BIC.tmp.tsv | perl -n -e 'chomp;
@a = split/\t/, $_;if($_ =~ /"/){$i = 1;next;}else{
print "$i\t$a[0]\tBICE\n"; print "$i\t$a[1]\tBICV\n";
$i++;}' | tsv-filter --str-ne 2:NA > BIC2.tmp.tsv
cat ICL.tmp.tsv | perl -n -e 'chomp;
@a = split/\t/, $_;if($_ =~ /"/){$i = 1;next;}else{
print "$i\t$a[0]\tICLE\n"; print "$i\t$a[1]\tICLV\n";
$i++;}' | tsv-filter --str-ne 2:NA > ICL2.tmp.tsv
cat ICL2.tmp.tsv BIC2.tmp.tsv | sed '1inum\tgrade\tmethod' > braZ_mclust_model.tsv
rm *.tmp.*
#使用BIC和ICL画图
library(readr)
library(ggplot2)
num <- read_tsv("braZ_mclust_model.tsv", show_col_types = FALSE)
p <- ggplot()+
geom_line(data = num, aes(x = num, y = grade, group = method, color = method)) +
scale_x_continuous(breaks = seq(0,9,1))
ggsave(p, file = "braZ_mclust.pdf", width = 5, height = 4)
```

## 7.3 mclust计算均值和平方差
```r
#检验是否符合正态分布
#1.绘制QQ图  ：qqnorm()可以绘制QQ图。通过绘制的图是否呈现一直线判断是否符合正态分布。另外还有一个qqline()函数，在QQ图中绘制一条直线，QQ图中的点越接近这条直线，表示数据越接近正态分布。
qqnorm(ratio$ratio)  qqline(ratio$ratio)
#2.夏皮罗检验shapiro.test()函数 p-value反应服从正态分布的概率，值越小越小的概率符合，通常0.05做标准，大于0.05则表示符合正态分布
nortest1<-shapiro.test(ratio$ratio) #缺失值应该在3至5000之间
#3.ks检验
ks.test(ratio$ratio,rnorm(5900,mean = mean(ratio$ratio), sd = sd(ratio$ratio)))
#ks >0.05，说明数据符合正太分布
#以data的平均值为平均值，以data的方差为方差模拟一个新的正态分布数据，和data做比较，看data符合不符合正态分布。

###braB选择1个峰
ratio<-read.table("braB_bins_mclust.tsv",header=T)
dens <- densityMclust(ratio,G=1)
summary(dens, parameters = TRUE)
br <- seq(min(ratio$ratio), max(ratio$ratio), length =339)
#hist(ratio$ratio,breaks=br)
x <- seq(min(ratio$ratio)-diff(range(ratio$ratio))/10,max(ratio$ratio)+diff(range(ratio$ratio))/10, length = 500)
cdens <-predict(dens, x, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*dens$parameters$pro))
plot(dens, what = "density", data = ratio$ratio, breaks = br)
###matplot(x, cdens, type = "l", lwd = 1, add = TRUE, lty = 1)
#平均值和方差
Means:
[1] 20121.66
Variances:
[1] 3509466
sd
[1] 1873.36
1873.36 x 3 +20121.66=25741.74

###braZ选择1个峰
ratio<-read.table("braZ_bins_mclust.tsv",header=T)
dens <- densityMclust(ratio,G=1)
summary(dens, parameters = TRUE)
br <- seq(min(ratio$ratio), max(ratio$ratio), length =339)
x <- seq(min(ratio$ratio)-diff(range(ratio$ratio))/10,max(ratio$ratio)+diff(range(ratio$ratio))/10, length = 500)
cdens <-predict(dens, x, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*dens$parameters$pro))
plot(dens, what = "density", data = ratio$ratio, breaks = br)
#平均值和方差
Means:
[1] 20483.46
Variances:
[1] 3320404
sd
[1] 1822.20
1822.20 x 3 +20483.46=25950.06
```

## 7.4ggplot正态分布直方密度图
```r
#绘制braB的正态密度分布直方图
cut -f 2 braB_whole_number.tsv >braB_density_ggplot.tsv
library(ggplot2)
setwd("D:/WGCNA/rma_WGCNA/condition_three_combine/result")
ratio<-read.table("braB_density_ggplot.tsv",header=T)
set.seed(1000)
df <- data.frame(ratio)
#geom_density会为数据集提供一条密度曲线，而不是正态分布,而stat_function(fun = dnorm, args = list(mean = mean(df$ratio), sd = sd(df$ratio)))提供正态分布密度曲线
#geom_histogram  stat="bin"直方图类型  position="stack" 位置  ...其他geom类函数的参数  binwidth #直方图的图距    bins=NULL直方个数
#...count...参数画计数直方图     ...density...频率分布直方图  alpha=0.9透明度
data_mean=mean(df$ratio) #20070.71
data_sd=sd(df$ratio) #1873.92
three_sd_plus=mean(df$ratio)+data_sd*3  #25692.47
three_sd_minus=mean(df$ratio)-data_sd*3 #14448.95

p1<-ggplot(df, aes(x = ratio)) +geom_histogram(aes(y =..density..),breaks = seq(min(df$ratio), max(df$ratio), length =100),colour ="#C6DBEF", fill ="#08519C", size = 0.1,alpha=0.8) + #使用density代替y轴
stat_function(fun = dnorm, args = list(mean = data_mean, sd = data_sd,color ="black", size = 1))+  #直方图上加密度曲线，也可以用geom_density()
geom_vline(aes(xintercept=data_mean), color="red", linetype="dashed", size=1)+ #添加均值线和两个标准差线
geom_vline(aes(xintercept=three_sd_plus), color="red", linetype="dashed", size=1)+
annotate(geom="text",fontface="bold",color="black",x=three_sd_plus,y=0.00015,label="25692.47",size=5)+   
 #通过annotate函数添加注释，将geom变量设置为"text"和"rect"分别代表文字与矩形，并且可以调整位置，大小颜色
annotate(geom="text",fontface="bold",color="black",x=40000,y=0.000220,label="μ+3σ=20070.71+5621.76",size=7)+
labs(x="frequency",y = "density")+
theme_bw()+
theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
theme(text=element_text(size=16,family="Arial",face="bold"))

#绘制braZ的正态密度分布直方图
cut -f 2 braZ_whole_number.tsv >braZ_density_ggplot.tsv
library(ggplot2)
library(showtext)#ggplot作图显示字体错误时需要用showtext包输入需要的字体
font_add('Arial','/Library/Fonts/Arial.ttf')
showtext_auto()
library(patchwork)
ratio1<-read.table("braZ_density_ggplot.tsv",header=T)
set.seed(1000)
df1 <- data.frame(ratio1)
#geom_density会为数据集提供一条密度曲线，而不是正态分布,而stat_function(fun = dnorm, args = list(mean = mean(df$ratio), sd = sd(df$ratio)))提供正态分布密度曲线
#添加均值线和三个标准差线
data_mean1=mean(df1$ratio) #20432.97
data_sd1=sd(df1$ratio) #1821.934
three_sd_plus1=mean(df1$ratio)+data_sd1*3  # 25898.77
three_sd_minus1=mean(df1$ratio)-data_sd1*3 #14967.16

p2<-ggplot(df1, aes(x = ratio)) +geom_histogram(aes(y =..density..),breaks = seq(min(df1$ratio), max(df1$ratio), length =100),colour ="#C6DBEF", fill ="#08519C", size = 0.1,alpha=0.8) + #使用density代替y轴
stat_function(fun = dnorm, args = list(mean = data_mean1, sd = data_sd1,color ="black", size = 1))+  #直方图上加密度曲线，也可以用geom_density()
geom_vline(aes(xintercept=data_mean1), color="red", linetype="dashed", size=1)+ #添加均值线和两个标准差线
geom_vline(aes(xintercept=three_sd_plus1), color="red", linetype="dashed", size=1)+
annotate(geom="text",fontface="bold",color="black",x=three_sd_plus,y=0.00015,label="25898.77",size=5)+   
 #通过annotate函数添加注释，将geom变量设置为"text"和"rect"分别代表文字与矩形，并且可以调整位置，大小颜色
annotate(geom="text",fontface="bold",color="black",x=40000,y=0.000225,label="μ+3σ=20432.97+5465.803",size=7)+
labs(x="frequency",y = "density")+
theme_bw()+
theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
theme(text=element_text(size=16,family="Arial",face="bold"))
p3<-p1+p2
ggsave(file="braB_braZ_density_picture.pdf",plot=p3,width=20,height=8)
```

# 8.eggnog在线注释
* 使用eggnog,即使用gene对应的序列去数据库里查功能，非PA号查询功能
## 8.1绘制三三组合的两个sigma和三个sigma的braB和braZ的韦恩图
```r
#three sigma   braB 25693  braZ 25899
#two sigma     braB  23819   braZ  24077
#分别提取和braB,braZ共表达的高频基因
cd /mnt/d/WGCNA/rma_WGCNA/condition_three_combine/result
cat PA1590.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | tsv-filter --ge 2:25693 | cut -f 1 >/mnt/d/braB_three_sigma_neighber_high_fre.tsv
cat PA1971.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | tsv-filter --ge 2:25899 | cut -f 1 >/mnt/d/braZ_three_sigma_neighber_high_fre.tsv

cat PA1590.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | tsv-filter --ge 2:23819 | cut -f 1 >/mnt/d/braB_two_sigma_neighber_high_fre.tsv
cat PA1971.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | tsv-filter --ge 2:24077 | cut -f 1 >/mnt/d/braZ_two_sigma_neighber_high_fre.tsv

#将braB_neighber_high_fre.tsv的gene id转换locus_tag id
cd /mnt/d/WGCNA/venn
cat  braB_three_sigma_neighber_high_fre.tsv | grep "PA" | perl -alne '$_=~/PA(.*?)\_/;$name=$1;print"PA$name"' >braB_three_Pfam_id.tsv
cat  braZ_three_sigma_neighber_high_fre.tsv | grep "PA" | perl -alne '$_=~/PA(.*?)\_/;$name=$1;print"PA$name"' >braZ_three_Pfam_id.tsv

cat  braB_two_sigma_neighber_high_fre.tsv | grep "PA" | perl -alne '$_=~/PA(.*?)\_/;$name=$1;print"PA$name"' >braB_two_Pfam_id.tsv
cat  braZ_two_sigma_neighber_high_fre.tsv | grep "PA" | perl -alne '$_=~/PA(.*?)\_/;$name=$1;print"PA$name"' >braZ_two_Pfam_id.tsv

#绘制3个sigma的braB和braZ的韦恩图
cd /mnt/d/WGCNA/venn
plotr venn braB_three_Pfam_id.tsv braZ_three_Pfam_id.tsv

```

## 8.2下载注释文件
```bash

cd /mnt/d/WGCNA/eggnog/
#braB_three_Pfam_id.tsv braZ_three_Pfam_id.tsv 不用更改
#提取2个sigma中共享的,braB独有的,braZ独有的
cat braB_two_Pfam_id.tsv | grep -f braZ_two_Pfam_id.tsv  >two_sigma_share_pfam_id.tsv
cat braB_two_Pfam_id.tsv | grep -v -f two_sigma_share_pfam_id.tsv >braB_two_special_pfam_id.tsv
cat braZ_two_Pfam_id.tsv | grep -v -f two_sigma_share_pfam_id.tsv >braZ_two_special_pfam_id.tsv
#提取PAO1的gbk文件中相应locus_tag id对应的氨基酸序列 
python fetch_faa.py PAO1.gbff PAO1.fa
perl -alne '$_=~s/\,(.*)//g;print"$_";' PAO1.fa >PAO1_replace.fa
#下载eggnog的注释文件
#删除注释文件带#的注释行
grep -v "##" out.emapper.annotations >PAO1_eggnog.tsv
sed -i 's/#//g' PAO1_eggnog.tsv
#利用R语言将注释结果整理成clusterProfiler包的enricher函数需要的输入格式
```

## 8.3改变注释文件的格式
* 读取注释文件并处理，得到背景文件。clusterprofiler要求的背景文件term2gene需要两列，第一列为goid,第二列为geneid。从eggnog下载处理后的文件，每个基因对应多个goid，要对其进行处理
```r
library(stringr)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
egg<-read.table("PAO1_eggnog.tsv",sep="\t",header=T,quote = "")
#####处理GO注释文件
#提取id列
gene_ids <- egg$query
#有的基因没有注释到会显示为"".需要使用逻辑值索引去除未注释到的
eggnog_lines_with_go <- egg$GOs!= ""
#将一个geneid对应多个goid的宽数据格式转换为长数据格式
eggnog_lines_with_go
eggnog_annoations_go <- str_split(egg[eggnog_lines_with_go,]$GOs, ",")
gene_to_go <- data.frame(gene = rep(gene_ids[eggnog_lines_with_go],
                                   times = sapply(eggnog_annoations_go, length)),
                         term = unlist(eggnog_annoations_go))
head(gene_to_go)
#交换geneid和goid
term2gene1<-gene_to_go[, c(2,1)]
#进一步处理，将间接注释补全，将GOid翻译为GOterm和GOontology
#将直接注释补充为间接注释
term2gene<-buildGOmap(term2gene1)
head(term2gene1)
#将GOid转换为GOterm
go2term<-go2term(term2gene$GO)
head(go2term)
#将GOid转换为GoOnt
go2ont<-go2ont(term2gene$GO)
head(go2ont)
#每个基因都会对应有一个或多个GO term（也就是GO功能）。Gene Ontology分为分子功能，生物过程和细胞组成三个部分。蛋白质或者基因可以通过ID对应或者序列注释的方法找到与之对应的GO号，而GO号可对应到Term，即功能类别或者细胞定位。这也是GO富集的一个基础。



#####处理KEGG注释文件
library(stringr)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
egg<-read.table("PAO1_eggnog.tsv",sep="\t",header=T,quote = "")
gene_ids <- egg$query
eggnog_lines_with_ko <- egg$KEGG_ko!= ""
#将一个geneid对应多个koid的宽数据格式转换为长数据格式
eggnog_annoations_ko <- str_split(egg[eggnog_lines_with_ko,]$KEGG_ko, ",")
gene_to_ko <- data.frame(gene = rep(gene_ids[eggnog_lines_with_ko],
                                   times = sapply(eggnog_annoations_ko, length)),
                         term = unlist(eggnog_annoations_ko))
#gene_to_ko重命名
colnames(gene_to_ko)<-c("GID","Ko")
gene_to_ko$Ko <-gsub("ko:","",gene_to_ko$Ko)
#kegg_info.RData这个文件里有pathway2name这个对象
#下载kegg和pathway相关信息
#下载json文件  https://www.genome.jp/kegg-bin/get_htext?ko00001
library(jsonlite)
library(tibble)
library(stringr)
library(dplyr)
setwd("D:/")
update_kegg <- function(json = "ko00001.json") {
pathway2name <- tibble(Pathway = character(), Name = character())
ko2pathway <- tibble(Ko = character(), Pathway = character())
kegg <- fromJSON(json)
for (a in seq_along(kegg[["children"]][["children"]])) {
A <- kegg[["children"]][["name"]][[a]]
for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
kos <- str_match(kos_info, "K[0-9]*")[,1]
ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))}}}
save(pathway2name, ko2pathway, file = "kegg_info.RData")}
#调用函数后在本地创建kegg_info.RData文件，以后只需要载入即可
update_kegg()
#载入文件
load(file = "kegg_info.RData")
pathway2gene1 <- gene_to_ko %>% left_join(ko2pathway, by = "Ko") %>% dplyr::select(GID, Pathway) %>%na.omit()
#交换geneid和koid
pathway2gene<-pathway2gene1[, c(2,1)]
pathway2name  
head(pathway2gene)    
head(pathway2name)                                    
```
## 8.4GO富集分析
```r
#pvalue(pval): 统计学差异显著性检验指标。
#qvalue(p-adjusted): 校正后的p值（qvalue=padj=FDR=Corrected p-Value=p-adjusted），是对p值进行了多重假设检验，能更好地控制假阳性率。校正后的p值不同的几种表现形式，都是基于BH的方法进行多重假设检验得到的。校正后的p值不同的展现形式是因为不同的分析软件产生的。

#读取基因列表，感兴趣的基因存为一列
gene_ins<-read.table("braZ_three_Pfam_id.tsv",sep="\t",header=F)
gene_list<-gene_ins$V1
##对BP MF CC分类
df<-go2term
df$ont<-go2ont$Ontology
df<-split(df,with(df,ont))
#在GOterm后添加ont信息
GO_BP<-enricher(gene=gene_list,pvalueCutoff =0.05,qvalueCutoff=1,pAdjustMethod = "none",TERM2GENE =term2gene,TERM2NAME=df$BP,minGSSize=0,maxGSSize=400)
GO_CC<-enricher(gene=gene_list,pvalueCutoff =0.05,qvalueCutoff=1,pAdjustMethod = "none",TERM2GENE =term2gene,TERM2NAME=df$CC,minGSSize=0,maxGSSize=400)
GO_MF<-enricher(gene=gene_list,pvalueCutoff =0.05,qvalueCutoff=1,pAdjustMethod = "none",TERM2GENE =term2gene,TERM2NAME=df$MF,minGSSize=0,maxGSSize=400)
p_BP<-dotplot(GO_BP)
P_CC<-dotplot(GO_CC)
p_MF<-dotplot(GO_MF)
plotc <- p_BP/P_CC/p_MF
pdf(file = "braB_BP_CC_MF.pdf",width =8,height = 16)
print(plotc)
dev.off()
#使用ggplot来画气泡图
#我们再用ggplot画气泡图，结果看起来和clusterprofiler画的差不多，不一样的地方在于clusterproliler的x轴默认的gene ratio,而ggplot我们用的是count
GO<-enricher(gene=gene_list,pvalueCutoff =0.05,qvalueCutoff=1,pAdjustMethod = "none",TERM2GENE =term2gene,TERM2NAME=df,minGSSize=0,maxGSSize=400)
head(GO@result)
GO_input<-GO@result[1:20,]
GO_input<-na.omit(GO_input)
df1<-go2ont(GO_input$ID)
GO_input$ont<-df1$Ontology
#使画出go term的顺序与输入的一致
GO_input$Description<-factor(GO_input$Description,levels = rev(GO_input$Description))
#reorder使纵轴按照go term 和count排序
library(ggplot2)
ggplot(GO_input,aes(x = Count ,y =reorder(Description,Count),shape=ont))+ 
  geom_point(aes(size=Count,color=p.adjust))+
  scale_colour_gradient(low="blue",high="red")+
  labs(
       color=expression(p.adjust),
       size=" Count Number",
       x="Gene Count"
      )+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = rel(1.8)),
    axis.title.x = element_text(size=rel(1.8)),
    axis.title.y = element_blank()
  )+ scale_size(range=c(5, 10))
​

##画图
GO_1<-barplot(GO,showCategory=10)
pdf(file = "braB_two_sigma_GO_1.pdf",width =12,height = 9)
print(GO_1)
dev.off()
GO_2<-dotplot(GO,showCategory=10)
pdf(file = "braB_two_sigma_GO_2.pdf",width =12,height = 9)
print(GO_2)
dev.off()
```

## 8.5KEGG富集分析
```r
gene_ins<-read.table("braB_Pfam_id.tsv",sep="\t",header=T)
gene_list<-gene_ins$id
KEGG<-enricher(gene=gene_list,TERM2GENE = pathway2gene,TERM2NAME = pathway2name,pvalueCutoff =1,qvalueCutoff=1,pAdjustMethod = "BH",minGSSize=0,maxGSSize=400)
KEGG_1<-barplot(KEGG,showCategory=20)
pdf(file = "braB_two_sigma_KEGG_1.pdf",width =12,height = 9)
print(KEGG_1)
dev.off()
KEGG_2<-dotplot(KEGG,showCategory=20)
pdf(file = "braB_two_sigma_KEGG_2.pdf",width =12,height = 9)
print(KEGG_2)
dev.off()
```


# 9.绘制三三组合的braZ和braB的基因表达谱
## 9.1绘制braB的基因表达谱
```bash
#数量是56663
cd /mnt/d/WGCNA/rma_WGCNA/condition_three_combine/result/whole
for filename in *.soft.nodes.txt
do
base=$(basename $filename .soft.nodes.txt)
echo $base >>three_condition_combine_whole_expression.tsv
done

#提取三三组合braB的表达量 56663
cd /mnt/d/WGCNA/rma_WGCNA/condition_three_combine/img
for f in $(cat three_condition_combine_whole_expression.tsv)
do
gene=$(cat $f.soft | grep "braB")
echo -e "$f\t$gene" >>three_condition_combine_braB_expression.tsv
done
#计算每一行的平均基因表达量,注意不同组合的样本量不同，所以会存在缺失值 (已经转置该文件并用excel计算前几个样本验证过计算无错误)
cat three_condition_combine_braB_expression.tsv | perl -MStatistics::Descriptive -alne '($gene,$exp)=(split/\t/,$_,3)[0,2];@data=(split/\t/,$exp);$stat = Statistics::Descriptive::Full->new();$stat->add_data(@data);$mean = $stat->mean();$mean=sprintf"%.2f",$mean;print"$gene\t$mean";' >three_condition_combine_braB_mean_expression.tsv

#查看56663个组合中braB的基因表达量最大值和最小值
cat three_condition_combine_braB_mean_expression.tsv |cut -f 2 | sort -n | uniq | head  #最小值3.82
cat three_condition_combine_braB_mean_expression.tsv |cut -f 2 | sort -n | uniq | tail  #最大值10.2

#使用mclust绘制braB基因表达量的mclust模型图
library(mclust)
ratio<-read.table("three_condition_combine_braB_mean_expression.tsv",header=F)
#mclust分析
BIC <- mclustBIC(ratio$V2)
ICL <- mclustICL(ratio$V2)
write.table(BIC, "BIC.tmp.tsv", row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(ICL, "ICL.tmp.tsv", row.names = FALSE, col.names = TRUE, sep = "\t")
#将BIC和ICL合并
cat BIC.tmp.tsv | perl -n -e 'chomp;
@a = split/\t/, $_;if($_ =~ /"/){$i = 1;next;}else{
print "$i\t$a[0]\tBICE\n"; print "$i\t$a[1]\tBICV\n";
$i++;}' | tsv-filter --str-ne 2:NA > BIC2.tmp.tsv
cat ICL.tmp.tsv | perl -n -e 'chomp;
@a = split/\t/, $_;if($_ =~ /"/){$i = 1;next;}else{
print "$i\t$a[0]\tICLE\n"; print "$i\t$a[1]\tICLV\n";
$i++;}' | tsv-filter --str-ne 2:NA > ICL2.tmp.tsv
cat ICL2.tmp.tsv BIC2.tmp.tsv | sed '1inum\tgrade\tmethod' > braB_mclust_expr_model.tsv
rm *.tmp.*
#使用BIC和ICL画图
library(readr)
library(ggplot2)
num <- read_tsv("braB_mclust_expr_model.tsv", show_col_types = FALSE)
p <- ggplot()+
geom_line(data = num, aes(x = num, y = grade, group = method, color = method)) +
scale_x_continuous(breaks = seq(0,9,1))
ggsave(p, file = "braB_expr_mclust.pdf", width = 5, height = 4)


#绘制braB的基因表达密度图
data<-read.table("three_condition_combine_braB_mean_expression.tsv",header=F)
braB<-data$V2
dens <- densityMclust(braB)
summary(dens, parameters = TRUE)  # model with 8 components
br <- seq(min(braB), max(braB), length =566)
x <- seq(min(braB)-diff(range(braB))/100,max(braB)+diff(range(braB))/100, length = 566)
cdens <-predict(dens, x, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*dens$parameters$pro))
plot(dens, what = "density", data = data$V2, breaks = br)
matplot(x, cdens, type = "l", lwd = 1, add = TRUE, lty = 1)
```

## 9.2绘制braZ的基因表达谱
```r
#提取三三组合braZ的表达量
cd /mnt/d/WGCNA/rma_WGCNA/condition_three_combine/img
for f in $(cat three_condition_combine_whole_expression.tsv)
do
gene=$(cat $f.soft | grep "braZ")
echo -e "$f\t$gene" >>three_condition_combine_braZ_expression.tsv
done
#计算每一行的平均基因表达量
cat three_condition_combine_braZ_expression.tsv | perl -MStatistics::Descriptive -alne '($gene,$exp)=(split/\t/,$_,3)[0,2];@data=(split/\t/,$exp);$stat = Statistics::Descriptive::Full->new();$stat->add_data(@data);$mean = $stat->mean();$mean=sprintf"%.2f",$mean;print"$gene\t$mean";' >three_condition_combine_braZ_mean_expression.tsv


#查看56663个组合中braB的基因表达量最大值和最小值
cat three_condition_combine_braZ_mean_expression.tsv |cut -f 2 | sort -n | uniq | head  #最小值1.65
cat three_condition_combine_braZ_mean_expression.tsv |cut -f 2 | sort -n | uniq | tail  #最大值7.08

#使用mclust绘制braB基因表达量的mclust模型图
library(mclust)
ratio<-read.table("three_condition_combine_braZ_mean_expression.tsv",header=F)
#mclust分析
BIC <- mclustBIC(ratio$V2)
ICL <- mclustICL(ratio$V2)
write.table(BIC, "BIC.tmp.tsv", row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(ICL, "ICL.tmp.tsv", row.names = FALSE, col.names = TRUE, sep = "\t")
#将BIC和ICL合并
cat BIC.tmp.tsv | perl -n -e 'chomp;
@a = split/\t/, $_;if($_ =~ /"/){$i = 1;next;}else{
print "$i\t$a[0]\tBICE\n"; print "$i\t$a[1]\tBICV\n";
$i++;}' | tsv-filter --str-ne 2:NA > BIC2.tmp.tsv
cat ICL.tmp.tsv | perl -n -e 'chomp;
@a = split/\t/, $_;if($_ =~ /"/){$i = 1;next;}else{
print "$i\t$a[0]\tICLE\n"; print "$i\t$a[1]\tICLV\n";
$i++;}' | tsv-filter --str-ne 2:NA > ICL2.tmp.tsv
cat ICL2.tmp.tsv BIC2.tmp.tsv | sed '1inum\tgrade\tmethod' > braZ_mclust_expr_model.tsv
rm *.tmp.*
#使用BIC和ICL画图
library(readr)
library(ggplot2)
num <- read_tsv("braZ_mclust_expr_model.tsv", show_col_types = FALSE)
p <- ggplot()+
geom_line(data = num, aes(x = num, y = grade, group = method, color = method)) +
scale_x_continuous(breaks = seq(0,9,1))
ggsave(p, file = "braZ_expr_mclust.pdf", width = 5, height = 4)


#绘制braB的基因表达密度图
data<-read.table("three_condition_combine_braZ_mean_expression.tsv",header=F)
braZ<-data$V2
dens <- densityMclust(braZ)
summary(dens, parameters = TRUE)  #model with 5 components
br <- seq(min(braZ), max(braZ), length =566)
x <- seq(min(braZ)-diff(range(braZ))/100,max(data$V2)+diff(range(braZ))/100, length = 566)
cdens <-predict(dens, x, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*dens$parameters$pro))
plot(dens, what = "density", data = braZ, breaks = br)
matplot(x, cdens, type = "l", lwd = 1, add = TRUE, lty = 1)
```

# 10.绘制40个条件的braZ和braB的基因表达谱
## 10.1绘制braB的基因表达谱
```bash
#提取四十个条件的braB的表达量 
cd /mnt/d/WGCNA/rma_WGCNA/condition_three_combine
for filename  in *.condition.txt
do
base=$(basename $filename .condition.txt)
gene=$(cat $base.condition.txt | grep "braB")
echo -e "$base\t$gene" >>braB_forty_condition_expression.tsv
done
#计算每一行的平均基因表达量,注意不同组合的样本量不同，所以会存在缺失值 (已经转置该文件并用excel计算前几个样本验证过计算无错误)
cat braB_forty_condition_expression.tsv | perl -MStatistics::Descriptive -alne '($gene,$exp)=(split/\t/,$_,3)[0,2];@data=(split/\t/,$exp);$stat = Statistics::Descriptive::Full->new();$stat->add_data(@data);$mean = $stat->mean();$mean=sprintf"%.2f",$mean;print"$gene\t$mean";' >braB_forty_condition_mean_expression.tsv

#查看56663个组合中braB的基因表达量最大值和最小值
cat braB_forty_condition_mean_expression.tsv |cut -f 2 | sort -n | uniq | head  #最小值3.82
cat braB_forty_condition_mean_expression.tsv |cut -f 2 | sort -n | uniq | tail  #最大值10.2

#使用mclust绘制braB基因表达量的mclust模型图
library(mclust)
ratio<-read.table("braB_forty_condition_mean_expression.tsv",header=F)
#mclust分析
BIC <- mclustBIC(ratio$V2)
ICL <- mclustICL(ratio$V2)
write.table(BIC, "BIC.tmp.tsv", row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(ICL, "ICL.tmp.tsv", row.names = FALSE, col.names = TRUE, sep = "\t")
#将BIC和ICL合并
cat BIC.tmp.tsv | perl -n -e 'chomp;
@a = split/\t/, $_;if($_ =~ /"/){$i = 1;next;}else{
print "$i\t$a[0]\tBICE\n"; print "$i\t$a[1]\tBICV\n";
$i++;}' | tsv-filter --str-ne 2:NA > BIC2.tmp.tsv
cat ICL.tmp.tsv | perl -n -e 'chomp;
@a = split/\t/, $_;if($_ =~ /"/){$i = 1;next;}else{
print "$i\t$a[0]\tICLE\n"; print "$i\t$a[1]\tICLV\n";
$i++;}' | tsv-filter --str-ne 2:NA > ICL2.tmp.tsv
cat ICL2.tmp.tsv BIC2.tmp.tsv | sed '1inum\tgrade\tmethod' > braB_mclust_Forty_model.tsv
rm *.tmp.*
#使用BIC和ICL画图
library(readr)
library(ggplot2)
num <- read_tsv("braB_mclust_Forty_model.tsv", show_col_types = FALSE)
p <- ggplot()+
geom_line(data = num, aes(x = num, y = grade, group = method, color = method)) +
scale_x_continuous(breaks = seq(0,9,1))
ggsave(p, file = "braB_mclust_Forty_model.pdf", width = 5, height = 4)


#绘制braB的基因表达密度图
data<-read.table("braB_forty_condition_mean_expression.tsv",header=F)
braB<-data$V2
dens <- densityMclust(braB)
summary(dens, parameters = TRUE)  # model with 1 component
br <- seq(min(braB), max(braB), length =40)
x <- seq(min(braB)-diff(range(braB))/1,max(braB)+diff(range(braB))/1, length = 40)
cdens <-predict(dens, x, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*dens$parameters$pro))
plot(dens, what = "density", data = data$V2, breaks = br)
matplot(x, cdens, type = "l", lwd = 1, add = TRUE, lty = 1)
```

## 10.2绘制braZ的基因表达谱
```bash
#提取四十个条件的braB的表达量 
cd /mnt/d/WGCNA/rma_WGCNA/condition_three_combine
for filename  in *.condition.txt
do
base=$(basename $filename .condition.txt)
gene=$(cat $base.condition.txt | grep "braZ")
echo -e "$base\t$gene" >>braZ_forty_condition_expression.tsv
done
#计算每一行的平均基因表达量,注意不同组合的样本量不同，所以会存在缺失值 (已经转置该文件并用excel计算前几个样本验证过计算无错误)
cat braZ_forty_condition_expression.tsv | perl -MStatistics::Descriptive -alne '($gene,$exp)=(split/\t/,$_,3)[0,2];@data=(split/\t/,$exp);$stat = Statistics::Descriptive::Full->new();$stat->add_data(@data);$mean = $stat->mean();$mean=sprintf"%.2f",$mean;print"$gene\t$mean";' >braZ_forty_condition_mean_expression.tsv

#查看56663个组合中braB的基因表达量最大值和最小值
cat braZ_forty_condition_mean_expression.tsv |cut -f 2 | sort -n | uniq | head  #最小值1.65
cat braZ_forty_condition_mean_expression.tsv |cut -f 2 | sort -n | uniq | tail  #最大值7.08

#使用mclust绘制braB基因表达量的mclust模型图
library(mclust)
ratio<-read.table("braZ_forty_condition_mean_expression.tsv",header=F)
#mclust分析
BIC <- mclustBIC(ratio$V2)
ICL <- mclustICL(ratio$V2)
write.table(BIC, "BIC.tmp.tsv", row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(ICL, "ICL.tmp.tsv", row.names = FALSE, col.names = TRUE, sep = "\t")
#将BIC和ICL合并
cat BIC.tmp.tsv | perl -n -e 'chomp;
@a = split/\t/, $_;if($_ =~ /"/){$i = 1;next;}else{
print "$i\t$a[0]\tBICE\n"; print "$i\t$a[1]\tBICV\n";
$i++;}' | tsv-filter --str-ne 2:NA > BIC2.tmp.tsv
cat ICL.tmp.tsv | perl -n -e 'chomp;
@a = split/\t/, $_;if($_ =~ /"/){$i = 1;next;}else{
print "$i\t$a[0]\tICLE\n"; print "$i\t$a[1]\tICLV\n";
$i++;}' | tsv-filter --str-ne 2:NA > ICL2.tmp.tsv
cat ICL2.tmp.tsv BIC2.tmp.tsv | sed '1inum\tgrade\tmethod' > braZ_mclust_Forty_model.tsv
rm *.tmp.*
#使用BIC和ICL画图
library(readr)
library(ggplot2)
num <- read_tsv("braZ_mclust_Forty_model.tsv", show_col_types = FALSE)
p <- ggplot()+
geom_line(data = num, aes(x = num, y = grade, group = method, color = method)) +
scale_x_continuous(breaks = seq(0,9,1))
ggsave(p, file = "braZ_mclust_Forty_model.pdf", width = 5, height = 4)

#绘制braZ的基因表达密度图
data<-read.table("braZ_forty_condition_mean_expression.tsv",header=F)
braZ<-data$V2
dens <- densityMclust(braZ)
summary(dens, parameters = TRUE)  # model with 2 component
br <- seq(min(braZ), max(braZ), length =40)
x <- seq(min(braZ)-diff(range(braZ))/1,max(braZ)+diff(range(braZ))/1, length = 40)
cdens <-predict(dens, x, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*dens$parameters$pro))
plot(dens, what = "density", data = data$V2, breaks = br)
matplot(x, cdens, type = "l", lwd = 1, add = TRUE, lty = 1)
```

## 10.3以14,13,13分为低中高组合
```r
cat braB_forty_condition_mean_expression.tsv |sort -k 2 -n | head -n 14 | tail -n 1  #6.90
cat braB_forty_condition_mean_expression.tsv |sort -k 2 -n | head -n 27 | tail -n 1  #7.56
cat braB_forty_condition_mean_expression.tsv |sort -k 2 -n | head -n 40 | tail -n 1  #10.20

cat braZ_forty_condition_mean_expression.tsv |sort -k 2 -n | head -n 14 | tail -n 1  #4.85
cat braZ_forty_condition_mean_expression.tsv |sort -k 2 -n | head -n 27 | tail -n 1  #6.12
cat braZ_forty_condition_mean_expression.tsv |sort -k 2 -n | head -n 40 | tail -n 1  #7.08

cat braB_forty_condition_mean_expression.tsv |sort -k 2 -n | head -n 14 | cut -f 1 >braB_low.tsv
cat braB_forty_condition_mean_expression.tsv |sort -k 2 -n | head -n 27 | tail -n 13 |cut -f 1 >braB_mid.tsv
cat braB_forty_condition_mean_expression.tsv |sort -k 2 -n | head -n 40 | tail -n 13 |cut -f 1 >braB_high.tsv

cat braZ_forty_condition_mean_expression.tsv |sort -k 2 -n | head -n 14 | cut -f 1 >braZ_low.tsv
cat braZ_forty_condition_mean_expression.tsv |sort -k 2 -n | head -n 27 | tail -n 13 |cut -f 1 >braZ_mid.tsv
cat braZ_forty_condition_mean_expression.tsv |sort -k 2 -n | head -n 40 | tail -n 13 |cut -f 1 >braZ_high.tsv

#绘制韦恩图
library(VennDiagram)
library(RColorBrewer)
set1<-read.table("braB_low.tsv",header=F)
set2<-read.table("braB_mid.tsv",header=F)
set3<-read.table("braB_high.tsv",header=F)
set4<-read.table("braZ_low.tsv",header=F)
set5<-read.table("braZ_mid.tsv",header=F)
set6<-read.table("braZ_high.tsv",header=F)
venn.diagram(x=list(braB_low=set1$V1,braZ_low=set4$V1),fill=c('red','blue'),filename='braB_low_braZ_low.tiff')
venn.diagram(x=list(braB_low=set1$V1,braZ_mid=set5$V1),fill=c('red','blue'),filename='braB_low_braZ_mid.tiff')
venn.diagram(x=list(braB_low=set1$V1,braZ_high=set6$V1),fill=c('red','blue'),filename='braB_low_braZ_high.tiff')

venn.diagram(x=list(braB_high=set3$V1,braZ_low=set4$V1),fill=c('red','blue'),filename='braB_high_braZ_low.tiff')
venn.diagram(x=list(braB_high=set3$V1,braZ_mid=set5$V1),fill=c('red','blue'),filename='braB_high_braZ_mid.tiff')
venn.diagram(x=list(braB_high=set3$V1,braZ_high=set6$V1),fill=c('red','blue'),filename='braB_high_braZ_high.tiff')
```

## 10.4plotr绘制40个条件的braZ和braB的密度分布图
```r
cd /mnt/d/WGCNA/expression/forty_condition_expr
#准备braB画图文件
cat braB_forty_condition_mean_expression.tsv | perl -alne '$exp=(split/\t/,$_,2)[1];print"A\t$exp"' | sed '1iconditon\texp' >braB_density_picture.tsv
#绘制频数图，非频率图
plotr hist braB_density_picture.tsv --xl braB --bins 20 --col 2


#准备braZ画图文件
cat braZ_forty_condition_mean_expression.tsv | perl -alne '$exp=(split/\t/,$_,2)[1];print"A\t$exp"' | sed '1iconditon\texp' >braZ_density_picture.tsv
#绘制频数图，非频率图
plotr hist braZ_density_picture.tsv --xl braZ --bins 20 --col 2
```

# 11.绘制40个条件中braB和braZ的基因表达量热图

## 11.1绘制40个条件中braB和braZ的平均基因表达量热图，不聚类
```r
#提取40个条件里braB和braZ的平均基因表达量
braB_forty_condition_mean_expression.tsv
braZ_forty_condition_mean_expression.tsv
cd /mnt/d/WGCNA/expression/raw_expr_heatmap
cat braB_forty_condition_mean_expression.tsv | tsv-join -d 1 -f  braZ_forty_condition_mean_expression.tsv -k 1 --append-fields 2 | tsv-join -d 1 -f condition_detail.tsv -k 2 --append-fields 1 >detail_braB_braZ_men_expr.tsv
tsv-select -f 4,2,3 detail_braB_braZ_men_expr.tsv | sed '1icondition\tbraB\tbraZ' >detail_braB_braZ_mean_heatmap.tsv
#绘制热图
library(pheatmap)
data<-read.table("detail_braB_braZ_mean_heatmap.tsv",header=T,row.names=1, sep="\t", quote="")
#将数据转换为数据框格式
data<-as.data.frame(data)
#pheatmap对数据进行z-score仅用参数scale="row"即对行进行标准化
#热图对行不进行聚类，对列不进行聚类
pheatmap(data,display_numbers = TRUE,scale="column",cluster_row = FALSE,cluster_cols = FALSE,border_color = "black",cellwidth = 30,cellheight = 10)

```

## 11.2绘制40个条件中braB和braZ和相关基因的平均基因表达量热图，对条件进行聚类，基因不聚类,对基因标准化
```r
#提取以下50个基因的基因表达量
for filename  in *.condition.txt
do
base=$(basename $filename .condition.txt)
cat $base.condition.txt | grep -f braB_braZ_cor_gene.tsv >$base.braB_braZ_cor.condition.tsv
done

#在上述文件中增添条件名,然后合并(50个基因，40个条件，共2000行)
for filename in *.braB_braZ_cor.condition.tsv
do
base=$(basename $filename .braB_braZ_cor.condition.tsv)
for f in {1..50}
do
gene=$(cat $base.braB_braZ_cor.condition.tsv | head -n $f | tail -n 1)
echo -e "$base\t$gene" >>braB_braZ_cor_condition_expr.tsv
done
done

#计算每一行的平均基因表达量
cat braB_braZ_cor_condition_expr.tsv | perl -MStatistics::Descriptive -alne '($condi,$gene,$exp)=(split/\t/,$_,3)[0,1,2];@data=(split/\t/,$exp);$stat = Statistics::Descriptive::Full->new();$stat->add_data(@data);$mean = $stat->mean();$mean=sprintf"%.2f",$mean;print"$condi\t$gene\t$mean";' >braB_braZ_cor_condition_mean_expr.tsv
#拼接条件详细信息
cat braB_braZ_cor_condition_mean_expr.tsv | tsv-join -d 1 -f condition_detail.tsv -k 2 --append-fields 1 | tsv-select -f 4,2,3 | sed '1icondi\tgene\texpr' >braB_braZ_cor_detail_mean_expr.tsv

#手动去除第一个列名
#读入数据并且将长数据变为短数据
setwd("D:/WGCNA/expression/raw_expr_heatmap")
library(pheatmap)
library(tidyr) #gather宽转长  spread长转宽
data1<-read.table("braB_braZ_cor_detail_mean_expr.tsv",header=T)
#关键的参数是key和value，指定后就可以根据相应的key和value进行长到宽的转换
data1<-spread(data1,key="gene",value="expr")
#将数据矩阵化，去除sample这个单词
rownames(data1)=data1[,1]
data1<-data1[,-1]
data1<-as.matrix(data1)
#载入分类信息
class<-read.table("brab_braZ_cor_class.tsv",header=T)
#type<-class$V2
#top = HeatmapAnnotation(df=data.frame(type))
#聚类方法为average,只对基因聚类，分成两个簇,不展示行名,顶部分组信息
annotation1 = HeatmapAnnotation(df = data.frame(type = c(rep("braB",21), rep("braZ", 29))))
pheatmap(data1,scale="column",show_rownames = FALSE,cluster_row = FALSE,cluster_cols = TRUE,border_color = "black",cellwidth = 10,cellheight = 10,clustering_method = "average",cutree_cols =2, annotation_col = class)
```

## 11.3使两个热图的条件顺序一致
```r
setwd("D:/WGCNA/expression/raw_expr_heatmap")
library(pheatmap)
library(tidyr)
 library(cowplot)
data1<-read.table("braB_braZ_cor_detail_mean_expr.tsv",header=T)
#关键的参数是key和value，指定后就可以根据相应的key和value进行长到宽的转换
data1<-spread(data1,key="gene",value="expr")
write.table(data1,file="condition_order.tsv",sep="\t")
#将数据矩阵化，去除sample这个单词
rownames(data1)=data1[,1]
data1<-data1[,-1]
data1<-as.matrix(data1)
#载入分类信息
class<-read.table("brab_braZ_cor_class.tsv",header=T)
#聚类方法为average,只对基因聚类，分成两个簇,不展示行名,顶部分组信息
#annotation1 = HeatmapAnnotation(df = data.frame(type = c(rep("braB",21), rep("braZ", 29))))  
pheatmap(data1,scale="column",show_rownames = T,show_colnames = T,cluster_row = F,cluster_cols = TRUE,border_color = "black",cellwidth = 10,cellheight = 10,clustering_method = "average",cutree_cols =2, annotation_col = class)

#更改braB和braZ条件的顺序
cat condition_order.tsv | cut -f 2 | grep -v "PA0378_at" | sed 's/"//g' >braB_braZ_condition_order.tsv
cat braB_braZ_condition_order.tsv | tsv-join -d 1 -f detail_braB_braZ_men_expr.tsv -k 4  --append-fields 2,3 |sed '1icondition\tbraB\tbraZ' >detail_braB_braZ_mean_heatmap.tsv
#读入文件
data<-read.table("detail_braB_braZ_mean_heatmap.tsv",header=T,row.names=1, sep="\t", quote="")
#将数据转换为数据框格式
data<-as.data.frame(data)
#pheatmap对数据进行z-score仅用参数scale="row"即对行进行标准化
#热图对行不进行聚类，对列不进行聚类
pheatmap(data,display_numbers = TRUE,cluster_row = FALSE,cluster_cols = FALSE,border_color = "black",cellwidth = 30,cellheight = 10)
```

## 11.3将热图的样本条件改成样本号码
```r
#将样本条件改为数据集号braB_braZ_condition_order.tsv
sed -i '1icondi\tgene\tmean' braB_braZ_cor_condition_mean_expr.tsv
setwd("D:/WGCNA/expression/raw_expr_heatmap")
library(pheatmap)
library(tidyr)
 library(cowplot)
data1<-read.table("braB_braZ_cor_condition_mean_expr.tsv",header=T)
#关键的参数是key和value，指定后就可以根据相应的key和value进行长到宽的转换
data1<-spread(data1,key="gene",value="mean")
write.table(data1,file="condition_order.tsv",sep="\t")
#将数据矩阵化，去除sample这个单词
rownames(data1)=data1[,1]
data1<-data1[,-1]
data1<-as.matrix(data1)
#载入分类信息
type<-read.table("brab_braZ_cor_class.tsv",header=T)
#聚类方法为average,只对基因聚类，分成两个簇,不展示行名,顶部分组信息
#annotation1 = HeatmapAnnotation(df = data.frame(type = c(rep("braB",21), rep("braZ", 29))))  
y<-pheatmap(data1,scale="column",show_rownames = T,show_colnames = T,cluster_row = F,cluster_cols = TRUE,border_color = "black",cellwidth = 10,cellheight = 10,clustering_method = "average",cutree_cols =2, annotation_col =type, fontfamily= "Arial",fontface="bold",fontsize=12,color = colorRampPalette(colors = c("#4DAF4A","#FFFFFF","#386CB0"))(100))
save_pheatmap_pdf <- function(y, filename, width=7, height=7) {
   stopifnot(!missing(y))
   stopifnot(!missing(filename))
   pdf(filename, width=10, height=10)
   grid::grid.newpage()
   grid::grid.draw(y$gtable)
   dev.off()
}
save_pheatmap_pdf(y, "braB_braZ_cor_condi.pdf")


#更改braB和braZ条件的顺序
cat condition_order.tsv | cut -f 2 | grep -v "PA0378_at" | sed 's/"//g' >braB_braZ_condition_order.tsv
cat braB_braZ_condition_order.tsv | tsv-join -d 1 -f detail_braB_braZ_men_expr.tsv  -k 1  --append-fields 2,3 |sed '1icondition\tbraB\tbraZ' >detail_braB_braZ_mean_heatmap.tsv
#读入文件
data<-read.table("detail_braB_braZ_mean_heatmap.tsv",header=T,row.names=1, sep="\t", quote="")
#将数据转换为数据框格式
data<-as.data.frame(data)
#pheatmap对数据进行z-score仅用参数scale="row"即对行进行标准化
#热图对行不进行聚类，对列不进行聚类
x<-pheatmap(data,display_numbers = TRUE,cluster_row = FALSE,cluster_cols = FALSE,border_color = "black",cellwidth = 30,cellheight = 10, fontfamily= "Arial",fontface="bold",fontsize=12,color = colorRampPalette(colors = c("#FF7F00","#FFFFFF","#F781BF"))(100))
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
save_pheatmap_pdf(x, "braB_braZ_condi.pdf")





```