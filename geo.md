
# 1.手动绘制地理信息表格
* 全国共34个省市，选取的样本没有甘肃省的，没有云南省的，没有澳门的，最终只有31个地点
```bash
cd /mnt/d/arg_blast/moran
#拼接样本，基因，有无，类型，地点信息
cat whole_sample_state.tsv | tsv-join  -d 1 -f whole_type_site_sample.tsv --k 3 --append-fields 1,2 >whole_sample_gene_state_type_site.tsv #90396
#仅保留sul2基因的样本，类型和地点信息
cat whole_sample_gene_state_type_site.tsv | grep "sul2" | tsv-select -f 1,4,5 | sort -n | uniq >sul2_sample_state_site.tsv #729 #31
#该样本包含地点和类型信息
wc -l  whole_site_type.tsv #310

#提取sul2基因比对上簇1的样本
perl AAA_提取contig和簇的id.pl sul2_clinical_identity_100.tsv >sul2_merge.tsv
sed -i 's/ids/cluster/g' sul2_group.tsv
tsv-join -H --filter-file sul2_group.tsv --key-fields cluster --append-fields group sul2_merge.tsv >sul2_merge_clusterid.tsv
tsv-select -f 4,2 sul2_merge_clusterid.tsv | sed 's/^/cluster&/g' | tsv-select -f 2,1 | sed 's/$/&cluster/g' | grep "cluster1cluster" | cut -f 1 | sort -n | uniq  >sul2_cluster1_sample.tsv #292

#提取比对上的样本对应的类型和地点信息
cat sul2_sample_state_site.tsv   | grep -f sul2_cluster1_sample.tsv | perl -alne '($sample,$type,$site)=(split/\s/,$_)[0,1,2];$gene=sul2;$state=1;print"$sample\t$type\t$site\t$gene\t$state";' | wc -l      #292
cat sul2_sample_state_site.tsv   | grep -v -f  sul2_cluster1_sample.tsv | perl -alne '($sample,$type,$site)=(split/\s/,$_)[0,1,2];$gene=sul2;$state=0;print"$sample\t$type\t$site\t$gene\t$state";' | wc -l #437

#sul2比对上簇1的样本对应类型和地点的去重
cat sul2_sample_state_site.tsv   | grep -f sul2_cluster1_sample.tsv | perl -alne '($sample,$type,$site)=(split/\s/,$_)[0,1,2];$gene=sul2;$state=1;print"$sample\t$type\t$site\t$gene\t$state";'|cut -f 2,3,5  | sort -n | uniq  >sul2_cluster1_type_site_exists.tsv #80
cat whole_site_type.tsv | grep -v -f <(cut -f 1,2 sul2_cluster1_type_site_exists.tsv| tsv-select -f 2,1) | tsv-select -f 2,1|perl -alne '($type,$site)=(split/\s/,$_)[0,1];$state=0;print"$type\t$site\t$state";' >sul2_cluster1_type_site_absence.tsv  #230
#合并sul2比对簇1的存在和不存在情况
cat sul2_cluster1_type_site_exists.tsv sul2_cluster1_type_site_absence.tsv | tsv-select -f 2,1,3 | sed '1isite\ttype\tstate' >sul2_cluster1_type_site_whole.tsv
 mlr --itsv --ocsv cat sul2_cluster1_type_site_whole.tsv  >sul2_cluster1_type_site_whole.csv
```

# 2.计算莫兰指数
```bash
install.packages(c("spdep",  "rgdal"))
library(spdep)
library(rgdal)
#读入属性数据
setwd("D:/arg_blast/moran")
mydata<-read.csv("sul2_cluster1_type_site_whole.csv",header=T,stringsAsFactors = F,encoding =  "UTF-8")
head(mydata)

#准备权重数据.有三个细节,一是原始的shape文件包含了甘肃，云南和澳门需要去掉,二是shape文件中的省份排列顺序和属性数据不一致，需要调整为和属性数据已知，三是海南没有邻居，指定广东为邻居，同时需要将两者的空间邻接关系调整为对称.台湾没有邻居，指定福建为邻居
shp <-  readOGR("China_province_shp\\provinces.shp", stringsAsFactors = F, encoding =  "UTF-8")
shp <- subset(shp, !NAME %in%  c("甘肃", "澳门", "云南"))
row.names(shp) <- shp$NAME
shp<-shp[mydata[mydata$type=="landfill","site"],]
plot(shp)
#首先生成每个区域的邻接list
mynb<-poly2nb(shp,row.names=shp$NAME)
#2 regions with no links:
#台湾 海南
number.hainan <- which(attr(mynb,  "region.id") == "海南")
number.guangdong <- which(attr(mynb,  "region.id") == "广东")
mynb[[number.hainan]] <-  number.guangdong
number.taiwan <- which(attr(mynb,  "region.id") == "台湾")
number.fujian <- which(attr(mynb,  "region.id") == "福建")
mynb[[number.taiwan]] <-  number.fujian
plot(mynb)
# make neighbour relation symmetric
mynb <- make.sym.nb(mynb)
# see if nb object is really symmetric
is.symmetric.nb(mynb)
#[1] TRUE
#再将邻接list转换为矩阵形式
mylistw <- nb2listw(mynb)
## 可视化邻接关系
map_crd <-coordinates(shp)
par(mar=rep(0,4))
plot(mylistw,coords=map_crd,pch=19, cex=0.1, col="gray")


#计算莫兰指数，莫兰指数的检验是基于标准化后的莫兰指数统计量的分布来进行的。计算这一分布的方差时，有两种做法：基于随机排列模拟计算和假设数据为正态分布根据理论公式推算。moran.test函数默认通过随机模拟来进行显著性检验。如参数设置为randomisation=FALSE,则为正态分布假设。由于空间自相关一般为正相关，因此默认的检验方向为"greater".
x<-mydata[mydata$type=="landfill","state"]
moran.test(x, mylistw, randomisation =  FALSE, alternative = "two.sided")

#莫兰指数模拟
#随机模拟
set.seed(12345)
morperm <- moran.mc(x, mylistw, 9999)
morperm$res[10000]
# plot simulated Moran's index
morp <- morperm$res[1:9999]
md <- density(morp)
plot(md, main="Moran's I  Permutation Test", xlab = "Reference  Distribution", xlim = c(-0.4, 0.55), ylim =  c(0, 4), lwd = 2, col = 2)
hist(morp, freq = F, add = T)
abline(v = morperm$statistic, lwd = 2,  col = 4)

#莫兰散点图
#将属性数据标准化后，利用lag.listw函数计算其空间滞后值，可画出莫兰散点图。同时还可发现，莫兰指数实际上等于空间滞后值对其本身值（需要事先标准化）的回归系数。
production <- x <- scale(x)[,1]
moran.plot(production,mylistw, labels=T)
wx <- lag.listw(mylistw, x)
lm(wx ~ x)

#批量计算八种类型的莫兰指数
# Moran's index equals standardized  regression coefficient
moran.all <- tapply(mydata$state, mydata$type, moran.test, mylistw)
moran <- t(sapply(moran.all,  "[[", "estimate"))[,1]
p.value <- sapply(moran.all,  "[[", "p.value")
moran.all <- data.frame(moran,  p.value)
write.table(moran.all,"sul2_cluster1_moran.tsv",sep="\t",row.names=TRUE,col.names=TRUE,quote=TRUE)
```



```bash
shp <-  readOGR("China_province_shp\\provinces.shp", stringsAsFactors = F, encoding =  "UTF-8")
shp <- subset(shp, !NAME %in%  c("甘肃", "澳门", "云南"))
row.names(shp) <- shp$NAME
shp<-shp[mydata[mydata$type=="landfill","site"],]
plot(shp)
#首先生成每个区域的邻接list
mynb<-poly2nb(shp,row.names=shp$NAME)
#2 regions with no links:
#台湾 海南
number.hainan <- which(attr(mynb,  "region.id") == "海南")
number.guangdong <- which(attr(mynb,  "region.id") == "广东")
mynb[[number.hainan]] <-  number.guangdong
number.taiwan <- which(attr(mynb,  "region.id") == "台湾")
number.fujian <- which(attr(mynb,  "region.id") == "福建")
mynb[[number.taiwan]] <-  number.fujian
plot(mynb)
# make neighbour relation symmetric
mynb <- make.sym.nb(mynb)
# see if nb object is really symmetric
is.symmetric.nb(mynb)
#[1] TRUE
#再将邻接list转换为矩阵形式
mylistw <- nb2listw(mynb)
## 可视化邻接关系
map_crd <-coordinates(shp)
par(mar=rep(0,4))
plot(mylistw,coords=map_crd,pch=19, cex=0.1, col="gray")


#计算莫兰指数，莫兰指数的检验是基于标准化后的莫兰指数统计量的分布来进行的。计算这一分布的方差时，有两种做法：基于随机排列模拟计算和假设数据为正态分布根据理论公式推算。moran.test函数默认通过随机模拟来进行显著性检验。如参数设置为randomisation=FALSE,则为正态分布假设。由于空间自相关一般为正相关，因此默认的检验方向为"greater".
x<-mydata[mydata$type=="landfill","state"]
moran.test(x, mylistw)


#莫兰散点图
#将属性数据标准化后，利用lag.listw函数计算其空间滞后值，可画出莫兰散点图。同时还可发现，莫兰指数实际上等于空间滞后值对其本身值（需要事先标准化）的回归系数。
production <- x <- scale(x)[,1]
moran.plot(production,mylistw, labels=T)
wx <- lag.listw(mylistw, x)
lm(wx ~ x)

#批量计算八种类型的莫兰指数
moran.all <- tapply(mydata$state, mydata$type, moran.test, mylistw)
moran_value <- t(sapply(moran.all,  "[[", "estimate"))[,1]
p.value <- sapply(moran.all,  "[[", "p.value")
moran.all <- data.frame(moran_value,  p.value)
write.table(moran.all,"sul2_cluster1_moran.tsv",sep="\t",row.names=TRUE,col.names=TRUE,quote=TRUE)
```