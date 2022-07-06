

# 1.使用diamond比对nr库搜索braz的来源(超算上运行)
* 准备数据和软件
```bash
#准备PAO1菌株的braz蛋白 Pseudom_aeru_PAO1_NP_250661
echo "Pseudom_aeru_PAO1_NP_250661" >branched-chain/branched-chain_PAO1_braz.tsv
faops some PROTEINS/all.replace.fa branched-chain/branched-chain_PAO1_braz.tsv  branched-chain/branched-chain_PAO1_braz.fa
#下载diamond
wget https://github.com/bbuchfink/diamond/releases/download/v2.0.5/diamond-linux64.tar.gz
cd data/blast/fasta
tar -zxvf diamond-linux64.tar.gz
#查看帮助信息
./diamond --help


cd data/blast/fasta
#本地建立nr数据库并比对
gunzip nr.gz
mv nr nr.fa
bsub -q mpi -n 24 -J "DIA" ./diamond makedb --in nr.fa -d nr --threads 20
bsub -q mpi -n 24-J "DIA" ./diamond blastp -d nr.dmnd -q branched-chain_PAO1_braz.fa -o braZ_nr_result.tsv --id 30 --threads 20 -e 1e-5 --more-sensitive
#结果只有25行,可能建库时对序列id有要求
#本地建立braz数据库并比对
./diamond makedb --in branched-chain_PAO1_braz.fa -d braz
bsub -q serial -n 20 -J "DIA" ./diamond blastp -d braz.dmnd -q nr.fa -o braz_nr_reverse.result.tsv --id 30 --subject-cover 30 --threads 16 -e 1e-5 --more-sensitive

#分割文件后使用diamond比对nr库搜索braz的来源(超算上运行)
#按序列名分割
faops size nr.fa | cut -f 1 >nr_name.tsv
split -l 5000000 nr_name.tsv -d -a 3 nr_
#-l 按行分割
#-d 添加数字后缀
#-a 3 表示用3位数来顺序命名 后缀长度
#url_ 分割后文件的前缀
#按序列分割
faops split-about nr.fa 3000000000 ./
##提取braZ抓取的序列对应的菌株蛋白名(37986)（第一次diamond）
cat nr.fa | grep ">" >nr_protein_name.tsv
cat nr_protein_name.tsv | grep -Ff <(cut -f 1 braz_nr_reverse.result.tsv | sed  's/^/>&/g') >braz_nr_strain_protein.tsv
#替换菌株名和蛋白名
perl -alne '/\>(.*?)\s.*(\[.*\])/;$name=$1;$strain=$2;$strain=~s/\[//g;$strain=~s/\]//g;$strain=~s/\s/\_/g;$id=$strain."\_".$name;print"$name\t$id";' braz_nr_strain_protein.tsv >braz_nr_strain_protein.replace.tsv
sed -i 's/(//g' braz_nr_strain_protein.replace.tsv
sed -i 's/)//g' braz_nr_strain_protein.replace.tsv
#替换序列蛋白名为菌株_蛋白名
faops replace braz_nr_reverse.fa  braz_nr_strain_protein.replace.tsv  braz_nr_strain_protein.replace.fa

#不聚类
cat braz_nr_strain_protein.replace.fa PAO1.fa >braz_nr_PAO1_whole.fa #37988
#建立树
bsub -q mpi -n 24 -J "mus" mafft --retree 1 --maxiterate 0 braz_nr_PAO1_whole.fa >braz_nr_PAO1_whole_mafft.fa
sed -i 's/://g' braz_nr_PAO1_whole_mafft.fa
sed -i 's/,//g' braz_nr_PAO1_whole_mafft.fa
FastTree braz_nr_PAO1_whole_mafft.fa >braz_nr_PAO1_whole_mafft.newick

#在树上标出PAO1中braz和braB的位置
NP_250661 braz 蓝色
NP_250281 braB 红色
#画图
setwd("D:/")
library(dplyr)
library(ggtree)
library(ape)
library(tidytree)
library(treeio)
tree<-read.newick("braz_nr_PAO1_whole_mafft.newick")
data<-fortify(tree)
#分别给braz和braB所在簇共同上色
node1<-grep("NP_250281",tree$tip.label)
node1 #35103
node2 <- grep("NP_250661", tree$tip.label)
node2 #36933
#对类群进行高亮显示
tree<-groupOTU(tree, c(node1,node2))
ggtree(tree, aes(colour = group)) + 
  scale_color_manual(values=c( "black","red")) +
  theme(legend.position = "none")
# 最近的父节点
nodes<-c(node1,node2)
clade <- MRCA(tree, nodes) 
sub_tree<- tree_subset(tree, clade, levels_back = 0)
ggtree(sub_tree) 
#输出树文件
write.tree(sub_tree,file = "sub_tree.nwk")
```

# 2.使用cd-hit对抓取出来的序列进行聚类分簇，然后将分簇的序列建树
```bash
#cd-hit使用
-i：输入文件，蛋白序列，fasta格式
-o：输出文件，有两个，一个代表性序列文件，一个聚类文件
-c：聚类阈值，1.0代表100%一致性，0.9代表90%一致性，以此类推。比对中相同氨基酸的数量除以较短序列的全长
-n：两两序列进行序列比对时选择的 word size，具体选值参考下面
-d：0表示使用 fasta 标题中第一个空格前的字段作为序列名字
-M：设置内存，16000，16G
-t：设置线程数
#cd-hit聚类共5876个簇(设置簇内相似性为90%)
cd cd_hit
cd-hit -i braz_nr_strain_protein.replace.fa -c 0.9 -n 5  -d 0 -M 16000 -T 8  -o braz_nr_cdhit_identity90.fa
cat braz_nr_cdhit_identity90.fa PAO1.fa >braz_nr_PAO1_cdhit_identity90.fa
#建立树
bsub -q mpi -n 24 -J "mus" mafft --retree 1 --maxiterate 0 braz_nr_PAO1_cdhit_identity90.fa >braz_nr_PAO1_cdhit_identity90_mafft.fa
sed -i 's/://g' braz_nr_PAO1_cdhit_identity90_mafft.fa
FastTree braz_nr_PAO1_cdhit_identity90_mafft.fa >braz_nr_PAO1_cdhit_identity90_aln.newick

#cd-hit聚类共2224个簇(设置簇内相似性为70%)
cd cd_hit
cd-hit -i braz_nr_strain_protein.replace.fa -c 0.7 -n 5  -d 0 -M 16000 -T 8 -o braz_nr_cdhit_identity70.fa 
cat braz_nr_cdhit_identity70.fa PAO1.fa >braz_nr_PAO1_cdhit_identity70.fa
#建立树
bsub -q mpi -n 24 -J "mus" mafft --retree 1 --maxiterate 0 braz_nr_PAO1_cdhit_identity70.fa >braz_nr_PAO1_cdhit_identity70_mafft.fa
sed -i 's/://g' braz_nr_PAO1_cdhit_identity70_mafft.fa
FastTree braz_nr_PAO1_cdhit_identity70_mafft.fa >braz_nr_PAO1_cdhit_identity70_aln.newick

#在树上标出PAO1中braz和braB的位置
NP_250661 braz 蓝色
NP_250281 braB 红色
#画图
setwd("D:/")
library(dplyr)
library(ggtree)
library(ape)
library(tidytree)
library(treeio)
tree<-read.newick("braz_nr_PAO1_cdhit_identity90_aln.newick")
data<-fortify(tree)
node1<-grep("NP_250281",tree$tip.label)
node1 #3709
node2 <- grep("NP_250661", tree$tip.label)
node2 #3778
node3<-grep("GFR65790",tree$tip.label)
node3 #4114
#对类群进行高亮显示
tree<-groupOTU(tree, c(node1,node3))
ggtree(tree, aes(colour = group)) + 
  scale_color_manual(values=c( "black","red")) +
  theme(legend.position = "none")
# 最近的父节点
nodes<-c(node1,node3)
clade <- MRCA(tree, nodes) 
sub_tree<- tree_subset(tree, clade, levels_back = 0)
grep("NP_250281", sub_tree$tip.label) #25,68
ggtree(sub_tree) + geom_tiplab() + xlim(0, 5)
#输出树文件
write.tree(sub_tree,file = "identity90_sub_tree.nwk")
#提取树文件中的gene_id
perl nwk_geneid.pl -i identity90_sub_tree.nwk  -o identity90_sub_tree_geneid.tsv

```

# 3.重新使用braB和braZ和nr库的序列比对
```bash
#diamond建库
sed  's/\s/\//g' nr.fa >nr.repalce
bsub -q mpi -n 24 -J "DIA"  ./diamond makedb --in PAO1_braB_braZ.fa  -d PAO1 --threads 20
#比对
bsub -q mpi -n 24 -J "DIA" ./diamond blastp --db PAO1.dmnd --query nr.repalce.fa --out nr_braB_braZ_diamond.tsv  --id 30 --subject-cover 30 --threads 20 -e 1e-5 
#提取braB和braZ抓取出来的序列
wc -l nr_braB_braZ_diamond.tsv #72684
cut -f 1 nr_braB_braZ_diamond.tsv | sort -n | uniq  >nr_braB_braZ_diamond_uniq.tsv   #37914
#只提取蛋白名
cut -d '/' -f 1 nr_braB_braZ_diamond_uniq.tsv | sort -n | uniq   >nr_braB_braZ_protein_name.tsv #37914
#根据蛋白名提取序列
faops some nr.fa  nr_braB_braZ_protein_name.tsv  nr_braB_braZ_protein_name.fa
#找到蛋白名对应的菌株蛋白名
cat nr_braB_braZ_diamond_uniq.tsv | grep -Ff <(cut -f 1 nr_braB_braZ_protein_name.tsv ) | sed  's/^/>&/g'| sed 's/\// /g'>nr_braB_braz_strain_protein.tsv #37914 
#替换菌株名和蛋白名
perl -alne '/\>(.*?)\s.*(\[.*\])/;$name=$1;$strain=$2;$strain=~s/\[//g;$strain=~s/\]//g;$strain=~s/\s/\_/g;$id=$strain."\_".$name;print"$name\t$id";' nr_braB_braz_strain_protein.tsv   > nr_braB_braz_strain_protein.repalce.tsv
sed -i 's/(//g'  nr_braB_braz_strain_protein.repalce.tsv
sed -i 's/)//g'   nr_braB_braz_strain_protein.repalce.tsv
#替换序列蛋白名为菌株_蛋白名
faops replace  nr_braB_braZ_protein_name.fa  nr_braB_braz_strain_protein.repalce.tsv nr_braB_braz_strain_protein.repalce.fa #37914
#cd-hit聚类共5876个簇(设置簇内相似性为90%)
cd-hit -i nr_braB_braz_strain_protein.repalce.fa -c 0.9 -n 5  -d 0 -M 16000 -T 8  -o nr_braB_braZ_cdhit_identity90.fa  #5994
cat  PAO1_braB_braZ.fa  >> nr_braB_braZ_cdhit_identity90.fa  ##5996
#多序列比对
mafft --retree 1 --maxiterate 0 nr_braB_braZ_cdhit_identity90.fa >nr_braB_braZ_cdhit_identity90_mafft.fa
#建树
bsub -q mpi -n 24 -J "iq" ./iqtree2 -s nr_braB_braZ_cdhit_identity90_mafft.fa -m MFP  --prefix nr_braB_braZ_cdhit_identity90_mafft  -T 20 -B 1000 -bnni 



```

























