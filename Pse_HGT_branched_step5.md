
<!-- TOC -->

- [1.分析braB和braZ的泛基因组频率](#1分析brab和braz的泛基因组频率)
  - [1.1roray分析泛基因组](#11roray分析泛基因组)
  - [1.2ppanggolin分析泛基因组](#12ppanggolin分析泛基因组)
  - [1.3分析braB和braZ的泛基因组频率](#13分析brab和braz的泛基因组频率)
  - [1.4对braB和braZ的泛基因组频率及泛基因组分布进行M-W检验和卡方检验](#14对brab和braz的泛基因组频率及泛基因组分布进行m-w检验和卡方检验)
- [2.使用braZ上下游50000bp在nt库里找braZ水平基因转移的证据](#2使用braz上下游50000bp在nt库里找braz水平基因转移的证据)
  - [2.1自己在NCBI下载的bacteria representive数据](#21自己在ncbi下载的bacteria-representive数据)
  - [2.2根据王老师的nwr/doc/assembly.md下载数据](#22根据王老师的nwrdocassemblymd下载数据)
- [3.统计二代测序宏基因组样本的抗性基因存在情况以及三代测序宏基因组样本的抗性基因存在情况](#3统计二代测序宏基因组样本的抗性基因存在情况以及三代测序宏基因组样本的抗性基因存在情况)
  - [3.1三代污泥宏基因组中的细菌物种存在情况(使用的是sludge_bin.fa)](#31三代污泥宏基因组中的细菌物种存在情况使用的是sludge_binfa)
  - [3.2三代chicken gut宏基因组中的细菌物种存在情况(使用的是chicken_bin.fa)](#32三代chicken-gut宏基因组中的细菌物种存在情况使用的是chicken_binfa)
  - [3.3二代wastewater和chicken中的ARG结果统计](#33二代wastewater和chicken中的arg结果统计)
  - [3.4三代污泥中的ARG结果统计(使用的是sludge.fa)](#34三代污泥中的arg结果统计使用的是sludgefa)
  - [3.5三代chicken gut 中的ARG结果统计(使用的是chicken.fa)](#35三代chicken-gut-中的arg结果统计使用的是chickenfa)
- [4.使用RColorBrewer](#4使用rcolorbrewer)

<!-- /TOC -->

# 1.分析braB和braZ的泛基因组频率
* (core genes)核心基因组  (accessory genes)附属基因组 (specific genes)特有基因。核心基因组即所有菌株共有的基因，这些基因参与基础生物学过程，如基因表达，能量产生，氨基酸代谢等，一般是house-keeping gene。附属基因指存在于部分菌株中的基因，这些基因与物种的多样性有关，赋予个体竞争优势。特有基因只存在于某一菌株中，这些基因一般是通过水平基因转移(HGT)而来，通常与该菌株的独特表型特征有关，比如对特定环境的适应性，或独特的致病性。
* 根据物种的泛基因组大小与菌株数目的关系，将物种的泛基因组分为开放型泛基因组（open）和闭合型泛基因组（close）。开放型的泛基因组是指，随着测序的基因组数目的增加，物种的泛基因组大小也不断增加。闭合性的泛基因组是指，随着测序的基因组数目增加，物种的泛基因组大小增加到一定的程度后收敛于某一值。
* 在附加基因组中，若仅有一个个体具有该基因，则可称之为独特基因（英语：unique gene）。为了允许注解及基因序列组装的错误，对核心基因组较为宽松的定义可称之为软核心基因（soft core gene），其定义为于 95％以上的个体具有此基因。核心基因：一般控制着生命体基本新陈代谢的功能，因为它们广泛存在所有个体中，是不可缺少的。可变基因：往往只存在于一部分个体中，可能就是导致个体产生特异性的性状（抗病性，抗寒性等重要农艺性状）的原因。
* 泛基因组分析软件roary和BPGA和PGAP和PanOCT和peppan等
* Roary要求输入文件必须是GFF3文件(.gff)，且文件的最后要有FASTA格式的全基因组信息。而在Genbank中的.gff文件只有注释信息，缺少基因组序列。因此需要下载.gbff格式的文件，利用bp_genbank2gff3.pl脚本将其转化成gff文件，即可作为Roary的输入文件。
## 1.1roray分析泛基因组
```r
#conda安装泛基因组分析软件roary
#Roary要求输入文件必须是GFF3文件(.gff)且文件的最后要有FASTA格式的全基因组信息，而在Genbank中的.gff文件只有注释信息且缺少基因组序列,因此需要下载.gbff格式的文件。
#方法一：利用bp_genbank2gff3.pl脚本将其转化成gff文件，即可作为Roary的输入文件。例如：
#方法二：使用Prokka软件将FASTA格式的基因组数据转化为GFF格式，也可以作为Roary的输入文件。例如：
conda install roary
roary -h
roary -e  -n -v -i 80 -p 4 -f ./result_roary   ./prokka_gff/*.gff
-f 输出目录
-e 使用prank创建一个核心基因的多fasta比对
-n 使用maff进行快速核心基因组比对
-v 输出的格式
-r 创建r 图，需要R和ggplot2
-i blastp的最小相似性
#结果文件
（1）summary_statistics.txt   包含了不同存在类型的gene个数及gene总数。
如果core genes数或者total genes数特别高，那么要小心你的输入文件是否混入了其他种的基因组，或者你的样品是否存在污染。
（2）gene_presence_absence.csv
该文件记录了每个输入样本中某个基因存在或缺失的信息，可以用Excel等软件打开，其中还包括基因名名称、功能注释、存在该cluster的菌株的数量等等。
（3）gene_presence_absence.Rtab
该文件第一行为每个输入样本的名称，第一列为每个基因名称。文件是由0/1构成的矩阵，0代表缺失，1代表存在。可以用R的read.table载入，进行后续的分析。
（4）pan_genome_reference.fa
  FASTA格式包含每个cluster的一个代表序列，组成pan-genome
（5）core_gene_alignment.aln
   core genes多重比对的输出结果，可以作为构建系统发生树的输入文件。

accessory_binary_genes.fa 非核心基因的二进制分布数据，以0/1表示携带或不携带
accessory_binary_genes.fa.newick 非核心基因的二进制分布数据的newick树图数据文件
accessory_graph.dot 非核心基因点图
accessory.header.embl 非核心基因数据头信息，以embl格式保存
accessory.tab 非核心基因信息
blast_identity_frequency.Rtab blast比对一致性结果的R语言工具
clustered_proteins 聚类的蛋白质
core_accessory_graph.dot 核心基因点图
core_accessory.header.embl embl 格式的文件显示各 accessory 基因
core_accessory.tab accessory 基因在所在的基因组
core_alignment_header.embl 核心序列比对结果的头信息，以embl格式保存
core_gene_alignment.aln 核心基因序列比对
core_gene_alignment.aln.reduced 核心基因序列比对，去除冗余数据
gene_presence_absence.csv csv 格式的基因在各个基因组中是否存在的数据文件
gene_presence_absence.Rtab Rtab 格式的基因在各个基因组中是否存在的数据文件
number_of_conserved_genes.Rtab Rtab 格式的不同数量基因组所共有基因数
number_of_genes_in_pan_genome.Rtab Rtab 格式的不同数量基因组的所有基因数
number_of_new_genes.Rtab Rtab 格式的不同数量基因组所新增的基因数
number_of_unique_genes.Rtab Rtab 格式的不同数量基因组所特有基因数
pan_genome_reference.fa 泛基因组参考序列
summary_statistics.txt pangenome 分析各种基因数量结果
```

## 1.2ppanggolin分析泛基因组
```r
#提取genbank文件
cat strains.taxon.tsv | grep "Pseudom_aeru" >Pse_aeru_strain_GCF.tsv
wc -l Pse_aeru_strain_GCF.tsv

for f in $(cat Pse_aeru_strain_GCF.tsv)
do
cp  ASSEMBLY/$f/*.gbff.gz ASSEMBLY/$f/$f.gbff.gz
done

for f in $(cat Pse_aeru_strain_GCF.tsv)
do
mv  ASSEMBLY/$f/$f.gbff.gz /mnt/d/pan
done

#将文件名和路径合并为一个列表，方便泛基因组分析
for filename in *.gbff
do
base=$(basename $filename .gbff)
echo -e "$base\t$base.gbff" >>pan_name_dir.tsv
done 

#进行pangeomone分析，默认覆盖度为0.8，默认相似性为0.8
ppanggolin workflow --anno ./pan_name_dir.tsv --cpu 8 -o ./result

#依次查看结果文件
##### 1.gene_presence_absence.Rtab是不同菌株中基因family的存在情况，存在记为1
cat gene_presence_absence.Rtab | grep -v "Gene" | wc -l  #25001
cat gene_presence_absence.Rtab | grep -o "Pse" | wc -l  #391
#共25001个基因，但是此处的基因名是泛基因组软件定义的gene family

##### 2.matrix.csv 查看矩阵文件
#共25001个基因family  分为persistent shell和cloud基因。还包括基因的平均长度，平均菌株中存在数量等

##### 3.mean_persistent_duplication.tsv 查看persistent_family的重复率和平均出现率或者是否是单拷贝marker
wc -l mean_persistent_duplication.tsv #5132

#### 4.organisms_statistics.tsv 391个菌株的核心基因组数目，附属基因组数目

###### 5.查看文件夹partitions里面的内容,共三种分类...此外将所有基因分为11个parititions
wc -l persistent.txt  #5132
wc -l  cloud.txt   #14316
wc -l  shell.txt   # 5552

wc -l exact_core.txt    #1437
wc -l exact_accessory.txt #23563

wc -l soft_accessory.txt  #19959
wc -l  soft_core.txt  #5041
#可以发现将shell基因组分为了11个级别
cat S1.txt   S11.txt  S3.txt  S5.txt  S7.txt  S9.txt S10.txt  S2.txt   S4.txt  S6.txt  S8.txt | wc -l

##### 6.查看文件夹project里面的内容，共391个菌株泛基因组的相似信息，包含基因名，基因起始位置，终止位置，蛋白family 菌株种的nb_copy 每个基因对应的partition  persistent_neighbors shell_nei cloud_nei
#目前存在问题：每个菌株里对应基因的名字都不一致
cat *.tsv | grep -v "gene"|cut -f 1 | sort -n | uniq | wc -l #2384834
#合并391个菌株的family,共#25001
cat *.tsv | grep -v "gene"| cut -f 6 | sort -n | uniq | wc -l  #25001
```

## 1.3分析braB和braZ的泛基因组频率
```r
#提取braB共表达的基因对应的family和braZ共表达的基因对应的famliy
cd pangenome 
mkdir braB_braZ
cp braB_three_Pfam_id.tsv braZ_three_Pfam_id.tsv Pseudom_aeru_PAO1.tsv ./braB_braZ
cd braB_braZ
cat Pseudom_aeru_PAO1.tsv | grep -f braB_three_Pfam_id.tsv |cut -f 6 >braB_family.tsv
cat Pseudom_aeru_PAO1.tsv | grep -f braZ_three_Pfam_id.tsv |cut -f 6 >braZ_family.tsv
#提取family对应的泛基因组分类以及在菌株里的平均序列
cat matrix.csv | tr "," "\t" | cut -f 1,2,6 | sed 's/"//g' >famliy_pan_class.tsv
cat famliy_pan_class.tsv | grep -f braB_family.tsv  #发现braB都在persistent上
cat famliy_pan_class.tsv | grep -f braZ_family.tsv  #发现braZ都在persistent上
#查看泛基因组的其他分类
#给partitions文件夹下的文件都加上一列文件名
for filename in *.txt
do
base=$(basename $filename .txt)
echo $base
perl name.pl $base.txt >$base.tsv
done
#行数和上面比对，没有问题
#分别合并三种分类
cat persistent.tsv cloud.tsv shell.tsv >pan_class1.tsv
cat exact_core.tsv exact_accessory.tsv >pan_class2.tsv
cat soft_accessory.tsv  soft_core.tsv  >pan_class3.tsv
#查看braB和braZ的信息
cd ../braB_braZ
cat pan_class1.tsv | grep -f braB_family.tsv >pan_class1_braB.tsv
cat pan_class1.tsv | grep -f braZ_family.tsv >pan_class1_braZ.tsv
cat pan_class2.tsv | grep -f braB_family.tsv >pan_class2_braB.tsv
cat pan_class2.tsv | grep -f braZ_family.tsv >pan_class2_braZ.tsv
cat pan_class3.tsv | grep -f braB_family.tsv >pan_class3_braB.tsv
cat pan_class3.tsv | grep -f braZ_family.tsv >pan_class3_braZ.tsv
#查看braB及braZ共表达的基因在391个菌株里面的频率
#awk '{i=1;while(i <= NF){col[i]=col[i] $i " ";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' gene_presence_absence.Rtab  | sed 's/[ \t]*$//g'
cat  gene_presence_absence.Rtab | grep -f braB_family.tsv  | perl -MStatistics::Descriptive -alne '($gene,$exp)=(split/\t/,$_,3)[0,2];@data=(split/\t/,$exp);$stat = Statistics::Descriptive::Full->new();$stat->add_data(@data);$mean = $stat->mean();$mean=sprintf"%.2f",$mean;print"$gene\t$mean";' >braB_mean_presence.tsv
cat  gene_presence_absence.Rtab | grep -f braZ_family.tsv  | perl -MStatistics::Descriptive -alne '($gene,$exp)=(split/\t/,$_,3)[0,2];@data=(split/\t/,$exp);$stat = Statistics::Descriptive::Full->new();$stat->add_data(@data);$mean = $stat->mean();$mean=sprintf"%.2f",$mean;print"$gene\t$mean";' >braZ_mean_presence.tsv

#将braB和braZ结果更改为柱状图格式
cat braB_mean_presence.tsv | tsv-join -d 1 -f Pseudom_aeru_PAO1.tsv -k 6 --append-fields 1 --allow-duplicate-keys | tsv-select -f 3,2 | sed '1iname\tnum' >braB_picture.tsv
cat braZ_mean_presence.tsv | tsv-join -d 1 -f Pseudom_aeru_PAO1.tsv -k 6 --append-fields 1 --allow-duplicate-keys | tsv-select -f 3,2 | sed '1iname\tnum' >braZ_picture.tsv

#绘制braB及cor频率分布柱状图
library(ggplot2)
setwd("D:/pangenome/braB_braZ")
data<-read.table("braB_picture.tsv",header=T)
library(showtext)#ggplot作图显示字体错误时需要用showtext包输入需要的字体
font_add('Arial','/Library/Fonts/Arial.ttf')
showtext_auto()
library(patchwork)
#设置柱状图的颜色,调整y轴的范围
#Default：stat="count" 表示从给定的数据里，统计每个类别出现的次数；此时aes()只需要给定x参数即可；
#stat="identity"表示直接指定每种类别的频数；此时aes()除了需要给定x参数交代类别，还需要指定y参数表示频数值。
p<-ggplot(data,aes(x=name,y=num))+geom_bar(stat="identity", position="dodge",fill= "#FC9272",width=0.9)+theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90,size=10,color='black', hjust = 1, vjust =1))+theme(text=element_text(size=16,family="Arial",face="bold"))+ylim(0,1)
ggsave('braB_cor_frequency.pdf', p,dpi = 480, width=8, height=6)   



#绘制braZ及cor在391个菌株的频率分布图
data1<-read.table("braZ_picture.tsv",header=T)
p1<-ggplot(data1,aes(x=name,y=num))+geom_bar(stat="identity", position="dodge",fill= "#6BAED6",width=0.9)+theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90,size=10,color='black', hjust = 1, vjust =1))+theme(text=element_text(size=16,family="Arial",face="bold"))+ylim(0,1)
ggsave('braZ_cor_frequency.pdf', p1,dpi = 480, width=8, height=6)   

```

## 1.4对braB和braZ的泛基因组频率及泛基因组分布进行M-W检验和卡方检验
```r
cd /mnt/d/pangenome/braB_braZ
cat <(perl -alne '$name=(split/\t/,$_)[1];print"braB_cor\t$name";' braB_braZ_cor_pan/pan_class2_braB.tsv)  <(perl -alne '$name=(split/\t/,$_)[1];print"braZ_cor\t$name";' braB_braZ_cor_pan/pan_class2_braZ.tsv) | sed '1igene\tpan' >braB_braZ_cor_aquare.tsv
setwd("D:/pangenome/braB_braZ")
set.seed(123)
data<-read.table("braB_braZ_cor_aquare.tsv",header=T)
braB_cor<-c(6,15)
braZ_cor<-c(14,15)
data_square<-data.frame(braB_cor,braZ_cor,row.names=c("exact_core","exact_accessory"))
chisq.test(t(data_square))
cat <(perl -alne '$name=(split/\t/,$_)[1];print"braB_cor\t$name";' braB_mean_presence.tsv)  <(perl -alne '$name=(split/\t/,$_)[1];print"braZ_cor\t$name";' braZ_mean_presence.tsv) | sed '1igene\tpan' >braB_braZ_cor_wc.tsv
data1<-read.table("braB_braZ_cor_wc.tsv",header=T)
wilcox.test(pan~gene,paired=FALSE)
```

# 2.使用braZ上下游50000bp在nt库里找braZ水平基因转移的证据

* nt数据只下载变形菌纲的，因此找到变形菌纲对应的taxonmy
* NR(Non-Redundant Protein Sequence Database)非冗余蛋白库，所有GenBank+EMBL+DDBJ+PDB中的非冗余蛋白序列，对于所有已知的或可能的编码序列，NR记录中都给出了相应的氨基酸序列（通过已知或可能的读码框推断而来）以及专门蛋白数据库中的序列号。NR库相当于一个以核酸序列为基础的交叉索引，将核酸数据和蛋白数据联系起来。NT(Nucleotide Sequence Database),核酸序列数据库，是NR库的子集。
* RefSeq(the reference sequence database,https://www.ncbi.nlm.nih.gov/refseq/ ).参考序列数据库，包含RefSeq_genomic(NCBI genomic reference sequences)，RefSeq_protein(NCBI protein reference sequences)和RefSeq transpans(NCBI transpans reference sequences)具有生物意义上的非冗余基因,转录本和蛋白质序列，是经过NCBI和其他组织校正的数据库，使用人类基因命名委员会定义的术语，并且包括了官方的基因符号和可选的符号。RefSeq记录有三种可以获得的状态：预测的、临时的和检查过的（reviewd）。预测的RefSeq记录是来自于那些未知功能的cDNA序列，它们有一个预测的蛋白编码区;临时的RefSeq记录还没有被检查过,它们是有自动的程序产生的；检查过的记录代表了目前关于一个基因和它的转录子的知识的汇编，它们很多都来自于GenBank记录、人类基因组命名委员会和OMIM，RefSeq标准为人类基因组的功能注解提供一个基础。
* RefSeq数据库和GenBank数据库的区别在于：GenBank是一个开放的数据库，对每个基因都含有许多序列。很多研究者或者公司都可以自己提交序列，另外这个数据库每天都要和EMBL和DDBJ交换数据。genbank的数据可能重复或者不准。而RefSeq数据库被设计成每个人类位点挑出一个代表序列来减少重复，是NCBI提供的校正的序列数据和相关的信息。数据库包括构建的基因组contig、mRNA、蛋白和整个染色体。refseq序列是NCBI筛选过的非冗余数据库，一般可信度比较高。
* 下载NR数据库,下载分类数据库,taxdump 目录中有两个重要文件：names.dmp：记录物种名及其分类编号；nodes.dmp：记录分类编号的分类节点信息）,下载accession与taxid的对应关系
wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz taxdump.tar.gz.md5
wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz 
wget -c https://github.com/shenwei356/taxonkit/releases/download/v0.6.0/taxonkit_linux_amd64.tar.gz
wget -c https://github.com/shenwei356/csvtk/releases/download/v0.20.0/csvtk_linux_amd64.tar.gz
* 用到的工具taxonkit和csvtk
* bacteria对应的taxid是2
* https://www.ncbi.nlm.nih.gov/assembly/?term=bacteria+representative

## 2.1自己在NCBI下载的bacteria representive数据
```bash

rsync常用参数
--copy-links 拷贝链接对应的文件或者目录而非链接本身(对于NCBI数据库很适用，比如Refseq数据库的数据其实是链接回genome/all数据库的)
--recursive -r 递归目录（有时候下载的是整个文件夹，可使用这个参数）
--progress 表示在同步的过程中可以看到同步的过程状态，比如统计要同步的文件数量、 同步的文件传输速度等。
--partial 保持部分传输的文件
--P 相当于 --partial --progress
--exclude=PATTERN 排除规则PATTERN。指定同步需要过滤掉的文件或子目录（即不需要同步过去的），后面直接跟不需要同步的单个文件或子目录（不需要跟路径），可以是通配符模式（如\ *.txt）。过滤多个文件或子目录，就使用多个--exclude
--exclude-from=FILE 从文件中读取排除规则。指定同步需要过滤掉的文件或子目录，后面跟文件（比如/root/exclude.txt）,不需要同步的文件和子目录放到/root/exclude.txt下。
--include=PATTERN 不要排除指定规则的文件
--include-from=FILE 从文件中读取包含的规则
--files-from=FILE 从文件中读取文件列表
--delete 删除 DEST 中 SRC 没有的文件
#选择id table 和taxonomy id并下载bacteria representitive的refer accession(直接在NCBI下载，名字是assembly_result.txt)(14912)
#下载refer数据库的bacteria的下载链接(文件名字是assembly_summary.txt)247607
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
#过滤得到bacteria的reprensentitive的genome链接(14912)
cat assembly_summary.txt | grep -f <(cut -f 1 assembly_result.txt) >bacteria_representitive.tsv
#准备下载信息,并且下载这14912个基因组文件
cut -f 20 bacteria_representitive.tsv | perl -alne '$name=$_;$id=(split/\//,$_)[9];print"$id\t$name"."\/"."$id"."_genomic.fna.gz";' >bacteria_repre_download.tsv
cat bacteria_repre_download.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 '
    wget {2}
    '
#解压这些基因组压缩文件，分割然后合并
gunzip *.gz
python floder_split.py 
#共有15条序列
faops size *.fna #836316
#更改序列名
cat bac_repre.1.fna  | sed 's/\s/\//g' >bac_repre.1.re.fna
faops size *.re.fna #836316
#建库
for filename in *.re.fna
do
base=$(basename $filename .re.fna)
echo $base
bsub -q mpi -n 24 -J "db" makeblastdb -in $base.re.fna -dbtype nucl -out $base
done
#比对
for filename in *.re.fna
do
base=$(basename $filename .re.fna)
echo $base
bsub -q mpi -n 24 -J "BL" blastn -db $base  -query braZ_50000bp.fa -out $base.braZ.tsv -outfmt 6 -num_threads 20  -evalue 1e-5
done 
#合并分割的文件
cat *.tsv >bac_repre_braZ_50000.tsv
```

## 2.2根据王老师的nwr/doc/assembly.md下载数据
* 查看细菌物种分类
```r
#安装nwr和初始化taxonomy,assembly数据库
#进入数据库
cd ~/.nwr
#查看文件assembly_summary_genbank.txt和assembly_summary_refseq.txt
wc -l assembly_summary_refseq.txt #253139
wc -l assembly_summary_genbank.txt #1223766
# assembly_summary_refseq.txt 该文件和自己从NCBI下载的assembly_summary.txt文件内容还是有区别的
#下载细菌基因组文件
#首先获得一个细菌的属的taxid并排序
GENUS=$(
    nwr member Bacteria -r genus |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        grep -v " sp." |
        sed '1d' |
        cut -f 1 |
        sort -n
)

#sqlite3创建数据库,所有的SQLite语句都是以关键字(如：SELECT，INSERT，UPDATE，DELETE，ALTER，DROP等)开始的。所有语句都以分号(;)结尾。
#SELECT column1, column2, columnN FROM table_name; 从某个表中获得某些列
#selite3 -tabs将数据库变成tab格式
#HAVING 子句允许指定条件来过滤将出现在最终结果中的分组结果
#最后输出文件第一列是taxid(species)，第二列是species，第三列是species在数据库中的数量
#遍历细菌属的taxid，从ar(ar_refseq)这个数据库中获得species_id和species这两列
#当属的genus_id等于tax_id(细菌的)时且组装水平为complete,chromosome时,对结果按照species_id分组
#只显示对于species_id或者species计数大于2的信息
for R in ${GENUS}; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = $R
            AND assembly_level IN ('Complete Genome', 'Chromosome')
        GROUP BY species_id
        HAVING count > 2
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |head

#nwr lineage 获取tax_id对应的分类种属科目纲门界信息
#parallel 中的{1}获取文件的第一列的内容
#grep -E "\sArchaea\s"匹配细菌 -A 1如果匹配成功，将匹配行及其后1行一起打印出来  -B匹配行前面,-C匹配行前后
#添加科信息nwr append stdin -r family 
for R in ${GENUS}; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = $R
            AND assembly_level IN ('Complete Genome', 'Chromosome')
        GROUP BY species_id
        HAVING count > 2
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
#sed '1d' |  这里不需要删除第一行因为并没有列名
    parallel --col-sep "\t" -j 1 '
        GROUP=$(
            nwr lineage {1} |
                grep -E "\sBacteria\s" -A 1 |
                sed "1d" |
                tsv-select -f 2
        )
        printf "%s\t%s\t%s\t%s\n" {1} "$GROUP" "{2}" {3}
    ' |
    sed "s/'//g" |
    nwr append stdin -r family |
    tsv-select -f 1,2,5,3,4 |
    tsv-sort -k2,2 -k3,3 -k4,4 |
    (echo -e '#tax_id\tgroup\tfamily\tspecies\tRS' && cat) \
    > species.count.tsv

    cat species.count.tsv |
    mlr --itsv --omd cat| wc -l
    #共929个物种
```
* 下载细菌reprensentive基因组
```r
#下载refer数据库所有的细菌assembly
GENUS=$(
    nwr member Bacteria -r genus |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        sed '1d' |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN ($GENUS)
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > assembly_one.tsv
    wc -l  assembly_one.tsv #237393


  #下载refer数据库的细菌的的representative和reference
  GENUS=$(
    nwr member Bacteria -r genus |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        sed '1d' |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)
echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN ($GENUS)
        AND refseq_category IN ('reference genome','representative genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > assembly_two.tsv
#查看基因组有多少个
wc -l assembly_two.tsv #14602
#查看组装水平有哪些 
tsv-summarize -g 6 --count  assembly_two.tsv
#Download all assemblies
cut -f 10 assembly_two.tsv | grep -v "ftp_path" |perl -alne '$name=$_;$id=(split/\//,$_)[9];print"$id\t$name"."\/"."$id"."_genomic.fna.gz";' >bacteria_repre_download.tsv
#bacteria_repre_download.tsv 14601
cat bacteria_repre_download.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 '
    wget {2}
    '
###比较自己下载的基因组和使用nwr下载的基因组名称重叠数量
#自己下载的14912   nwr下载的14602  重叠的14524
#自己下载的基因组可能还是存在部分问题，后续有空查找原因


#解压这些基因组压缩文件，分割然后合并
gunzip *.gz
python floder_split.py 
#共有15条序列
faops size *.fna #836316
#更改序列名
cat bac_repre.1.fna  | sed 's/\s/\//g' >bac_repre.1.re.fna
faops size *.re.fna #836316
#建库
for filename in *.re.fna
do
base=$(basename $filename .re.fna)
echo $base
bsub -q mpi -n 24 -J "db" makeblastdb -in $base.re.fna -dbtype nucl -out $base
done
#比对
for filename in *.re.fna
do
base=$(basename $filename .re.fna)
echo $base
bsub -q mpi -n 24 -J "BL" blastn -db $base  -query braZ_50000bp.fa -out $base.braZ.tsv -outfmt 6 -num_threads 20  -evalue 1e-5
done 
#合并分割的文件
cat *.tsv >bac_repre_braZ_50000.tsv
```






# 3.统计二代测序宏基因组样本的抗性基因存在情况以及三代测序宏基因组样本的抗性基因存在情况
## 3.1三代污泥宏基因组中的细菌物种存在情况(使用的是sludge_bin.fa)
```bash
#提取sludge组装的序列
cd hifiasm/sludge
cat sludge.hifiasm-meta.p_ctg.gfa | grep "S" | perl -alne '($id,$seq)=(split/\t/,$_)[1,2];print">$id\n$seq";' >sludge.fa
#将sludge组装的序列和sludge.hifiasm-meta.binning.tsv结果比对
faops size sludge.fa >sludge.abstract.tsv  #25408
cat sludge.hifiasm-meta.binning.tsv | grep "CT" >sludge_hifisam_bin.tsv #8905
#两者重叠的(序列名字相同，序列大小也相同)
cat sludge.abstract.tsv | grep -f <(cut -f 2,3 sludge_hifisam_bin.tsv | grep -v "CT") | wc -l  #8905
#两者多余的
cat sludge.abstract.tsv | grep -v -f <(cut -f 2,3 sludge_hifisam_bin.tsv | grep -v "CC") | wc -l #16503
#多余的为什么没有划分到bin里面，不清楚？？？
cat sludge.abstract.tsv | grep -v -f <(cut -f 2,3 sludge_hifisam_bin.tsv | grep -v "CC") | sort -n | uniq | tail  -n 1
#s9999.ctg020666l        44864
#猜测原因：可能是contig比较短以及评估质量较低，被舍弃了
#因此提取bin文件中的序列名字和序列然后进行比对
faops some sludge.fa <( cut -f 2,3 sludge_hifisam_bin.tsv | grep -v "CT" | cut -f 1)  sludge_bin.fa
#将nt的bacteria的repsentitive和hifi组装的污泥基因组进行比对
for filename in *.re.fna
do
base=$(basename $filename .re.fna)
echo $base
bsub -q mpi -n 24 -J "BL" blastn -db $base  -query sludge_bin.fa -out $base.sludge.bin.tsv -outfmt 6 -num_threads 20  -evalue 1e-5 
done
#统计三代宏基因组污泥的比对结果
#细菌基因组特点： 基因组：一般在0.16- 13Mb，大部分在5M左右500,0000
cat  bac_repre_sludge.tsv | tsv-filter --gt 4:50000 
#共61条比对结果，且相似性在99%左右,查看几条比对上的结果，有的是肠道内发现的新物种，有的是发酵菌
```

## 3.2三代chicken gut宏基因组中的细菌物种存在情况(使用的是chicken_bin.fa)
```bash
#提取chicken gut 组装的序列
cd hifiasm/chicken
cat chicken.hifiasm-meta.p_ctg.gfa | grep "S" | perl -alne '($id,$seq)=(split/\t/,$_)[1,2];print">$id\n$seq";' >chicken.fa
#将chicken组装的序列和chicken.hifiasm-meta.binning.tsv结果比对
faops size chicken.fa >chicken.abstract.tsv        #17073
cat chicken.hifiasm-meta.binning.tsv | grep "CT" >chicken_hifisam_bin.tsv #8365
#两者重叠的(序列名字相同，序列大小也相同)
cat chicken.abstract.tsv | grep -f <(cut -f 2,3 chicken_hifisam_bin.tsv | grep -v "CT") | wc -l  #8365
#两者多余的
cat chicken.abstract.tsv | grep -v -f <(cut -f 2,3 chicken_hifisam_bin.tsv | grep -v "CT") | wc -l #8708
#多余的为什么没有划分到bin里面，不清楚？？？
cat chicken.abstract.tsv | grep -v -f <(cut -f 2,3 chicken_hifisam_bin.tsv | grep -v "CT") | sort -n | uniq | tail  -n 1
#s9999.ctg020666l        14111
#猜测原因：可能是contig比较短以及评估质量较低，被舍弃了
#因此提取bin文件中的序列名字和序列然后进行比对
faops some chicken.fa <( cut -f 2,3 chicken_hifisam_bin.tsv | grep -v "CT" | cut -f 1)  chicken_bin.fa
#将nt的bacteria的repsentitive和hifi组装的chicken gut基因组进行比对
for filename in *.re.fna
do
base=$(basename $filename .re.fna)
echo $base
bsub -q mpi -n 24 -J "BL" blastn -db $base  -query chicken_bin.fa -out $base.chicken.bin.tsv -outfmt 6 -num_threads 20  -evalue 1e-5 
done
#统计三代宏基因组chicken gut 的比对结果
cat bac_repre_chicken.tsv | tsv-filter --gt 4:50000 
#共393条比对结果，且相似性在99%左右,查看几条比对上的结果，有很多是临床分离的病原菌，如大肠杆菌，志贺氏菌，粪肠球菌等
cat bac_repre_chicken.tsv | tsv-filter --gt 4:50000 | cut -f 2 | cut -d "/" -f 2,3 | sort -n | uniq  #45
```

## 3.3二代wastewater和chicken中的ARG结果统计
```bash
#统计wastewater平均每个样本比对上的基因的数量，基因长度（大于300bp）
#有56个样本可以比对上基因
cat AAA_waste_300.tsv | sed 's/,/\t/g'| sed 's/sample=//g' | sed 's/gene=//g' | sed 's/contig=//g' | tsv-summarize -g 1 --count | wc -l 
#有66个基因可以比对上样本 
cat AAA_waste_300.tsv | sed 's/,/\t/g'| sed 's/sample=//g' | sed 's/gene=//g' | sed 's/contig=//g' | tsv-summarize -g 2 --count 
#所有样本比对上的基因数目是2316  
cat AAA_waste_300.tsv | sed 's/,/\t/g'| sed 's/sample=//g' | sed 's/gene=//g' | sed 's/contig=//g' | tsv-summarize -g 2 --count | tsv-summarize --sum 2  


#统计chicken平均每个样本比对上的基因的数量，基因长度（大于300bp）
#有74个样本可以比对上基因
cat AAA_chicken_300.tsv | sed 's/,/\t/g'| sed 's/sample=//g' | sed 's/gene=//g' | sed 's/contig=//g' | tsv-summarize -g 1 --count | wc -l 
#有69个基因可以比对上样本 
cat AAA_chicken_300.tsv | sed 's/,/\t/g'| sed 's/sample=//g' | sed 's/gene=//g' | sed 's/contig=//g' | tsv-summarize -g 2 --count  | wc -l 
#所有样本比对上的基因数目是3010 
cat AAA_chicken_300.tsv | sed 's/,/\t/g'| sed 's/sample=//g' | sed 's/gene=//g' | sed 's/contig=//g' | tsv-summarize -g 2 --count | tsv-summarize --sum 2 
#平均每个样本的基因
3010/74=40

```
## 3.4三代污泥中的ARG结果统计(使用的是sludge.fa)
```bash
#建库
for filename in *.fa
do
base=$(basename $filename .fa)
echo $base
bsub -q mpi -n 24 -J "db" makeblastdb -in $base.fa -dbtype nucl -out $base
done
#比对
for filename in  *.fasta
do 
base=$(basename $filename .fasta)
echo $base
bsub -q mpi -n 24 -J "BL" blastn -db  sludge  -query $base.fasta -out $base.sludge_hifi.tsv -outfmt 6 -num_threads 20  -evalue 1e-5 
done
#有10个基因可以比对上一个样本(>300bp)
#比对上一个样本的基因数目是10条(>300bp)
cat  *.tsv | tsv-filter --ge 4:300 | wc -l #23

```

## 3.5三代chicken gut 中的ARG结果统计(使用的是chicken.fa)
```bash
#比对
for filename in  *.fasta
do 
base=$(basename $filename .fasta)
echo $base
bsub -q mpi -n 24 -J "BL" blastn -db  chicken  -query $base.fasta -out $base.chicken_hifi.tsv -outfmt 6 -num_threads 20  -evalue 1e-5 
done
#比对上一个样本的基因数目是289条(>300bp)
cat  *.tsv | tsv-filter --ge 4:300 | wc -l #415
#有21个基因可以比对上一个样本(>300bp)
cat  *.tsv | tsv-filter --ge 4:300 | tsv-summarize -g 1 --count | wc -l #22
```

# 4.使用RColorBrewer
```r
library(ggplot2)
library(RColorBrewer)
#创建数据集
x<-matrix(1:200,nrow=20,ncol=10)
x1<-x[,1:2]
colnames(x1)<-c("sample","num")
data<-data.frame(x1) 
data$sample<-as.factor(data$sample)
#将第一列定为样本类别，第二列定义为样本数量
#展示所有色板
display.brewer.all()
#展示"Set3"色板中的9个颜色
display.brewer.pal(9,'Set3')
#统计sample的数量 20
colorcount = length(data$sample)
samplecolor=colorRampPalette(brewer.pal(9,'Set3'))
#使用palette控制scale_fill_brewer()中的颜色选择
ggplot(data,aes(x=sample,y=num,fill=sample))+geom_bar(stat = "identity")
ggplot(data,aes(x=sample,y=num,fill=sample))+geom_bar(stat = "identity")+scale_fill_manual(values=samplecolor(colorcount))
```