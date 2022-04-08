- [1.使用diamond比对nr库搜索braz的来源](#1.使用diamond比对nr库搜索braz的来源)

* 1.上一步主要是分别统计铜绿假单胞菌等中的两个拷贝的ka/ks，比较两个拷贝的进化趋势，由结果可知两个拷贝braZ和braB的进化趋势不一致，功能发生了分化  
* 2.查看铜绿假单胞菌中braZ和braB的motif，发现在所有菌株中存在两个保守motif  
* 3.查看铜绿PAO1菌株的基因岛，并未发现有braz和braB  
* 4.查看假单胞属内铜绿假单胞菌中两个拷贝braz和braB在物种内都具有共线性,在物种间只有braB具有共线性  
* 5.下一步主要用braZ比对nr库，查看braZ从何转移而来  

# 1.使用diamond比对nr库搜索braz的来源(超算上运行)
* 准备数据和软件
```bash
#准备PAO1菌株的braz蛋白 Pseudom_aeru_PAO1_NP_250661
echo "Pseudom_aeru_PAO1_NP_250661" >branched-chain/branched-chain_PAO1_braz.tsv
faops some PROTEINS/all.replace.fa branched-chain/branched-chain_PAO1_braz.tsv  branched-chain/branched-chain_PAO1_braz.fa
#下载diamond
wget https://github.com/bbuchfink/diamond/releases/download/v2.0.5/diamond-linux64.tar.gz
#使用compare4上传diamond和branched-chain_PAO1_braz.fa到超算目录data/blast/fasta
cd data/blast/fasta
tar -zxvf diamond-linux64.tar.gz
#查看帮助信息
./diamond help

数据库参数：
--in <string>    default: STDIN
    输入FASTA格式的蛋白序列数据库文件。
--db | -d <string>
    设置数据库文件路径和前缀。创建数据库时，会生成一个后缀为.dmnd的数据库文件。比对时，则是输入相应的数据库文件。

输入参数：
--db | -d <string>
    设置数据库文件路径和前缀。创建数据库时，会生成一个后缀为.dmnd的数据库文件。比对时，则是输入相应的数据库文件。
--query | -q <string>    default: STDIN
    输入需要注释的FASTA或FASTQ格式的序列文件。可以是带.gz后缀的压缩文件。

比对参数：
--sensitive
    添加该参数，则能得到更多比对结果。该模式适合比对较长的序列。默认模式主要适用于比对short reads序列（Illumina reads），搜寻比对长度为30~40aa且bit得分大于50的匹配结果。
--more-sensitive
    相比于sensitive，能得到更全的比对结果。


输出参数：
--out | -o <string>    default：STDOUT
    设置输出文件。
--outfmt | -f <int>    default: 6
    设置输出格式。支持的格式有：0，两两比对格式；5，XML格式；6，BLAST表格格式；100，DIAMOND匹配存档（DDA）格式，该格式可以使用diamond view命令转换成其它格式；101，SAM格式；102，分类鉴定结果，该模式下结果文件分三列，QueryID、NCBI物种分类ID、最佳匹配evaule；103，PAF格式。
--evalue | -e <float>    default: 0.001
    设置比对的evalue阈值。
--id <int>
    设置identity阈值。
--query-cover <int>
    设置对query序列的覆盖度阈值。该阈值是HSP的阈值。
--subject-cover <int>
    设置对subject序列的覆盖度阈值。该阈值是HSP的阈值。

性能参数：
--threads | -p <int>    default: Max
    设置程序运行所使用的CPU线程数。默认是服务器可用的最大CPU线程数。
```

```bash
cd data/blast/fasta
#本地建立nr数据库并比对
gunzip nr.gz
mv nr nr.fa
bsub -q mpi -n 24 -J "DIA" ./diamond makedb --in nr.fa -d nr --threads 20
bsub -q mpi -n 24-J "DIA" ./diamond blastp -d nr.dmnd -q branched-chain_PAO1_braz.fa -o braZ_nr_result.tsv --id 50 --threads 20 -e 1e-5 --more-sensitive
#结果只有25行,可能建库时对序列id有要求

#本地建立braz数据库并比对
./diamond makedb --in branched-chain_PAO1_braz.fa -d braz
bsub -q serial -n 20 -J "DIA" ./diamond blastp -d braz.dmnd -q nr.fa -o braz_nr_reverse.result.tsv --id 30 --subject-cover 30 --threads 16 -e 1e-5 --more-sensitive
```

# 2.分割文件后使用diamond比对nr库搜索braz的来源(超算上运行)
```bash
#按序列名分割
faops size nr.fa | cut -f 1 >nr_name.tsv
split -l 5000000 nr_name.tsv -d -a 3 nr_
-l 按行分割
-d 添加数字后缀
-a 3 表示用3位数来顺序命名 后缀长度
url_ 分割后文件的前缀

#按序列分割
faops split-about nr.fa 3000000000 ./
#建库
for name in {000..059} 
do
echo $name
./diamond makedb --in $name.fa -d aaa_$name --threads 20
done
#比对
for name in {042..053} 
do
echo $name
bsub -q mpi -n 24 -J "DIA" ./diamond blastp -d braz.dmnd -q $name.fa -o aaa.$name.result.tsv --id 30  --subject-cover 30 --threads 20 -e 1e-5 --more-sensitive
done
cat *.tsv >braz_nr_whole.tsv
#分开比对的braz_nr_whole.tsv与未分开比对的braz_nr_reverse.result.tsv结果一致

#提取序列
faops some nr.fa <(cut -f 1 braz_nr_reverse.result.tsv) braz_nr_reverse.fa 
#筛选e值显著的序列（7050）
tsv-filter --lt 11:1.0e-100 braz_nr_reverse.result.tsv >braz_nr_evalue.tsv
faops some braz_nr_reverse.fa <(cut -f 1 braz_nr_evalue.tsv) braz_nr_evalue.fa
#筛选e值显著的序列对应的菌株蛋白名
cat nr_protein_name.tsv | grep -Ff <(cut -f 1 braz_nr_evalue.tsv) >braz_nr_strain_protein.tsv
#替换菌株名和蛋白名
perl -alne '/\>(.*?)\s.*(\[.*\])/;$name=$1;$strain=$2;$strain=~s/\[//g;$strain=~s/\]//g;$strain=~s/\s/\_/g;$id=$strain."\_".$name;print"$name\t$id";' braz_nr_strain_protein.tsv >braz_nr_strain_protein.replace.tsv
sed -i 's/(//g' braz_nr_strain_protein.replace.tsv
sed -i 's/)//g' braz_nr_strain_protein.replace.tsv
#替换序列蛋白名为菌株_蛋白名
faops replace braz_nr_evalue.fa  braz_nr_strain_protein.replace.tsv  braz_nr_strain_protein.replace.fa

```

# 3.使用cd-hit对抓取出来的序列进行聚类分簇，然后将分簇的序列建树
```bash
#cd-hit聚类共469个簇(设置簇内相似性为90%)
mkdir -p cd_hit
cd cd_hit
cd-hit -i braz_nr_strain_protein.replace.fa -o braz_nr_cdhit_identity90.fa -c 0.9
cat braz_nr_cdhit_identity90.fa PAO1.fa >braz_nr_PAO1_cdhit_identity90.fa
#建立树
muscle -in braz_nr_PAO1_cdhit_identity90.fa -out braz_nr_PAO1_cdhit_identity90_aln.fa
FastTree braz_nr_PAO1_cdhit_identity90_aln.fa >braz_nr_PAO1_cdhit_identity90_aln.newick

#cd-hit聚类共91个簇(设置簇内相似性为90%)
cd cd_hit
cd-hit -i braz_nr_strain_protein.replace.fa -o braz_nr_cdhit_identity70.fa -c 0.7
cat braz_nr_cdhit_identity70.fa PAO1.fa >braz_nr_PAO1_cdhit_identity70.fa
#建立树
muscle -in braz_nr_PAO1_cdhit_identity70.fa -out braz_nr_PAO1_cdhit_identity70_aln.fa
FastTree braz_nr_PAO1_cdhit_identity70_aln.fa >braz_nr_PAO1_cdhit_identity70_aln.newick

```










```







