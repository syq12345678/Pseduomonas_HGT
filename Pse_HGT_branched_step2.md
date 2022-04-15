<!-- TOC -->

- [1. ka/ks的文件准备](#1-kaks的文件准备)
  - [1.1 分别提取branched在铜绿假单胞菌中的两个拷贝蛋白序列](#11-分别提取branched在铜绿假单胞菌中的两个拷贝蛋白序列)
  - [1.2提取branched在铜绿假单胞菌中的两个CDS          branched-chain/branched-chain.CDS.fa](#12提取branched在铜绿假单胞菌中的两个cds----------branched-chainbranched-chaincdsfa)
  - [1.3 铜绿假单胞菌copy1准备列表文件，cds和protein](#13-铜绿假单胞菌copy1准备列表文件cds和protein)
  - [1.4 铜绿假单胞菌copy2准备列表文件，cds和protein](#14-铜绿假单胞菌copy2准备列表文件cds和protein)
- [2. 使用paraAT2和caculator2计算铜绿假单胞菌的kaks](#2-使用paraat2和caculator2计算铜绿假单胞菌的kaks)
  - [2.1安装paraAT2](#21安装paraat2)
  - [2.2 计算铜绿假单胞菌copy1的kaks](#22-计算铜绿假单胞菌copy1的kaks)
  - [2.3 计算铜绿假单胞菌copy2的kaks](#23-计算铜绿假单胞菌copy2的kaks)
  - [2.4 铜绿假单胞菌kaks结果画图(不能去除kaks等于0的点)](#24-铜绿假单胞菌kaks结果画图不能去除kaks等于0的点)
- [3. 使用paraAT2和caculator2计算其他假单胞菌的kaks(待定)](#3-使用paraat2和caculator2计算其他假单胞菌的kaks待定)
  - [4. 使用mega计算dn/ds](#4-使用mega计算dnds)
  - [5. 使用bp_pairwise_kaks计算kaks(bp_pairwise的结果与mega和kaks_calculator不一致，可能是算法的原因)](#5-使用bp_pairwise_kaks计算kaksbp_pairwise的结果与mega和kaks_calculator不一致可能是算法的原因)
- [5. 基因共线性](#5-基因共线性)
  - [5.1 准备假单胞菌属的典型菌株的genbank文件，共9种12个](#51-准备假单胞菌属的典型菌株的genbank文件共9种12个)
  - [2.2使用clinker查看基因共线性](#22使用clinker查看基因共线性)
- [3.MEME查找motif](#3meme查找motif)
- [4.islandviewer4查看基因岛](#4islandviewer4查看基因岛)
- [5.gephi查看pangenome](#5gephi查看pangenome)

<!-- /TOC -->

* 1.上一步主要是统计branched在不同菌株基因组中的拷贝数，可知铜绿假单胞菌等中具有两个拷贝    
* 2.分别统计铜绿假单胞菌等中的两个拷贝的ka/ks


# 1. ka/ks的文件准备
* 1.在遗传学中，Ka/Ks表示的是两个蛋白编码基因的非同义替换率(Ka)和同义替换率(Ks)之间的比例。这个比例可以判断是否有选择压力作用于这个蛋白编码基因。
* 2.如果有两个不同物种的同一个基因的序列，比如人和小鼠的p53基因，然后把这两个基因的序列进行比对，你会发现这两段序列有差异（进化！）。  
再仔细观察，你会发现有些碱基的变化导致了编码氨基酸的变化（非同义替换），有些没有导致编码氨基酸的变化（同义替换）。  
这是由密码子的简并性造成的，因为3个碱基决定1个氨基酸，所以64种碱基组合决定20种氨基酸，会有冗余出现。  
* 3.一般情况下，第三个碱基变化会造成同义替换，而第一二个碱基的变化会造成非同义替换。  
* 4.Ka和Ks的计算公式：  
Ka=发生非同义替换的SNP数/非同义替换位点数  
Ks=发生同义替换的SNP数/同义替换位点数
其中非同义替换位点数就是会造成氨基酸变化的位点数的总和，比如编码丝氨酸（ser）的第一二位碱基。  
而同义替换位点数就是不会造成氨基酸变化的位点数的总和，比如编码丝氨酸的第三位碱基。  
* 5.计算Ka/Ks时，不考虑start codon和stop codon  
* 6.转换（transition，嘌呤变嘌呤，嘧啶变嘧啶）发生的概率要高于颠换（transversion，嘌呤变嘧啶，嘧啶变嘌呤）发生的概率  
* 7.ka/ks与进化的关系：如果一个基因没有受到自然选择的压力，那么非同义替换率和同义替换率是相同的。而一般情况下，非同义替换会造成氨基酸变化  
改变蛋白质的构象和功能因此造成适应性的变化，从而带来自然选择的优势和劣势。而同义替换没有改变蛋白质的组成，不受自然选择的影响，因此ks反映  
过程中的背景碱基替换率。ka/ks的比值说明了这个基因受到了何种选择
* 8.Ka>>Ks或者Ka/Ks >> 1，基因受正选择(positive selection)  
Ka＝Ks或者Ka/Ks=1，基因中性进化(neutral evolution)  
Ka << Ks或者Ka/Ks << 1，基因受纯化选择(purify selection)  
* 9.实际情况下Ka/Ks << 1，因为一般非同义替换带来的都是有害的性状，只有极少数情况下会造成进化上的优势  
* 10.当Ka/Ks>>1时，基因受到强烈正选择，这样的基因即为近期正在快速进化的基因，对于物种的进化有着非常重要的意义
* 11.计算ka和ks的步骤
步骤一：假设比较的两条DNA序列之间的长度数为n,它们之间的替换数为m。计算同义(S)和非同义(N)位点的数量(S+N=n)以及同义(Sd)和非同义替换的数量  
(Sd+Nd=m)
步骤二:校正多个替换后,(Nd/N)和(Sd/S)分别代表ka和ks


## 1.1 分别提取branched在铜绿假单胞菌中的两个拷贝蛋白序列
```BASH
#提取铜绿假单胞菌的两个拷贝的蛋白名称(382*2=764)
cd ~/data/Pseudomonas
mkdir -p ~/data/Pseudomonas/branched-chain/kaks
cat branched-chain/branched-chain_hmmscan_copy.pfam.tsv | tsv-filter --ge 3:2 |
tsv-filter --str-in-fld 1:"Pseudom_aeru" | cut -f 1 >branched-chain/kaks/b
ranched-chain.Pseudom_aeru.two.copy.tsv
cat branched-chain/branched-chain_minevalue.tsv | grep -f branched-chain/kaks/branched-chain.Pseudom_aeru.two.copy.tsv |
cut -f 1 >branched-chain/kaks/branched-chain.Pseudom_aeru.protein.tsv
#提取铜绿假单胞菌中branced的两个拷贝蛋白序列(764)
faops some PROTEINS/all.replace.fa branched-chain/kaks/branched-chain.Pseudom_aeru.protein.tsv branched-chain/kaks/branched-chain.Pseudom_aeru.protein.fa

#使用cd-hit提取出branched-chain/branched-chain.Pseudom_aeru.protein.fa中差异最大的两条序列最为参考序列
cd-hit -i branched-chain/kaks/branched-chain.Pseudom_aeru.protein.fa -o branched-chain/kaks/branched-chain_refer.fa -c 0.8
#提取参考序列1
faops size branched-chain/kaks/branched-chain_refer.fa | cut -f 1 | head -n 1  >branched-chain/kaks/branched-chain_refer1.tsv
faops some branched-chain/kaks/branched-chain.Pseudom_aeru.protein.fa  branched-chain/kaks/branched-chain_refer1.tsv branched-chain/kaks/branched-chain_refer1.fa
#根据参考序列1提取簇1的序列(382)
diamond makedb --in branched-chain/kaks/branched-chain_refer1.fa -d branched-chain/kaks/branched-chain_refer1
diamond blastp -d branched-chain/kaks/branched-chain_refer1.dmnd -q branched-chain/kaks/branched-chain.Pseudom_aeru.protein.fa -o branched-chain/kaks/branched-chain_refer1_result
tsv-filter --gt 3:90 branched-chain/kaks/branched-chain_refer1_result  | cut -f 1 | sort -n | uniq >branched-chain/kaks/branched-chain_cluster1.tsv
#验证簇1,共382正确
cat branched-chain/kaks/branched-chain_cluster1.tsv | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 | cut -f 2 | sort -n | uniq | wc -l 

#蛋白簇1的序列名称和序列)382)
branched-chain/kaks/branched-chain_cluster1.tsv
faops some branched-chain/kaks/branched-chain.Pseudom_aeru.protein.fa  branched-chain/kaks/branched-chain_cluster1.tsv  branched-chain/kaks/branched-chain_cluster1.fa
#蛋白簇2的序列名称和序列（382）
cat branched-chain/kaks/branched-chain.Pseudom_aeru.protein.tsv | grep -v -f branched-chain/kaks/branched-chain_cluster1.tsv \
>branched-chain/kaks/branched-chain_cluster2.tsv
faops some branched-chain/kaks/branched-chain.Pseudom_aeru.protein.fa  branched-chain/kaks/branched-chain_cluster2.tsv branched-chain/kaks/branched-chain_cluster2.fa

#蛋白簇1进行排列组合(145542)
cat branched-chain/kaks/branched-chain_cluster1.fa | mash sketch -k 21 -s 1000 -i -p 8 - -o branched-chain/kaks/branched-chain_cluster1.msh
mash dist branched-chain/kaks/branched-chain_cluster1.msh branched-chain/kaks/branched-chain_cluster1.msh >branched-chain/kaks/branched-chain_cluster1.mash
cut -f 1,2  branched-chain/kaks/branched-chain_cluster1.mash | tsv-filter --ff-str-ne 1:2 >branched-chain/kaks/branched-chain_cluster1.arrange.tsv
rm -rf branched-chain/kaks/branched-chain_cluster1.msh branched-chain/kaks/branched-chain_cluster1.mash
#蛋白簇2进行排列组合(145542)
cat branched-chain/kaks/branched-chain_cluster2.fa | mash sketch -k 21 -s 1000 -i -p 8 - -o branched-chain/kaks/branched-chain_cluster2.msh
mash dist branched-chain/kaks/branched-chain_cluster2.msh branched-chain/kaks/branched-chain_cluster2.msh >branched-chain/kaks/branched-chain_cluster2.mash
cut -f 1,2  branched-chain/kaks/branched-chain_cluster2.mash | tsv-filter --ff-str-ne 1:2 >branched-chain/kaks/branched-chain_cluster2.arrange.tsv
rm -rf branched-chain/kaks/branched-chain_cluster2.msh branched-chain/kaks/branched-chain_cluster2.mash
```

## 1.2提取branched在铜绿假单胞菌中的两个CDS          branched-chain/branched-chain.CDS.fa
```BASH
#将branched-chain_minevalue.tsv结果改成replace文件格式(2140)
cut -f 1 branched-chain/branched-chain_minevalue.tsv | tsv-join -d 1 -f branched-chain/branched.tigerfam.replace.tsv -k 2 --append-fields 1 \
>branched-chain/kaks/branched.WP.tsv
#提取含有branched蛋白的所有菌株中的CDS名称(2142)
faops size CDS/all.cds.fa | grep -f <(cat branched-chain/kaks/branched.WP.tsv | cut -f 2 ) >branched-chain/kaks/branched-chain.CDS.tsv
#提取含有branched蛋白的所有菌株中的CDS序列(2142)
faops some CDS/all.cds.fa <(cut -f 1 branched-chain/kaks/branched-chain.CDS.tsv) stdout | sed -f CDS/sed.script >branched-chain/kaks/branched-chain.CDS.fa

#提取铜绿假单胞菌中branced的cds簇1名称（382）
faops size branched-chain/kaks/branched-chain.CDS.fa | cut -f 1  | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/kaks/branched-chain_cluster1.tsv) >branched-chain/kaks/branched-chain_cluster1.cds.tsv
#提取铜绿假单胞菌中branced的cds簇2名称（383）,多了一个cds名称，具体情况如下图
faops size branched-chain/kaks/branched-chain.CDS.fa | cut -f 1 | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/kaks/branched-chain_cluster2.tsv) >branched-chain/kaks/branched-chain_cluster2.cds.tsv
sed -i 's/Pseudom_aeru_GCF_001516005_1_cds_WP_003114263.1_10//g;'  branched-chain/kaks/branched-chain_cluster2.cds.tsv
sed -i '/^$/d' branched-chain/kaks/branched-chain_cluster2.cds.tsv

#更改铜绿假单胞菌中branced的cds名称
perl -alne 'm/(.*)\_cds(.*)\.(.*)$/;$cds=$_;$wp=$1.$2; print"$cds\t$wp";' branched-chain/kaks/branched-chain_cluster1.cds.tsv |
tsv-join -d 2 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 >branched-chain/kaks/branched-chain_cluster1_del.cds.tsv
perl -alne 'm/(.*)\_cds(.*)\.(.*)$/;$cds=$_;$wp=$1.$2; print"$cds\t$wp";' branched-chain/kaks/branched-chain_cluster2.cds.tsv |
tsv-join -d 2 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 >branched-chain/kaks/branched-chain_cluster2_del.cds.tsv
#将铜绿假单胞菌中branced的cds的名称替换为wp的名称
cat branched-chain/kaks/branched-chain_cluster1_del.cds.tsv  branched-chain/kaks/branched-chain_cluster2_del.cds.tsv | tsv-select -f 1,2 |
    perl -nla -e '
        print q{s/^>} . quotemeta($F[0]) . q{/>} . quotemeta($F[1]) . q{/g;};
    ' \
    > branched-chain/kaks/sed.script
rm -rf branched-chain/kaks/branched-chain_cluster1_del.cds.tsv  branched-chain/kaks/branched-chain_cluster2_del.cds.tsv

#提取铜绿假单胞菌中branced的cds簇1序列(名称替换为wp)/(382)
faops some branched-chain/kaks/branched-chain.CDS.fa branched-chain/kaks/branched-chain_cluster1.cds.tsv stdout | sed -f branched-chain/kaks/sed.script \
>branched-chain/kaks/branched-chain_cluster1.cds.fa
#提取铜绿假单胞菌中branced的cds簇2序列(名称替换为wp)/(382)
faops some branched-chain/kaks/branched-chain.CDS.fa branched-chain/kaks/branched-chain_cluster2.cds.tsv stdout | sed -f branched-chain/kaks/sed.script \
>branched-chain/kaks/branched-chain_cluster2.cds.fa

```
## 1.3 铜绿假单胞菌copy1准备列表文件，cds和protein
```BASH
#列表文件 (145542*2=291084)
perl -alne 'print"$F[0]\n$F[1]" ' branched-chain/kaks/branched-chain_cluster1.arrange.tsv | 
perl -alne 'BEGIN{%seen;$h;} $h=$_;$seen{$h}++;$name=$h."\_".$seen{$h};print"$name";' >branched-chain/kaks/branched-chain_cluster1.arrange.rename.tsv
sed -i 's/_1$//g' branched-chain/kaks/branched-chain_cluster1.arrange.rename.tsv
#protein序列  branched-chain/kaks/branched-chain_cluster1.fa(382*381*2) 
seqkit duplicate -n 762  branched-chain/kaks/branched-chain_cluster1.fa >branched-chain/kaks/branched-chain_cluster1.seqkit.fa
seqkit rename branched-chain/kaks/branched-chain_cluster1.seqkit.fa >branched-chain/kaks/branched-chain_cluster1.seqkit.rename.fa
faops order branched-chain/kaks/branched-chain_cluster1.seqkit.rename.fa branched-chain/kaks/branched-chain_cluster1.arrange.rename.tsv  branched-chain/kaks/branched-chain_copy1.protein.fa
#cds序列  branched-chain/kaks/branched-chain_cluster1.cds.fa
seqkit duplicate -n 762 branched-chain/kaks/branched-chain_cluster1.cds.fa  >branched-chain/kaks/branched-chain_cluster1.seqkit.cds.fa
seqkit rename branched-chain/kaks/branched-chain_cluster1.seqkit.cds.fa >branched-chain/kaks/branched-chain_cluster1.seqkit.rename.cds.fa
faops order branched-chain/kaks/branched-chain_cluster1.seqkit.rename.cds.fa  branched-chain/kaks/branched-chain_cluster1.arrange.rename.tsv branched-chain/kaks/branched-chain_copy1.cds.fa
```
## 1.4 铜绿假单胞菌copy2准备列表文件，cds和protein
```BASH
#列表文件
perl -alne 'print"$F[0]\n$F[1]" ' branched-chain/kaks/branched-chain_cluster2.arrange.tsv | 
perl -alne 'BEGIN{%seen;$h;} $h=$_;$seen{$h}++;$name=$h."\_".$seen{$h};print"$name";' >branched-chain/kaks/branched-chain_cluster2.arrange.rename.tsv
sed -i 's/_1$//g' branched-chain/kaks/branched-chain_cluster2.arrange.rename.tsv
#protein序列  branched-chain/kaks/branched-chain_cluster2.fa
seqkit duplicate -n 762  branched-chain/kaks/branched-chain_cluster2.fa >branched-chain/kaks/branched-chain_cluster2.seqkit.fa
seqkit rename branched-chain/kaks/branched-chain_cluster2.seqkit.fa >branched-chain/kaks/branched-chain_cluster2.seqkit.rename.fa
faops order branched-chain/kaks/branched-chain_cluster2.seqkit.rename.fa  branched-chain/kaks/branched-chain_cluster2.arrange.rename.tsv  branched-chain/kaks/branched-chain_copy2.protein.fa
#cds序列  branched-chain/kaks/branched-chain_cluster2.cds.fa
seqkit duplicate -n 762 branched-chain/branched-chain_cluster2.cds.fa  >branched-chain/kaks/branched-chain_cluster2.seqkit.cds.fa
seqkit rename branched-chain/kaks/branched-chain_cluster2.seqkit.cds.fa >branched-chain/kaks/branched-chain_cluster2.seqkit.rename.cds.fa
faops order branched-chain/kaks/branched-chain_cluster2.seqkit.rename.cds.fa branched-chain/kaks/branched-chain_cluster2.arrange.rename.tsv  branched-chain/kaks/branched-chain_copy2.cds.fa

```

# 2. 使用paraAT2和caculator2计算铜绿假单胞菌的kaks
## 2.1安装paraAT2
```BASH
#安装paraat2和caculator2
cd ~
wget ftp://download.big.ac.cn/bigd/tools/ParaAT2.0.tar.gz
tar -xzvf ParaAT2.0.tar.gz
vim ~/.bashrc
export PATH="$PATH:/home/syq/ParaAT2.0"
source ~/.bashrc
conda install kaks_caculator
#使用paraAT2计算kaks #proc一定要在文件目录下
-h, 同源基因名称文件
-n, 指定核酸序列文件
-a, 指定蛋白序列文件
-p, 指定多线程文件
-m, 指定比对工具
-g, 去除比对有gap的密码子
-k, 用KaKs_Calculator 计算kaks值
-o, 输出结果的目录
-f, 输出比对文件的格式
-t：移除mismatched codons；
-k：用KaKs_Calculator计算(需要输出axt格式)Ka和Ks，获得axt文件后自动计算kaks值，使用MA模型，比YN模型慢，推荐输出axt后自己用KaKs_Calculator计算并用YN模型
```
## 2.2 计算铜绿假单胞菌copy1的kaks
```bash
#copy1输出axt格式
cd ~/data/Pseudomonas
ParaAT.pl -h branched-chain/kaks/branched-chain_cluster1.arrange.tsv -n branched-chain/kaks/branched-chain_copy1.cds.fa \
-a branched-chain/kaks/branched-chain_copy1.protein.fa -p proc -m muscle -f axt  -o branched-chain/kaks/branched_copy1
#计算copy1的kaks,使用GNG模型(该模型MGEA中有，易于mega中结果比对,且该模型更新时间较近)
cd branched-chain/kaks/branched_copy1
for filename in *.axt
do
base=$(basename $filename .axt)
echo $base
KaKs_Calculator -i ${base}.axt -o ${base}.kaks -m GNG 
done
###由于比对完后的kaks文件太多，共291084，无法查看，需要将大文件夹分成多个小文件夹
python floder_split.py
#分别合并kaks结果
for f in img_{1..30}
do
echo $f
cat $f/*.kaks >$f.kaks
done
#上述结果可直接得到每一对同源基因的ka，ks值，可通过如下命令将其整合(145542)
cat *.kaks | cut -f 1,3,4,5 |grep -v 'Sequence' >../branched_copy1.kaks
```

## 2.3 计算铜绿假单胞菌copy2的kaks
```bash
#copy2输出axt格式
cd ~/data/Pseudomonas
ParaAT.pl -h branched-chain/kaks/branched-chain_cluster2.arrange.tsv -n branched-chain/kaks/branched-chain_copy2.cds.fa \
-a branched-chain/kaks/branched-chain_copy2.protein.fa -p proc -m muscle -f axt  -o branched-chain/kaks/branched_copy2
#计算copy1的kaks,使用GNG模型(该模型MGEA中有，易于mega中结果比对,且该模型更新时间较近)
cd branched-chain/kaks/branched_copy2
for filename in *.axt
do
base=$(basename $filename .axt)
echo $base
KaKs_Calculator -i ${base}.axt -o ${base}.kaks -m GNG 
done
###由于比对完后的kaks文件太多，共291084，无法查看，需要将大文件夹分成多个小文件夹
python floder_split.py
#分别合并kaks结果
for f in img_{1..30}
do
echo $f
cat $f/*.kaks >$f.kaks
done
#上述结果可直接得到每一对同源基因的ka，ks值，可通过如下命令将其整合(145542)
cat *.kaks | cut -f 1,3,4,5 |grep -v 'Sequence' >../branched_copy2.kaks
``` 
## 2.4 铜绿假单胞菌kaks结果画图(不能去除kaks等于0的点)
```bash
cd ~/data/Pseudomonas
#统计kaks等于0的数量(58410和127070)
cat branched-chain/kaks/branched_copy1.kaks | tsv-filter --str-ne 4:"-0" | wc -l 
cat branched-chain/kaks/branched_copy2.kaks | tsv-filter --str-ne 4:"-0" | wc -l 
#统计kaks等于NA的数量(4516和11400)
cat branched-chain/kaks/branched_copy1.kaks | tsv-filter --str-eq 4:"NA" | wc -l
cat branched-chain/kaks/branched_copy2.kaks | tsv-filter --str-eq 4:"NA" | wc -l
#准备绘图的格式(141026和134142)
cat branched-chain/kaks/branched_copy1.kaks | tsv-filter --str-ne 4:"NA" | cut -f 4 >branched-chain/kaks/branched_copy1.kaks.tsv
cat branched-chain/kaks/branched_copy2.kaks | tsv-filter --str-ne 4:"NA" | cut -f 4 >branched-chain/kaks/branched_copy2.kaks.tsv
perl -alne '$num=$F[0];$name=A;print"$num\t$name";' branched-chain/kaks/branched_copy1.kaks.tsv | sed '1ikaks\tfrequency' >branched-chain/kaks/branched_copy1.kaks.picture.tsv
perl -alne '$num=$F[0];$name=B;print"$num\t$name";' branched-chain/kaks/branched_copy2.kaks.tsv >branched-chain/kaks/branched_copy2.kaks.picture.tsv
cat branched-chain/kaks/branched_copy1.kaks.picture.tsv  branched-chain/kaks/branched_copy2.kaks.picture.tsv >branched-chain/kaks/branched.kaks.picture.tsv
sed -i 's/-//g' branched-chain/kaks/branched.kaks.picture.tsv
#使用plotr绘图
plotr hist --xl Ka/Ks --yl Frequency -g 2 --bins 20 --xmm -0.1,1 --ymm 0,0.5 -p branched-chain/kaks/branched.kaks.picture.tsv
#添加平均值+-标准误SEM信息   SEM是standard error of mean是平均数的抽样误差，反应平均数的抽样准确性 使用excel计算
#copy1的ka(83560)  0.002215855+-0.003147
cat branched-chain/kaks/branched_copy1.kaks | tsv-filter --str-ne 2:"NA" | cut -f 2 >branched-chain/kaks/branched_copy1_ka.tsv
#copy1的ks(141026)  0.024234602+-0.028889
cat branched-chain/kaks/branched_copy1.kaks | tsv-filter --str-ne 3:"NA" | cut -f 3 >branched-chain/kaks/branched_copy1_ks.tsv
#copy2的ka(120258)  0.03124638+-0.096849
cat branched-chain/kaks/branched_copy2.kaks | tsv-filter --str-ne 2:"NA" | cut -f 2 >branched-chain/kaks/branched_copy2_ka.tsv
#copy2的ks(134142) 0.084684557+-0.230414
cat branched-chain/kaks/branched_copy2.kaks | tsv-filter --str-ne 3:"NA" | cut -f 3 >branched-chain/kaks/branched_copy2_ks.tsv
```

# 3. 使用paraAT2和caculator2计算其他假单胞菌的kaks(待定)
* Pseudomonas chlororaphis 绿刺假单胞菌，在农业和园艺中用作土壤接种剂的细菌,可以作为生物防治剂,通过生产吩嗪类抗生素来对抗某些真菌植物病原体
* Pseudomonas putida 恶臭假单胞菌 生活在土壤中，多样化代谢可以用于生物修复，生物降解，生物防治
* Pseudomonas syringae 丁香假单胞菌是一种具有极鞭毛的杆状革兰氏阴性细菌。作为一种植物病原体，可以感染多种物种
* Pseudomonas fluorescens 荧光假单胞菌是一种常见的革兰氏阴性杆状细菌，代谢灵活，在土壤和水中发现，是一种专性需氧菌
* Pseudomonas stutzeri 施氏假单胞菌,在人类脊髓液中分离出来
  

## 4. 使用mega计算dn/ds
```BASH
#dnds结果画图
sed -i '1d' mega.copy1.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' mega.copy1.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d'   >mega.copy1.dnds.tsv
perl -alne '$num=$F[0];$name=A;print"$num\t$name";' mega.copy1.dnds.tsv | sed '1idnds\tfrequency' >mega.copy1.dnds.picture.tsv

sed -i '1d' mega.copy2.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' mega.copy2.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d'   >mega.copy2.dnds.tsv
perl -alne '$num=$F[0];$name=B;print"$num\t$name";' mega.copy2.dnds.tsv >mega.copy2.dnds.picture.tsv
cat mega.copy1.dnds.picture.tsv mega.copy2.dnds.picture.tsv >mega.dnds.picture.tsv
#plotr画图
plotr hist --xl dN/dS --yl Frequency -g 2 --bins 20 --xmm -0.1,1 --ymm 0,0.5 -p mega.dnds.picture.tsv
```

## 5. 使用bp_pairwise_kaks计算kaks(bp_pairwise的结果与mega和kaks_calculator不一致，可能是算法的原因)
```BASH
bp_pairwise_kaks -i branched-chain/kaks/branched-chain_copy1.cds.fa -o branched-chain/kaks/bp_pairwise_branched-chain_cluster1.kaks
bp_pairwise_kaks -i branched-chain/kaks/branched-chain_copy1.cds.fa -o branched-chain/kaks/bp_pairwise_branched-chain_cluster2.kaks
cat  branched-chain/kaks/bp_pairwise_branched-chain_cluster1.kaks  | grep -v "Ks" | cut -f 5 | perl -alne '$num=$F[0];$name=A;print"$num\t$name";' | sed '1ikaks\tfrequency' >branched-chain/kaks/bp_pairwise.kaks.copy1.tsv
cat  branched-chain/kaks/bp_pairwise_branched-chain_cluster2.kaks  | grep -v "Ks" | cut -f 5 | perl -alne '$num=$F[0];$name=B;print"$num\t$name";' >branched-chain/kaks/bp_pairwise.kaks.copy2.tsv
cat branched-chain/kaks/bp_pairwise.kaks.copy1.tsv bp_pairwise.kaks.copy2.tsv >branched-chain/kaks/bp_pairwise.kaks.picuture.tsv
#plotr画图
plotr hist --xl Ka/Ks --yl Frequency -g 2 --bins 20 --xmm -0.1,1 --ymm 0,1 -p branched-chain/kaks/bp_pairwise.kaks.picuture.tsv
```

# 5. 基因共线性
## 5.1 准备假单胞菌属的典型菌株的genbank文件，共9种12个
```BASH
cd ~/data/Pseudomonas
for S in \
    Pseudom_aeru_PAO1 \
    Pseudom_puti_KT2440_GCF_000007565_2 \
    Pseudom_chl_aureofaciens_30_84_GCF_000281915_1 \
    Pseudom_entomophi_L48_GCF_000026105_1 \
    Pseudom_fluo_SBW25_GCF_000009225_2 \
    Pseudom_prot_Pf_5_GCF_000012265_1 \
    Pseudom_sav_pv_phaseolicola_1448A_GCF_000012205_1 \
    Pseudom_stu_A1501_GCF_000013785_1 \
    Pseudom_syr_pv_syringae_B728a_GCF_000012245_1 \
    Pseudom_aeru_UCBPP_PA14_GCF_000014625_1 \
    Pseudom_aeru_PA7_GCF_000017205_1 \
    Pseudom_aeru_LESB58_GCF_000026645_1 \
    ; do
    echo ${S}
done \
    > typical.lst

mkdir -p ~/data/Pseudomonas/branched-chain/collinearity/genebank
cd ~/data/Pseudomonas
#提取上述菌株的gbk文件
for name in $(cat typical.lst)
do
echo $name
cp ASSEMBLY/$name/*_genomic.gbff.gz ASSEMBLY/$name/$name.gbff.gz
done

for name in $(cat typical.lst)
do
echo $name
mv ASSEMBLY/$name/$name.gbff.gz branched-chain/collinearity/genebank
done
#解压上述文件
gunzip branched-chain/collinearity/genebank/*.gz

```
## 2.2使用clinker查看基因共线性
```bash
#目前查阅文献可得方法
#1.直接提取目标基因上下游各10000bp的gbk文件
#2.使用clinker查看基因共线性
python  fetch_gbk.py -i Pseudom_aeru_PAO1.gbff -e 10000 -l braB -o PAO1_braB
python  fetch_gbk.py -i Pseudom_aeru_PAO1.gbff -e 10000 -l braZ -o PAO1_braZ
python  fetch_gbk.py -i Pseudom_aeru_MTB_1_GCF_000504045.gbff -e 10000 -l brnQ -o MTB_brnQ
python  fetch_gbk.py -i Pseudom_fluo_SBW25_GCF_000009225_2.gbff -e 10000 -l brnQ -o fluo_brnQ
python  fetch_gbk.py -i Pseudom_psychrop_GCF_011040435_1.gbff -e 10000 -l brnQ -o psychrop_brnQ
#将截取的gbk文件一起画图
clinker *.gbk -p 
#在clinker链接的网页上将group名调整为对应的基因
```

# 3.MEME查找motif

# 4.islandviewer4查看基因岛

# 5.gephi查看pangenome


