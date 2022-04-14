<!-- TOC -->

- [1.计算ka/ks](#1计算kaks)
  - [1.1提取branched在假单胞菌中的两个WP](#11提取branched在假单胞菌中的两个wp)
  - [1.2提取branched在假单胞菌中的两个CDS](#12提取branched在假单胞菌中的两个cds)
  - [1.3 copy1准备列表文件，cds和protein](#13-copy1准备列表文件cds和protein)
  - [1.4 copy2准备列表文件，cds和protein](#14-copy2准备列表文件cds和protein)
  - [1.5使用paraAT2和caculator2计算kaks(以copy1为例)](#15使用paraat2和caculator2计算kaks以copy1为例)
  - [1.6使用mega计算dn/ds](#16使用mega计算dnds)
  - [1.7使用bp_pairwise_kaks计算kaks(以copy1为例)](#17使用bp_pairwise_kaks计算kaks以copy1为例)
- [2.基因共线性](#2基因共线性)
  - [2.1使用easyfig查看基因共线性](#21使用easyfig查看基因共线性)
  - [2.2使用clinker查看基因共线性](#22使用clinker查看基因共线性)
- [3.MEME查找motif](#3meme查找motif)
- [4.islandviewer4查看基因岛](#4islandviewer4查看基因岛)
- [5.gephi查看pangenome](#5gephi查看pangenome)

<!-- /TOC -->

* 1.上一步主要是统计branched在不同菌株基因组中的拷贝数，可知铜绿假单胞菌等中具有两个拷贝    
* 2.分别统计铜绿假单胞菌等中的两个拷贝的ka/ks


# 1.计算ka/ks
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

## 1.1提取branched在假单胞菌中的两个WP
```BASH
#提取假单胞菌的两个拷贝名称
cd ~/data/Pseudomonas
cat branched-chain/branched-chain_hmmscan_copy.pfam.tsv | tsv-filter --ge 3:2 |
tsv-filter --str-in-fld 1:"Pseudom_aeru" | cut -f 1 >branched-chain/branched-chain.Pseudom_aeru.two.copy.tsv
cat branched-chain/branched-chain_minevalue.tsv | grep -f branched-chain/branched-chain.Pseudom_aeru.two.copy.tsv |
cut -f 1 >branched-chain/branched-chain.Pseudom_aeru.protein.tsv
#提取假单胞菌中branced的两个拷贝WP序列
faops some PROTEINS/all.replace.fa branched-chain/branched-chain.Pseudom_aeru.protein.tsv branched-chain/branched-chain.Pseudom_aeru.protein.fa

#使用cd-hit提取出branched-chain/branched-chain.Pseudom_aeru.protein.fa中差异最大的两条序列最为参考序列
cd-hit -i branched-chain/branched-chain.Pseudom_aeru.protein.fa -o branched-chain/branched-chain_refer.fa -c 0.8
#提取参考序列1
echo "Pseudom_aeru_B136_33_GCF_000359505_1_WP_015503046" >branched-chain/branched-chain_refer1.tsv
faops some branched-chain/branched-chain.Pseudom_aeru.protein.fa branched-chain/branched-chain_refer1.tsv branched-chain/branched-chain_refer1.fa
#根据参考序列1提取簇1的序列
diamond makedb --in branched-chain/branched-chain_refer1.fa -d branched-chain/branched-chain_refer1
diamond blastp -d branched-chain/branched-chain_refer1.dmnd -q branched-chain/branched-chain.Pseudom_aeru.protein.fa -o branched-chain/branched-chain_refer1_result
tsv-filter --gt 3:90 branched-chain/branched-chain_refer1_result  | cut -f 1 | sort -n | uniq >branched-chain/branched-chain_cluster1.tsv
#验证簇1,共382正确
cat branched-chain/branched-chain_cluster1.tsv | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 |
cut -f 2 | sort -n | uniq | wc -l 

#wp簇1的序列名称
branched-chain/branched-chain_cluster1.tsv
faops some branched-chain/branched-chain.Pseudom_aeru.protein.fa branched-chain/branched-chain_cluster1.tsv branched-chain/branched-chain_cluster1.fa
#wp簇2的序列名称
cat branched-chain/branched-chain.Pseudom_aeru.protein.tsv | grep -v -f branched-chain/branched-chain_cluster1.tsv \
>branched-chain/branched-chain_cluster2.tsv
faops some branched-chain/branched-chain.Pseudom_aeru.protein.fa branched-chain/branched-chain_cluster2.tsv branched-chain/branched-chain_cluster2.fa

#wp簇1进行排列组合
cat branched-chain/branched-chain_cluster1.fa | mash sketch -k 21 -s 1000 -i -p 8 - -o branched-chain/branched-chain_cluster1.msh
mash dist branched-chain/branched-chain_cluster1.msh branched-chain/branched-chain_cluster1.msh >branched-chain/branched-chain_cluster1.mash
cut -f 1,2  branched-chain/branched-chain_cluster1.mash | tsv-filter --ff-str-ne 1:2 >branched-chain/branched-chain_cluster1.arrange.tsv
wc -l branched-chain/branched-chain_cluster1.arrange.tsv #145542
rm -rf branched-chain/branched-chain_cluster1.msh branched-chain/branched-chain_cluster1.mash
#wp簇2进行排列组合
cat branched-chain/branched-chain_cluster2.fa | mash sketch -k 21 -s 1000 -i -p 8 - -o branched-chain/branched-chain_cluster2.msh
mash dist branched-chain/branched-chain_cluster2.msh branched-chain/branched-chain_cluster2.msh >branched-chain/branched-chain_cluster2.mash
cut -f 1,2  branched-chain/branched-chain_cluster2.mash | tsv-filter --ff-str-ne 1:2 >branched-chain/branched-chain_cluster2.arrange.tsv
wc -l branched-chain/branched-chain_cluster2.arrange.tsv #145542
rm -rf branched-chain/branched-chain_cluster2.msh branched-chain/branched-chain_cluster2.mash
```

## 1.2提取branched在假单胞菌中的两个CDS
```BASH 
#查看branched在菌株中的蛋白名称
cat branched-chain/branched.tigerfam.replace.tsv | grep -f <(cat branched-chain/branched-chain_minevalue.tsv | cut -f 1 ) \
 >branched-chain/branched.WP.tsv
#提取branched在菌株中的CDS名称
faops size CDS/all.cds.fa | grep -f <(cat branched-chain/branched.WP.tsv | cut -f 1 ) >branched-chain/branched-chain.CDS.tsv
#提取branched的CDS序列
faops some CDS/all.cds.fa <(cut -f 1 branched-chain/branched-chain.CDS.tsv) stdout | sed -f CDS/sed.script >branched-chain/branched-chain.CDS.fa
#提取假单胞菌中branced的两个拷贝CDS序列名称
faops size branched-chain/branched-chain.CDS.fa | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/branched-chain.Pseudom_aeru.protein.tsv) |
cut -f 1 >branched-chain/branched-chain.Pseudom_aeru.CDS.tsv
#提取假单胞菌中branced的两个拷贝CDS
faops some branched-chain/branched-chain.CDS.fa branched-chain/branched-chain.Pseudom_aeru.CDS.tsv  branched-chain/branched-chain.Pseudom_aeru.CDS.fa

#cds簇1的序列名称382个
cat branched-chain/branched-chain.Pseudom_aeru.CDS.tsv | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/branched-chain_cluster1.tsv) \
>branched-chain/branched-chain_cluster1.cds.tsv
#cds簇2的序列名称383个,多了一个cds名称，具体情况如下图
cat branched-chain/branched-chain.Pseudom_aeru.CDS.tsv | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/branched-chain_cluster2.tsv) \
>branched-chain/branched-chain_cluster2.cds.tsv
sed -i 's/Pseudom_aeru_GCF_001516005_1_cds_WP_003114263.1_10//g;'  branched-chain/branched-chain_cluster2.cds.tsv
sed -i '/^$/d' branched-chain/branched-chain_cluster2.cds.tsv

#更改cds名称
perl -alne 'm/(.*)\_cds(.*)\.(.*)$/;$cds=$_;$wp=$1.$2; print"$cds\t$wp";' branched-chain/branched-chain_cluster1.cds.tsv |
tsv-join -d 2 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 >branched-chain/branched-chain_cluster1_del.cds.tsv
perl -alne 'm/(.*)\_cds(.*)\.(.*)$/;$cds=$_;$wp=$1.$2; print"$cds\t$wp";' branched-chain/branched-chain_cluster2.cds.tsv |
tsv-join -d 2 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 >branched-chain/branched-chain_cluster2_del.cds.tsv
#将cds的名称替换为wp的名称
cat branched-chain/branched-chain_cluster1_del.cds.tsv  branched-chain/branched-chain_cluster2_del.cds.tsv | tsv-select -f 1,2 |
    perl -nla -e '
        print q{s/^>} . quotemeta($F[0]) . q{/>} . quotemeta($F[1]) . q{/g;};
    ' \
    > branched-chain/sed.script
cat branched-chain/branched-chain.Pseudom_aeru.CDS.fa | sed -f branched-chain/sed.script >branched-chain/branched-chain.Pseudom_aeru.CDS.replace.fa
rm -rf branched-chain/branched-chain_cluster1_del.cds.tsv  branched-chain/branched-chain_cluster2_del.cds.tsv

#提取cds簇1和簇2的序列
faops some branched-chain/branched-chain.Pseudom_aeru.CDS.replace.fa \
<(cat branched-chain/branched-chain_cluster1.tsv  branched-chain/branched-chain_cluster2.tsv) \
branched-chain/branched-chain.Pseudom_aeru.CDS.delete.fa
faops some branched-chain/branched-chain.Pseudom_aeru.CDS.delete.fa  branched-chain/branched-chain_cluster1.tsv branched-chain/branched-chain_cluster1.cds.fa
faops some branched-chain/branched-chain.Pseudom_aeru.CDS.delete.fa  branched-chain/branched-chain_cluster2.tsv branched-chain/branched-chain_cluster2.cds.fa
```

## 1.3 copy1准备列表文件，cds和protein
```BASH
#列表文件 branched-chain/branched-chain_cluster1.arrange.tsv 
perl -alne 'print"$F[0]\n$F[1]" ' branched-chain/branched-chain_cluster1.arrange.tsv | 
perl -alne 'BEGIN{%seen;$h;} $h=$_;$seen{$h}++;$name=$h."\_".$seen{$h};print"$name";' >branched-chain/branched-chain_cluster1.arrange.rename.tsv
sed -i 's/_1$//g' branched-chain/branched-chain_cluster1.arrange.rename.tsv
#protein序列  branched-chain/branched-chain_copy1.protein.fa
seqkit duplicate -n 762  branched-chain/branched-chain_cluster1.fa >branched-chain/branched-chain_cluster1.seqkit.fa
seqkit rename branched-chain/branched-chain_cluster1.seqkit.fa >branched-chain/branched-chain_cluster1.seqkit.rename.fa
faops order branched-chain/branched-chain_cluster1.seqkit.rename.fa branched-chain/branched-chain_cluster1.arrange.rename.tsv  branched-chain/branched-chain_copy1.protein.fa
#cds序列  branched-chain/branched-chain_copy1.cds.fa
seqkit duplicate -n 762 branched-chain/branched-chain_cluster1.cds.fa  >branched-chain/branched-chain_cluster1.seqkit.cds.fa
seqkit rename branched-chain/branched-chain_cluster1.seqkit.cds.fa >branched-chain/branched-chain_cluster1.seqkit.rename.cds.fa
faops order branched-chain/branched-chain_cluster1.seqkit.rename.cds.fa branched-chain/branched-chain_cluster1.arrange.rename.tsv  branched-chain/branched-chain_copy1.cds.fa

```


## 1.4 copy2准备列表文件，cds和protein
```BASH
#列表文件 branched-chain/branched-chain_cluster1.arrange.tsv 
perl -alne 'print"$F[0]\n$F[1]" ' branched-chain/branched-chain_cluster2.arrange.tsv | 
perl -alne 'BEGIN{%seen;$h;} $h=$_;$seen{$h}++;$name=$h."\_".$seen{$h};print"$name";' >branched-chain/branched-chain_cluster2.arrange.rename.tsv
sed -i 's/_1$//g' branched-chain/branched-chain_cluster2.arrange.rename.tsv
#protein序列  branched-chain/branched-chain_copy2.protein.fa
seqkit duplicate -n 762  branched-chain/branched-chain_cluster2.fa >branched-chain/branched-chain_cluster2.seqkit.fa
seqkit rename branched-chain/branched-chain_cluster1.seqkit.fa >branched-chain/branched-chain_cluster1.seqkit.rename.fa
faops order branched-chain/branched-chain_cluster2.seqkit.rename.fa branched-chain/branched-chain_cluster2.arrange.rename.tsv  branched-chain/branched-chain_copy2.protein.fa
#cds序列  branched-chain/branched-chain_copy1.cds.fa
seqkit duplicate -n 762 branched-chain/branched-chain_cluster2.cds.fa  >branched-chain/branched-chain_cluster2.seqkit.cds.fa
seqkit rename branched-chain/branched-chain_cluster2.seqkit.cds.fa >branched-chain/branched-chain_cluster2.seqkit.rename.cds.fa
faops order branched-chain/branched-chain_cluster2.seqkit.rename.cds.fa branched-chain/branched-chain_cluster2.arrange.rename.tsv  branched-chain/branched-chain_copy2.cds.fa

```

```BASH
#多出了一条序列名称
wc -l branched-chain/branched-chain.Pseudom_aeru.CDS.tsv
Pseudom_aeru_GCF_001516005_1_cds_WP_003114263.1_10
Pseudom_aeru_GCF_001516005_1_cds_WP_003113469.1_5330
Pseudom_aeru_GCF_001516005_1_cds_WP_003114263.1_5716
wc -l branched-chain/branched-chain.Pseudom_aeru.protein.tsv 764
Pseudom_aeru_GCF_001516005_1_WP_003113469
Pseudom_aeru_GCF_001516005_1_WP_003114263
#使用mega查看多出来的一条序列名称，发现是Pseudom_aeru_GCF_001516005_1_cds_WP_003114263蛋白出现了重复序列
echo "Pseudom_aeru_PAO1_cds_NP_250281.1_1591
Pseudom_aeru_PAO1_cds_NP_250661.1_1974
Pseudom_aeru_GCF_001516005_1_cds_WP_003114263.1_10
Pseudom_aeru_GCF_001516005_1_cds_WP_003113469.1_5330
Pseudom_aeru_GCF_001516005_1_cds_WP_003114263.1_5716" \
>test.tsv
faops some branched-chain/branched-chain.CDS.fa test.tsv test.fas
```
![pse_aeru](https://github.com/syq12345678/Pseduomonas_HGT/blob/main/branched_three_tree/branech_pse_auer.png)


## 1.5使用paraAT2和caculator2计算kaks(以copy1为例)

```BASH
#安装paraat2
cd ~
wget ftp://download.big.ac.cn/bigd/tools/ParaAT2.0.tar.gz
tar -xzvf ParaAT2.0.tar.gz
vim ~/.bashrc
export PATH="$PATH:/home/syq/ParaAT2.0"
source ~/.bashrc
conda install kaks_caculator
#使用paraAT2计算kaks #proc一定要在文件目录下
ParaAT.pl -h branched-chain/branched-chain_cluster1.arrange.tsv -n branched-chain/branched-chain_copy1.cds.fa \
-a branched-chain/branched-chain_copy1.protein.fa -p proc -m muscle -f axt  -o branched-chain/branched_copy1
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
#比对完后计算kaks
cd branched-chain/branched_copy1
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
#上述结果可直接得到每一对同源基因的ka，ks值，可通过如下命令将其整合
cat *.kaks | cut -f 1,3,4,5 |grep -v 'Sequence' >branched_copy1.kaks
#同上得branched_copy2.kaks

#kaks结果画图
cat branched_copy1.kaks | tsv-filter --str-ne 4:"-0" | tsv-filter --str-ne 4:"NA" | cut -f 4 >branched_copy1.kaks.tsv
cat branched_copy2.kaks | tsv-filter --str-ne 4:"-0" | tsv-filter --str-ne 4:"NA" | cut -f 4 >branched_copy2.kaks.tsv
perl -alne '$num=$F[0];$name=A;print"$num\t$name";' branched_copy1.kaks.tsv | sed '1ikaks\tfrequency' >branched_copy1.kaks.picture.tsv
perl -alne '$num=$F[0];$name=B;print"$num\t$name";' branched_copy2.kaks.tsv >branched_copy2.kaks.picture.tsv
cat branched_copy1.kaks.picture.tsv branched_copy2.kaks.picture.tsv >branched.kaks.picture.tsv

plotr hist --xl Ka/Ks --yl Frequency -g 2 --bins 20 --xmm 0,1 --ymm 0,1 -p branched.kaks.picture.tsv
```

## 1.6使用mega计算dn/ds
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

#plotr画图注意事项
plotr hist --xl dN/dS --yl Frequency -g 2 --bins 20 --xmm 0,1 --ymm 0,1 -p mega.dnds.picture.tsv
```

## 1.7使用bp_pairwise_kaks计算kaks(以copy1为例)
```BASH
bp_pairwise_kaks -i branched-chain_copy1.cds.fa -o branched_copy1.kaks
#去除序列的回车符
cat  branched-chain_cluster1.cds.fa | perl -p -e 's/\r?\n//;s/^>(.+)$/>$1\n/;s/^>/\n>/'

#kaks结果画图
cat branched_copy1.kaks | tsv-filter --str-ne 4:"-0" | tsv-filter --str-ne 4:"NA" | cut -f 4 >branched_copy1.kaks.tsv
cat branched_copy2.kaks | tsv-filter --str-ne 4:"-0" | tsv-filter --str-ne 4:"NA" | cut -f 4 >branched_copy2.kaks.tsv
perl -alne '$num=$F[0];$name=A;print"$num\t$name";' branched_copy1.kaks.tsv | sed '1ikaks\tfrequency' >branched_copy1.kaks.picture.tsv
perl -alne '$num=$F[0];$name=B;print"$num\t$name";' branched_copy2.kaks.tsv >branched_copy2.kaks.picture.tsv
cat branched_copy1.kaks.picture.tsv branched_copy2.kaks.picture.tsv >branched.kaks.picture.tsv

plotr hist --xl Ka/Ks --yl Frequency -g 2 --bins 20 --xmm 0,1 --ymm 0,1 -p branched.kaks.picture.tsv
```
# 2.基因共线性
## 2.1使用easyfig查看基因共线性
* 选择模式菌株和代表菌株中的铜绿假单胞菌和非铜绿假单胞菌进行查看
```BASH
cat branched-chain/branched-chain_minevalue.tsv | grep -f easyfig/typical.lst  | cut -f 1 >easyfig/branched.typical.pro.tsv

mkdir ~/data/Pseudomonas/easyfig
cd ~/data/Pseudomonas/easyfig

for S in \
    Pseudom_aeru_PAO1 \
    Pseudom_aeru_DK2_GCF_000271365_1\
    Pseudom_aeru_MTB_1_GCF_000504045_1\
    Pseudom_aeru_SCV20265_GCF_000510305_1\
    Pseudom_aeru_UCBPP_PA14_GCF_000014625_1 \
    ; do
    echo ${S}
done \
    > test.lst

mkdir test
cat test.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        compgen -G "../ASSEMBLY/{}/*_genomic.gbff.gz"
    ' |
    paste - - \
    > test/genbank.list
for G in $(cat test/genbank.list);do
    cp $G test
done
gzip -d test/*.gz


cat test.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        compgen -G "../ASSEMBLY/{}/*_genomic.fna.gz" |
        grep -v "from"
    ' |
    paste - - \
    > test/genome.list
for G in $(cat test/genome.list);do
    cp $G test
done
gzip -d test/*.gz
```
## 2.2使用clinker查看基因共线性
```bash
#######!!!!!! 该方法目前行不通
#目前查阅文献可得方法
#1.使用seqkit或者bedtools或者faops size 截取braZ的上下游各1000bp
#2.然后使用prokka注释braz及其上下游序列(prokka使用的是hmmscan和blast)
#3.使用clinker查看基因共线性
perl -alne 'next until(/gene/);$_=~s/\(//g;$_=~s/\)//g;print"$_";' Pseudom_aeru_PAO1.gbff | uniq >PAO1.gbff
#查看PAO1的braz 
grep -B 1 braZ 1.gbff | perl -alne '/(\d*)\.\.(\d*)/;$sta=$1;$end=$2;print"$sta:$end";' | uniq 
#结果
2151755:2153068
#提取PAO1的braz
seqkit subseq  Pseudom_aeru_PAO1.fna -r 2141755:2163068 >PAO1_braz.fa
#使用prokka或者dfast注释的结果与gbff文件中的注释结果不符合,舍弃
```
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


