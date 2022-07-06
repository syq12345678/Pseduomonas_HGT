<!-- TOC -->

- [1. kaks的文件准备 （铜绿假单胞菌的拷贝信息需要修正，不能用)!!!!!](#1-kaks的文件准备-铜绿假单胞菌的拷贝信息需要修正不能用)
  - [1.1 分别提取branched在铜绿假单胞菌中的两个拷贝蛋白序列](#11-分别提取branched在铜绿假单胞菌中的两个拷贝蛋白序列)
  - [1.2提取branched在铜绿假单胞菌中的两个CDS](#12提取branched在铜绿假单胞菌中的两个cds)
  - [1.3 铜绿假单胞菌copy1准备列表文件，cds和protein](#13-铜绿假单胞菌copy1准备列表文件cds和protein)
  - [1.4 铜绿假单胞菌copy2准备列表文件，cds和protein](#14-铜绿假单胞菌copy2准备列表文件cds和protein)
- [2. 使用paraAT2和caculator2计算铜绿假单胞菌的kaks （铜绿假单胞菌的拷贝信息需要修正，不能用)](#2-使用paraat2和caculator2计算铜绿假单胞菌的kaks-铜绿假单胞菌的拷贝信息需要修正不能用)
  - [2.1 安装paraAT2](#21-安装paraat2)
  - [2.2 计算铜绿假单胞菌copy1的kaks](#22-计算铜绿假单胞菌copy1的kaks)
  - [2.3 计算铜绿假单胞菌copy2的kaks](#23-计算铜绿假单胞菌copy2的kaks)
  - [2.4 铜绿假单胞菌kaks结果画图(不能去除kaks等于0的点)](#24-铜绿假单胞菌kaks结果画图不能去除kaks等于0的点)
- [3. 使用mega计算铜绿假单胞菌的dn/ds （铜绿假单胞菌的拷贝信息需要修正，不能用)](#3-使用mega计算铜绿假单胞菌的dnds-铜绿假单胞菌的拷贝信息需要修正不能用)
- [4. 使用bp_pairwise_kaks计算铜绿假单胞菌的kaks(bp_pairwise的结果与mega和kaks_calculator不一致，可能是算法的原因)  （铜绿假单胞菌的拷贝信息需要修正，不能用)](#4-使用bp_pairwise_kaks计算铜绿假单胞菌的kaksbp_pairwise的结果与mega和kaks_calculator不一致可能是算法的原因--铜绿假单胞菌的拷贝信息需要修正不能用)
- [5.重新修正的铜绿假单胞菌拷贝文件](#5重新修正的铜绿假单胞菌拷贝文件)
  - [5.1重新修正铜绿假单胞菌的拷贝蛋白序列](#51重新修正铜绿假单胞菌的拷贝蛋白序列)
  - [5.2重新修正铜绿假单胞菌的拷贝CDS序列](#52重新修正铜绿假单胞菌的拷贝cds序列)
  - [5.3使用mega计算铜绿假单胞菌的braz和braB的dn/ds](#53使用mega计算铜绿假单胞菌的braz和brab的dnds)
  - [5.4提取临床环境的序列，非临床环境的序列](#54提取临床环境的序列非临床环境的序列)
  - [5.5 使用mega计算铜绿假单胞菌临床环境的braB和braZ，非临床环境的braB和braZ的dN和dS](#55-使用mega计算铜绿假单胞菌临床环境的brab和braz非临床环境的brab和braz的dn和ds)
- [6. 使用mega计算其他假单胞菌的kaks](#6-使用mega计算其他假单胞菌的kaks)
  - [6.1分别提取branched在绿刺假单胞菌中的CDS序列](#61分别提取branched在绿刺假单胞菌中的cds序列)
  - [6.2分别提取branched在丁香假单胞菌中的CDS序列](#62分别提取branched在丁香假单胞菌中的cds序列)
  - [6.3分别提取branched在荧光假单胞菌中的CDS序列](#63分别提取branched在荧光假单胞菌中的cds序列)
  - [6.4分别提取branched在恶臭假单胞菌中的CDS序列](#64分别提取branched在恶臭假单胞菌中的cds序列)
  - [6.5分别提取branched在施式假单胞菌中的CDS序列](#65分别提取branched在施式假单胞菌中的cds序列)
- [7. 基因共线性](#7-基因共线性)
  - [7.1 准备假单胞菌属的典型菌株的genbank文件，共9种12个](#71-准备假单胞菌属的典型菌株的genbank文件共9种12个)
  - [7.2 查看铜绿假单胞菌物种内的基因共线性](#72-查看铜绿假单胞菌物种内的基因共线性)
  - [7.3 查看假单胞菌属内物种间的基因共线性](#73-查看假单胞菌属内物种间的基因共线性)
  - [7.4 13个两个braZ拷贝的共线性分析（gbk转gff或者gff转gbk）](#74-13个两个braz拷贝的共线性分析gbk转gff或者gff转gbk)
  - [7.5 查看9个单拷贝铜绿假单胞菌的共线性](#75-查看9个单拷贝铜绿假单胞菌的共线性)
  - [7.6查看假单胞菌属的基因树两个clade的物种分别与braB和braZ的共线性](#76查看假单胞菌属的基因树两个clade的物种分别与brab和braz的共线性)
  - [7.6查看假单胞菌属的基因树两个clade的物种分别与braZ上游基因PA1064的共线性](#76查看假单胞菌属的基因树两个clade的物种分别与braz上游基因pa1064的共线性)

<!-- /TOC -->


# 1. kaks的文件准备 （铜绿假单胞菌的拷贝信息需要修正，不能用)!!!!!
##  1.1 分别提取branched在铜绿假单胞菌中的两个拷贝蛋白序列
```BASH
#提取铜绿假单胞菌的两个拷贝的蛋白名称(382*2=764)  （注意还有9个铜绿假单胞菌是单拷贝)
cd ~/data/Pseudomonas
mkdir -p ~/data/Pseudomonas/branched-chain/kaks
cat branched-chain/branched-chain_hmmscan_copy.pfam.tsv | tsv-filter --ge 3:2 | grep "Pseudom_aeru" | cut -f 1 >branched-chain/kaks/branched-chain.Pseudom_aeru.two.copy.tsv
cat branched-chain/branched-chain_minevalue.pfam.tsv | grep -f branched-chain/kaks/branched-chain.Pseudom_aeru.two.copy.tsv |
cut -f 1 >branched-chain/kaks/branched-chain.Pseudom_aeru.protein.tsv
#提取铜绿假单胞菌中branced的两个拷贝蛋白序列(764)
faops some PROTEINS/all.replace.fa branched-chain/kaks/branched-chain.Pseudom_aeru.protein.tsv branched-chain/kaks/branched-chain.Pseudom_aeru.protein.fa

#使用cd-hit提取出branched-chain/branched-chain.Pseudom_aeru.protein.fa中差异最大的两条序列最为参考序列
cd-hit -i branched-chain/kaks/branched-chain.Pseudom_aeru.protein.fa -o branched-chain/kaks/branched-chain_refer.fa -c 0.8
#提取参考序列1
faops some branched-chain/kaks/branched-chain.Pseudom_aeru.protein.fa <(faops size branched-chain/kaks/branched-chain_refer.fa | cut -f 1 | head  -n 1 ) branched-chain/kaks/branched-chain_refer1.fa
#根据参考序列1提取簇1的序列(382)
diamond makedb --in branched-chain/kaks/branched-chain_refer1.fa -d branched-chain/kaks/branched-chain_refer1
diamond blastp -d branched-chain/kaks/branched-chain_refer1.dmnd -q branched-chain/kaks/branched-chain.Pseudom_aeru.protein.fa -o branched-chain/kaks/branched-chain_refer1_result
tsv-filter --gt 3:90 branched-chain/kaks/branched-chain_refer1_result  | cut -f 1 | sort -n | uniq >branched-chain/kaks/branched-chain_cluster1.tsv
#验证簇1,共382正确
cat branched-chain/kaks/branched-chain_cluster1.tsv | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 | cut -f 2 | sort -n | uniq | wc -l 

#蛋白簇1的序列名称和序列)(382)
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

## 1.2提取branched在铜绿假单胞菌中的两个CDS    
```BASH
#将branched-chain_minevalue.tsv结果改成replace文件格式(2140)
cut -f 1 branched-chain/branched-chain_minevalue.pfam.tsv | tsv-join -d 1 -f branched-chain/branched.tigerfam.replace.tsv -k 2 --append-fields 1 \
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

# 2. 使用paraAT2和caculator2计算铜绿假单胞菌的kaks （铜绿假单胞菌的拷贝信息需要修正，不能用)
## 2.1 安装paraAT2
```bash
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
```
# 3. 使用mega计算铜绿假单胞菌的dn/ds （铜绿假单胞菌的拷贝信息需要修正，不能用)
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

# 4. 使用bp_pairwise_kaks计算铜绿假单胞菌的kaks(bp_pairwise的结果与mega和kaks_calculator不一致，可能是算法的原因)  （铜绿假单胞菌的拷贝信息需要修正，不能用)
```BASH
bp_pairwise_kaks -i branched-chain/kaks/branched-chain_copy1.cds.fa -o branched-chain/kaks/bp_pairwise_branched-chain_cluster1.kaks
bp_pairwise_kaks -i branched-chain/kaks/branched-chain_copy1.cds.fa -o branched-chain/kaks/bp_pairwise_branched-chain_cluster2.kaks
cat  branched-chain/kaks/bp_pairwise_branched-chain_cluster1.kaks  | grep -v "Ks" | cut -f 5 | perl -alne '$num=$F[0];$name=A;print"$num\t$name";' | sed '1ikaks\tfrequency' >branched-chain/kaks/bp_pairwise.kaks.copy1.tsv
cat  branched-chain/kaks/bp_pairwise_branched-chain_cluster2.kaks  | grep -v "Ks" | cut -f 5 | perl -alne '$num=$F[0];$name=B;print"$num\t$name";' >branched-chain/kaks/bp_pairwise.kaks.copy2.tsv
cat branched-chain/kaks/bp_pairwise.kaks.copy1.tsv bp_pairwise.kaks.copy2.tsv >branched-chain/kaks/bp_pairwise.kaks.picuture.tsv
#plotr画图
plotr hist --xl Ka/Ks --yl Frequency -g 2 --bins 20 --xmm -0.1,1 --ymm 0,1 -p branched-chain/kaks/bp_pairwise.kaks.picuture.tsv
```

# 5.重新修正的铜绿假单胞菌拷贝文件
## 5.1重新修正铜绿假单胞菌的拷贝蛋白序列
```bash
cd ~/data/Pseudomonas
#提取所有的铜绿假单胞菌的单双拷贝蛋白序列(773)
mkdir -p ~/data/Pseudomonas/branched-chain/kaks/aeru
cat branched-chain/branched-chain_hmmscan_copy.pfam.tsv | grep "Pseudom_aeru" | cut -f 1 >branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.two.copy.tsv
faops some PROTEINS/all.replace.fa <(cat branched-chain/branched-chain_minevalue.pfam.tsv |grep -f branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.two.copy.tsv | cut -f 1 ) branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.protein.fa
#将PAO1中的braB和braZ分别作为簇1和簇2的参考序列
faops some branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.protein.fa <(faops size branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.protein.fa | grep "PAO1" | grep -v "GCF" | cut -f 1 | head -n 1) branched-chain/kaks/aeru/branched-chain_refer1.fa
faops some branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.protein.fa <(faops size branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.protein.fa | grep "PAO1" | grep -v "GCF" | cut -f 1 | tail -n 1) branched-chain/kaks/aeru/branched-chain_refer2.fa
#根据参考序列1提取簇1的序列(375)
diamond makedb --in branched-chain/kaks/aeru/branched-chain_refer1.fa -d branched-chain/kaks/aeru/branched-chain_refer1
diamond blastp -d branched-chain/kaks/aeru/branched-chain_refer1.dmnd -q branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.protein.fa -o branched-chain/kaks/aeru/branched-chain_refer1_result
tsv-filter --gt 3:90 branched-chain/kaks/aeru/branched-chain_refer1_result | cut -f 1 | sort -n | uniq | wc -l #375
tsv-filter --gt 3:80 branched-chain/kaks/aeru/branched-chain_refer1_result | cut -f 1 | sort -n | uniq | wc -l #375
tsv-filter --gt 3:70 branched-chain/kaks/aeru/branched-chain_refer1_result | cut -f 1 | sort -n | uniq | wc -l #375
tsv-filter --gt 3:70 branched-chain/kaks/aeru/branched-chain_refer1_result | cut -f 1 | sort -n | uniq >branched-chain/kaks/aeru/branched-chain_cluster1.tsv
faops some branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.protein.fa branched-chain/kaks/aeru/branched-chain_cluster1.tsv  branched-chain/kaks/aeru/branched-chain_cluster1.fa

#根据参考序列2提取簇2的序列(398)
diamond makedb --in branched-chain/kaks/aeru/branched-chain_refer2.fa -d branched-chain/kaks/aeru/branched-chain_refer2
diamond blastp -d branched-chain/kaks/aeru/branched-chain_refer2.dmnd -q branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.protein.fa -o branched-chain/kaks/aeru/branched-chain_refer2_result
tsv-filter --gt 3:90 branched-chain/kaks/aeru/branched-chain_refer2_result | cut -f 1 | sort -n | uniq | wc -l #(384)
tsv-filter --gt 3:80 branched-chain/kaks/aeru/branched-chain_refer2_result | cut -f 1 | sort -n | uniq | wc -l #(398)
tsv-filter --gt 3:70 branched-chain/kaks/aeru/branched-chain_refer2_result | cut -f 1 | sort -n | uniq | wc -l #(398)
tsv-filter --gt 3:70 branched-chain/kaks/aeru/branched-chain_refer2_result | cut -f 1 | sort -n | uniq >branched-chain/kaks/aeru/branched-chain_cluster2.tsv
faops some branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.protein.fa branched-chain/kaks/aeru/branched-chain_cluster2.tsv  branched-chain/kaks/aeru/branched-chain_cluster2.fa

#统计braB与braB(簇1相似性)（140625)
mkdir -p branched-chain/blast
cp branched-chain/kaks/aeru/branched-chain_cluster1.fa branched-chain/blast
mv branched-chain/blast/branched-chain_cluster1.fa  branched-chain/blast/cluster1.fa
makeblastdb -in branched-chain/blast/cluster1.fa -dbtype prot -out branched-chain/blast/cluster1
blastp -query branched-chain/blast/cluster1.fa -out branched-chain/blast/cluster1_cluster1.tsv -db branched-chain/blast/cluster1 -outfmt 6 -evalue 1e-5

#统计braZ与braz(簇2相似性) (158404)
cp branched-chain/kaks/aeru/branched-chain_cluster2.fa branched-chain/blast
mv branched-chain/blast/branched-chain_cluster2.fa  branched-chain/blast/cluster2.fa
makeblastdb -in branched-chain/blast/cluster2.fa -dbtype prot -out branched-chain/blast/cluster2
blastp -query branched-chain/blast/cluster2.fa -out branched-chain/blast/cluster2_cluster2.tsv -db branched-chain/blast/cluster2 -outfmt 6 -evalue 1e-5

#统计braZ与braB之间相似性 (149250)
blastp -query branched-chain/blast/cluster1.fa -out branched-chain/blast/cluster1_cluster2.tsv -db branched-chain/blast/cluster2 -outfmt 6 -evalue 1e-5

cat branched-chain/blast/cluster1_cluster1.tsv | perl -alne '$identity=(split/\t/,$_)[2];$name="braB_braB";print"$name\t$identity";' | sed '1itype\tidentity' >branched-chain/blast/braB_braB.tsv
cat branched-chain/blast/cluster2_cluster2.tsv | perl -alne '$identity=(split/\t/,$_)[2];$name="braZ_braZ";print"$name\t$identity";' >branched-chain/blast/braZ_braZ.tsv
cat branched-chain/blast/cluster1_cluster2.tsv | perl -alne '$identity=(split/\t/,$_)[2];$name="braB_braZ";print"$name\t$identity";' >branched-chain/blast/braB_braZ.tsv
#合并文件
cat branched-chain/blast/braB_braB.tsv branched-chain/blast/braB_braZ.tsv   branched-chain/blast/braZ_braZ.tsv >branched-chain/blast/whole.identity.picture.tsv
#绘制相似性的箱线图
library(ggplot2)
data<-read.table("whole.identity.picture.tsv",header=T)
p<-ggplot(data,aes(x=type,y=identity,col=type),show.legend = F)+ geom_violin()+ geom_boxplot()+theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+theme(text=element_text(size=16,family="Arial",face="bold"))
ggsave('pse_aeru_identity.png', p,dpi = 480, width=8, height=6)   
```
## 5.2重新修正铜绿假单胞菌的拷贝CDS序列    
```BASH
cd ~/data/Pseudomonas
#将branched-chain_minevalue.tsv结果改成replace文件格式(2140)
cut -f 1 branched-chain/branched-chain_minevalue.pfam.tsv | tsv-join -d 1 -f branched-chain/branched.tigerfam.replace.tsv -k 2 --append-fields 1 \
>branched-chain/kaks/aeru/branched.WP.tsv
#提取含有branched蛋白的所有菌株中的CDS名称(2142)
faops size CDS/all.cds.fa | grep -f <(cat branched-chain/kaks/aeru/branched.WP.tsv | cut -f 2 ) >branched-chain/kaks/aeru/branched-chain.CDS.tsv
#提取含有branched蛋白的所有菌株中的CDS序列(2142)
faops some CDS/all.cds.fa <(cut -f 1 branched-chain/kaks/aeru/branched-chain.CDS.tsv) stdout | sed -f CDS/sed.script >branched-chain/kaks/aeru/branched-chain.CDS.fa

#提取铜绿假单胞菌中branced的cds簇1名称和序列（376）
faops size branched-chain/kaks/aeru/branched-chain.CDS.fa | cut -f 1  | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/kaks/aeru/branched-chain_cluster1.tsv) >branched-chain/kaks/aeru/branched-chain_cluster1.cds.tsv
#sed -i 's/Pseudom_aeru_GCF_001516005_1_cds_WP_003114263.1_10//g;'  branched-chain/kaks/aeru/branched-chain_cluster1.cds.tsv
#sed -i '/^$/d' branched-chain/kaks/aeru/branched-chain_cluster1.cds.tsv
faops some branched-chain/kaks/aeru/branched-chain.CDS.fa branched-chain/kaks/aeru/branched-chain_cluster1.cds.tsv branched-chain/kaks/aeru/branched-chain_cluster1.pse_aeru.cds.fa
#提取铜绿假单胞菌中branced的cds簇2名称（398）
faops size branched-chain/kaks/aeru/branched-chain.CDS.fa | cut -f 1 | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/kaks/aeru/branched-chain_cluster2.tsv) >branched-chain/kaks/aeru/branched-chain_cluster2.cds.tsv
faops some branched-chain/kaks/aeru/branched-chain.CDS.fa branched-chain/kaks/aeru/branched-chain_cluster2.cds.tsv branched-chain/kaks/aeru/branched-chain_cluster2.pse_aeru.cds.fa
```    
## 5.3使用mega计算铜绿假单胞菌的braz和braB的dn/ds
```bash
cd ~/data/Pseudomonas
mkdir -p ~/data/Pseudomonas/branched-chain/kaks/aeru/mega
#计算braB的dn/ds(64668) #64303
sed -i '1d' branched-chain/kaks/aeru/mega/pseudom_auer-braB_dN-dS.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";'  branched-chain/kaks/aeru/mega/pseudom_auer-braB_dN-dS.tsv | sed 's/NA//g' | sed  '/^\s*$/d' | perl -alne '$num=$F[0];$name=braB;print"$num\t$name";' | sed '1idnds\tfrequency' >branched-chain/kaks/aeru/mega/mega_braB_dnds.tsv
#计算braB的dn(70500) #70125
sed -i '1d' branched-chain/kaks/aeru/mega/pseudom_auer-braB_dN.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/aeru/mega/pseudom_auer-braB_dN.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/aeru/mega/mega_braB_dn.tsv
#计算braB的ds(70500) #70125
sed -i '1d' branched-chain/kaks/aeru/mega/pseudom_auer-braB_dS.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/aeru/mega/pseudom_auer-braB_dS.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/aeru/mega/mega_braB_ds.tsv
#查看braB的dn的平均值和平均值标准误差
library(plotrix)
setwd("~/data/Pseudomonas")
data<-read.table("branched-chain/kaks/aeru/mega/mega_braB_dn.tsv")
data1<-data$V1
mean(data1)  # 0.00224695
std.error(data1) # 0.00001016217
#查看braB的ds的平均值和平均值标准误差
data<-read.table("branched-chain/kaks/aeru/mega/mega_braB_ds.tsv")
data1<-data$V1
mean(data1)  # 0.02046081
std.error(data1) # 0.0002066775


#计算braZ的dn/ds(76720)
sed -i '1d' branched-chain/kaks/aeru/mega/pseudom_auer-braZ_dN-dS.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/aeru/mega/pseudom_auer-braZ_dN-dS.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d' | perl -alne '$num=$F[0];$name=braZ;print"$num\t$name";'  >branched-chain/kaks/aeru/mega/mega_braZ_dnds.tsv
#计算braZ的dn(79003)
sed -i '1d' branched-chain/kaks/aeru/mega/pseudom_auer-braZ_dN.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/aeru/mega/pseudom_auer-braZ_dN.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/aeru/mega/mega_braZ_dn.tsv
#计算braZ的ds(79003)
sed -i '1d' branched-chain/kaks/aeru/mega/pseudom_auer-braZ_dS.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/aeru/mega/pseudom_auer-braZ_dS.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/aeru/mega/mega_braZ_ds.tsv
#查看braZ的dn的平均值和平均值标准误差
library(plotrix)
setwd("~/data/Pseudomonas")
data<-read.table("~/data/Pseudomonas/branched-chain/kaks/aeru/mega/mega_braZ_dn.tsv")
data1<-data$V1
mean(data1)  # 0.01099577
std.error(data1) # 0.0001283826
#查看braZ的ds的平均值和平均值标准误差
data<-read.table("~/data/Pseudomonas/branched-chain/kaks/aeru/mega/mega_braZ_ds.tsv")
data1<-data$V1
mean(data1)  # 0.05826713
std.error(data1) # 0.0004722907

#合并braB和braZ 141389
cat branched-chain/kaks/aeru/mega/mega_braB_dnds.tsv branched-chain/kaks/aeru/mega/mega_braZ_dnds.tsv >branched-chain/kaks/aeru/mega/mega_picture.tsv
#plotr画图
plotr hist --xl dN/dS --yl Frequency -g 2 --bins 20 --xmm -0.1,1.0 --ymm 0,0.5 -p branched-chain/kaks/aeru/mega/mega_picture.tsv
#查看braB和braZ的dnds为0的数量，方便判断矩形颜色
 cat  branched-chain/kaks/aeru/mega/mega_picture.tsv | grep -v "dnds" | tsv-filter --le 1:0.1 | grep "braB" | wc -l #19659
 cat branched-chain/kaks/aeru/mega/mega_picture.tsv | grep -v "dnds" | tsv-filter --le 1:0.1 | grep "braZ" | wc -l #59071

```


## 5.4提取临床环境的序列，非临床环境的序列
```bash
cd ~/data/Pseudomonas
mkdir -p branched-chain/kaks/biosample
wc -l biosample/biosample_env_anno.tsv  #391
#查看生境为临床的铜绿菌株数
cat biosample/biosample_env_anno.tsv | grep "clinical" | wc -l #334
#生境为临床的铜绿菌株蛋白 660
cat branched-chain/branched-chain_minevalue.pfam.tsv | cut -f 1 | grep -f <(cat biosample/biosample_env_anno.tsv | grep "clinical" | cut -f 2) >branched-chain/kaks/biosample/biosample_clincal_pro.tsv
#提取对应的cds序列名字661
faops size branched-chain/kaks/aeru/branched-chain.CDS.fa | cut -f 1  | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/kaks/biosample/biosample_clincal_pro.tsv) >branched-chain/kaks/biosample/biosample_clincal_cds.tsv
#簇1braB的序列数376
wc -l branched-chain/kaks/aeru/branched-chain_cluster1.cds.tsv
#簇2braZ的序列数398
wc -l branched-chain/kaks/aeru/branched-chain_cluster2.cds.tsv
#提取生境为临床的braB的cds序列 321
cat branched-chain/kaks/aeru/branched-chain_cluster1.cds.tsv | grep -f branched-chain/kaks/biosample/biosample_clincal_cds.tsv | wc -l 
faops some branched-chain/kaks/aeru/branched-chain.CDS.fa <(cat branched-chain/kaks/aeru/branched-chain_cluster1.cds.tsv | grep -f branched-chain/kaks/biosample/biosample_clincal_cds.tsv ) branched-chain/kaks/biosample/biosample_clincal_braB.fa #321
#提取生境为临床的braZ的cds序列 340
cat branched-chain/kaks/aeru/branched-chain_cluster2.cds.tsv | grep -f branched-chain/kaks/biosample/biosample_clincal_cds.tsv | wc -l 
faops some branched-chain/kaks/aeru/branched-chain.CDS.fa <(cat branched-chain/kaks/aeru/branched-chain_cluster2.cds.tsv | grep -f branched-chain/kaks/biosample/biosample_clincal_cds.tsv ) branched-chain/kaks/biosample/biosample_clincal_braZ.fa #340

#生境为非临床的铜绿菌株数 
cat biosample/biosample_env_anno.tsv | grep -v "unknow" | grep -v "clinical" | wc -l #42
#生境为非临床的铜绿蛋白
cat branched-chain/branched-chain_minevalue.pfam.tsv | cut -f 1| grep -f <(cat biosample/biosample_env_anno.tsv | grep -v "unknow" | grep -v "clinical" | cut -f 2) >branched-chain/kaks/biosample/biosample_other_pro.tsv   #84
#提取对应的cds序列名字84
faops size branched-chain/kaks/aeru/branched-chain.CDS.fa | cut -f 1  | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/kaks/biosample/biosample_other_pro.tsv) >branched-chain/kaks/biosample/biosample_other_cds.tsv
#簇1braB的序列数376
wc -l branched-chain/kaks/aeru/branched-chain_cluster1.cds.tsv
#簇2braZ的序列数398
wc -l branched-chain/kaks/aeru/branched-chain_cluster2.cds.tsv
#提取生境为非临床的braB的cds序列 41
cat branched-chain/kaks/aeru/branched-chain_cluster1.cds.tsv | grep -f branched-chain/kaks/biosample/biosample_other_cds.tsv | wc -l 
faops some branched-chain/kaks/aeru/branched-chain.CDS.fa <(cat branched-chain/kaks/aeru/branched-chain_cluster1.cds.tsv | grep -f branched-chain/kaks/biosample/biosample_other_cds.tsv ) branched-chain/kaks/biosample/biosample_other_braB.fa #41
#提取生境为非临床的braZ的cds序列 43
cat branched-chain/kaks/aeru/branched-chain_cluster2.cds.tsv | grep -f branched-chain/kaks/biosample/biosample_other_cds.tsv | wc -l 
faops some branched-chain/kaks/aeru/branched-chain.CDS.fa <(cat branched-chain/kaks/aeru/branched-chain_cluster2.cds.tsv | grep -f branched-chain/kaks/biosample/biosample_other_cds.tsv ) branched-chain/kaks/biosample/biosample_other_braZ.fa #43
```

## 5.5 使用mega计算铜绿假单胞菌临床环境的braB和braZ，非临床环境的braB和braZ的dN和dS
```bash
cd ~/data/Pseudomonas
mkdir -p branched-chain/kaks/biosample
#生境为临床的
#计算braB的dn/ds(47433)
sed -i '1d' branched-chain/kaks/biosample/clinical_braB_dnds.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";'  branched-chain/kaks/biosample/clinical_braB_dnds.tsv | sed 's/NA//g' | sed  '/^\s*$/d' | perl -alne '$num=$F[0];$name=braB;print"$num\t$name";' | sed '1idnds\tfrequency' >branched-chain/kaks/biosample/mega_clinical_braB_dnds.tsv

#计算braZ的dn/ds(55652)
sed -i '1d' branched-chain/kaks/biosample/clinical_braZ_dnds.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";'  branched-chain/kaks/biosample/clinical_braZ_dnds.tsv | sed 's/NA//g' | sed  '/^\s*$/d' | perl -alne '$num=$F[0];$name=braZ;print"$num\t$name";' | sed '1idnds\tfrequency' >branched-chain/kaks/biosample/mega_clinical_braZ_dnds.tsv

#合并braB和braZ 103085
cat branched-chain/kaks/biosample/mega_clinical_braB_dnds.tsv  branched-chain/kaks/biosample/mega_clinical_braZ_dnds.tsv >branched-chain/kaks/biosample/mega_clinical.picture.tsv
#plotr画图
plotr hist --xl dN/dS --yl Frequency -g 2 --bins 20 --xmm -0.1,1.0 --ymm 0,0.5 -p branched-chain/kaks/biosample/mega_clinical.picture.tsv
#查看braB和braZ的dnds为0的数量，方便判断矩形颜色
 cat  branched-chain/kaks/biosample/mega_clinical.picture.tsv | grep -v "dnds" | tsv-filter --le 1:0.1 | grep "braB" | wc -l #14260
 cat  branched-chain/kaks/biosample/mega_clinical.picture.tsv | grep -v "dnds" | tsv-filter --le 1:0.1 | grep "braZ" | wc -l #42609


#计算braB的dn(51360)
sed -i '1d' branched-chain/kaks/biosample/clinical_braB_dn.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/biosample/clinical_braB_dn.tsv | sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/biosample/mega_clinical_braB_dn.tsv
#计算braB的ds(51360)
sed -i '1d' branched-chain/kaks/biosample/clinical_braB_ds.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/biosample/clinical_braB_ds.tsv | sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/biosample/mega_clinical_braB_ds.tsv
#查看braB的dn的平均值和平均值标准误差
library(plotrix)
setwd("~/data/Pseudomonas")
data<-read.table("branched-chain/kaks/biosample/mega_clinical_braB_dn.tsv")
data1<-data$V1
mean(data1)  # 0.002160744
std.error(data1) # 0.00001095158
#查看braB的ds的平均值和平均值标准误差
data<-read.table("branched-chain/kaks/biosample/mega_clinical_braB_ds.tsv")
data1<-data$V1
mean(data1)  # 0.01856531
std.error(data1) #0.0002155396

#计算braZ的dn(57630)
sed -i '1d' branched-chain/kaks/biosample/clinical_braZ_dn.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/biosample/clinical_braZ_dn.tsv | sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/biosample/mega_clinical_braZ_dn.tsv
#计算braZ的ds(57630)
sed -i '1d' branched-chain/kaks/biosample/clinical_braZ_ds.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/biosample/clinical_braZ_ds.tsv | sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/biosample/mega_clinical_braZ_ds.tsv
#查看braZ的dn的平均值和平均值标准误差
library(plotrix)
setwd("~/data/Pseudomonas")
data<-read.table("branched-chain/kaks/biosample/mega_clinical_braZ_dn.tsv")
data1<-data$V1
mean(data1)  # 0.01098205
std.error(data1) # 0.0001506622
#查看braZ的ds的平均值和平均值标准误差
data<-read.table("branched-chain/kaks/biosample/mega_clinical_braZ_ds.tsv")
data1<-data$V1
mean(data1)  # 0.05722712
std.error(data1) #0.0005527106

```


```bash
#生境为非临床的
#计算braB的dn/ds(696)
sed -i '1d' branched-chain/kaks/biosample/other_braB_dnds.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";'  branched-chain/kaks/biosample/other_braB_dnds.tsv | sed 's/NA//g' | sed  '/^\s*$/d' | perl -alne '$num=$F[0];$name=braB;print"$num\t$name";' | sed '1idnds\tfrequency' >branched-chain/kaks/biosample/mega_other_braB_dnds.tsv

#计算braZ的dn/ds(877)
sed -i '1d' branched-chain/kaks/biosample/other_braZ_dnds.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";'  branched-chain/kaks/biosample/other_braZ_dnds.tsv | sed 's/NA//g' | sed  '/^\s*$/d' | perl -alne '$num=$F[0];$name=braZ;print"$num\t$name";' >branched-chain/kaks/biosample/mega_other_braZ_dnds.tsv

#合并braB和braZ 1574
cat branched-chain/kaks/biosample/mega_other_braB_dnds.tsv  branched-chain/kaks/biosample/mega_other_braZ_dnds.tsv >branched-chain/kaks/biosample/mega_other.picture.tsv
#查看braB和braZ的dnds为0的数量，方便判断矩形颜色
 cat  branched-chain/kaks/biosample/mega_other.picture.tsv | grep -v "dnds" | tsv-filter --le 1:0.1 | grep "braB" | wc -l #212
 cat  branched-chain/kaks/biosample/mega_other.picture.tsv | grep -v "dnds" | tsv-filter --le 1:0.1 | grep "braZ" | wc -l #710
#plotr画图
plotr hist --xl dN/dS --yl Frequency  -g 2 --bins 20 --xmm -0.1,1.0 --ymm 0,0.5 -p branched-chain/kaks/biosample/mega_other.picture.tsv


#计算braB的dn(820)
sed -i '1d' branched-chain/kaks/biosample/other_braB_dn.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/biosample/other_braB_dn.tsv | sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/biosample/mega_other_braB_dn.tsv
#计算braB的ds(820)
sed -i '1d' branched-chain/kaks/biosample/other_braB_ds.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/biosample/other_braB_ds.tsv | sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/biosample/mega_other_braB_ds.tsv
#查看braB的dn的平均值和平均值标准误差
library(plotrix)
setwd("~/data/Pseudomonas")
data<-read.table("branched-chain/kaks/biosample/mega_other_braB_dn.tsv")
data1<-data$V1
mean(data1)  # 0.003035082
std.error(data1) # 0.0001411965
#查看braB的ds的平均值和平均值标准误差
data<-read.table("branched-chain/kaks/biosample/mega_other_braB_ds.tsv")
data1<-data$V1
mean(data1)  # 0.03932062
std.error(data1) #0.003167266

#计算braZ的dn(903)
sed -i '1d' branched-chain/kaks/biosample/other_braZ_dn.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/biosample/other_braZ_dn.tsv | sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/biosample/mega_other_braZ_dn.tsv
#计算braZ的ds(903)
sed -i '1d' branched-chain/kaks/biosample/other_braZ_ds.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/biosample/other_braZ_ds.tsv | sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/biosample/mega_other_braZ_ds.tsv
#查看braZ的dn的平均值和平均值标准误差
library(plotrix)
setwd("~/data/Pseudomonas")
data<-read.table("branched-chain/kaks/biosample/mega_other_braZ_dn.tsv")
data1<-data$V1
mean(data1)  # 0.008649314
std.error(data1) # 0.00100089
#查看braZ的ds的平均值和平均值标准误差
data<-read.table("branched-chain/kaks/biosample/mega_other_braZ_ds.tsv")
data1<-data$V1
mean(data1)  # 0.05772698
std.error(data1) #0.003851761
```

# 6. 使用mega计算其他假单胞菌的kaks
* Pseudomonas chlororaphis 绿刺假单胞菌，在农业和园艺中用作土壤接种剂的细菌,可以作为生物防治剂,通过生产吩嗪类抗生素来对抗某些真菌植物病原体
* Pseudomonas putida 恶臭假单胞菌 生活在土壤中，多样化代谢可以用于生物修复，生物降解，生物防治
* Pseudomonas syringae 丁香假单胞菌是一种具有极鞭毛的杆状革兰氏阴性细菌。作为一种植物病原体，可以感染多种物种
* Pseudomonas fluorescens 荧光假单胞菌是一种常见的革兰氏阴性杆状细菌，代谢灵活，在土壤和水中发现，是一种专性需氧菌
* Pseudomonas stutzeri 施氏假单胞菌,在人类脊髓液中分离出来

## 6.1分别提取branched在绿刺假单胞菌中的CDS序列
```bash           
#提取绿刺假单胞菌的单个拷贝的蛋白名称(59)
cd ~/data/Pseudomonas
mkdir -p ~/data/Pseudomonas/branched-chain/kaks/chl
cat branched-chain/branched-chain_minevalue.pfam.tsv | grep -f <(cat branched-chain/branched-chain_hmmscan_copy.pfam.tsv | grep "Pseudom_chl" | cut -f 1 ) |  cut -f 1 >branched-chain/kaks/chl/branched-chain.Pseudom_chl.protein.tsv
#提取绿刺假单胞菌中branced的cds名称和序列（59）
faops size branched-chain/kaks/aeru/branched-chain.CDS.fa | cut -f 1  | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/kaks/chl/branched-chain.Pseudom_chl.protein.tsv) >branched-chain/kaks/chl/branched-chain_pseudom_chl.cds.tsv
faops some branched-chain/kaks/aeru/branched-chain.CDS.fa branched-chain/kaks/chl/branched-chain_pseudom_chl.cds.tsv  branched-chain/kaks/chl/branched-chain_pseudom_chl.cds.fa

#绿刺假单胞菌的dnds结果画图(用mega计算生成三角矩阵文件)
#计算chl的brnQ的dn/ds
cd ~/data/Pseudomonas
#计算brnQ的dn/ds(1641)
sed -i '1d' branched-chain/kaks/chl/pseudom_chl-dN-dS.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/chl/pseudom_chl-dN-dS.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d' | perl -alne '$num=$F[0];$name=A;print"$num\t$name";' | sed '1idnds\tfrequency' >branched-chain/kaks/chl/mega_chl_dnds.tsv
#计算brnQ的dn(1711)
sed -i '1d' branched-chain/kaks/chl/pseudom_chl-dN.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/chl/pseudom_chl-dN.tsv |
sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/chl/mega_chl-dN.tsv
#计算brnQ的ds(1711)
sed -i '1d' branched-chain/kaks/chl/pseudom_chl-dS.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/chl/pseudom_chl-dS.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/chl/mega_chl-dS.tsv
#查看brnQ的dn的平均值和平均值标准误差
library(plotrix)
setwd("~/data/Pseudomonas")
data<-read.table("~/data/Pseudomonas/branched-chain/kaks/chl/mega_chl-dN.tsv")
data1<-data$V1
mean(data1)  # 0.003620104
std.error(data1) # 0.00007738851
#查看brnQ的ds的平均值和平均值标准误差
data<-read.table("~/data/Pseudomonas/branched-chain/kaks/chl/mega_chl-dS.tsv")
data1<-data$V1
mean(data1)  # 0.09491241
std.error(data1) # 0.001306388
#plotr画图
plotr hist --xl dN/dS --yl Frequency  --bins 20 --xmm -0.1,1.0 --ymm 0,1 -p branched-chain/kaks/chl/mega_chl_dnds.tsv
```

## 6.2分别提取branched在丁香假单胞菌中的CDS序列
```bash            
#提取丁香假单胞菌的单个拷贝的蛋白名称(40)
cd ~/data/Pseudomonas
mkdir -p ~/data/Pseudomonas/branched-chain/kaks/syr
cat branched-chain/branched-chain_minevalue.pfam.tsv | grep -f <(cat branched-chain/branched-chain_hmmscan_copy.pfam.tsv | grep "Pseudom_syr" | cut -f 1 ) |  cut -f 1 >branched-chain/kaks/syr/branched-chain.Pseudom_syr.protein.tsv
#提取丁香假单胞菌中branced的cds名称和序列（40）
faops size branched-chain/kaks/aeru/branched-chain.CDS.fa | cut -f 1  | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/kaks/syr/branched-chain.Pseudom_syr.protein.tsv) >branched-chain/kaks/syr/branched-chain_pseudom_syr.cds.tsv
faops some branched-chain/kaks/aeru/branched-chain.CDS.fa branched-chain/kaks/syr/branched-chain_pseudom_syr.cds.tsv  branched-chain/kaks/syr/branched-chain_pseudom_syr.cds.fa

#丁香假单胞菌的dnds结果画图(用mega计算生成三角矩阵文件)
#计算syr的brnQ的dn/ds
cd ~/data/Pseudomonas
#计算brnQ的dn/ds(721)
sed -i '1d' branched-chain/kaks/syr/pseudom_syr-dN-dS.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/syr/pseudom_syr-dN-dS.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d' | perl -alne '$num=$F[0];$name=A;print"$num\t$name";' | sed '1idnds\tfrequency' >branched-chain/kaks/syr/mega_syr_dnds.tsv
#计算brnQ的dn(780)
sed -i '1d' branched-chain/kaks/syr/pseudom_syr-dN.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/syr/pseudom_syr-dN.tsv |
sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/syr/mega_syr-dN.tsv
#计算brnQ的ds(780)
sed -i '1d' branched-chain/kaks/syr/pseudom_syr-dS.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/syr/pseudom_syr-dS.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/syr/mega_syr-dS.tsv
#查看brnQ的dn的平均值和平均值标准误差
library(plotrix)
setwd("~/data/Pseudomonas")
data<-read.table("~/data/Pseudomonas/branched-chain/kaks/syr/mega_syr-dN.tsv")
data1<-data$V1
mean(data1)  # 0.01634968
std.error(data1) # 0.000422107
#查看brnQ的ds的平均值和平均值标准误差
data<-read.table("~/data/Pseudomonas/branched-chain/kaks/syr/mega_syr-dS.tsv")
data1<-data$V1
mean(data1)  # 0.4171773
std.error(data1) # 0.009293657
#plotr画图
plotr hist --xl dN/dS --yl Frequency  --bins 20 --xmm -0.1,1.0 --ymm 0,1 -p branched-chain/kaks/syr/mega_syr_dnds.tsv
```

## 6.3分别提取branched在荧光假单胞菌中的CDS序列
```bash            
#提取荧光假单胞菌的单个拷贝的蛋白名称(38)
cd ~/data/Pseudomonas
mkdir -p ~/data/Pseudomonas/branched-chain/kaks/flu
cat branched-chain/branched-chain_minevalue.pfam.tsv | grep -f <(cat branched-chain/branched-chain_hmmscan_copy.pfam.tsv | grep "Pseudom_flu" | cut -f 1 ) |  cut -f 1 >branched-chain/kaks/flu/branched-chain.Pseudom_flu.protein.tsv
#提取荧光假单胞菌中branced的cds名称和序列（38）
faops size branched-chain/kaks/aeru/branched-chain.CDS.fa | cut -f 1  | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/kaks/flu/branched-chain.Pseudom_flu.protein.tsv) >branched-chain/kaks/flu/branched-chain_pseudom_flu.cds.tsv
faops some branched-chain/kaks/aeru/branched-chain.CDS.fa branched-chain/kaks/flu/branched-chain_pseudom_flu.cds.tsv  branched-chain/kaks/flu/branched-chain_pseudom_flu.cds.fa

#计算flu的brnQ的dn/ds
cd ~/data/Pseudomonas
#计算brnQ的dn/ds(698)
sed -i '1d' branched-chain/kaks/flu/pseudom_flu-dN-dS.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/flu/pseudom_flu-dN-dS.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d' | perl -alne '$num=$F[0];$name=A;print"$num\t$name";' | sed '1idnds\tfrequency' >branched-chain/kaks/flu/mega_flu_dnds.tsv
#计算brnQ的dn(703)
sed -i '1d' branched-chain/kaks/flu/pseudom_flu-dN.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/flu/pseudom_flu-dN.tsv |
sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/flu/mega_flu-dN.tsv
#计算brnQ的ds(703)
sed -i '1d' branched-chain/kaks/flu/pseudom_flu-dS.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/flu/pseudom_flu-dS.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/flu/mega_flu-dS.tsv
#查看brnQ的dn的平均值和平均值标准误差
library(plotrix)
setwd("~/data/Pseudomonas")
data<-read.table("~/data/Pseudomonas/branched-chain/kaks/flu/mega_flu-dN.tsv")
data1<-data$V1
mean(data1)  # 0.03741564
std.error(data1) # 0.0006310837
#查看brnQ的ds的平均值和平均值标准误差
data<-read.table("~/data/Pseudomonas/branched-chain/kaks/flu/mega_flu-dS.tsv")
data1<-data$V1
mean(data1)  # 0.5835503
std.error(data1) # 0.005500969
#plotr画图
plotr hist --xl dN/dS --yl Frequency  --bins 20 --xmm -0.1,1.0 --ymm 0,1 -p branched-chain/kaks/flu/mega_flu_dnds.tsv
```

## 6.4分别提取branched在恶臭假单胞菌中的CDS序列
```bash            
#提取恶臭假单胞菌的单个拷贝的蛋白名称(49)
cd ~/data/Pseudomonas
mkdir -p ~/data/Pseudomonas/branched-chain/kaks/puti
cat branched-chain/branched-chain_minevalue.pfam.tsv | grep -f <(cat branched-chain/branched-chain_hmmscan_copy.pfam.tsv | grep "Pseudom_puti" | cut -f 1 ) |  cut -f 1 >branched-chain/kaks/puti/branched-chain.Pseudom_puti.protein.tsv
#提取恶臭假单胞菌中branced的cds名称和序列（49）
faops size branched-chain/kaks/aeru/branched-chain.CDS.fa | cut -f 1  | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/kaks/puti/branched-chain.Pseudom_puti.protein.tsv) >branched-chain/kaks/puti/branched-chain_pseudom_puti.cds.tsv
faops some branched-chain/kaks/aeru/branched-chain.CDS.fa branched-chain/kaks/puti/branched-chain_pseudom_puti.cds.tsv  branched-chain/kaks/puti/branched-chain_pseudom_puti.cds.fa

#恶臭假单胞菌的dnds结果画图(用mega计算生成三角矩阵文件)
#计算puti的brnQ的dn/ds
cd ~/data/Pseudomonas
#计算brnQ的dn/ds(1118)
sed -i '1d' branched-chain/kaks/puti/pseudom_puti-dN-dS.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/puti/pseudom_puti-dN-dS.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d' | perl -alne '$num=$F[0];$name=A;print"$num\t$name";' | sed '1idnds\tfrequency' >branched-chain/kaks/puti/mega_puti_dnds.tsv
#计算brnQ的dn(1128)
sed -i '1d' branched-chain/kaks/puti/pseudom_puti-dN.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/puti/pseudom_puti-dN.tsv |
sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/puti/mega_puti-dN.tsv
#计算brnQ的ds(1128)
sed -i '1d' branched-chain/kaks/puti/pseudom_puti-dS.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/puti/pseudom_puti-dS.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/puti/mega_puti-dS.tsv
#查看brnQ的dn的平均值和平均值标准误差
library(plotrix)
setwd("~/data/Pseudomonas")
data<-read.table("~/data/Pseudomonas/branched-chain/kaks/puti/mega_puti-dN.tsv")
data1<-data$V1
mean(data1)  # 0.0168802
std.error(data1) # 0.0003544275
#查看brnQ的ds的平均值和平均值标准误差
data<-read.table("~/data/Pseudomonas/branched-chain/kaks/puti/mega_puti-dS.tsv")
data1<-data$V1
mean(data1)  #  0.3353134
std.error(data1) # 0.004536949
#plotr画图
plotr hist --xl dN/dS --yl Frequency  --bins 20 --xmm -0.1,1.0 --ymm 0,1 -p branched-chain/kaks/puti/mega_puti_dnds.tsv
```

## 6.5分别提取branched在施式假单胞菌中的CDS序列
```bash            
#提取施式假单胞菌的单个拷贝的蛋白名称(30)
cd ~/data/Pseudomonas
mkdir -p ~/data/Pseudomonas/branched-chain/kaks/stu
cat branched-chain/branched-chain_minevalue.pfam.tsv | grep -f <(cat branched-chain/branched-chain_hmmscan_copy.pfam.tsv | grep "Pseudom_stu" | cut -f 1 ) |  cut -f 1 >branched-chain/kaks/stu/branched-chain.Pseudom_stu.protein.tsv
#提取恶臭假单胞菌中branced的cds名称和序列（30）
faops size branched-chain/kaks/aeru/branched-chain.CDS.fa | cut -f 1  | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/kaks/stu/branched-chain.Pseudom_stu.protein.tsv) >branched-chain/kaks/stu/branched-chain_pseudom_stu.cds.tsv
faops some branched-chain/kaks/aeru/branched-chain.CDS.fa branched-chain/kaks/stu/branched-chain_pseudom_stu.cds.tsv  branched-chain/kaks/stu/branched-chain_pseudom_stu.cds.fa


#施式假单胞菌的dnds结果画图(用mega计算生成三角矩阵文件)
#计算stu的brnQ的dn/ds
cd ~/data/Pseudomonas
#计算brnQ的dn/ds(433)
sed -i '1d' branched-chain/kaks/stu/pseudom_stu-dN-dS.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/stu/pseudom_stu-dN-dS.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d' | perl -alne '$num=$F[0];$name=A;print"$num\t$name";' | sed '1idnds\tfrequency' >branched-chain/kaks/stu/mega_stu_dnds.tsv
#计算brnQ的dn(435)
sed -i '1d' branched-chain/kaks/stu/pseudom_stu-dN.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/stu/pseudom_stu-dN.tsv |
sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/stu/mega_stu-dN.tsv
#计算brnQ的ds(435)
sed -i '1d' branched-chain/kaks/stu/pseudom_stu-dS.tsv
perl -ne 'chomp;($id,$name)=(split/\t/,$_,2)[0,1];@data=(split/\t/,$name);$num=join"\n",@data;print"$num\n";' branched-chain/kaks/stu/pseudom_stu-dS.tsv  |
sed 's/NA//g' | sed  '/^\s*$/d' >branched-chain/kaks/stu/mega_stu-dS.tsv
#查看brnQ的dn的平均值和平均值标准误差
library(plotrix)
setwd("~/data/Pseudomonas")
data<-read.table("~/data/Pseudomonas/branched-chain/kaks/stu/mega_stu-dN.tsv")
data1<-data$V1
mean(data1)  # 0.06714633
std.error(data1) # 0.002394998
#查看brnQ的ds的平均值和平均值标准误差
data<-read.table("~/data/Pseudomonas/branched-chain/kaks/stu/mega_stu-dS.tsv")
data1<-data$V1
mean(data1)  #   0.8034892
std.error(data1) #  0.01867506
#plotr画图
plotr hist --xl dN/dS --yl Frequency  --bins 20 --xmm -0.1,1.0 --ymm 0,1 -p branched-chain/kaks/stu/mega_stu_dnds.tsv
```

# 7. 基因共线性
## 7.1 准备假单胞菌属的典型菌株的genbank文件，共9种12个
```BASH
cd ~/data/Pseudomonas
for S in \
    Pseudom_aeru_PAO1 \
    Pseudom_puti_KT2440_GCF_000007565_2 \
    Pseudom_chl_aureofaciens_GCF_003851365_1 \
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

## 7.2 查看铜绿假单胞菌物种内的基因共线性
```bash
cd ~/data/Pseudomonas
#截取目标gbk
mkdir -p  ~/data/Pseudomonas/branched-chain/collinearity/pse_auer
python  fetch_gbk.py -i branched-chain/collinearity/genebank/Pseudom_aeru_PAO1.gbff -e 10000 -l braB -o branched-chain/collinearity/pse_auer/PAO1_braB
python  fetch_gbk.py -i branched-chain/collinearity/genebank/Pseudom_aeru_PAO1.gbff -e 10000 -l braZ -o branched-chain/collinearity/pse_auer/PAO1_braZ
python  fetch_gbk.py -i branched-chain/collinearity/genebank/Pseudom_aeru_PA7_GCF_000017205_1.gbff -e 10000 -l brnQ -o branched-chain/collinearity/pse_auer/PA7_brnQ
python  fetch_gbk.py -i branched-chain/collinearity/genebank/Pseudom_aeru_LESB58_GCF_000026645_1.gbff -e 10000 -l brnQ -o branched-chain/collinearity/pse_auer/LESB58_brnQ
python  fetch_gbk.py -i branched-chain/collinearity/genebank/Pseudom_aeru_UCBPP_PA14_GCF_000014625_1.gbff -e 10000 -l brnQ -o branched-chain/collinearity/pse_auer/UCBPP_brnQ
#绘制braZ和brnQ1的图
mkdir -p branched-chain/collinearity/pse_auer/braZ_brnQ1
for f in $(find branched-chain/collinearity/pse_auer -maxdepth 2 -type f -name "*.gbk")
do
cp $f  branched-chain/collinearity/pse_auer/braZ_brnQ1
done
rm -rf branched-chain/collinearity/pse_auer/braZ_brnQ1/*_brnQ_2_rc.gbk
rm -rf branched-chain/collinearity/pse_auer/braZ_brnQ1/*_braB_1.gbk
#在clinker链接的网页上将group名调整为对应的基因
#调整菌株名字
clinker -p branched-chain/collinearity/pse_auer/braZ_brnQ1/*.gbk -p

#绘制braB和brnQ2的图
mkdir branched-chain/collinearity/pse_auer/braB_brnQ2
for f in $(find branched-chain/collinearity/pse_auer -maxdepth 2 -type f -name "*.gbk")
do
cp $f  branched-chain/collinearity/pse_auer/braB_brnQ2
done
rm -rf branched-chain/collinearity/pse_auer/braB_brnQ2/*_1_brnQ_1.gbk
rm -rf branched-chain/collinearity/pse_auer/braB_brnQ2/*_braZ_1_rc.gbk
clinker -p branched-chain/collinearity/pse_auer/braB_brnQ2/*.gbk -p
```
## 7.3 查看假单胞菌属内物种间的基因共线性
```bash
cd ~/data/Pseudomonas
#截取目标gbk   
mkdir -p  ~/data/Pseudomonas/branched-chain/collinearity/pse_other
python  fetch_gbk.py -i branched-chain/collinearity/genebank/Pseudom_aeru_PAO1.gbff -e 10000 -l braB -o branched-chain/collinearity/pse_other/PAO1_braB
python  fetch_gbk.py -i branched-chain/collinearity/genebank/Pseudom_aeru_PAO1.gbff -e 10000 -l braZ -o branched-chain/collinearity/pse_other/PAO1_braZ
python  fetch_gbk.py -i branched-chain/collinearity/genebank/Pseudom_puti_KT2440_GCF_000007565_2.gbff -e 10000 -l brnQ -o branched-chain/collinearity/pse_other/puti_brnQ
python  fetch_gbk.py -i branched-chain/collinearity/genebank/Pseudom_entomophi_L48_GCF_000026105_1.gbff -e 10000 -l brnQ -o branched-chain/collinearity/pse_other/entomophi_brnQ
python  fetch_gbk.py -i branched-chain/collinearity/genebank/Pseudom_fluo_SBW25_GCF_000009225_2.gbff -e 10000 -l brnQ -o branched-chain/collinearity/pse_other/fluo_brnQ
python  fetch_gbk.py -i branched-chain/collinearity/genebank/Pseudom_prot_Pf_5_GCF_000012265_1.gbff -e 10000 -l brnQ -o branched-chain/collinearity/pse_other/prot_brnQ
python  fetch_gbk.py -i branched-chain/collinearity/genebank/Pseudom_sav_pv_phaseolicola_1448A_GCF_000012205_1.gbff  -e 10000 -l brnQ -o branched-chain/collinearity/pse_other/sav_brnQ
python  fetch_gbk.py -i branched-chain/collinearity/genebank/Pseudom_syr_pv_syringae_B728a_GCF_000012245_1.gbff  -e 10000 -l brnQ -o branched-chain/collinearity/pse_other/syr_brnQ
python  fetch_gbk.py -i branched-chain/collinearity/genebank/Pseudom_stu_A1501_GCF_000013785_1.gbff  -e 10000 -l brnQ -o branched-chain/collinearity/pse_other/stu_brnQ
python  fetch_gbk.py -i branched-chain/collinearity/genebank/Pseudom_chl_aureofaciens_GCF_003851365_1.gbff -e 10000 -l brnQ -o branched-chain/collinearity/pse_other/chl_brnQ

#绘制braZ和brnQ的图
mkdir -p branched-chain/collinearity/pse_other/braZ_brnQ
for f in $(find branched-chain/collinearity/pse_other -maxdepth 2 -type f -name "*.gbk")
do
cp $f  branched-chain/collinearity/pse_other/braZ_brnQ
done
rm -rf branched-chain/collinearity/pse_other/braZ_brnQ/*_braB_1.gbk
#在clinker链接的网页上将group名调整为对应的基因
#调整菌株名字
clinker -p branched-chain/collinearity/pse_other/braZ_brnQ/*.gbk -p

#绘制braB和brnQ的图
mkdir -p branched-chain/collinearity/pse_other/braB_brnQ
for f in $(find branched-chain/collinearity/pse_other -maxdepth 2 -type f -name "*.gbk")
do
cp $f  branched-chain/collinearity/pse_other/braB_brnQ
done
rm -rf branched-chain/collinearity/pse_other/braB_brnQ/*_braZ_1_rc.gbk
clinker -p branched-chain/collinearity/pse_other/braB_brnQ/*.gbk -p

#查看菌株对应菌株名
cat typical.lst | tsv-join -d 1 -f strains.taxon.tsv -k 1 --append-fields 2
Pseudom_aeru_PAO1       Pseudomonas aeruginosa PAO1
Pseudom_puti_KT2440_GCF_000007565_2     Pseudomonas putida KT2440
Pseudom_chl_aureofaciens_GCF_003851365_1        Pseudomonas chlororaphis subsp. aureofaciens
Pseudom_entomophi_L48_GCF_000026105_1   Pseudomonas entomophila L48
Pseudom_fluo_SBW25_GCF_000009225_2      Pseudomonas fluorescens SBW25
Pseudom_prot_Pf_5_GCF_000012265_1       Pseudomonas protegens Pf-5
Pseudom_sav_pv_phaseolicola_1448A_GCF_000012205_1       Pseudomonas savastanoi pv. phaseolicola 1448A
Pseudom_stu_A1501_GCF_000013785_1       Pseudomonas stutzeri A1501
Pseudom_syr_pv_syringae_B728a_GCF_000012245_1   Pseudomonas syringae pv. syringae B728a
Pseudom_aeru_UCBPP_PA14_GCF_000014625_1 Pseudomonas aeruginosa UCBPP-PA14
Pseudom_aeru_PA7_GCF_000017205_1        Pseudomonas aeruginosa PA7
Pseudom_aeru_LESB58_GCF_000026645_1     Pseudomonas aeruginosa LESB58
```

## 7.4 13个两个braZ拷贝的共线性分析（gbk转gff或者gff转gbk）
```bash
cd ~/data/Pseudomonas
#提取所有的含有两个braZ的菌株(13)
mkdir -p  ~/data/Pseudomonas/branched-chain/collinearity/two_braz
cat branched-chain/kaks/aeru/branched-chain_cluster2.tsv  | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 | tsv-summarize -g 2 --count| tsv-filter --eq 2:2 | cut -f 1  >branched-chain/collinearity/two_braz/pse_auer_two_braz.tsv
#提取上述菌株的gbk文件
for name in $(cat branched-chain/collinearity/two_braz/pse_auer_two_braz.tsv)
do
echo $name
cp ASSEMBLY/$name/*.gbff.gz ASSEMBLY/$name/$name.gbff.gz
done

for name in $(cat branched-chain/collinearity/two_braz/pse_auer_two_braz.tsv)
do
echo $name
mv ASSEMBLY/$name/$name.gbff.gz branched-chain/collinearity/two_braz
done
#解压上述文件
gunzip branched-chain/collinearity/two_braz/*.gz
#截取目标gbk  
for name in $(cat branched-chain/collinearity/two_braz/pse_auer_two_braz.tsv)
do
echo $name
python  fetch_gbk.py -i branched-chain/collinearity/two_braz/$name.gbff -e 10000 -l brnQ -o branched-chain/collinearity/two_braz/$name
done

#有三个特殊菌株的gbk文件无法提取，可以使用gff文件和genomic文件转换为gbk文件再提取
Pseudom_aeru_GCF_015697465_1
Pseudom_aeru_GCF_015697605_1
Pseudom_aeru_GCF_018141645_1
#注意有个模块需要改biopython的版本
pip3 unintall biopython
pip3 install biopython==1.76
python  gff_to_gbk.py Pseudom_aeru_PA1R_GCF_000496645_1.gff  Pseudom_aeru_PA1R_GCF_000496645_1.fna
```
## 7.5 查看9个单拷贝铜绿假单胞菌的共线性
```bash
cd ~/data/Pseudomonas
#提取所有的铜绿假单胞菌的单拷贝蛋白序列(9)
cat branched-chain/branched-chain_hmmscan_copy.pfam.tsv | tsv-filter --lt 3:2 | grep "Pseudom_aeru" | cut -f 1 >branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.single.copy.tsv
faops some PROTEINS/all.replace.fa <(cat branched-chain/branched-chain_minevalue.pfam.tsv|grep -f branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.single.copy.tsv | cut -f 1 ) branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.single.copy.fa
#提取上述菌株的gbk文件
mkdir -p branched-chain/collinearity/single_pse_aeru
for name in $(cat branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.single.copy.tsv)
do
echo $name
cp ASSEMBLY/$name/*.gbff.gz ASSEMBLY/$name/$name.gbff.gz
done

for name in $(cat branched-chain/kaks/aeru/branched-chain.Pseudom_aeru.single.copy.tsv)
do
echo $name
mv ASSEMBLY/$name/$name.gbff.gz 、branched-chain/collinearity/single_pse_aeru
done
#解压上述文件
gunzip branched-chain/collinearity/single_pse_aeru/*.gz

#截取上述文件的目标序列上下游各10000bp
for f in $(find branched-chain/collinearity/single_pse_aeru -maxdepth 2 -type f -name "*.gbff")
do
python  fetch_gbk.py -i $f  -e 10000 -l brnQ -o branched-chain/collinearity/single_pse_aeru/brnQ
done

#删除注释文件中没有氨基酸序列但是有DNA序列的文件
rm -rf branched-chain/collinearity/single_pse_aeru/PAO1/Pseudom_aeru_GCF_001516205_2_brnQ_1.gbk
rm -rf branched-chain/collinearity/single_pse_aeru/PAO1/Pseudom_aeru_GCF_002205355_1_brnQ_1_rc.gbk
rm -rf branched-chain/collinearity/single_pse_aeru/PAO1/Pseudom_aeru_GCF_002205355_1_brnQ_2.gbk
rm -rf branched-chain/collinearity/single_pse_aeru/PAO1/Pseudom_aeru_PAK_GCF_000568855_2_brnQ_1.gbk
rm -rf branched-chain/collinearity/single_pse_aeru/PAO1/Pseudom_aeru_PAK_GCF_902172305_2_brnQ_2_rc.gbk

#有一个菌株的gbk文件无法提取，可以使用gff文件和genomic文件转换为gbk文件再提取
Pseudom_aeru_PA1R_GCF_000496645_1
#转换gbk
python  gff_to_gbk.py Pseudom_aeru_PA1R_GCF_000496645_1.gff  Pseudom_aeru_PA1R_GCF_000496645_1.fna
#提取gbk
python  fetch_gbk.py -i Pseudom_aeru_PA1R_GCF_000496645_1.gb -e 10000 -l brnQ -o Pseudom_aeru_PA1R_GCF_000496645_1
```
## 7.6查看假单胞菌属的基因树两个clade的物种分别与braB和braZ的共线性
```bash
cd ~/data/Pseudomonas
mkdir -p branched-chain/collinearity/pseudomons_two_clade
mkdir -p branched-chain/collinearity/pseudomons_two_clade/clade1
#提取clade1的菌株蛋白名
echo -e "Pseudom_oti_GCF_011397855_1_WP_074969323\nPseudom_alcalig_GCF_001597285_1_WP_061904350\nPseudom_men_S5_2_GCF_000733715_2_WP_017360782\nPseudom_stu_GCF_019704535_1_WP_014596526\nPseudom_aeru_PAO1_NP_250281" >branched-chain/collinearity/pseudomons_two_clade/clade1/clade1.tsv
#提取clade1的菌株名
cat branched-chain/collinearity/pseudomons_two_clade/clade1/clade1.tsv | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 | cut -f 2 >branched-chain/collinearity/pseudomons_two_clade/clade1/clade1_species.tsv
#提取clade1菌株的gbff文件
for name in $(cat branched-chain/collinearity/pseudomons_two_clade/clade1/clade1_species.tsv)
do
echo $name
cp ASSEMBLY/$name/*.gbff.gz ASSEMBLY/$name/$name.gbff.gz
done

for name in $(cat branched-chain/collinearity/pseudomons_two_clade/clade1/clade1_species.tsv)
do
echo $name
mv ASSEMBLY/$name/$name.gbff.gz branched-chain/collinearity/pseudomons_two_clade/clade1
done
#解压上述文件
gunzip branched-chain/collinearity/pseudomons_two_clade/clade1/*.gz
#截取其他菌株的gbk  
for name in $(cat branched-chain/collinearity/pseudomons_two_clade/clade1/clade1_species.tsv)
do
echo $name
python  fetch_gbk.py -i branched-chain/collinearity/pseudomons_two_clade/clade1/$name.gbff -e 10000 -l brnQ -o branched-chain/collinearity/pseudomons_two_clade/clade1/$name
done
#截取PAO1的gbk
python  fetch_gbk.py -i branched-chain/collinearity/pseudomons_two_clade/clade1/Pseudom_aeru_PAO1.gbff -e 10000 -l braB -o branched-chain/collinearity/pseudomons_two_clade/clade1/Pseudom_aeru_PAO1
python  fetch_gbk.py -i branched-chain/collinearity/pseudomons_two_clade/clade1/Pseudom_aeru_PAO1.gbff -e 10000 -l braZ -o branched-chain/collinearity/pseudomons_two_clade/clade1/Pseudom_aeru_PAO1
###注意clade1中的Pseudom_alcalig_GCF_001597285_1和Pseudom_oti_GCF_011397855_1没有gbff文件
for name in $(cat branched-chain/collinearity/pseudomons_two_clade/clade1/clade1_species.tsv)
do
echo $name
cp ASSEMBLY/$name/*.gff.gz ASSEMBLY/$name/$name.gff.gz
done

for name in $(cat branched-chain/collinearity/pseudomons_two_clade/clade1/clade1_species.tsv)
do
echo $name
mv ASSEMBLY/$name/$name.gff.gz branched-chain/collinearity/pseudomons_two_clade/clade1
done

cp ASSEMBLY/Pseudom_alcalig_GCF_001597285_1/GCF_001597285.1_ASM159728v1_genomic.fna.gz ASSEMBLY/Pseudom_alcalig_GCF_001597285_1/Pseudom_alcalig_GCF_001597285_1.fna.gz
mv ASSEMBLY/Pseudom_alcalig_GCF_001597285_1/Pseudom_alcalig_GCF_001597285_1.na.gz branched-chain/collinearity/pseudomons_two_clade/clade1
cp ASSEMBLY/Pseudom_oti_GCF_011397855_1/GCF_011397855.1_ASM1139785v1_genomic.fna.gz ASSEMBLY/Pseudom_oti_GCF_011397855_1/Pseudom_oti_GCF_011397855_1.fna.gz
mv ASSEMBLY/Pseudom_oti_GCF_011397855_1/Pseudom_oti_GCF_011397855_1.fna.gz branched-chain/collinearity/pseudomons_two_clade/clade1
gunzip branched-chain/collinearity/pseudomons_two_clade/clade1/*.gz
pip3 unintall biopython
pip3 install biopython==1.76
python  gff_to_gbk.py branched-chain/collinearity/pseudomons_two_clade/clade1/Pseudom_alcalig_GCF_001597285_1.gff  branched-chain/collinearity/pseudomons_two_clade/clade1/Pseudom_alcalig_GCF_001597285_1.fna
python  gff_to_gbk.py branched-chain/collinearity/pseudomons_two_clade/clade1/Pseudom_oti_GCF_011397855_1.gff  branched-chain/collinearity/pseudomons_two_clade/clade1/Pseudom_oti_GCF_011397855_1.fna
###Pseudom_alcalig_GCF_001597285_1和Pseudom_oti_GCF_011397855_1的gbk提取
python  fetch_gbk.py -i branched-chain/collinearity/pseudomons_two_clade/clade1/Pseudom_alcalig_GCF_001597285_1.gb -e 10000 -l brnQ -o branched-chain/collinearity/pseudomons_two_clade/clade1/Pseudom_alcalig_GCF_001597285_1
python  fetch_gbk.py -i branched-chain/collinearity/pseudomons_two_clade/clade1/Pseudom_oti_GCF_011397855_1.gb -e 10000 -l brnQ -o branched-chain/collinearity/pseudomons_two_clade/clade1/Pseudom_oti_GCF_011397855_1


#提取clade2的菌株蛋白名
mkdir -p branched-chain/collinearity/pseudomons_two_clade/clade2
echo -e "Pseudom_lal_GCF_008807375_1_WP_151133060\nPseudom_kna_GCF_009911755_1_WP_160286596\nPseudom_citro_GCF_001586155_1_WP_009619386\nPseudom_puti_NBRC_14164_GCF_000412675_1_WP_016498696\nPseudom_syr_GCF_016694755_2_WP_103688145\nPseudom_fluo_GCF_900215245_1_WP_053255234\nPseudom_aeru_PAO1_NP_250281" >branched-chain/collinearity/pseudomons_two_clade/clade2/clade2.tsv
#提取clade2的菌株名
cat branched-chain/collinearity/pseudomons_two_clade/clade2/clade2.tsv | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 | cut -f 2 >branched-chain/collinearity/pseudomons_two_clade/clade2/clade2_species.tsv
#提取clade2菌株的gbff文件
for name in $(cat branched-chain/collinearity/pseudomons_two_clade/clade2/clade2_species.tsv)
do
echo $name
cp ASSEMBLY/$name/*.gbff.gz ASSEMBLY/$name/$name.gbff.gz
done

for name in $(cat branched-chain/collinearity/pseudomons_two_clade/clade2/clade2_species.tsv)
do
echo $name
mv ASSEMBLY/$name/$name.gbff.gz branched-chain/collinearity/pseudomons_two_clade/clade2
done
#解压上述文件
gunzip branched-chain/collinearity/pseudomons_two_clade/clade2/*.gz
#截取其他菌株的gbk  
for name in $(cat branched-chain/collinearity/pseudomons_two_clade/clade2/clade2_species.tsv)
do
echo $name
python  fetch_gbk.py -i branched-chain/collinearity/pseudomons_two_clade/clade2/$name.gbff -e 10000 -l brnQ -o branched-chain/collinearity/pseudomons_two_clade/clade2/$name
done
#截取PAO1的gbk
python  fetch_gbk.py -i branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_aeru_PAO1.gbff -e 10000 -l braB -o branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_aeru_PAO1
python  fetch_gbk.py -i branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_aeru_PAO1.gbff -e 10000 -l braZ -o branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_aeru_PAO1

###注意clade2中的Pseudom_citro_GCF_001586155_1和Pseudom_kna_GCF_009911755_1和Pseudom_lal_GCF_008807375_1没有gbff文件
cp ASSEMBLY/Pseudom_citro_GCF_001586155_1/GCF_001586155.1_ASM158615v1_genomic.gff.gz ASSEMBLY/Pseudom_citro_GCF_001586155_1/Pseudom_citro_GCF_001586155_1.gff.gz 
cp ASSEMBLY/Pseudom_citro_GCF_001586155_1/GCF_001586155.1_ASM158615v1_genomic.fna.gz ASSEMBLY/Pseudom_citro_GCF_001586155_1/Pseudom_citro_GCF_001586155_1.fna.gz
mv ASSEMBLY/Pseudom_citro_GCF_001586155_1/Pseudom_citro_GCF_001586155_1.gff.gz   branched-chain/collinearity/pseudomons_two_clade/clade2
mv ASSEMBLY/Pseudom_citro_GCF_001586155_1/Pseudom_citro_GCF_001586155_1.fna.gz   branched-chain/collinearity/pseudomons_two_clade/clade2

cp ASSEMBLY/Pseudom_kna_GCF_009911755_1/GCF_009911755.1_ASM991175v1_genomic.gff.gz ASSEMBLY/Pseudom_kna_GCF_009911755_1/Pseudom_kna_GCF_009911755_1.gff.gz
cp ASSEMBLY/Pseudom_kna_GCF_009911755_1/GCF_009911755.1_ASM991175v1_genomic.fna.gz ASSEMBLY/Pseudom_kna_GCF_009911755_1/Pseudom_kna_GCF_009911755_1.fna.gz 
mv ASSEMBLY/Pseudom_kna_GCF_009911755_1/Pseudom_kna_GCF_009911755_1.gff.gz   branched-chain/collinearity/pseudomons_two_clade/clade2
mv ASSEMBLY/Pseudom_kna_GCF_009911755_1/Pseudom_kna_GCF_009911755_1.fna.gz   branched-chain/collinearity/pseudomons_two_clade/clade2

cp ASSEMBLY/Pseudom_lal_GCF_008807375_1/GCF_008807375.1_ASM880737v1_genomic.gff.gz ASSEMBLY/Pseudom_lal_GCF_008807375_1/Pseudom_lal_GCF_008807375_1.gff.gz
cp ASSEMBLY/Pseudom_lal_GCF_008807375_1/GCF_008807375.1_ASM880737v1_genomic.fna.gz ASSEMBLY/Pseudom_lal_GCF_008807375_1/Pseudom_lal_GCF_008807375_1.fna.gz 
mv ASSEMBLY/Pseudom_lal_GCF_008807375_1/Pseudom_lal_GCF_008807375_1.gff.gz  branched-chain/collinearity/pseudomons_two_clade/clade2
mv ASSEMBLY/Pseudom_lal_GCF_008807375_1/Pseudom_lal_GCF_008807375_1.fna.gz   branched-chain/collinearity/pseudomons_two_clade/clade2
gunzip branched-chain/collinearity/pseudomons_two_clade/clade2/*.gz
#转换格式
python  gff_to_gbk.py branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_citro_GCF_001586155_1.gff  branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_citro_GCF_001586155_1.fna
python  gff_to_gbk.py branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_kna_GCF_009911755_1.gff branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_kna_GCF_009911755_1.fna
python  gff_to_gbk.py branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_lal_GCF_008807375_1.gff branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_lal_GCF_008807375_1.fna
###Pseudom_citro_GCF_001586155_1和Pseudom_kna_GCF_009911755_1和Pseudom_lal_GCF_008807375_1的gbk提取
python  fetch_gbk.py -i branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_citro_GCF_001586155_1.gb -e 10000 -l brnQ -o branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_citro_GCF_001586155_1
python  fetch_gbk.py -i branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_kna_GCF_009911755_1.gb -e 10000 -l brnQ -o branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_kna_GCF_009911755_1
python  fetch_gbk.py -i branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_lal_GCF_008807375_1.gb -e 10000 -l brnQ -o branched-chain/collinearity/pseudomons_two_clade/clade2/Pseudom_lal_GCF_008807375_1
```


## 7.6查看假单胞菌属的基因树两个clade的物种分别与braZ上游基因PA1064的共线性
```bash
cd ~/data/Pseudomonas
mkdir -p branched-chain/collinearity/pseudomons_two_clade/PA1964
#clade1的菌株名
branched-chain/collinearity/pseudomons_two_clade/clade1/clade1_species.tsv
#clade1菌株的gbff文件
ls -lh branched-chain/collinearity/pseudomons_two_clade/clade1/*.gbff
###注意clade1中的Pseudom_alcalig_GCF_001597285_1和Pseudom_oti_GCF_011397855_1没有gbff文件
ls -lh branched-chain/collinearity/pseudomons_two_clade/clade1/*.gb

#clade1的菌株名
branched-chain/collinearity/pseudomons_two_clade/clade2/clade2_species.tsv
#clade1菌株的gbff文件
ls -lh branched-chain/collinearity/pseudomons_two_clade/clade2/*.gbff
###注意clade1中的Pseudom_alcalig_GCF_001597285_1和Pseudom_oti_GCF_011397855_1没有gbff文件
ls -lh branched-chain/collinearity/pseudomons_two_clade/clade2/*.gb

#查看10个菌株中PA1964对应的名字
Pseudom_aeru_PAO1                      PA1964        
Pseudom_men_S5_2_GCF_000733715_2       DW68_RS10725   WP_016486587.1
Pseudom_stu_GCF_019704535_1            PszF2a_RS10410   WP_003449774.1
Pseudom_fluo_GCF_900215245_1           CPH89_RS24310   WP_073637695.1
Pseudom_puti_NBRC_14164_GCF_000412675_1  PP4_RS14490   WP_016486587.1
Pseudom_syr_GCF_016694755_2             JJQ97_RS13305  WP_003395611.1

Pseudom_alcalig_GCF_001597285_1        A0T30_RS12715  WP_016486587.1
Pseudom_citro_GCF_001586155_1          PcP3B5_RS11985  WP_061561347.1
Pseudom_kna_GCF_009911755_1          GVI47_RS15020      WP_043252681.1
Pseudom_lal_GCF_008807375_1          FXN65_RS10545    WP_003449774.1
Pseudom_oti_GCF_011397855_1          PotMrB4_RS09930  WP_044409009.1

#截取PAO1的gbk
cd /mnt/d/pseudomonas_HGT/branched_easyfig/two_clade_50k
python  fetch_gbk.py -i Pseudom_aeru_PAO1.gbff -e 50000 -l PA1964 -o PA1094

#截取其他菌株gbff文件的gbk  
python  fetch_gbk.py -i Pseudom_men_S5_2_GCF_000733715_2.gbff -e 50000 -l DW68_RS10725 -o PA1964
python  fetch_gbk.py -i Pseudom_stu_GCF_019704535_1.gbff -e 50000 -l PszF2a_RS10410 -o PA1964
python  fetch_gbk.py -i Pseudom_fluo_GCF_900215245_1.gbff -e 50000 -l  CPH89_RS24310 -o PA1964
python  fetch_gbk.py -i Pseudom_puti_NBRC_14164_GCF_000412675_1.gbff -e 50000 -l PP4_RS14490 -o PA1964
python  fetch_gbk.py -i Pseudom_syr_GCF_016694755_2.gbff -e 50000 -l JJQ97_RS13305 -o PA1964

#截取其他菌株gb文件的gbk  
python  fetch_gbk.py -i Pseudom_alcalig_GCF_001597285_1.gb -e 50000 -l A0T30_RS12715 -o PA1964
python  fetch_gbk.py -i Pseudom_citro_GCF_001586155_1.gb -e 50000 -l PcP3B5_RS11985 -o PA1964
python  fetch_gbk.py -i Pseudom_kna_GCF_009911755_1.gb -e 50000 -l  GVI47_RS15020 -o PA1964
python  fetch_gbk.py -i Pseudom_lal_GCF_008807375_1.gb -e 50000 -l FXN65_RS10545 -o PA1964
python  fetch_gbk.py -i Pseudom_oti_GCF_011397855_1.gb -e 50000 -l PotMrB4_RS09930 -o PA1964


```








* other
```bash
#提取PAO1的braz和braB上下游各10000bp
perl -alne 'next until(/gene/);$_=~s/\(//g;$_=~s/\)//g;print"$_";' Pseudom_aeru_PAO1.gbff | uniq >PAO1.gbff
#查看PAO1的braz 
grep -B 1 braZ PAO1.gbff | perl -alne '/(\d*)\.\.(\d*)/;$sta=$1;$end=$2;print"$sta:$end";' | uniq 
#结果
2151755:2153068
#提取PAO1的braz
seqkit subseq  Pseudom_aeru_PAO1.fna -r 2141755:2151755 >PAO1_braz_front.fa
seqkit subseq  Pseudom_aeru_PAO1.fna -r 2153068:2163068 >PAO1_braz_after.fa
#查看PAO1的braB
grep -B 1 braB PAO1.gbff | perl -alne '/(\d*)\.\.(\d*)/;$sta=$1;$end=$2;print"$sta:$end";' | uniq 
#结果
1732545:1733858
#提取PAO1的braB
seqkit subseq  Pseudom_aeru_PAO1.fna -r 1722545:1732545 >PAO1_braB_front.fa
seqkit subseq  Pseudom_aeru_PAO1.fna -r 1733858:1743858 >PAO1_braB_after.fa
```

