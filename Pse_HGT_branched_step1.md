<!-- TOC -->

- [1. codes from Professor Wang](#1-codes-from-professor-wang)
  - [1.1 外类群：使用一些模型生物作为外类群，共有15种.后期会去掉放线菌门，衣原体和厚壁菌门仅进行gamma变形菌纲的分析](#11-外类群使用一些模型生物作为外类群共有15种后期会去掉放线菌门衣原体和厚壁菌门仅进行gamma变形菌纲的分析)
  - [1.2 species group](#12-species-group)
  - [1.3 菌株基因组蛋白信息文件](#13-菌株基因组蛋白信息文件)
  - [1.4 refer genome和representative genome和假单胞菌属的typical genome](#14-refer-genome和representative-genome和假单胞菌属的typical-genome)
- [2. 统计branched-chain在菌株中拷贝分布](#2-统计branched-chain在菌株中拷贝分布)
  - [2.1 使用hmmsearch抓取相同domain的基因或者基因家族](#21-使用hmmsearch抓取相同domain的基因或者基因家族)
  - [2.2 将branched-chain提取的蛋白序列与tigerfam数据库比对](#22-将branched-chain提取的蛋白序列与tigerfam数据库比对)
  - [2.3 将branched-chain提取的蛋白序列与pfam数据库比对](#23-将branched-chain提取的蛋白序列与pfam数据库比对)
  - [2.4 多轮diamond（将hmmer匹配的蛋白作为query序列与所有蛋白进行比对，再次进行筛选）all.replace.fa中有两条序列存在大量空字符(null)导致无法建库，需要删除](#24-多轮diamond将hmmer匹配的蛋白作为query序列与所有蛋白进行比对再次进行筛选allreplacefa中有两条序列存在大量空字符null导致无法建库需要删除)
  - [2.5 统计diamond比对和hmmer比对时不同species中BCAA个数](#25-统计diamond比对和hmmer比对时不同species中bcaa个数)
  - [2.6 菌株中丢失了BCAA记为0](#26-菌株中丢失了bcaa记为0)
  - [2.7 绘制统计不同物种中BCAA平均拷贝数的表格](#27-绘制统计不同物种中bcaa平均拷贝数的表格)
  - [2.8 重新统计去掉epsilon，alpha变形菌和放线菌门，支原体门,芽孢杆菌纲的砂眼衣原体等后的基因组数，物种数，BCAA拷贝数](#28-重新统计去掉epsilonalpha变形菌和放线菌门支原体门芽孢杆菌纲的砂眼衣原体等后的基因组数物种数bcaa拷贝数)
- [3. 两种树](#3-两种树)
  - [3.1铜绿假单胞菌物种内所有的braB和braZ序列建立基因树（不需要外类群,因为蛋白序列有一半是braB,有一半是braZ,无法区分) (OK)](#31铜绿假单胞菌物种内所有的brab和braz序列建立基因树不需要外类群因为蛋白序列有一半是brab有一半是braz无法区分-ok)
  - [3.2铜绿假单胞菌物种内所有的braB和braZ序列建立物种树（需要外类群)](#32铜绿假单胞菌物种内所有的brab和braz序列建立物种树需要外类群)
  - [3.3假单胞菌属内(模式菌株15+代表菌株533 共有548)建立基因树(需要外类群)](#33假单胞菌属内模式菌株15代表菌株533-共有548建立基因树需要外类群)
  - [3.4假单胞菌属内(模式菌株15+代表菌株533 共有548)建立物种树(OK)](#34假单胞菌属内模式菌株15代表菌株533-共有548建立物种树ok)
  - [3.5gamma变形菌纲内(模式菌株15+代表菌株533 共有548)建立基因树 (需要外类群) (OK)](#35gamma变形菌纲内模式菌株15代表菌株533-共有548建立基因树-需要外类群-ok)
  - [3.6gamma变形菌纲(模式菌株15+代表菌株533 共有548)建立物种树，alpha变形菌纲作为外类群，厚壁菌门作为外类群(金黄色葡萄球菌和枯草芽孢杆菌)](#36gamma变形菌纲模式菌株15代表菌株533-共有548建立物种树alpha变形菌纲作为外类群厚壁菌门作为外类群金黄色葡萄球菌和枯草芽孢杆菌)
- [4.生物环境](#4生物环境)
  - [4.1铜绿样本生物环境注释信息](#41铜绿样本生物环境注释信息)
  - [4.2注释信息自动生成配色](#42注释信息自动生成配色)
- [5.铜绿基因树的注释信息](#5铜绿基因树的注释信息)
- [6.共线性的注意事项](#6共线性的注意事项)
- [7.变形菌纲基因树和物种树的注释信息](#7变形菌纲基因树和物种树的注释信息)

<!-- /TOC -->
* 1.由于外类群不够，新添了一些基因组文件，由1514增加至1953个基因组文件，由此需要重新操作  
* 2.hmmsearch搜索结构域，hmmscan过滤,diamond重复计算拷贝数，绘制拷贝数的表格，绘制Reference Genome的树和Representative Genome与Reference Genome的树  
* 3.目前以branched为例,branched在pfam panthern和tigrfam中均存在    
# 1. codes from Professor Wang
## 1.1 外类群：使用一些模型生物作为外类群，共有15种.后期会去掉放线菌门，衣原体和厚壁菌门仅进行gamma变形菌纲的分析
| #tax\_id | organism\_name | phylum |
| :--- | :--- | :--- |
| 565050 | Caulobacter vibrioides NA1000 | Alphaproteobacteria |
| 192222 | Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819 | Epsilonproteobacteria |
| 208964 | Pseudomonas aeruginosa PAO1 | Gammaproteobacteria |
| 871585 | Acinetobacter pittii PHEA-2 | Gammaproteobacteria |
| 511145 | Escherichia coli str. K-12 substr. MG1655 | Gammaproteobacteria |
| 386585 | Escherichia coli O157:H7 str. Sakai | Gammaproteobacteria |
| 1125630 | Klebsiella pneumoniae subsp. pneumoniae HS11286 | Gammaproteobacteria |
| 99287 | Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 | Gammaproteobacteria |
| 198214 | Shigella flexneri 2a str. 301 | Gammaproteobacteria |
| 227377 | Coxiella burnetii RSA 493 | Gammaproteobacteria |
| 272561 | Chlamydia trachomatis D/UW-3/CX | Chlamydiae |
| 93061 | Staphylococcus aureus subsp. aureus NCTC 8325 | Firmicutes |
| 224308 | Bacillus subtilis subsp. subtilis str. 168 | Firmicutes |
| 169963 | Listeria monocytogenes EGD-e | Firmicutes |
| 83332 | Mycobacterium tuberculosis H37Rv | Actinobacteria |

## 1.2 species group
```bash
mkdir -p ~/data/Pseudomonas
cd ~/data/Pseudomonas

nwr member Pseudomonas |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat
#关于假单胞菌属，共有405个物种，6个species group
rank	count
genus	1
species	405
strain	747
subspecies	12
no rank	120
species group	6
species subgroup	5
isolate	1
#假单胞菌属共有405个物种，6个species group

nwr member Pseudomonas -r "species group"
#共有6个
#tax_id sci_name        rank    division
136849  Pseudomonas syringae group      species group   Bacteria
136846  Pseudomonas stutzeri group      species group   Bacteria
136845  Pseudomonas putida group        species group   Bacteria
136843  Pseudomonas fluorescens group   species group   Bacteria
136842  Pseudomonas chlororaphis group  species group   Bacteria
136841  Pseudomonas aeruginosa group    species group   Bacteria

```

## 1.3 菌株基因组蛋白信息文件
```shell
#包含NCBI的样本信息且去除了错误的菌株基因组文件信息，共1952个  
#assembly level:Complete Genome(1810)或Chromosome(142)
ASSEMBLY/Pseudomonas.assembly.pass.csv  

#包含菌株基因组名(group1)，属名(group4),共1952个  
strains.taxon.tsv  
#包含菌株名，共1952个  
strains.lst  
order.lst

#菌株基因组蛋白名(group1),菌株基因组名(group2)  
PROTEINS/all.strain.tsv  
#菌株基因组蛋白名和响应的蛋白序列  
PROTEINS/all.replace.fa  
#菌株基因组蛋白名，菌株基因组名，蛋白序列长度，蛋白注释  
PROTEINS/all.info.tsv  
#使用fasops concat将bac120相同菌株基因组名字的蛋白序列合并，然后去除低质量的序列，最后剩1952条基因组序列    
PROTEINS/bac120.trim.fa  
bac120.reroot.newick  

```

## 1.4 refer genome和representative genome和假单胞菌属的typical genome
```bash
#共有15个Reference Genome参考菌株，分别属于变形菌纲，厚壁菌纲，放线菌纲和衣原体
#共有533个Representative Genome代表菌株都属于Gamma变形菌纲
cut -d , -f 18 ASSEMBLY/Pseudomonas.assembly.pass.csv | tsv-summarize -g 1 --count
#提取参考菌株的基因组名称
cat ASSEMBLY/Pseudomonas.assembly.pass.csv  | grep "Reference Genome" |  
cut -d , -f 1 >reference.tsv
#提取代表菌株和模式菌株的基因组名名称
cut -d , -f 1,18 ASSEMBLY/Pseudomonas.assembly.pass.csv | tsv-filter --d , --str-in-fld 2:"Genome" |
cut -d , -f 1 >representative.tsv
#提取代表菌株的基因组名称并添加物种分类等信息
cat ASSEMBLY/Pseudomonas.assembly.pass.csv | grep "Representative Genome" | cut -d , -f 1 | tsv-join -d 1 -f strains.taxon.tsv -k 1  --append-fields 4 |  tsv-select -f 2,1 | nwr append stdin -r class -r order -r family -r genus  -r "species group" >representative_add_rank.tsv

#提取假单胞菌属的典型菌株名
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

```


# 2. 统计branched-chain在菌株中拷贝分布
## 2.1 使用hmmsearch抓取相同domain的基因或者基因家族
```bash
#查看branched-chain的hmm文件所在数据库
cd ~/data/Pseudomonas
cat STRAINS/Pseudom_aer_PAO1/*.tsv |
    grep "IPR004685"
#可以发现数据库中的登录号是
PANTHER PTHR30588
TIGRFAM TIGR00796
Pfam    PF05525

#下载hmm文件
mkdir -p branched-chain/HMM
curl -L http://www.pantherdb.org/panther/exportHmm.jsp?acc=PTHR30588 >branched-chain/HMM/branched.panther.hmm
curl -L https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR00796.1.HMM >branched-chain/HMM/branched.tigerfam.hmm
curl -L http://pfam.xfam.org/family/PF05525/hmm >branched-chain/HMM/branched.pfam.hmm

E_VALUE=1e-20
for domain in branched.panther branched.tigerfam branched.pfam  ; do
    >&2 echo "==> domain [${domain}]"

    if [ -e branched-chain/${domain}.replace.tsv ]; then
        continue;
    fi

    for ORDER in $(cat order.lst); do
        >&2 echo "==> ORDER [${ORDER}]"

        cat taxon/${ORDER} |
            parallel --no-run-if-empty --linebuffer -k -j 8 "
                gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                    hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw branched-chain/HMM/${domain}.hmm - |
                    grep '>>' |
                    perl -nl -e '
                        m{>>\s+(\S+)} or next;
                        \$n = \$1;
                        \$s = \$n;
                        \$s =~ s/\.\d+//;
                        printf qq{%s\t%s_%s\n}, \$n, {}, \$s;
                    '
            "
    done \
        > branched-chain/${domain}.replace.tsv

    >&2 echo
done

#统计不同数据库根据相同domain抓取出的基因或基因家族序列数目
wc -l branched-chain/*.tsv
       4280 branched-chain/branched.panther.replace.tsv
       4280 branched-chain/branched.pfam.replace.tsv
       4280 branched-chain/branched.tigerfam.replace.tsv
```
## 2.2 将branched-chain提取的蛋白序列与tigerfam数据库比对
```bash
#下载tigerfam数据库
mkdir -p ~/data/HMM/TIGRFAM
cd ~/data/HMM/TIGRFAM
wget -N --content-disposition ftp://ftp.jcvi.org/data/TIGRFAMs/14.0_Release/TIGRFAMs_14.0_HMM.tar.gz
tar -vzvf TIGRFAMs_14.0_HMM.tar.gz
cat  *.HMM >tigerfam.hmm

#格式化tigerfam数据库
cd ~/data/Pseudomonas
hmmpress ~/data/HMM/TIGRFAM/tigerfam.hmm

#提取抓出来的序列2140
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 branched-chain/branched.tigerfam.replace.tsv) branched-chain/branched-chain.fa

E_VALUE=1e-10
NAME=branched-chain
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.tigerfam.txt --tblout ${NAME}/${NAME}.tigerfam.tbl  \
   ~/data/HMM/TIGRFAM/tigerfam.hmm  ${NAME}/${NAME}.fa

perl abstract1.pl branched-chain/branched-chain.tigerfam.tbl >branched-chain/branched-chain.tigerfam.abstract.tsv
tsv-filter --le 4:1e-50 branched-chain/branched-chain.tigerfam.abstract.tsv >branched-chain/branched-chain.tigerfam.cutoff.tsv
tsv-summarize  -g 8 --count branched-chain/branched-chain.tigerfam.cutoff.tsv
#查看以e值过滤后的结果
livcs:_branched-chain_amino_acid_transport_system_II_carrier_protein    2140

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名 2140
perl compare1.pl  branched-chain/branched-chain.tigerfam.cutoff.tsv > branched-chain/branched-chain_tigerfam.minevalue.tsv

#拼接属名等信息并统计拷贝数
cat branched-chain/branched-chain_tigerfam.minevalue.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"branched-chain" |
tsv-join -d 1 \
-f PROTEINS/all.strain.tsv -k 1 \
--append-fields 2 |
tsv-join -d 3 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-summarize -g 3,4 --count |
keep-header -- tsv-sort -k3,3n >branched-chain/branched-chain_hmmscan_tigerfam.copy.tsv

#统计拷贝数的分布390
tsv-summarize -g 3,2 --count  branched-chain/branched-chain_hmmscan_tigerfam.copy.tsv > branched-chain/branched-chain_hmmscan_tigerfam.GCF_copy.tsv
sed -i '1icopy\tgenus\tGCF'  branched-chain/branched-chain_hmmscan_tigerfam.GCF_copy.tsv

```
## 2.3 将branched-chain提取的蛋白序列与pfam数据库比对
```bash
#下载pfam数据库
mkdir -p ~/data/HMM/PFAM
cd  ~/data/HMM/PFAM
for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
done

for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    echo "==> ${basename}"
    gzip -dcf ${basename}.gz > ${basename}
done

#格式化pfam数据库
cd ~/data/Pseudomonas
hmmpress ~/data/HMM/PFAM/Pfam-A.hmm

#提取抓出来的序列2140
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 branched-chain/branched.pfam.replace.tsv) branched-chain/branched-chain.pfam.fa
E_VALUE=1e-10
NAME=branched-chain
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.pfam.txt --tblout ${NAME}/${NAME}.pfam.tbl  \
   ~/data/HMM/PFAM/Pfam-A.hmm  ${NAME}/${NAME}.pfam.fa

perl abstract1.pl branched-chain/branched-chain.pfam.tbl >branched-chain/branched-chain.abstract.pfam.tsv
tsv-filter --le 4:1e-50 branched-chain/branched-chain.abstract.pfam.tsv >branched-chain/branched-chain.cutoff.pfam.tsv
tsv-summarize  -g 8 --count branched-chain/branched-chain.cutoff.pfam.tsv
#查看以e值过滤后的结果
Branched-chain_amino_acid_transport_protein     2140
#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare1.pl branched-chain/branched-chain.cutoff.pfam.tsv > branched-chain/branched-chain_minevalue.pfam.tsv

#拼接属名等信息并统计拷贝数
cat branched-chain/branched-chain_minevalue.pfam.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"Branched-chain" |
tsv-join -d 1 \
-f PROTEINS/all.strain.tsv -k 1 \
--append-fields 2 |
tsv-join -d 3 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-summarize -g 3,4 --count |
keep-header -- tsv-sort -k3,3n >branched-chain/branched-chain_hmmscan_copy.pfam.tsv

#统计拷贝数的分布390
tsv-summarize -g 3,2 --count  branched-chain/branched-chain_hmmscan_copy.pfam.tsv > branched-chain/branched-chain_hmmscan_GCF_copy.pfam.tsv
sed -i '1icopy\tgenus\tGCF'  branched-chain/branched-chain_hmmscan_GCF_copy.pfam.tsv

#查看pfam和tigerfam结果的交集2140
cat branched-chain/branched-chain_minevalue.pfam.tsv | grep -f <(cut -f 1 branched-chain/branched-chain_tigerfam.minevalue.tsv ) | wc -l
```
## 2.4 多轮diamond（将hmmer匹配的蛋白作为query序列与所有蛋白进行比对，再次进行筛选）all.replace.fa中有两条序列存在大量空字符(null)导致无法建库，需要删除
```bash
cd ~/data/Pseudomonas
#删除all.replace.fa的空字符
sed -i 's/\x0//g' PROTEINS/all.replace.fa
#提取hmmsearch和hmmscan后的结果
cat branched-chain/branched-chain_minevalue.pfam.tsv
#第一轮blastp 
mkdir  -p branched-chain/blastp
faops some PROTEINS/all.replace.fa  <(cut -f 1 branched-chain/branched-chain_minevalue.pfam.tsv) branched-chain/blastp/branched-chain_blastp1.fa
faops size branched-chain/blastp/branched-chain_blastp1.fa | wc -l  #2140
makeblastdb -in PROTEINS/all.replace.fa -dbtype prot -out branched-chain/blastp/whole
blastp -db branched-chain/blastp/whole -query branched-chain/blastp/branched-chain_blastp1.fa -out branched-chain/blastp/branched-chain_result1.tsv -outfmt 6 -evalue 1e-5
#第二轮blastp 
faops some PROTEINS/all.replace.fa  <(cut -f 2 branched-chain/blastp/branched-chain_result1.tsv | sort -n | uniq) branched-chain/blastp/branched-chain_blastp2.fa
blastp -db branched-chain/blastp/whole -query branched-chain/blastp/branched-chain_blastp2.fa -out branched-chain/blastp/branched-chain_result2.tsv -outfmt 6 -evalue 1e-5


#第一轮diamond(将hmmsearhc和hmmscan的结果放在所有的蛋白序列库里进行blast)
#相似性选择和覆盖度选择（结构阈330aa,braB437，330/437=0.755, 因此覆盖度选择50，相似性最低为55.606 ，因此相似性选50 ）
mkdir -p branched-chain/diamond
faops some PROTEINS/all.replace.fa  <(cut -f 1 branched-chain/branched-chain_minevalue.pfam.tsv)  branched-chain/diamond/branched-chain_diamond1.fa
faops size branched-chain/diamond/branched-chain_diamond1.fa | wc -l #2140
diamond makedb --in PROTEINS/all.replace.fa --db branched-chain/diamond/whole
diamond blastp --db branched-chain/diamond/whole.dmnd --query  branched-chain/diamond/branched-chain_diamond1.fa -e 1e-5 --outfmt 6 --threads 4 --out branched-chain/diamond/branched-chain_result1.tsv  
--id 50 --subject-cover 50
#第二轮diamond(将第一轮从蛋白数据库中抓取出的序列作为query序列继续进行blast)
faops some PROTEINS/all.replace.fa  <(cut -f 2 branched-chain/diamond/branched-chain_result1.tsv | sort -n | uniq) branched-chain/diamond/branched-chain_diamond2.fa
faops size branched-chain/diamond/branched-chain_diamond2.fa | wc -l #1492
diamond blastp --db branched-chain/diamond/whole.dmnd --query  branched-chain/diamond/branched-chain_diamond2.fa -e 1e-5 --outfmt 6 --threads 4 --out branched-chain/diamond/branched-chain_result2.tsv  --id 50 --subject-cover 50

#第一轮diamond的query
cut -f 1 branched-chain/diamond/branched-chain_result1.tsv | sort -n | uniq | wc -l #2140
#第一轮diamond的target
cut -f 2 branched-chain/diamond/branched-chain_result1.tsv | sort -n | uniq | wc -l #1492
#第二轮diamond的query
cut -f 1 branched-chain/diamond/branched-chain_result2.tsv | sort -n | uniq | wc -l #1492
#第二轮diamond的target
cut -f 2 branched-chain/diamond/branched-chain_result2.tsv | sort -n | uniq | wc -l #1494





#三轮diamond结果一致
#查看hmmer结果和diamond结果的相同之处2108
cat branched-chain/branched-chain_minevalue.tsv | cut -f 1 | grep -f <(cut -f 1 branched-chain/diamond/branched-chain_result4.tsv | sort -n | uniq)



```

## 2.5 统计diamond比对和hmmer比对时不同species中BCAA个数
```bash
#统计diamond抓取的braB个数
cat branched-chain/diamond/branched-chain_result3.tsv | tsv-filter --eq 3:100 |cut -f 1 | sort -n | uniq |
tsv-join -d 1 \
-f PROTEINS/all.strain.tsv -k 1 \
--append-fields 2 |
tsv-join -d 2 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 |
tsv-summarize -g 3 --count |
keep-header -- tsv-sort -k2,2n >branched-chain/diamond/branched_chain_diamond_speices_num.tsv

#统计hmmer抓取的braB个数
cat branched-chain/branched-chain_minevalue.pfam.tsv | tsv-select -f 1,3 | 
tsv-filter --str-in-fld 2:"branched-chain" | cut -f 1 | 
tsv-join -d 1 \
-f PROTEINS/all.strain.tsv -k 1 \
--append-fields 2 |
tsv-join -d 2 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 |
tsv-summarize -g 3 --count |
keep-header -- tsv-sort -k2,2n >branched-chain/diamond/branched_chain_hmmer_species_num.tsv
```

## 2.6 菌株中丢失了BCAA记为0
```bash
#统计assembly总个数
cat strains.taxon.tsv | cut -f 1 | sort -n | uniq | wc -l #1952 
#统计不同species的所有assembly #572
cat strains.taxon.tsv | tsv-summarize -g 4 --count | keep-header -- tsv-sort -k2,2n >branched-chain/diamond/assembly_species_num.tsv
#统计不同species的含有copy的assembly #386
wc -l branched-chain/diamond/branched_chain_hmmer_species_num.tsv  
##根据不同物种的菌株数量和含有copy的菌株数量差异，发现有的菌株存在基因丢失情况
##因此基因存在丢失的菌株需要将拷贝记为0
#差异数目为186
cat branched-chain/diamond/assembly_species_num.tsv | grep -v -f <(cut -f 1 branched-chain/diamond/branched_chain_hmmer_species_num.tsv) | tr "\t" "," | perl -F, -alne '$name=$F[0];$num=0;print"$name\t$num";'  >branched-chain/diamond/branched_chain_hmmer_species_num_addinformation.tsv
#合并未丢失和丢失的共有572个
cat branched-chain/diamond/branched_chain_hmmer_species_num.tsv  branched-chain/diamond/branched_chain_hmmer_species_num_addinformation.tsv >branched-chain/diamond/branched_chain_hmmer_species_copy_num.tsv
```

## 2.7 绘制统计不同物种中BCAA平均拷贝数的表格
* 表格内容
```bash
#species  | number of assemblies  | numbers of BCAA | average per genome 
#合并不同菌株的组装体数目和不同菌株中的BCAA拷贝数并且相除计算copy per genome(572)
cat branched-chain/diamond/branched_chain_hmmer_species_copy_num.tsv | 
tsv-join -d 1 \
-f branched-chain/diamond/assembly_species_num.tsv -k 1 \
--append-fields 2 |
tsv-select -f 1,3,2  | tr "\t" "," |
perl -F, -alne '$per=$F[2]/$F[1];$per=sprintf "%.1f",$per;print"$F[0]\t$F[1]\t$F[2]\t$per";' >branched-chain/diamond/branched_chain_species_assembly_copy.tsv
#添加属，科,目，纲和group species信息(572)
cat branched-chain/diamond/branched_chain_species_assembly_copy.tsv | nwr append stdin -r class -r order -r family -r genus -r "species group" \
 >branched-chain/diamond/branched_chain_species_mean_copy.tsv
#将不同物种的assembly个数和braz总个数及平均个数进行排序统计
sed -i '1ispecies\tnumber of assemblies\tnumber of braZ\taverage per genome\tclass\torder\tfamily\tgenus\tspecies group' \
branched-chain/diamond/branched_chain_species_mean_copy.tsv
mlr --itsv --ocsv cat  <(cat branched-chain/diamond/branched_chain_species_mean_copy.tsv |  keep-header -- tsv-sort -k3,3rn -k4,4rn  )\
 > branched-chain/diamond/branched_chain_species_mean_copy.csv
```
* 表格形式
```bash
#查看假单胞菌属的typical菌株的拷贝数目(9）
cut -f 1,2,3 branched-chain/diamond/branched_chain_species_mean_copy.tsv | grep -v "number" | grep -f <(cat typical.lst | tsv-join -d 1 -f strains.taxon.tsv -k 1 --append-fields 4 | cut -f 2 | sort -n | uniq) | tr "\t" "," | perl -F, -alne '$per=$F[2]/$F[1];$per=sprintf "%.1f",$per;print"$F[0]\t$F[1]\t$F[2]\t$per";' |  
tsv-sort -k1,1rn >branched-chain/diamond/branched_chain_copy_format1.tsv
#查看假单胞菌属的非typical菌株的拷贝数目(57) 共190
cut -f 1,2,3 branched-chain/diamond/branched_chain_species_mean_copy.tsv | grep "Pseudomonas" | grep -v -f <(cut -f 1 branched-chain/diamond/branched_chain_copy_format1.tsv) | nwr append stdin -r class | tsv-summarize -g 4 --sum 2,3 | sed 's/Gammaproteobacteria/Other Pseudomonas/g'  | tr "\t" "," | 
perl -F, -alne '$per=$F[2]/$F[1];$per=sprintf "%.1f",$per;print"$F[0]\t$F[1]\t$F[2]\t$per";'  >branched-chain/diamond/branched_chain_copy_format2.tsv
#提取gamma变形菌纲的非假单胞菌属的其他科的结果(506/499) 去除拷贝为0的科
cut -f 1,2,3,4,5,6,7 branched-chain/diamond/branched_chain_species_mean_copy.tsv |  grep "Gammaproteobacteria" | grep -v "Pseudomonas"  | tsv-select -f 7,2,3 |  
tsv-summarize -g 1 --sum 2,3 | tr "\t" "," | perl -F, -alne '$per=$F[2]/$F[1];$per=sprintf "%.1f",$per;print"$F[0]\t$F[1]\t$F[2]\t$per";' |  
tsv-sort -k1,1rn  >branched-chain/diamond/branched_chain_copy_format3.tsv
#合并三者的结果(由于只关注gamma变形菌纲，去除厚壁菌门的芽孢杆菌目和衣原体目)

cat branched-chain/diamond/branched_chain_copy_format1.tsv  branched-chain/diamond/branched_chain_copy_format2.tsv branched-chain/diamond/branched_chain_copy_format3.tsv |  sed '1iTaxonomy\tNumber of assemblies\tNumber of braB\braZ\tAverage per genome'  >branched-chain/diamond/branched_chain_copy_whole.tsv
mlr --itsv --ocsv cat branched-chain/diamond/branched_chain_copy_whole.tsv >branched-chain/diamond/branched_chain_copy_whole.csv
```
  
## 2.8 重新统计去掉epsilon，alpha变形菌和放线菌门，支原体门,芽孢杆菌纲的砂眼衣原体等后的基因组数，物种数，BCAA拷贝数
```bash
#只保留变形菌纲和芽孢杆菌纲的(金葡和枯草芽孢杆菌和李斯特菌)三个物种作为外类群
#统计剩余的基因组数量1947
cat strains.taxon.tsv | cut -f 4  | nwr append stdin -r class | grep -E "Gammaproteobacteria|Bacilli" | grep -v "Listeria monocytogenes" | wc -l
#统计剩余的基因组数量对应的物种数和family数量和order数量567和45
cat strains.taxon.tsv | cut -f 4  | nwr append stdin -r class | grep -E "Gammaproteobacteria|Bacilli" |grep -v "Listeria monocytogenes" | tsv-summarize -g 1 --count | wc -l
cat strains.taxon.tsv | cut -f 4  | nwr append stdin -r class | grep -E "Gammaproteobacteria|Bacilli" | grep -v "Listeria monocytogenes" |nwr append stdin -r family | tsv-summarize -g 3 --count | wc -l

#统计变形菌纲对应的基因组数目和物种数和family数量1945和565和43
cat strains.taxon.tsv | cut -f 4  | nwr append stdin -r class | grep  "Gammaproteobacteria" | wc -l
cat strains.taxon.tsv | cut -f 4  | nwr append stdin -r class | grep "Gammaproteobacteria" | tsv-summarize -g 1 --count | wc -l
cat strains.taxon.tsv | cut -f 4  | nwr append stdin -r class | grep "Gammaproteobacteria" | nwr append stdin -r family | tsv-summarize -g 3 --count | wc -l
#统计只有变形菌纲的BCAA拷贝数2134(-6)
cat branched-chain/branched-chain_minevalue.pfam.tsv | cut -f 1 |
tsv-join -d 1 \
-f PROTEINS/all.strain.tsv -k 1 \
--append-fields 2 |
tsv-join -d 2 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | tsv-select -f 3,1,2 | nwr append stdin -r class | grep  "Gammaproteobacteria" | wc -l 
#统计只有变形菌纲的BCAA拷贝数对应的基因组数量1662(-3)
cat branched-chain/branched-chain_minevalue.pfam.tsv | cut -f 1 |
tsv-join -d 1 \
-f PROTEINS/all.strain.tsv -k 1 \
--append-fields 2 |
tsv-join -d 2 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | tsv-select -f 3,1,2 | nwr append stdin -r class | grep  "Gammaproteobacteria" | tsv-summarize -g 3 --count | wc -l
#统计只有变形菌纲的BCAA拷贝数对应的物种数量 383(-3)
cat branched-chain/branched-chain_minevalue.pfam.tsv | cut -f 1 |
tsv-join -d 1 \
-f PROTEINS/all.strain.tsv -k 1 \
--append-fields 2 |
tsv-join -d 2 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | tsv-select -f 3,1,2 | nwr append stdin -r class | grep  "Gammaproteobacteria" | tsv-summarize -g 1 --count | wc -l
#统计只有变形菌纲的BCAA拷贝数对应的科的数量 18
cat branched-chain/branched-chain_minevalue.pfam.tsv | cut -f 1 |
tsv-join -d 1 \
-f PROTEINS/all.strain.tsv -k 1 \
--append-fields 2 |
tsv-join -d 2 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | tsv-select -f 3,1,2 | nwr append stdin -r class | grep  "Gammaproteobacteria" |nwr append stdin -r family | tsv-summarize -g 5 --count | wc -l
#统计变形菌纲没有鉴定到BCAA的基因组数目283
cat branched-chain/branched-chain_minevalue.pfam.tsv | cut -f 1 |
tsv-join -d 1 \
-f PROTEINS/all.strain.tsv -k 1 \
--append-fields 2 |
tsv-join -d 2 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | tsv-select -f 3,1,2 | nwr append stdin -r class | grep  "Gammaproteobacteria" | cut -f 3 >1.tsv
cat strains.taxon.tsv | grep -v -f 1.tsv | tsv-select -f 4,1 | nwr append stdin -r class | grep "Gammaproteobacteria" | tsv-summarize -g 2 --count | wc -l
#统计变形菌纲没有鉴定到BCAA的物种数目185
cat strains.taxon.tsv | grep -v -f 1.tsv | tsv-select -f 4,1 | nwr append stdin -r class | grep "Gammaproteobacteria" | tsv-summarize -g 1 --count | wc -l
#统计变形菌纲没有鉴定到BCAA的科数目85
cat strains.taxon.tsv | grep -v -f 1.tsv | tsv-select -f 4,1 | nwr append stdin -r class -r family | grep "Gammaproteobacteria" | tsv-summarize -g 4 --count | wc -l
```




# 3. 两种树
* [系统发育和系统发育生态学标记PhyEco marker](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0077033)
* [系统识别细菌和古细菌不同分类水平的系统发育标记](https://figshare.com/articles/dataset/Systematically_identify_phylogenetic_markers_at_different_taxonomic_levels_for_bacteria_and_archaea/722713/1)
* 在不同分类水平上自动识别系统发育标记的协议，该协议使用快速搜索和聚类算法为选择的基因组生成蛋白质家族,然后建立系统发育树，并自动对来自树的进化枝进行采样和评估，以评估使它们对系统发育分析(例如感兴趣的基因组的普遍性）和生态学研究（例如复制均匀性），然后使用多重比较和系统发育分析进一步评估潜在的标记家族.使用该种方法，作者已经为所有细菌确定了114个系统发育标记,其中包括40个同时涵盖细菌和古细菌的标记
* [举例：Izimaplasma的系统发育分析](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5113845/)
* 基因组树构建:基因组树是由38个蛋白质编码基因的串联比对构建的，这些细菌在古细菌和细菌中普遍分布并且以单拷贝形式存在.该树的基础包括来自NCBI和IMG数据库的草图和完成的基因组，使用hmmer在每个基因组中鉴定标记基因.使用默认设置的FastTree创建了最大似然树.
* 基因树构建:使用来自Tigrfams的HMM模型和IMG基因组下载的基因组中鉴定RpoB或EF-tu基因序列，使用muscle比对，使用FastTree构建系统发育树
* 建树(使用raxml)  brew install raxml   raxml版本的选择：Sequential 版本适合于中小型的数据； PThreads 版本适合于长序列或多条序列；MPI 版本适合于较大(100~1000) bootstraps 的运行。  RAxML 是用极大似然法建立进化树的软件之一，可以处理超大规模的序列数据，包括上千至上万个物种，几百至上万个已经比对好的碱基序列。
* 建树法包括(raxml iq-tree fasttree phyml四个最大似然法maximum likehood方法) .FastTtree采用的是SH检验来判断每个节点的可信度。该值的范围在0~1之间，与一般用的bootstrap值高度相关。 
* 比对速度（Muscle>MAFFT>ClustalW>T-Coffee） 比对准确性（MAFFT>Muscle>T-Coffee>ClustalW） 但是实际运行并非如上述文献引用，反而是mafft的运行速度远高于muscle
* 比对后处理:PhyDE 序列编辑软件，可以多开，方便多个比对结果之间进行比较。FastGap 可以对比对后的序列的gap进行重编码，提高序列信息使用效率。Gblocks 在线工具，选择序列保守区，使得后续系统发育分析免受变异过大的比对区域的不良影响。DAMBE 一个低调但全能的系统发育软件，定位与Mega类似，包括饱和度检测功能点。DNAsp 序列分析软件，计算各种序列多样性数据，如核苷酸多样性、序列信息位点含量、单倍型等
* 目前常用的构建系统发育树的方法有：邻位归并法(Neighbor joining, NJ)、最大似然法(Maximum likelihood method, ML) 以及贝叶斯法（BI）。综合速度和准确度，ML用得较多。ML对替代模型非常敏感，因此利用ML法构建系统发育树之前，选择合适的替代模型是必不可少的过程。(如果序列的相似度较高，每种方法和模型构建的系统发育树差别不大)

## 3.1铜绿假单胞菌物种内所有的braB和braZ序列建立基因树（不需要外类群,因为蛋白序列有一半是braB,有一半是braZ,无法区分) (OK)
```bash
cd ~/data/Pseudomonas
#提取所有的铜绿假单胞菌的单双拷贝蛋白序列(773)
mkdir -p ~/data/Pseudomonas/branched-chain/tree
faops some PROTEINS/all.replace.fa <(cat branched-chain/branched-chain_minevalue.tsv |grep -f <(cat branched-chain/branched-chain_hmmscan_copy.pfam.tsv | grep "Pseudom_aeru" | cut -f 1) | cut -f 1 ) branched-chain/tree/branched-chain.Pseudom_aeru.protein.fa
#查看773
faops size branched-chain/tree/branched-chain.Pseudom_aeru.protein.fa | wc -l
#比对
mafft --auto   branched-chain/tree/branched-chain.Pseudom_aeru.protein.fa > branched-chain/tree/branched-chain.Pseudom_aeru.aln.mafft.fa
#使用iqtree2建树（超算上)773
bsub -q mpi -n 24 -J "iq" ./iqtree2 -s branched-chain.Pseudom_aeru.aln.mafft.fa  -m MFP  --prefix branched-chain.Pseudom_aeru -T 20 -b 100
#改名
mv branched-chain.Pseudom_aeru.treefile branched-chain.Pseudom_aeru.newick


#modeltest-ng查看最佳建树模型
modeltest-ng -d nt -i 输入fas -o 输出名称 -T mrbayes/raxml --force -p 2
#-d 表示数据类型，nt表示核酸，aa表示蛋白质
#-i 表示输入文件，可以是fas或者phylip格式
#-o 表示输出名称
#-T 表示模版，有raxml, phyml, mrbayes, paup 可以选择
#-p 表示处理器使用个数
#可以在 log 文件中查看输出的最佳模型。

#raxml-ng建树速度提升，得到更准确的结果，支持22个核酸进化模型，多线程效率提升
brew install  raxml-ng 
#第一步：检查比对后的序列MSA是否可以读取（MSA可以使用FASTA、PHYLIP格式）。注：这一步骤还会给出哪些序列是相同序列（推荐执行）
raxml-ng --check --msa prim.phy --model GTR --prefix T1
#第二步：构建核苷酸树（模型GTR，1000次自检值抽样）。 注：蛋白树可用LG模型 
raxml-ng --all --msa prim.phy --model LG --prefix T15 --threads 10  --bs-trees 1000
#-msa 比对好的序列文件
#--model 模型，可以参照 modeltest-ng 的输出
#--prefix 前缀
#--bs-trees BootStrap 树的个数

#重定根
#nw_stats显示树的属性
#nw_clade： 提取由节点标签指定的子树
#nw_labels：打印节点标签
#nw_order：排序
#nw_reroot：重置根
#nw_rename：标签重命名

#使用iqtree2建树
wget -c https://github.com/iqtree/iqtree2/releases/download/v2.1.3/iqtree-2.1.3-Linux.tar.gz
tar -zxvf iqtree-2.1.3-Linux.tar.gz
echo 'export PATH="~/iqtree-2.1.3-Linux/bin:$PATH"'>>~/.bashrc
source ~/.bashrc
cd iqtree-2.1.3-Linux/bin
#根据ModelFinder自动选择的最适合模型，从序列比对(example.phy)中推断出最大似然树：
#-s指定 IQ-TREE 需要的比对文件的名称,可选PHYLIP、FASTA、NEXUS、CLUSTAL、MSF
#-m是指定要在分析期间使用的模型名称的选项。特殊MFP关键字代表ModelFinder Plus，它告诉IQ-TREE使用所选模型执行ModelFinder和剩余分析.ModelFinder计算许多不同模型和Akaike信息准则(AIC)、校正Akaike信息准则(AICc)和贝叶斯信息准则(BIC)的初始简约树的对数似然。然后ModelFinder选择最小化BIC分数的模型（您也可以通过添加选项-AIC或-AICc分别更改为AIC或AICc）。当使用了ModelFinder模型时，IQ-TREE将写入一个附加文件：example.phy.model：所有测试模型的对数似然。它用作检查点文件以恢复中断的模型选择。如果您现在查看，example.phy.iqtree您将看到IQ-TREE被选TIM2+I+G4为此示例数据的最佳拟合模型
#-o：外类群列表，不同物种之间逗号隔开
#--prefix 更改输出前缀
#--seed：随机数种子，主要出于调试目的
#-T：程序运行使用的核数，可设置具体数字或者AUTO（推荐），默认为1
#-b：非参数bootstrap次数，大于等于100
#-B：超快速bootstrap次数，大于等于1000
#--bnni：使用NNI优化超快速bootstrap的树，搭配-B使用
#--alrt：SH近似似然比检验重复次数
#-m：模型选择，设置MF自动选择最佳模型但不建树；设置MFP自动检测最佳模型并建树。此外还可以设置具体的模型，或者多个可选模型，例如-m LG,WAG
#--ancestral：基于经验贝叶斯的祖先状态重建
#正常情况下，在运行结束时IQ-TREE将写入几个输出文件，包括：
#example.phy.iqtree：自读的主要报告文件。您可以查看此文件以查看计算结果。它还包含最终树的文本表示（见下文）。
#example.phy.treefile：NEWICK格式的ML树，可以通过任何支持的树查看器程序（如FigTree或iTOL）进行可视化。
#example.phy.log：整个运行的日志文件（也打印在屏幕上）。要报告错误，请将此日志文件和原始比对文件发送给作者。

#-m 设置MFP自动检测最佳模型并建树
iqtree2 -s example.phy -m MFP  --prefix myprefix -T 20
#使用bootstrap自助法计算节点支持率（类似于RaxML）：
iqtree -s example.phy -m MFP -b 100 -T AUTO
#使用超快速bootstrap自助法计算节点支持率：
iqtree -s example.phy -m MFP -B 1000 --bnni -T AUTO
#使用上述设置构建500个基因组的120个串联蛋白树需要两天左右。超快自助法ultrafast bootstrap1000次比普通自助法100次要快10倍左右，是该软件的特有算法，所以一般使用Iqtree的超快自助法建树。
```


## 3.2铜绿假单胞菌物种内所有的braB和braZ序列建立物种树（需要外类群)
```bash
cd ~/data/Pseudomonas
#外类群菌株名(恶臭，丁香，荧光)
cat branched-chain/branched-chain_minevalue.tsv | cut -f 1 | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 | cut -f 2 | sort -n | uniq | grep "Pseudom_puti"            #恶臭Pseudom_puti_GCF_000691565_1
cat branched-chain/branched-chain_minevalue.tsv | cut -f 1 | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 | cut -f 2 | sort -n | uniq |  grep "Pseudom_syr"            #丁香 Pseudom_syr_GCF_004323795_1  
cat branched-chain/branched-chain_minevalue.tsv | cut -f 1 | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 | cut -f 2 | sort -n | uniq |  grep "Pseudom_flu"            #荧光 Pseudom_fluo_GCF_000730425_1 
echo -e "Pseudom_puti_GCF_000691565_1\nPseudom_syr_GCF_004323795_1\nPseudom_fluo_GCF_000730425_1" >branched-chain/tree/branched-chain.Pseudom_aeru.bac120.outgroup.tsv
faops some PROTEINS/bac120.trim.fa branched-chain/tree/branched-chain.Pseudom_aeru.bac120.outgroup.tsv branched-chain/tree/branched-chain.Pseudom_aeru.bac120.outgroup.fa
#提取铜绿假单胞菌的物种的菌株名字 #391
cat branched-chain/branched-chain_minevalue.tsv | grep "Pseudom_aeru" | cut -f 1 | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 | cut -f 2 | sort -n | uniq >branched-chain/tree/branched-chain.Pseudom_aeru.bac120.species.tsv
#提取上述菌株名相应的bac120序列 #391
faops some PROTEINS/bac120.trim.fa branched-chain/tree/branched-chain.Pseudom_aeru.bac120.species.tsv  branched-chain/tree/branched-chain.Pseudom_aeru.bac120.species.fa
#合并铜绿假单胞菌序列和外类群 #394
cat branched-chain/tree/branched-chain.Pseudom_aeru.bac120.outgroup.fa >>branched-chain/tree/branched-chain.Pseudom_aeru.bac120.species.fa
#mafft比对 3 个参数都设置为最不消耗时间的类型，适合于 ~10,000 到 ~50,000 条序列的比对。 mafft --retree 1 --maxiterate 0 --nofft
mafft --auto branched-chain/tree/branched-chain.Pseudom_aeru.bac120.species.fa >branched-chain/tree/branched-chain.Pseudom_aeru.bac120.aln.mafft.fa
#使用iqtree2建树（超算上)#394
bsub -q mpi -n 24 -J "iq" ./iqtree2 -s branched-chain.Pseudom_aeru.bac120.aln.mafft.fa  -m MFP  --prefix branched-chain.Pseudom_aeru.120.aln -T 20 -B 1000 --bnni  -o Pseudom_puti_GCF_000691565_1,Pseudom_syr_GCF_004323795_1,Pseudom_fluo_GCF_000730425_1
#改名
mv branched-chain.Pseudom_aeru.bac120.aln.treefile branched-chain.Pseudom_aeru.bac120.newick
```

## 3.3假单胞菌属内(模式菌株15+代表菌株533 共有548)建立基因树(需要外类群)
```bash
cd ~/data/Pseudomonas
#外类群序列名(只选择2个，和braB和braZ的关系都较远)
Shewanella baltica
Vibrio cholerae
#提取外类群物种菌株对应的GCF名并且每个外类群仅保留一个
cat  strains.taxon.tsv | grep -f <(echo -e "Shewanella baltica\nVibrio cholerae") | cut -f 1 
echo -e 'She_balt_GCF_003030925_1\nVi_cho_GCF_008369605_1' >branched-chain/tree/branched-chain.Pseudomonas.bac120.outgroup.tsv
#提取外类群菌株对应的GCF_WP名字
cat branched-chain/branched-chain_minevalue.tsv | grep -f branched-chain/tree/branched-chain.Pseudomonas.bac120.outgroup.tsv | cut -f 1 >branched-chain/tree/branched-chain.Pseudomonas.protein.outgroup.tsv
She_balt_GCF_003030925_1_WP_006084009
Vi_cho_GCF_008369605_1_WP_000815020
#提取外类群序列
faops some PROTEINS/all.replace.fa branched-chain/tree/branched-chain.Pseudomonas.protein.outgroup.tsv branched-chain/tree/branched-chain.Pseudomonas.protein.outgroup.fa

#提取假单胞菌属的的（模式菌株+代表菌株)的菌株名字#49
cat  representative.tsv | tsv-join -d 1 -f strains.taxon.tsv -k 1 --append-fields 4 | tsv-select -f 2,1 | nwr append stdin -r genus | tsv-filter --str-in-fld 3:"Pseudomonas" | sort -n | uniq >branched-chain/tree/branched-chain.Pseudomonas.protein.tsv
#提取上述菌株名相应的braz或braB蛋白(50) #多出一个蛋白是因为PAO1里面有两个拷贝braB和braZ   #多出了一个Pseudom_aeru_PAO1_GCF_013001005_1舍弃
faops some PROTEINS/all.replace.fa <(cat branched-chain/branched-chain_minevalue.tsv | grep -f <(cut -f 2 branched-chain/tree/branched-chain.Pseudomonas.protein.tsv) | cut -f 1| grep -v "Pseudom_aeru_PAO1_GCF_013001005_1") branched-chain/tree/branched-chain.Pseudomonas.protein.fa
#验证蛋白对应的菌株数(49)  
cat branched-chain/branched-chain_minevalue.tsv |  grep -v "Pseudom_aeru_PAO1_GCF_013001005_1" | grep -f <(cut -f 2 branched-chain/tree/branched-chain.Pseudomonas.protein.tsv) | cut -f 1 | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 | cut -f 2 | sort -n | uniq | wc -l
#将外类群合并到铜绿序列中(52)
cat branched-chain/tree/branched-chain.Pseudomonas.protein.outgroup.fa >>branched-chain/tree/branched-chain.Pseudomonas.protein.fa
#比对(52)
mafft --auto  branched-chain/tree/branched-chain.Pseudomonas.protein.fa >branched-chain/tree/branched-chain.Pseudomonas.protein.aln.fa
#建树
bsub -q mpi -n 24 -J "iq" ./iqtree2 -s branched-chain.Pseudomonas.protein.aln.fa -m MFP  --prefix branched-chain.Pseudomonas.protein -T 20 -B 1000 --bnni   -o  She_balt_GCF_003030925_1_WP_006084009,Vi_cho_GCF_008369605_1_WP_000815020
#改名
mv  branched-chain.Pseudomonas.protein.treefile   branched-chain.Pseudomonas.protein.aln.newick
```

## 3.4假单胞菌属内(模式菌株15+代表菌株533 共有548)建立物种树(OK)
```bash
cd ~/data/Pseudomonas
#提取假单胞菌属的的（模式菌株+代表菌株)的菌株名字 #49
cat branched-chain/branched-chain_minevalue.tsv| grep -v "Pseudom_aeru_PAO1_GCF_013001005_1" | grep -f <(cut -f 2 branched-chain/tree/branched-chain.Pseudomonas.protein.tsv) | cut -f 1 | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 | cut -f 2 | sort -n | uniq >branched-chain/tree/branched-chain.Pseudomonas.bac120.species.tsv
#在假单胞菌属的物种中加上外类群 #51
cat branched-chain/tree/branched-chain.Pseudomonas.bac120.outgroup.tsv >>branched-chain/tree/branched-chain.Pseudomonas.bac120.species.tsv
#提取上述菌株名相应的bac120序列 #51
faops some PROTEINS/bac120.trim.fa branched-chain/tree/branched-chain.Pseudomonas.bac120.species.tsv branched-chain/tree/branched-chain.Pseudomonas.bac120.species.fa
#比对(51)
mafft --auto  branched-chain/tree/branched-chain.Pseudomonas.bac120.species.fa >branched-chain/tree/branched-chain.Pseudomonas.bac120.aln.mafft.fa
#建树
bsub -q mpi -n 24 -J "iq" ./iqtree2 -s branched-chain.Pseudomonas.bac120.aln.mafft.fa -m MFP  --prefix branched-chain.Pseudomonas.bac120 -T 20 -B 1000 --bnni  -o  She_balt_GCF_003030925_1,Vi_cho_GCF_008369605_1
#改名
mv branched-chain.Pseudomonas.bac120.treefile branched-chain.Pseudomonas.bac120.newick
```


## 3.5gamma变形菌纲内(模式菌株15+代表菌株533 共有548)建立基因树 (需要外类群) (OK)
```bash
cd ~/data/Pseudomonas
#选择非gamma变形菌纲的作为外类群
cat   representative.tsv  | tsv-join -d 1 -f strains.taxon.tsv -k 1 --append-fields 4 | tsv-select -f 2,1 | nwr append stdin -r phylum -r class |  grep -v "Gammaproteobacteria"
#结果(去除衣原体门Chlamydiia)(共6个外类群)
Bacillus subtilis       Bac_subti_subtilis_168  Firmicutes      Bacilli
Campylobacter jejuni    Cam_jej_jejuni_NCTC_11168_ATCC_700819   Proteobacteria  Epsilonproteobacteria
Caulobacter vibrioides  Cau_vib_NA1000  Proteobacteria  Alphaproteobacteria
Chlamydia trachomatis   Chl_tracho_D_UW_3_CX    Chlamydiae      Chlamydiia
Listeria monocytogenes  Lis_mono_EGD_e  Firmicutes      Bacilli
Mycobacterium tuberculosis      My_tube_H37Rv   Actinobacteria  Actinomycetia
Staphylococcus aureus   Sta_aure_aureus_NCTC_8325       Firmicutes      Bacilli
#提取外类群菌株的名字
cat branched-chain/branched-chain_minevalue.tsv | grep -f <(cat   representative.tsv  | tsv-join -d 1 -f strains.taxon.tsv -k 1 --append-fields 4 | tsv-select -f 2,1 | nwr append stdin -r phylum -r class | grep -v "Gammaproteobacteria" | cut -f 2) | cut -f 1 |tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 | cut -f 2 >branched-chain/tree/branched-chain.Gammaproteobacteria.outgroup.tsv
#共有6个3种菌株名字,需要去除支原体，最后只剩下两个外类群菌株
Bac_subti_subtilis_168
Sta_aure_aureus_NCTC_8325
##添加外类群蛋白序列名称
cat branched-chain/branched-chain_minevalue.tsv | grep -f <(cat branched-chain/tree/branched-chain.Gammaproteobacteria.outgroup.tsv |  grep -v "Chl_tracho_D_UW_3_CX" ) |cut -f 1 |sort -n | uniq  
Bac_subti_subtilis_168_NP_390546
Sta_aure_aureus_NCTC_8325_YP_498750
#提取变形菌纲的的（模式菌株+代表菌株)的菌株名字(#541个)
cat  representative.tsv | tsv-join -d 1 -f strains.taxon.tsv -k 1 --append-fields 4 | tsv-select -f 2,1 | nwr append stdin -r class | tsv-filter --str-in-fld 3:"Gammaproteobacteria" | sort -n | uniq >branched-chain/tree/branched-chain.Gammaproteobacteria.strain.tsv
#提取上述菌株名相应的braz或braB蛋白名字(435)
cat branched-chain/branched-chain_minevalue.tsv | grep -f <(cut -f 2 branched-chain/tree/branched-chain.Gammaproteobacteria.strain.tsv) | cut -f 1 >branched-chain/tree/branched-chain.Gammaproteobacteria.protein.tsv
#加上外类群437
echo -e 'Bac_subti_subtilis_168_NP_390546\nSta_aure_aureus_NCTC_8325_YP_498750' >>branched-chain/tree/branched-chain.Gammaproteobacteria.protein.tsv
faops some PROTEINS/all.replace.fa <(cat branched-chain/branched-chain_minevalue.tsv | grep -f <(cut -f 2 branched-chain/tree/branched-chain.Gammaproteobacteria.protein.tsv) | cut -f 1) branched-chain/tree/branched-chain.Gammaproteobacteria.protein.fa
#比对(437)
mafft --auto  branched-chain/tree/branched-chain.Gammaproteobacteria.protein.fa > branched-chain/tree/branched-chain.Gammaproteobacteria.protein.aln.mafft.fa
#建树
bsub -q mpi -n 24 -J "iq" ./iqtree2 -s branched-chain.Gammaproteobacteria.protein.aln.mafft.fa -m MFP  --prefix branched-chain.Gammaproteobacteria.protein  -T 20 -B 1000 -bnni -o Bac_subti_subtilis_168_NP_390546,Sta_aure_aureus_NCTC_8325_YP_498750
#改名
mv branched-chain.Gammaproteobacteria.protein.treefile branched-chain.Gammaproteobacteria.protein.newick
```


## 3.6gamma变形菌纲(模式菌株15+代表菌株533 共有548)建立物种树，alpha变形菌纲作为外类群，厚壁菌门作为外类群(金黄色葡萄球菌和枯草芽孢杆菌)
```bash
cd ~/data/Pseudomonas
#含有braB和braZ的变形菌纲的（模式菌株+代表菌株)的菌株名字 #364
cat branched-chain/branched-chain_minevalue.tsv | grep -f <(cut -f 2 branched-chain/tree/branched-chain.Gammaproteobacteria.strain.tsv) | cut -f 1 | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 | cut -f 2 | sort -n | uniq >branched-chain/tree/branched-chain.Gammaproteobacteria.bac120.tsv
#验证braB和braZ在变形菌纲的数目#364
cat branched-chain/branched-chain_minevalue.tsv | grep -f representative.tsv | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 | cut -f 4 | tsv-join -d 1 -f strains.taxon.tsv -k 1 --append-fields 4 | tsv-select -f 2,1 | nwr append stdin -r class | tsv-filter --str-in-fld 3:"Gammaproteobacteria" | sort -n | uniq | wc -l
#加上外类群名字366
echo -e 'Bac_subti_subtilis_168_NP_390546\nSta_aure_aureus_NCTC_8325_YP_498750' >>branched-chain/tree/branched-chain.Gammaproteobacteria.bac120.tsv
#提取上述菌株名相应的bac120序列 #366
faops some PROTEINS/bac120.trim.fa branched-chain/tree/branched-chain.Gammaproteobacteria.bac120.tsv branched-chain/tree/branched-chain.Gammaproteobacteria.bac120.fa
mafft --auto branched-chain/tree/branched-chain.Gammaproteobacteria.bac120.fa >branched-chain/tree/branched-chain.Gammaproteobacteria.bac120.mafft.aln.fa
#建树
bsub -q mpi -n 24 -J "iq" ./iqtree2 -s branched-chain.Gammaproteobacteria.bac120.mafft.aln.fa -m MFP  --prefix branched-chain.Gammaproteobacteria.bac120 -T 20  -B 1000 --bnni -o  Bac_subti_subtilis_168,Sta_aure_aureus_NCTC_8325
#改名
mv branched-chain.Gammaproteobacteria.bac120.treefile branched-chain.Gammaproteobacteria.bac120.newick

```

# 4.生物环境
## 4.1铜绿样本生物环境注释信息
```BASH
#共有391个铜绿组装体(有390个有样本信息)(有373个有生境信息)
for filename in *.txt
do
base=$(basename $filename .txt)
echo $base
cat $base.txt | grep -A 20 "Pseudomonas aeruginosa" | grep "/" |grep -E "isolation source|environmental medium|description|host">$base.json
done  

#提取样本环境信息
for filename in *.json
do
base=$(basename $filename .json)
sample=$(cat $base.json | cut -d "\"" -f 2)
echo -e "$base\t$sample"
done >env_sample.tsv

#选择样本信息和菌株信息391
 cat Pseudomonas.assembly.pass.csv| grep "Pseudom_aeru"| cut -d "," -f 1,2,6 | tr "," "\t" >biosample_strain.tsv
#拼接样本环境信息和菌株信息 347
cat env_sample.tsv | tsv-join -d 1 -f biosample_strain.tsv -k 3 --append-fields 1 | tsv-select -f 1,3,2 >biosample_strain_env.tsv
#最终手工分类注释的信息表!!!!!!
biosamle_env_anno.tsv
cut -f 2,4 biosample_env_anno.tsv |sed '1itip\tclass1' >biosample_env_class.tsv 

#添加拷贝数信息
cat branched-chain/branched-chain_minevalue.pfam.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"branched-chain" |
tsv-join -d 1 \
-f PROTEINS/all.strain.tsv -k 1 \
--append-fields 2 | grep "Pseudom_aeru" |tsv-summarize -g 3 --count >branched-chain/branched-chain_Pseudom.aeru.copy.num.tsv

#cat branched-chain_Pseudom.aeru.copy.num.tsv  | perl -alne '$_=~s/\t1/\tA/g;$_=~s/\t2/\tB/g;print"$_";' | sed '1itip\tclass2' >biosample_env_class2.tsv
cut -f 2,4  biosample_env_copy_special.tsv | sed 's/one braB-like copy and one braZ-like copy/braB_braZ/g' | sed 's/one braB-like copy without braZ-like copies/one_braB/g' | sed 's/one braZ-like copy without braB-like copies/one_braZ/g'  | sed 's/two braZ-like copies without braB-like copies/two_braZ/g'  | sed '1itip\tclass2' >biosample_env_class2.tsv
```
## 4.2注释信息自动生成配色
* table2itol是在GitHup上公开的R语言包，其作用是专门为iTOL生成所需的注释文件，只需准备表格形式的数据，包含配色方案的注释文件就会自动生成，极大提高了准备注释文件的效率。
* [table2itol地址](https://github.com/mgoeker/table2itol)
```r
#linux安装
wget https://github.com/mgoeker/table2itol/archive/master.zip 
unzip master.zip 
mv table2itol-master table2itol
## 测试 
chmod +x table2itol/table2itol.R  
Rscript ./table2itol/table2itol.R 
#optparse--强大的命令行参数处理包 
#optparse使得我们很方便地给R脚本设置命令行参数，从而使得R脚本的复用和流程化使用
install.packages("optparse")
library("optparse")
Rscript ./table2itol/table2itol.R  -D plan1 -i tip  -b class  -w 0.5 biosample_env_class.tsv biosample_env_class2.tsv
-D 输出目录
-b设置背景色
--na-strings X 颜色圈放在名称外面
 -w 指定颜色带和宽度区域
 -a 找不到输入列将终止运行（默认不执行）
 -c 将整数列转换为factor或具有小数点的数字，
 -t 偏离提示标签时转换ID列，
 -w 颜色带，区域宽度等， 
 -D输出目录，
 -i OTU列名，
 -l OTU显示名称如种/属/科名，
#第一个文件第一行是表头，第一列是进化树中的叶标签tip，第二列为分类信息class1
#第二个文件第一行是表头，第一列是进化树中的叶标签tip，第二列为分类信息class2

#手工替换颜色
sed -i 's/#1f78b4/#e040fb/g' iTOL_colorstrip-class1.txt
sed -i 's/#fb9a99/#ffb74d/g' iTOL_colorstrip-class1.txt
sed -i 's/#33a02c/#b2ff59/g' iTOL_colorstrip-class1.txt
sed -i 's/#a6cee3/#8d6e63/g' iTOL_colorstrip-class1.txt
sed -i 's/#b2df8a/#b2ebf2/g' iTOL_colorstrip-class1.txt
sed -i 's/#e31a1c/#eeeeee/g' iTOL_colorstrip-class1.txt

sed -i 's/#a6cee3/#fafafa/g' iTOL_colorstrip-class2.txt
sed -i 's/#33a02c/#304ffe/g' iTOL_colorstrip-class2.txt
sed -i 's/#b2df8a/#bbdefb/g' iTOL_colorstrip-class2.txt
sed -i 's/#1f78b4/#ef9a9a/g' iTOL_colorstrip-class2.txt


#window安装
setwd("D:/braB_Z_HGT/biosample")
#table2itol使用需要依赖很多包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN" 
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题 
package_list = c("grid","ggplot2","gridExtra","vegan","reshape2","readODS") 
# 判断R包加载是否成功来决定是否安装后再加载 
for(p in package_list){  
    if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){   
        install.packages(p, repos=site)   
        suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))   
    } 
}  
#### 安装缺少的R包 
source("http://bioconductor.org/biocLite.R") 
biocLite(c("optparse", "plotrix", "readODS", "readxl", "yaml"))
#输入文件准备

source("table2itol.R") 
create_itol_files(infiles ="2.tsv",
   identifier = "tip", label = "class", na.strings = "X")
```


# 5.铜绿基因树的注释信息
```R
#参考菌株的branched蛋白树
setwd("D:/")
library(dplyr)
library(ggtree)
library(ape)
library(tidytree)
library(treeio)
tree=read.newick("branched-chain.Pseudom_aeru.reroot.newick")
data<-fortify(tree)
#查看PAO1的编号
nodes <- grep("PAO1", tree$tip.label)
node1<-grep("GCF_000017205_1_WP_012076059", tree$tip.label) #374
node2<-grep("GCF_012276675_1_WP_024917320", tree$tip.label) #383
nodes<-c(node1,node2)
tree <- groupOTU(tree, nodes)
#查看以上四个节点在树的哪个位置
ggtree(tree, aes(colour = group)) + 
  scale_color_manual(values=c("black", "red")) +
  theme(legend.position = "none")
# 最近的父节点794
clade <- MRCA(tree, nodes) 
sub_tree<- tree_subset(tree, clade, levels_back = 0)
#导出树文件
write.tree(sub_tree,file = "branched_representative_sub.nwk")
#导出geneid
data1<-fortify(sub_tree)
write.table(data1,file="classA.txt")
write.table(data,file="whole.txt")
cat classA.txt | grep "Pse" | cut -d " " -f 5  | sed 's/"//g' | perl -alne '$id=(split/\s+/,$_)[0];print"$id\tbraZ_like";' >classA.tsv
cat whole.txt | grep "Pse" | cut -d " " -f 5  | sed 's/"//g' | grep -v -f <(cat classA.txt | grep "Pse" | cut -d " " -f 5  | sed 's/"//g' ) | perl -alne '$id=(split/\s+/,$_)[0];print"$id\tbraB_like";' >classB.tsv
cat classA.tsv classB.tsv >pse_aeru_structure_class.tsv
cat  pse_aeru_structure_class.tsv  |  sed '1itip\tclass1' >pse_aeru_structure_class1.tsv
cat  pse_aeru_structure_class.tsv  |  sed '1itip\tclass2' | sed 's/braB_//g' | sed 's/braZ_//g' >pse_aeru_structure_class2.tsv
#输出注释文件
Rscript ./table2itol/table2itol.R  -D plan2 -i tip  -b class  -w 0.5 pse_aeru_structure_class1.tsv  pse_aeru_structure_class2.tsv

###铜绿基因树的注意事项
#1.区分PAO1-braB和PAO1-braZ  树枝    深红和深蓝 
#2.区分braZ-like和braB-like  外圈   浅红和浅蓝
#3.去掉菌株蛋白名字   树枝加粗width=6 px
#4.树的角度 190°  外圈宽度  datasets: width 90px  border style 1px  complete border yes

#修改plan2文件夹中的 iTOL_colorstrip-class1.txt文件
#class1替换为class 
#颜色#d95f02替换为#FB6A4A     颜色#1b9e77替换为#4292C6
cat iTOL_colorstrip-class1.txt | sed 's/#d95f02/#FCC5C0/g' | sed 's/#1b9e77/#C6DBEF/g' | sed 's/class1/class/g' >pse_aeru_structure_class_color.txt
cat iTOL_colorstrip-class2.txt | sed 's/#d95f02/#ffffff/g' | sed 's/#1b9e77/#7A0177/g' | sed 's/class1/class/g' >pse_aeru_structure_class_color1.txt

#颜色选择
 library(RColorBrewer)
display.brewer.all() #展示所有的色板
display.brewer.pal(11，'RdBu')#展示色板中的颜色
brewer.pal(11,"RdBu")  #查看颜色码
#深红 #CB181D 
#深蓝 #08519C  
#浅红 #FCC5C0   braB 
#浅蓝 #C6DBEF  braZ
```

# 6.共线性的注意事项
```bash
#pse_aeru
python  fetch_gbk.py -i Pseudom_aeru_PAO1.gbff -u 4000 -d 6000 -l braB -o PAO1_braB
python  fetch_gbk.py -i Pseudom_aeru_PAO1.gbff -u 7000 -d 4000 -l braZ -o PAO1_braZ
#将共线性的基因调整到10个以内，铜绿假单胞菌其他菌株都是上下游各5000bp
#共线性的页面结果调整scale factor为25 vertical spacing为30

#two_braZ
python  fetch_gbk.py -i Pseudom_aeru_PAO1.gbff -u 7000 -d 7000 -l braB -o PAO1_braB
python  fetch_gbk.py -i Pseudom_aeru_PAO1.gbff -u 7000 -d 7000 -l braZ -o PAO1_braZ

for filename in *.gbff 
do
base=$(basename $filename .gbff)
python fetch_gbk.py -i $base.gbff -u 7000 -d 7000 -l brnQ -o braB
done


```

# 7.变形菌纲基因树和物种树的注释信息
```r
#提取变形菌纲的科分类信息，方便后续分类,不加外类群448
#共有18个对应的科
cat branched-chain/tree/branched-chain.Gammaproteobacteria.protein.tsv | grep -v "Bac_subti_subtilis" | grep -v "Sta_aure_aureus" | tsv-join -d 1 -f PROTEINS/all.strain.tsv  -k 1 --append-fields 2 | tsv-join -d 2 -f strains.taxon.tsv -k 1 --append-fields 4 | tsv-select -f 3,1,2 | nwr append stdin -r family | tsv-summarize -g 4 --count | wc -l
#基因树对应的注释信息448
cat branched-chain/tree/branched-chain.Gammaproteobacteria.protein.tsv | grep -v "Bac_subti_subtilis" | grep -v "Sta_aure_aureus" | tsv-join -d 1 -f PROTEINS/all.strain.tsv  -k 1 --append-fields 2 | tsv-join -d 2 -f strains.taxon.tsv -k 1 --append-fields 4 | tsv-select -f 3,1,2 | nwr append stdin -r family | cut -f 2,4 >branched-chain.Gammaproteobacteria.protein.anno.tsv
#物种树对应的注释信息374
cat branched-chain/tree/branched-chain.Gammaproteobacteria.protein.tsv | grep -v "Bac_subti_subtilis" | grep -v "Sta_aure_aureus" | tsv-join -d 1 -f PROTEINS/all.strain.tsv  -k 1 --append-fields 2 | tsv-join -d 2 -f strains.taxon.tsv -k 1 --append-fields 4 | tsv-select -f 3,1,2 | nwr append stdin -r family | cut -f 3,4 |sort -n | uniq  >branched-chain.Gammaproteobacteria.bac120.anno.tsv
sed -i '1ispe\tfamily' branched-chain.Gammaproteobacteria.bac120.anno.tsv
#提取注释信息(COLORSTRIP)
Rscript ./table2itol/table2itol.R -D plan3 -i spe  branched-chain.Gammaproteobacteria.bac120.anno.tsv
#提取注释信息(TREE_COLORS，包括range)
Rscript ./table2itol/table2itol.R -D plan3 -i spe -b family -w 0.5  branched-chain.Gammaproteobacteria.bac120.anno.tsv
cat branched-chain.Gammaproteobacteria.bac120.anno.tsv | cut -f 2 | sort -n | uniq | perl -e '$num=0;while(<>){chomp;$name=(split/\t/,$_)[0];$num+=1;print"$name\t$num\n";}' >family_num.tsv
Rscript ./table2itol/table2itol.R -D plan3 -i spe -b family -w 0   branched-chain.Gammaproteobacteria.protein.anno.tsv
```

