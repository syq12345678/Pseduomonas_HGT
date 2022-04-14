<!-- TOC -->

- [1.codes from Professor Wang](#1codes-from-professor-wang)
  - [1.1外类群：使用一些模型生物作为外类群，具体信息在reference.tsv中，共有15种](#11外类群使用一些模型生物作为外类群具体信息在referencetsv中共有15种)
  - [1.2菌株基因组蛋白信息文件](#12菌株基因组蛋白信息文件)
- [2.统计branched-chain在菌株中拷贝分布](#2统计branched-chain在菌株中拷贝分布)
  - [2.1使用hmmsearch抓取相同domain的基因或者基因家族](#21使用hmmsearch抓取相同domain的基因或者基因家族)
  - [2.2将branched-chain提取的蛋白序列与tigerfam数据库比对](#22将branched-chain提取的蛋白序列与tigerfam数据库比对)
  - [2.3将branched-chain提取的蛋白序列与pfam数据库比对](#23将branched-chain提取的蛋白序列与pfam数据库比对)
  - [2.4多轮diamond（将hmmer匹配的蛋白作为query序列与所有蛋白进行比对，再次进行筛选）！！！！！all.replace.fa中有两条序列存在大量空字符(null)导致无法建库，需要删除](#24多轮diamond将hmmer匹配的蛋白作为query序列与所有蛋白进行比对再次进行筛选allreplacefa中有两条序列存在大量空字符null导致无法建库需要删除)
  - [2.5统计diamond比对和hmmer比对时不同species中braB个数](#25统计diamond比对和hmmer比对时不同species中brab个数)
  - [2.5绘制统计不同物种中braB和braz拷贝数的表格](#25绘制统计不同物种中brab和braz拷贝数的表格)
- [3.两种树](#3两种树)
  - [3.1模式细菌的bac120蛋白树](#31模式细菌的bac120蛋白树)
  - [3.2模式生物的branched蛋白树](#32模式生物的branched蛋白树)
  - [3.3参考菌株的branced蛋白树](#33参考菌株的branced蛋白树)
- [4.计算ka/ks](#4计算kaks)
- [以下方法错误!!!](#以下方法错误)
- [以下方法错误!!!](#以下方法错误-1)
- [以下方法错误!!!](#以下方法错误-2)
  - [4.1生成两个copy列表](#41生成两个copy列表)
  - [4.2使用paraAT2和kakscalculator2计算kaks](#42使用paraat2和kakscalculator2计算kaks)

<!-- /TOC -->
* 1.由于外类群不够，新添了一些基因组文件，由1514增加至1953个基因组文件，由此需要重新操作  
* 2.hmmsearch搜索结构域，hmmscan过滤,diamond重复计算拷贝数，绘制拷贝数的表格，绘制Reference Genome的树和Representative Genome与Reference Genome的树  
* 3.目前以branched为例,branched在pfam panthern和tigrfam中均存在    
# 1.codes from Professor Wang
## 1.1外类群：使用一些模型生物作为外类群，具体信息在reference.tsv中，共有15种
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

## 1.2菌株基因组蛋白信息文件  
```shell
#包含NCBI的样本信息且去除了错误的菌株基因组文件信息，共1952个  
ASSEMBLY/Pseudomonas.assembly.pass.csv  
#assembly level:Complete Genome(1810)或Chromosome(142)
#共有15个Reference Genome参考菌株，分别属于变形菌纲，厚壁菌纲，放线菌纲和衣原体
#共有533个Representative Genome代表菌株都属于Gamma变形菌纲
cat ASSEMBLY/Pseudomonas.assembly.pass.csv | grep "Representative Genome" | cut -d , -f 1 | tsv-join -d 1 -f strains.taxon.tsv -k 1  --append-fields 4 |  
tsv-select -f 2,1 | nwr append stdin -r class -r order -r family -r genus >representative_add_rank.tsv

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
# 2.统计branched-chain在菌株中拷贝分布
## 2.1使用hmmsearch抓取相同domain的基因或者基因家族
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
## 2.2将branched-chain提取的蛋白序列与tigerfam数据库比对
```bash
#下载tigerfam数据库
mkdir -p ~/data/HMM/TIGERFAM
cd ~/data/HMM/TIGERFAM
wget -N --content-disposition https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz
tar -xzvf hmm_PGAP.HMM.tgz
cat *.HMM >tigrfams.hmm

#格式化tigerfam数据库
cd ~/data/Pseudomonas
hmmpress ~/data/HMM/TIGERFAM/tigerfam.hmm

#提取抓出来的序列
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 branched-chain/branched.tigerfam.replace.tsv) branched-chain/branched-chain.fa

E_VALUE=1e-10
NAME=branched-chain
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.txt --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/TIGERFAM/tigrfams.hmm  ${NAME}/${NAME}.fa

perl abstract1.pl branched-chain/branched-chain.tbl >branched-chain/branched-chain.abstract.tsv
tsv-filter --le 4:1e-50 branched-chain/branched-chain.abstract.tsv >branched-chain/branched-chain.cutoff.tsv
tsv-summarize  -g 8 --count branched-chain/branched-chain.cutoff.tsv
#查看以e值过滤后的结果
livcs:_branched-chain_amino_acid_transport_system_II_carrier_protein    2140

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare1.pl branched-chain/branched-chain.cutoff.tsv > branched-chain/branched-chain_minevalue.tsv

#拼接属名等信息并统计拷贝数
cat branched-chain/branched-chain_minevalue.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"branched-chain" |
tsv-join -d 1 \
-f PROTEINS/all.strain.tsv -k 1 \
--append-fields 2 |
tsv-join -d 3 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-summarize -g 3,4 --count |
keep-header -- tsv-sort -k3,3n >branched-chain/branched-chain_hmmscan_copy.tsv

#统计拷贝数的分布
tsv-summarize -g 3,2 --count  branched-chain/branched-chain_hmmscan_copy.tsv > branched-chain/branched-chain_hmmscan_GCF_copy.tsv
sed -i '1icopy\tgenus\tGCF'  branched-chain/branched-chain_hmmscan_GCF_copy.tsv

```
## 2.3将branched-chain提取的蛋白序列与pfam数据库比对
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

#提取抓出来的序列
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 branched-chain/branched.pfam.replace.tsv) branched-chain/branched-chain.pfam.fa

E_VALUE=1e-10
NAME=branched-chain
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.pfam.txt --tblout ${NAME}/${NAME}.pfam.tbl  \
   ~/data/HMM/TIGERFAM/tigrfams.hmm  ${NAME}/${NAME}.pfam.fa

perl abstract1.pl branched-chain/branched-chain.pfam.tbl >branched-chain/branched-chain.abstract.pfam.tsv
tsv-filter --le 4:1e-50 branched-chain/branched-chain.abstract.pfam.tsv >branched-chain/branched-chain.cutoff.pfam.tsv
tsv-summarize  -g 8 --count branched-chain/branched-chain.cutoff.pfam.tsv
#查看以e值过滤后的结果
livcs:_branched-chain_amino_acid_transport_system_II_carrier_protein    2140

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare1.pl branched-chain/branched-chain.cutoff.pfam.tsv > branched-chain/branched-chain_minevalue.pfam.tsv

#拼接属名等信息并统计拷贝数
cat branched-chain/branched-chain_minevalue.pfam.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"branched-chain" |
tsv-join -d 1 \
-f PROTEINS/all.strain.tsv -k 1 \
--append-fields 2 |
tsv-join -d 3 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-summarize -g 3,4 --count |
keep-header -- tsv-sort -k3,3n >branched-chain/branched-chain_hmmscan_copy.pfam.tsv

#统计拷贝数的分布
tsv-summarize -g 3,2 --count  branched-chain/branched-chain_hmmscan_copy.pfam.tsv > branched-chain/branched-chain_hmmscan_GCF_copy.pfam.tsv
sed -i '1icopy\tgenus\tGCF'  branched-chain/branched-chain_hmmscan_GCF_copy.pfam.tsv
```
## 2.4多轮diamond（将hmmer匹配的蛋白作为query序列与所有蛋白进行比对，再次进行筛选）！！！！！all.replace.fa中有两条序列存在大量空字符(null)导致无法建库，需要删除
```bash
cd ~/data/Pseudomonas
mkdir -p branched-chain/diamond
#删除all.replacd.fa的空字符
sed -i 's/\x0//g' PROTEINS/all.replace.fa
#提取hmmsearch和hmmscan结果
cat branched-chain/branched-chain_minevalue.pfam.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"branched-chain" | cut -f 1 >branched-chain/diamond/branched-chain_diamond1.tsv
#第一轮diamond(任取一条菌株里的braB序列进行blast)
faops some PROTEINS/all.replace.fa  branched-chain/diamond/branched-chain_diamond1.tsv branched-chain/diamond/branched-chain_diamond1.fa
diamond makedb --in branched-chain/diamond/branched-chain_diamond1.fa --db branched-chain/diamond/branched-chain_diamond1
diamond blastp --db branched-chain/diamond/branched-chain_diamond1.dmnd --query  PROTEINS/all.replace.fa -e 1e-5 --outfmt 6 --threads 4 --out branched-chain/diamond/branched-chain_result1.tsv --more-sensitive

#第二轮diamond(将第一轮从蛋白数据库中抓取出的序列进行blast)
faops some PROTEINS/all.replace.fa <(cat branched-chain/diamond/branched-chain_result1.tsv | cut -f 2 | sort -n | uniq)   branched-chain/diamond/branched-chain_diamond2.fa 
diamond makedb --in branched-chain/diamond/branched-chain_diamond2.fa --db branched-chain/diamond/branched-chain_diamond2
diamond blastp --db branched-chain/diamond/branched-chain_diamond2.dmnd --query  PROTEINS/all.replace.fa -e 1e-5 --outfmt 6 --threads 4 --out branched-chain/diamond/branched-chain_result2.tsv  --more-sensitive 

#第三轮diamond
faops some PROTEINS/all.replace.fa <(cat branched-chain/diamond/branched-chain_result2.tsv | cut -f 2 | sort -n | uniq)   branched-chain/diamond/branched-chain_diamond3.fa 
diamond makedb --in branched-chain/diamond/branched-chain_diamond3.fa --db branched-chain/diamond/branched-chain_diamond3
diamond blastp --db branched-chain/diamond/branched-chain_diamond3.dmnd --query  PROTEINS/all.replace.fa -e 1e-5 --outfmt 6 --threads 4 --out branched-chain/diamond/branched-chain_result3.tsv  --more-sensitive 



#hmmer结果
cat branched-chain/diamond/branched-chain_diamond1.tsv | wc -l  #2140
#第一轮diamond的query
cut -f 1 branched-chain/diamond/branched-chain_result1.tsv | sort -n | uniq | wc -l #2144
#第一轮diamond的target
cut -f 2 branched-chain/diamond/branched-chain_result1.tsv | sort -n | uniq | wc -l #1492

#第二轮diamond的query
cut -f 1 branched-chain/diamond/branched-chain_result2.tsv | sort -n | uniq | wc -l #2144
#第二轮diamond的target
cut -f 2 branched-chain/diamond/branched-chain_result2.tsv | sort -n | uniq | wc -l #1492

#第三轮diamond的query
cut -f 1 branched-chain/diamond/branched-chain_result3.tsv | sort -n | uniq | wc -l #2144
#第三轮diamond的target
cut -f 2 branched-chain/diamond/branched-chain_result3.tsv | sort -n | uniq | wc -l #1492

#三轮diamond结果一致
```


## 2.5统计diamond比对和hmmer比对时不同species中braB个数
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
keep-header -- tsv-sort -k2,2n >branched-chain/diamond/branched_chain_diamond_genus_num.tsv

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
keep-header -- tsv-sort -k2,2n >branched-chain/diamond/branched_chain_hmmer_genus_num.tsv
```

## 2.5绘制统计不同物种中braB和braz拷贝数的表格
```bash
#最终的表格形式如下
#species  | number of assemblies  | numbers of mltB | average per genome 

#统计assembly总个数
cat strains.taxon.tsv | cut -f 1 | sort -n | uniq | wc -l #1952 
#统计不同species的所有assembly
cat strains.taxon.tsv | tsv-summarize -g 4 --count | keep-header -- tsv-sort -k2,2n >branched-chain/diamond/assembly_genus_num.tsv
cat branched-chain/diamond/assembly_genus_num.tsv | wc -l  #572
#统计不同species的含有copy的assembly
cat branched-chain/diamond/branched_chain_hmmer_genus_num.tsv | wc -l  #386
##根据不同物种的菌株数量和含有copy的菌株数量差异，发现有的菌株存在基因丢失情况
##因此基因存在丢失的菌株需要将拷贝记为0
#差异数目为186
cat branched-chain/diamond/assembly_genus_num.tsv | grep -v -f <(cut -f 1 branched-chain/diamond/branched_chain_hmmer_genus_num.tsv) | tr "\t" "," | 
perl -F, -alne '$name=$F[0];$num=0;print"$name\t$num";'  >branched-chain/diamond/branched_chain_hmmer_genus_num_addinformation.tsv
cat branched-chain/diamond/branched_chain_hmmer_genus_num.tsv  branched-chain/diamond/branched_chain_hmmer_genus_num_addinformation.tsv >branched-chain/diamond/branched_chain_hmmer_species_copy_num.tsv

#两者合并并且相除计算copy per genome
cat branched-chain/diamond/branched_chain_hmmer_species_copy_num.tsv | 
tsv-join -d 1 \
-f branched-chain/diamond/assembly_genus_num.tsv -k 1 \
--append-fields 2 |
tsv-select -f 1,3,2  | tr "\t" "," |
perl -F, -alne '$per=$F[2]/$F[1];$per=sprintf "%.1f",$per;print"$F[0]\t$F[1]\t$F[2]\t$per";' >branched-chain/diamond/branched_chain_genus_assembly_copy.tsv

#添加门信息，科信息
cat branched-chain/diamond/branched_chain_genus_assembly_copy.tsv | nwr append stdin -r class -r order -r family -r genus  >branched-chain/diamond/branched_chain_class_order_family_genus.tsv

#将不同物种的assembly个数和braz总个数及平均个数进行排序统计
sed -i '1ispecies\tnumber of assemblies\tnumber of braZ\taverage per genome\tclass\torder\t\tfamily\tgenus' branched-chain/diamond/branched_chain_class_order_family_genus.tsv
cat branched-chain/diamond/branched_chain_class_order_family_genus.tsv |  keep-header -- tsv-sort -k3,3rn -k4,4rn   >branched-chain/diamond/branched_chain_mean_copy.tsv
plotr tsv branched-chain/diamond/branched_chain_mean_copy.tsv --header

```


# 3.两种树
## 3.1模式细菌的bac120蛋白树
```bash
cat strains.taxon.tsv |
    grep -v "GCF" | 
    cut -f 1 > model.lst 
#提取序列
faops some PROTEINS/bac120.trim.fa model.lst  PROTEINS/bac120.model.fa
muscle -in PROTEINS/bac120.model.fa  -out PROTEINS/bac120.model.aln.fa
#建树
FastTree PROTEINS/bac120.model.aln.fa > PROTEINS/bac120.model.aln.newick
nw_reroot  PROTEINS/bac120.model.aln.newick $(nw_labels PROTEINS/bac120.model.aln.newick | grep -E "Bac_subti|Sta_aure" ) |
    nw_order -c n - \
    > PROTEINS/bac120.model.reroot.newick
 ```

## 3.2模式生物的branched蛋白树
```bash
cat branched-chain/branched-chain_minevalue.tsv | grep -f model.lst | grep -v 'GCF' | cut -f 1 >branched-chain/branched-chain.model.tsv
#提取序列
faops some PROTEINS/all.replace.fa branched-chain/branched-chain.model.tsv  branched-chain/branched-chain.model.fa
muscle -in branched-chain/branched-chain.model.fa -out branched-chain/branched-chain.model.aln.fa
#建树
FastTree branched-chain/branched-chain.model.aln.fa > branched-chain/branched-chain.model.aln.newick
nw_reroot branched-chain/branched-chain.model.aln.newick $(nw_labels branched-chain/branched-chain.model.aln.newick | grep -E "POA1") |
    nw_order -c n - \
    > branched-chain/branched-chain.model.reroot.newick
```

```R
## 模式生物的bac120蛋白树和模式生物的branched蛋白树
setwd("D:/")
library(dplyr)
library(ggtree)
library(ape)
tree1=read.tree("bac120.model.reroot.newick")
tree2=read.tree("branched-chain.model.reroot.newick")
p1 <- ggtree(tree1)
p2 <- ggtree(tree2)
d1 <- p1$data
d2 <- p2$data
d2$x <- max(d2$x) - d2$x + max(d1$x)+1 #翻转第二棵树
d2$y <- d2$y+1 #将第二棵树向上移动对齐
p3 <- ggtree::rotate(p1,29) #旋转node，使两棵树拓扑结构一致
p1+geom_tiplab(offset=0.05,size=3)+geom_text2(aes(subset=!isTip, label=node), hjust=-.3) #注释每个节点的编号
p3 + geom_tiplab(offset=0.05,size=3) + geom_treescale()+ geom_highlight(node=25,fill="red")+ geom_tree(data=d2) + geom_tiplab(data = d2, hjust=1, offset =-0.05,size=3)
```

## 3.3参考菌株的branced蛋白树
```bash
#共有533个代表菌株，15个模式菌株
cut -d , -f 18 ASSEMBLY/Pseudomonas.assembly.pass.csv | tsv-summarize -g 1 --count
RefSeq_category 1
Representative Genome   533
        1404
Reference Genome        15
#筛选保留了genome的行
cut -d , -f 1,18 ASSEMBLY/Pseudomonas.assembly.pass.csv | tsv-filter --d , --str-in-fld 2:"Genome" |
cut -d , -f 1 >representative.tsv
#提取代表菌株的菌株基因组蛋白号
cat branched-chain/branched-chain_minevalue.tsv | grep -f representative.tsv | cut -f 1 >branched-chain/branched-chain_repre.tsv
faops some PROTEINS/all.replace.fa branched-chain/branched-chain_repre.tsv branched-chain/branched-chain_repre.fa
muscle -in branched-chain/branched-chain_repre.fa  -out branched-chain/branched-chain_repre.aln.fa
#建树
FastTree branched-chain/branched-chain_repre.aln.fa > branched-chain/branched-chain_repre.aln.newick
nw_reroot branched-chain/branched-chain_repre.aln.newick $(nw_labels branched-chain/branched-chain_repre.aln.newick | grep -E "POA1") |
    nw_order -c n - \
    > branched-chain/branched-chain_repre.reroot.newick
```
```R
#参考菌株的branched蛋白树
setwd("D:/")
library(dplyr)
library(ggtree)
library(ape)
library(tidytree)
library(treeio)
tree=read.newick("branched-chain_repre.reroot.newick")
data<-fortify(tree)
#查看PAO1的编号
nodes <- grep("Pseudom_aeru_PAO1", tree$tip.label)
tree <- groupOTU(tree, nodes)
#查看以上四个节点在树的哪个位置
ggtree(tree, aes(colour = group)) + 
  scale_color_manual(values=c("black", "red")) +
  theme(legend.position = "none")
# 最近的父节点621
clade <- MRCA(tree, nodes) 
sub_tree<- tree_subset(tree, clade, levels_back = 0)
#导出树文件
write.tree(sub_tree,file = "branched_repre_sub.nwk")
```

# 4.计算ka/ks
1.在遗传学中，Ka/Ks表示的是两个蛋白编码基因的非同义替换率(Ka)和同义替换率(Ks)之间的比例。这个比例可以判断是否有选择压力作用于这个蛋白编码基因。
2.如果有两个不同物种的同一个基因的序列，比如人和小鼠的p53基因，然后把这两个基因的序列进行比对，你会发现这两段序列有差异（进化！）。  
再仔细观察，你会发现有些碱基的变化导致了编码氨基酸的变化（非同义替换），有些没有导致编码氨基酸的变化（同义替换）。  
这是由密码子的简并性造成的，因为3个碱基决定1个氨基酸，所以64种碱基组合决定20种氨基酸，会有冗余出现。  
3.一般情况下，第三个碱基变化会造成同义替换，而第一二个碱基的变化会造成非同义替换。  
4.Ka和Ks的计算公式：  
Ka=发生非同义替换的SNP数/非同义替换位点数  
Ks=发生同义替换的SNP数/同义替换位点数
其中非同义替换位点数就是会造成氨基酸变化的位点数的总和，比如编码丝氨酸（ser）的第一二位碱基。  
而同义替换位点数就是不会造成氨基酸变化的位点数的总和，比如编码丝氨酸的第三位碱基。  
5.计算Ka/Ks时，不考虑start codon和stop codon  
6.转换（transition，嘌呤变嘌呤，嘧啶变嘧啶）发生的概率要高于颠换（transversion，嘌呤变嘧啶，嘧啶变嘌呤）发生的概率  
7.ka/ks与进化的关系：如果一个基因没有受到自然选择的压力，那么非同义替换率和同义替换率是相同的。而一般情况下，非同义替换会造成氨基酸变化  
改变蛋白质的构象和功能因此造成适应性的变化，从而带来自然选择的优势和劣势。而同义替换没有改变蛋白质的组成，不受自然选择的影响，因此ks反映  
过程中的背景碱基替换率。ka/ks的比值说明了这个基因受到了何种选择
8.Ka>>Ks或者Ka/Ks >> 1，基因受正选择(positive selection)  
Ka＝Ks或者Ka/Ks=1，基因中性进化(neutral evolution)  
Ka << Ks或者Ka/Ks << 1，基因受纯化选择(purify selection)  
9.实际情况下Ka/Ks << 1，因为一般非同义替换带来的都是有害的性状，只有极少数情况下会造成进化上的优势  
10.当Ka/Ks>>1时，基因受到强烈正选择，这样的基因即为近期正在快速进化的基因，对于物种的进化有着非常重要的意义
11.计算ka和ks的步骤
步骤一：假设比较的两条DNA序列之间的长度数为n,它们之间的替换数为m。计算同义(S)和非同义(N)位点的数量(S+N=n)以及同义(Sd)和非同义替换的数量  
(Sd+Nd=m)
步骤二:校正多个替换后,(Nd/N)和(Sd/S)分别代表ka和ks

# 以下方法错误!!!
# 以下方法错误!!!
# 以下方法错误!!!

## 4.1生成两个copy列表
```bash
#提取假单胞菌的两个拷贝名称
cd ~/data/Pseudomonas
cat branched-chain/branched-chain_hmmscan_copy.pfam.tsv | tsv-filter --ge 3:2 |
tsv-filter --str-in-fld 1:"Pseudom_aeru" | cut -f 1 >branched-chain/branched-chain.Pseudom_aeru.two.copy.tsv
cat branched-chain/branched-chain_minevalue.tsv | grep -f branched-chain/branched-chain.Pseudom_aeru.two.copy.tsv |
cut -f 1 >branched-chain/branched-chain.Pseudom_aeru.protein.tsv
#提取假单胞菌中branced的两个拷贝WP序列
faops some PROTEINS/all.replace.fa branched-chain/branched-chain.Pseudom_aeru.protein.tsv branched-chain/branched-chain.Pseudom_aeru.protein.fa


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


###使用cd-hit提取出branched-chain/branched-chain.Pseudom_aeru.protein.fa中差异最大的两条序列最为参考序列
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
#wp簇2的序列名称
cat branched-chain/branched-chain.Pseudom_aeru.protein.tsv | grep -v -f branched-chain/branched-chain_cluster1.tsv \
>branched-chain/branched-chain_cluster2.tsv
#wp合并簇一和簇二
cat branched-chain/branched-chain_cluster1.tsv | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 \
>branched-chain/branched-chain_cluster1_del.tsv
cat branched-chain/branched-chain_cluster2.tsv | tsv-join -d 1 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 \
>branched-chain/branched-chain_cluster2.del.tsv
cat branched-chain/branched-chain_cluster1_del.tsv | tsv-join -d 2 -f branched-chain/branched-chain_cluster2.del.tsv -k 2 \
--append-fields 1 | cut -f 1,3 >branched-chain/branched-chain_cluster.tsv


#cds簇1的序列名称
cat branched-chain/branched-chain.Pseudom_aeru.CDS.tsv | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/branched-chain_cluster1.tsv) \
>branched-chain/branched-chain_cluster1.cds.tsv
#cds簇2的序列名称
cat branched-chain/branched-chain.Pseudom_aeru.CDS.tsv | grep -f <(perl -alne 's/\_([W|N]P)/\_cds\_$1/;print"$_";' branched-chain/branched-chain_cluster2.tsv) \
>branched-chain/branched-chain_cluster2.cds.tsv
sed -i 's/Pseudom_aeru_GCF_001516005_1_cds_WP_003114263.1_10//g;'  branched-chain/branched-chain_cluster2.cds.tsv
sed -i '/^$/d' branched-chain/branched-chain_cluster2.cds.tsv
#cds合并簇一和簇二
perl -alne 'm/(.*)\_cds(.*)\.(.*)$/;$cds=$_;$wp=$1.$2; print"$cds\t$wp";' branched-chain/branched-chain_cluster1.cds.tsv |
tsv-join -d 2 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 >branched-chain/branched-chain_cluster1_del.cds.tsv
perl -alne 'm/(.*)\_cds(.*)\.(.*)$/;$cds=$_;$wp=$1.$2; print"$cds\t$wp";' branched-chain/branched-chain_cluster2.cds.tsv |
tsv-join -d 2 -f PROTEINS/all.strain.tsv -k 1 --append-fields 2 >branched-chain/branched-chain_cluster2_del.cds.tsv
cat branched-chain/branched-chain_cluster1_del.cds.tsv | tsv-join -d 3 -f branched-chain/branched-chain_cluster2_del.cds.tsv -k 3 \
--append-fields 1 | cut -f 1,4 >branched-chain/branched-chain_cluster.cds.tsv


#准备列表文件
branched-chain/branched-chain_cluster.tsv
#准备protein序列
faops order branched-chain/branched-chain.Pseudom_aeru.protein.fa \
<(perl -alne 'print"$F[0]\n$F[1]" ' branched-chain/branched-chain_cluster.tsv) branched-chain/branched-chain.Pseudom_aeru.PRO.fa
#准备cds序列
faops some branched-chain/branched-chain.Pseudom_aeru.CDS.fa \
<(cat branched-chain/branched-chain_cluster1_del.tsv  branched-chain/branched-chain_cluster2_del.tsv | cut -f 1) \
branched-chain/branched-chain.Pseudom_aeru.CDS.delete.fa
#将cds的名称替换为wp的名称
cat branched-chain/branched-chain_cluster1_del.tsv  branched-chain/branched-chain_cluster2_del.tsv | tsv-select -f 1,2 |
    perl -nla -e '
        print q{s/^>} . quotemeta($F[0]) . q{/>} . quotemeta($F[1]) . q{/g;};
    ' \
    > branched-chain/sed.script
cat branched-chain/branched-chain.Pseudom_aeru.CDS.delete.fa | sed -f branched-chain/sed.script \
>branched-chain/branched-chain.Pseudom_aeru.CDS.replace.fa
faops order branched-chain/branched-chain.Pseudom_aeru.CDS.replace.fa \
<(perl -alne 'print"$F[0]\n$F[1]" ' branched-chain/branched-chain_cluster.tsv) branched-chain/branched-chain.Pseudom_aeru.CDS.WP.fa

```
```
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


## 4.2使用paraAT2和kakscalculator2计算kaks
```bash
#安装paraat2
cd ~
wget ftp://download.big.ac.cn/bigd/tools/ParaAT2.0.tar.gz
tar -xzvf ParaAT2.0.tar.gz
vim ~/.bashrc
export PATH="$PATH:/home/syq/ParaAT2.0"
source ~/.bashrc
#安装kakscalculator2
conda install kakscalculator2 
#使用paraAT2计算kaks
ParaAT.pl -h branched-chain/branched-chain_cluster.tsv -n branched-chain/branched-chain.Pseudom_aeru.CDS.WP.fa \
-a branched-chain/branched-chain.Pseudom_aeru.PRO.fa -p proc -m muscle -f axt -k -o branched-chain/branched-chain_para

-h, 同源基因名称文件
-n, 指定核酸序列文件
-a, 指定蛋白序列文件
-p, 指定多线程文件
-m, 指定比对工具
-g, 去除比对有gap的密码子
-k, 用KaKs_Calculator 计算kaks值
-o, 输出结果的目录
-f, 输出比对文件的格式

#上述结果可直接得到每一对同源基因的ka，ks值，可通过如下命令将其整合
cat  branched-chain/branched-chain_para/*.kaks | cut -f 1,3,4,5 |grep -v 'Sequence' | head

```