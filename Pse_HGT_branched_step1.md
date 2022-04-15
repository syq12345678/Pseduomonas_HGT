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
- [3. 两种树](#3-两种树)
  - [3.1 参考细菌的bac120蛋白树](#31-参考细菌的bac120蛋白树)
  - [3.2 参考细菌的branched蛋白树](#32-参考细菌的branched蛋白树)
  - [3.3 代表细菌的bac120蛋白树](#33-代表细菌的bac120蛋白树)
  - [3.4 代表细菌的branched蛋白树](#34-代表细菌的branched蛋白树)

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
## 2.4 多轮diamond（将hmmer匹配的蛋白作为query序列与所有蛋白进行比对，再次进行筛选）all.replace.fa中有两条序列存在大量空字符(null)导致无法建库，需要删除
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
#统计不同species的所有assembly
cat strains.taxon.tsv | tsv-summarize -g 4 --count | keep-header -- tsv-sort -k2,2n >branched-chain/diamond/assembly_species_num.tsv
cat branched-chain/diamond/assembly_species_num.tsv | wc -l  #572
#统计不同species的含有copy的assembly
cat branched-chain/diamond/branched_chain_hmmer_speices_num.tsv | wc -l  #386
##根据不同物种的菌株数量和含有copy的菌株数量差异，发现有的菌株存在基因丢失情况
##因此基因存在丢失的菌株需要将拷贝记为0
#差异数目为186
cat branched-chain/diamond/assembly_species_num.tsv | grep -v -f <(cut -f 1 branched-chain/diamond/branched_chain_hmmer_species_num.tsv) | tr "\t" "," | 
perl -F, -alne '$name=$F[0];$num=0;print"$name\t$num";'  >branched-chain/diamond/branched_chain_hmmer_species_num_addinformation.tsv
cat branched-chain/diamond/branched_chain_hmmer_species_num.tsv  branched-chain/diamond/branched_chain_hmmer_species_num_addinformation.tsv >branched-chain/diamond/branched_chain_hmmer_species_copy_num.tsv
```

## 2.7 绘制统计不同物种中BCAA平均拷贝数的表格
* 表格内容
```bash
#species  | number of assemblies  | numbers of BCAA | average per genome 
#合并不同菌株的组装体数目和不同菌株中的BCAA拷贝数并且相除计算copy per genome
cat branched-chain/diamond/branched_chain_hmmer_species_copy_num.tsv | 
tsv-join -d 1 \
-f branched-chain/diamond/assembly_species_num.tsv -k 1 \
--append-fields 2 |
tsv-select -f 1,3,2  | tr "\t" "," |
perl -F, -alne '$per=$F[2]/$F[1];$per=sprintf "%.1f",$per;print"$F[0]\t$F[1]\t$F[2]\t$per";' >branched-chain/diamond/branched_chain_species_assembly_copy.tsv
#添加属，科,目，纲和group species信息
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
#查看假单胞菌属的非typical菌株的拷贝数目(57)
cut -f 1,2,3 branched-chain/diamond/branched_chain_species_mean_copy.tsv | grep "Pseudomonas" | grep -v -f <(cut -f 1 branched-chain/diamond/branched_chain_copy_format1.tsv) | nwr append stdin -r class | tsv-summarize -g 4 --sum 2,3 | sed 's/Gammaproteobacteria/Other Pseudomonas/g'  | tr "\t" "," | 
perl -F, -alne '$per=$F[2]/$F[1];$per=sprintf "%.1f",$per;print"$F[0]\t$F[1]\t$F[2]\t$per";'  >branched-chain/diamond/branched_chain_copy_format2.tsv
#提取非假单胞菌属的其他目的结果(506/499) 去除拷贝为0的目
cut -f 1,2,3,4,5,6 branched-chain/diamond/branched_chain_species_mean_copy.tsv |  grep -v "Pseudomonas" | grep -v "order" | grep "Gammaproteobacteria" | tsv-select -f 6,2,3 |  
tsv-summarize -g 1 --sum 2,3 | tr "\t" "," | perl -F, -alne '$per=$F[2]/$F[1];$per=sprintf "%.1f",$per;print"$F[0]\t$F[1]\t$F[2]\t$per";' |  
tsv-sort -k1,1rn  >branched-chain/diamond/branched_chain_copy_format3.tsv
#合并三者的结果(由于只关注gamma变形菌纲，去除厚壁菌门的芽孢杆菌目和衣原体目)
cat branched-chain/diamond/branched_chain_copy_format1.tsv  branched-chain/diamond/branched_chain_copy_format2.tsv branched-chain/diamond/branched_chain_copy_format3.tsv |  sed '1iTaxonomy\tNumber of assemblies\tNumber of braB\braZ\tAverage per genome'  >branched-chain/diamond/branched_chain_copy_whole.tsv
mlr --itsv --ocsv cat branched-chain/diamond/branched_chain_copy_whole.tsv >branched-chain/diamond/branched_chain_copy_whole.csv
```
  


# 3. 两种树
* [系统发育和系统发育生态学标记PhyEco marker](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0077033)
* [系统识别细菌和古细菌不同分类水平的系统发育标记](https://figshare.com/articles/dataset/Systematically_identify_phylogenetic_markers_at_different_taxonomic_levels_for_bacteria_and_archaea/722713/1)
* 在不同分类水平上自动识别系统发育标记的协议，该协议使用快速搜索和聚类算法为选择的基因组生成蛋白质家族,然后建立系统发育树，并自动对来自树的进化枝进行采样和评估，以评估使它们对系统发育分析(例如感兴趣的基因组的普遍性）和生态学研究（例如复制均匀性），然后使用多重比较和系统发育分析进一步评估潜在的标记家族.使用该种方法，作者已经为所有细菌确定了114个系统发育标记,其中包括40个同时涵盖细菌和古细菌的标记
* [举例：Izimaplasma的系统发育分析](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5113845/)
* 基因组树构建:基因组树是由38个蛋白质编码基因的串联比对构建的，这些细菌在古细菌和细菌中普遍分布并且以单拷贝形式存在.该树的基础包括来自NCBI和IMG数据库的草图和完成的基因组，使用hmmer在每个基因组中鉴定标记基因.使用默认设置的FastTree创建了最大似然树.
* 基因树构建:使用来自Tigrfams的HMM模型和IMG基因组下载的基因组中鉴定RpoB或EF-tu基因序列，使用muscle比对，使用FastTree构建系统发育树


## 3.1 参考细菌的bac120蛋白树
```bash
cd ~/data/Pseudomonas/
mkdir -p ~/data/Pseudomonas/branched-chain/tree
#提取参考菌株的bac120名字(去除非gamma变形菌纲的)
cat strains.taxon.tsv | grep -f reference.tsv | grep -v "GCF" | cut -f 1 > PROTEINS/bac120_refer.lst 
#提取序列
faops some PROTEINS/bac120.trim.fa PROTEINS/bac120_refer.lst   PROTEINS/bac120.refer.fa
muscle -in PROTEINS/bac120.refer.fa  -out PROTEINS/bac120.refer.aln.fa
#建树(外类群选择gamma变形菌纲的)
FastTree PROTEINS/bac120.refer.aln.fa > PROTEINS/bac120.refer.aln.newick
nw_reroot  PROTEINS/bac120.refer.aln.newick $(nw_labels PROTEINS/bac120.refer.aln.newick | grep -E "Bac_subti|Sta_aure" ) |
    nw_order -c n - \
    > PROTEINS/bac120.refer.reroot.newick
 ```

## 3.2 参考细菌的branched蛋白树
```bash
#提取参考菌株的braz名字(去除非gamma变形菌纲的)
cat branched-chain/branched-chain_minevalue.tsv | grep -f reference.tsv | grep -v 'GCF' | cut -f 1 >branched-chain/tree/branched-chain_refer.tsv
#提取序列
faops some PROTEINS/all.replace.fa branched-chain/tree/branched-chain_refer.tsv  branched-chain/tree/branched-chain_refer.fa
muscle -in branched-chain/tree/branched-chain_refer.fa -out branched-chain/tree/branched-chain_refer.aln.fa
#建树(外类群选择gamma变形菌纲的)
FastTree branched-chain/tree/branched-chain_refer.aln.fa > branched-chain/tree/branched-chain_refer.aln.newick
nw_reroot branched-chain/tree/branched-chain_refer.aln.newick $(nw_labels branched-chain/branched-chain_refer.aln.newick | grep -E "POA1") |nw_order -c n - > branched-chain/tree/branched-chain_refer.reroot.newick
```

```R
## 参考细菌的bac120蛋白树和参考生物的branched蛋白树
setwd("D:/")
library(dplyr)
library(ggtree)
library(ape)
tree1=read.tree("bac120.refer.reroot.newick")
tree2=read.tree("branched-chain.refer.reroot.newick")
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

## 3.3 代表细菌的bac120蛋白树
```bash
cd ~/data/Pseudomonas
#提取参考菌株的bac120名字(去除非gamma变形菌纲的)
cat strains.taxon.tsv | grep -f representative.tsv | grep -v "GCF" | cut -f 1 > PROTEINS/bac120_representative.lst 
#提取序列
faops some PROTEINS/bac120.trim.fa PROTEINS/bac120_representative.lst    PROTEINS/bac120.representative.fa
muscle -in PROTEINS/bac120.representative.fa  -out PROTEINS/bac120.representative.aln.fa
#建树(外类群选择gamma变形菌纲的)
FastTree PROTEINS/bac120.representative.aln.fa > PROTEINS/bac120.representative.aln.newick
nw_reroot  PROTEINS/bac120.refer.representative.newick $(nw_labels PROTEINS/bac120.representative.aln.newick | grep -E "Bac_subti|Sta_aure" ) |
nw_order -c n -  > PROTEINS/bac120.representative.reroot.newick
```

## 3.4 代表细菌的branched蛋白树
```bash
cd ~/data/Pseudomonas
#提取参考菌株的braz名字(去除非gamma变形菌纲的)
cat branched-chain/branched-chain_minevalue.tsv | grep -f representative.tsv | grep -v 'GCF' | cut -f 1 >branched-chain/tree/branched-chain_representative.tsv
#提取序列
faops some PROTEINS/all.replace.fa branched-chain/tree/branched-chain_representative.tsv  branched-chain/tree/branched-chain_representative.fa
muscle -in branched-chain/tree/branched-chain_representative.fa -out branched-chain/tree/branched-chain_representative.aln.fa
#建树(外类群选择gamma变形菌纲的)
FastTree branched-chain/tree/branched-chain_representative.aln.fa > branched-chain/tree/branched-chain_representative.aln.newick
nw_reroot branched-chain/tree/branched-chain_representative.aln.newick $(nw_labels branched-chain/branched-chain_representative.aln.newick | grep -E "POA1") |nw_order -c n - > branched-chain/tree/branched-chain_representative.reroot.newick
```
```R
#参考菌株的branched蛋白树
setwd("D:/")
library(dplyr)
library(ggtree)
library(ape)
library(tidytree)
library(treeio)
tree=read.newick("branched-chain_representative.reroot.newick")
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
write.tree(sub_tree,file = "branched_representative_sub.nwk")
```
