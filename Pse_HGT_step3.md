<!-- TOC -->

- [1.以MltB基因为例，将hmmer匹配的蛋白作为query序列与所有的蛋白进行比对，再次进行筛选](#1以mltb基因为例将hmmer匹配的蛋白作为query序列与所有的蛋白进行比对再次进行筛选)
- [2.统计不同species中的copy结果](#2统计不同species中的copy结果)
- [3.两种树(模式生物的)](#3两种树模式生物的)
- [4.提取所有的DNA序列](#4提取所有的dna序列)

<!-- /TOC -->
* 1.使用hmmsearch最大范围内统计不同菌株中蛋白的copy数    
* 2.使用hmmscan保留最显著的domain为给定famliy的蛋白，并统计不同菌株中蛋白的copy数据,e值为1e-50    
* 3.使用blastp将种子序列在所有蛋白里面检索直至没有新的结果出现

# 1.以MltB基因为例，将hmmer匹配的蛋白作为query序列与所有的蛋白进行比对，再次进行筛选
```BASH
#diamond比对结果
1. qseqid |  query  sequence id 
2. sseqid |  subject sequence id 
3. pident |  percentage of identical matches 
4. length |  alignment length 
5. mismatch |  number of mismatches 
6. gapopen |  number of gap openings 
7. q start |  start of alignment in query 
8. q end |  end of alignment in query 
9. s start |  start of alignment in subject 
10. s end |  end of alignment in subject 
11. evalue
12. bitscore

cd ~/data/Pseudomonas
#提取hmmsearch和hmmscan结果
cat MltB/MltB_minevalue.tsv |  tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"MltB:_lytic_murein_transglycosylase_B" |
cut -f 1 >MltB/MltB_diamond1.tsv

#第一轮diamond
faops some PROTEINS/all.replace.fa  MltB/MltB_diamond1.tsv MltB/MltB_diamond1.fa
diamond makedb --in MltB/MltB_diamond1.fa --db MltB/MltB1
diamond blastp --db MltB/MltB1.dmnd --query PROTEINS/all.replace.fa -e 1e-10 --outfmt 6 --threads 8 --out MltB/MltB_result1.tsv

#第二轮diamond
faops some PROTEINS/all.replace.fa <(cat MltB/MltB_result1.tsv | cut -f 2 | sort -n | uniq) MltB/MltB_diamond2.fa
diamond makedb --in MltB/MltB_diamond2.fa --db MltB/MltB2
diamond blastp --db MltB/MltB2.dmnd --query PROTEINS/all.replace.fa -e 1e-10 --outfmt 6 --threads 8 --out MltB/MltB_result2.tsv


#hmmer结果
cat MltB/MltB_diamond1.tsv | wc -l  #1790

#第一轮diamond的query
cut -f 1 MltB/MltB_result1.tsv | sort -n | uniq | wc -l #3672
#第一轮diamond的target
cut -f 2 MltB/MltB_result1.tsv | sort -n | uniq | wc -l #1047

#第二轮diamond的query
cut -f 1 MltB/MltB_result2.tsv | sort -n | uniq | wc -l #3762
#第二轮diamond的target
cut -f 2 MltB/MltB_result2.tsv | sort -n | uniq | wc -l #1047

两轮结果一致
```
# 2.统计不同species中的copy结果
```BASH
#species  | number of assemblies  | numbers of mltB | average per genome 

#统计diamond不同species中MltB个数（percentage=90)
cat MltB/MltB_result2.tsv | tsv-filter --eq 3:100 |
cut  -f 1 |
sort -n | uniq |
tsv-join -d 1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 2 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 |
tsv-summarize -g 3 --count |
keep-header -- tsv-sort -k2,2n >MltB/MltB_diamond_genus_num.tsv


#统计hmmer不同species中MltB个数
cat MltB/MltB_minevalue.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"MltB:_lytic_murein_transglycosylase_B" |
tsv-join -d 1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 3 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-summarize -g 4 --count | keep-header -- tsv-sort -k2,2n >MltB/MltB_hmmer_genus_num.tsv

###hmmer和diamond结果一致

#统计不同species的所有assembly
cat strains.taxon.tsv | tsv-summarize -g 4 --count | keep-header -- tsv-sort -k2,2n >MltB/assembly_genus_num.tsv

#两者合并并且相除计算copy per genome
cat MltB/MltB_hmmer_genus_num.tsv | 
tsv-join -d 1 \
-f MltB/assembly_genus_num.tsv -k 1 \
--append-fields 2 |
tsv-select -f 1,3,2  | tr "\t" "," |
perl -F, -alne '$per=$F[2]/$F[1];$per=sprintf "%.1f",$per;print"$F[0]\t$F[1]\t$F[2]\t$per";' >MltB/MltB_assem_copy.tsv

sed -i '1ispecies\tnumber of assemblies\tnumber of MltB\taverage per genome' MltB/MltB_assem_copy.tsv 
cat MltB/MltB_assem_copy.tsv |  keep-header -- tsv-sort -k4,4rn >MltB/MltB_assembly_copy.tsv
plotr tsv MltB/MltB_assembly_copy.tsv --header

```
# 3.两种树(模式生物的)
```BASH
#模式细菌的bac120蛋白树
cat strains.taxon.tsv |
    grep -v "GCF" | 
    cut -f 1 > model.lst 
提取序列
faops some PROTEINS/bac120.trim.fa model.lst  PROTEINS/bac120.model.fa
muscle -in PROTEINS/bac120.model.fa  -out PROTEINS/bac120.model.aln.fa
#建树
FastTree PROTEINS/bac120.model.aln.fa > PROTEINS/bac120.model.aln.newick
nw_reroot  PROTEINS/bac120.model.aln.newick $(nw_labels PROTEINS/bac120.model.aln.newick | grep -E "B_sub|St_aur" ) |
    nw_order -c n - \
    > PROTEINS/bac120.model.reroot.newick
 
#模式生物的MltB蛋白树
cat  MltB/MltB_minevalue.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"MltB:_lytic_murein_transglycosylase_B" |
grep -f model.lst | grep -v 'GCF' | cut -f 1 >MltB/MltB.model.tsv
#提取序列
faops some PROTEINS/all.replace.fa MltB/MltB.model.tsv MltB/MltB.model.fa
muscle -in MltB/MltB.model.fa -out MltB/MltB.model.aln.fa
#建树
FastTree MltB/MltB.model.aln.fa > MltB/MltB.model.aln.newick
nw_reroot MltB/MltB.model.aln.newick $(nw_labels MltB/MltB.model.aln.newick | grep "PAO1" ) |
    nw_order -c n - \
    > MltB/MltB.model.reroot.newick


###模式生物的两种树
setwd("D:/")
library(dplyr)
library(ggtree)
library(ape)
tree1=read.tree("bac120.model.reroot.newick")
tree2=read.tree("MltB.model.reroot.newick")
p1 <- ggtree(tree1)
p2 <- ggtree(tree2)
d1 <- p1$data
d2 <- p2$data
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1 #翻转第二棵树
d2$y <- d2$y + 8 #将第二棵树向上移动对齐
p3 <- ggtree::rotate(p1,29) #旋转node，使两棵树拓扑结构一致
pp <- p3 + geom_tiplab(offset=0.05) + geom_treescale()+ geom_highlight(node=25,fill="red")+ geom_tree(data=d2) + geom_tiplab(data = d2, hjust=1, offset =-0.05)+geom_treescale(x=2)
```

# 4.提取所有的DNA序列
```BASH
for GENUS in $(cat genus.lst); do
    echo 1>&2 "==> GENUS [${GENUS}]"

    for STRAIN in $(cat taxon/${GENUS}); do
        find ASSEMBLY/${STRAIN} -type f -name "*_genomic.fna.gz" |
            grep -v "_from_" |
            xargs gzip -dcf 
done
done \
>DNA/all.DNA.fa

cat DNA/all.DNA.fa |
    perl -nl -e '
        BEGIN { our %seen; our $h; }

        if (/^>/) {
            $h = (split(" ", $_))[0];
            $seen{$h}++;
            $_ = $h;
        }
        print if $seen{$h} == 1;
    ' \
    > DNA/all.uniq.DNA.fa
```



