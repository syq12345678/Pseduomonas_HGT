1.由于外类群不够，新添了一些基因组文件，由1514增加至1953个基因组文件，由此需要重新操作  
2.hmmsearch搜索结构域，hmmscan过滤,diamond重复计算拷贝数，绘制拷贝数的表格，绘制Reference Genome的树和Representative Genome与  
Reference Genome的树  
3.目前以branched为例,branched在pfam panthern和tigrfam中均存在  

# 1.codes from Professor Wang
## 1.外类群：使用一些模型生物作为外类群，具体信息在reference.tsv中，共有15种
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
```
#包含NCBI的样本信息且去除了错误的菌株基因组文件信息，共1953个  
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

# 统计branched-chain在菌株中拷贝分布
## 2.1使用hmmsearch抓取相同domain的基因或者基因家族
```
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
```
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
```
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
## 2.4 多轮diamond
```



```


# 3.两种树
```
#模式细菌的bac120蛋白树
cat strains.taxon.tsv |
    grep -v "GCF" | 
    cut -f 1 > model.lst 
提取序列
faops some PROTEINS/bac120.trim.fa model.lst  PROTEINS/bac120.model.fa
muscle -in PROTEINS/bac120.model.fa  -out PROTEINS/bac120.model.aln.fa
#建树
FastTree PROTEINS/bac120.model.aln.fa > PROTEINS/bac120.model.aln.newick
nw_reroot  PROTEINS/bac120.model.aln.newick $(nw_labels PROTEINS/bac120.model.aln.newick | grep -E "Bac_subti|Sta_aure" ) |
    nw_order -c n - \
    > PROTEINS/bac120.model.reroot.newick
 

#模式生物的branched蛋白树
cat branched-chain/branched-chain_minevalue.tsv | grep -f model.lst | grep -v 'GCF' | cut -f 1 >branched-chain/branched-chain.model.tsv
#提取序列
faops some PROTEINS/all.replace.fa branched-chain/branched-chain.model.tsv  branched-chain/branched-chain.model.fa
muscle -in branched-chain/branched-chain.model.fa -out branched-chain/branched-chain.model.aln.fa
#建树
FastTree branched-chain/branched-chain.model.aln.fa > branched-chain/branched-chain.model.aln.newick
nw_reroot branched-chain/branched-chain.model.aln.newick $(nw_labels branched-chain/branched-chain.model.aln.newick | grep -E "POA1") |
    nw_order -c n - \
    > branched-chain/branched-chain.model.reroot.newick


#参考菌株的branced蛋白树
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

#模式生物的bac120蛋白树和模式生物的branched蛋白树
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
Ka<<Ks或者Ka/Ks << 1，基因受纯化选择(purify selection)  
9.实际情况下Ka/Ks << 1，因为一般非同义替换带来的都是有害的性状，只有极少数情况下会造成进化上的优势  
10.当Ka/Ks>>1时，基因受到强烈正选择，这样的基因即为近期正在快速进化的基因，对于物种的进化有着非常重要的意义
11.计算ka和ks的步骤
步骤一：假设比较的两条DNA序列之间的长度数为n,它们之间的替换数为m。计算同义(S)和非同义(N)位点的数量(S+N=n)以及同义(Sd)和非同义替换的数量  
(Sd+Nd=m)
步骤二:校正多个替换后,(Nd/N)和(Sd/S)分别代表ka和ks


  
  
 



