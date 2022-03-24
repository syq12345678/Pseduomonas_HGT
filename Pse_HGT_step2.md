
[TOC]


# 1.使用hmmscan搜索蛋白序列相应的domian
## 1.1使用tigrfams数据库
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

```

### 1.1.1将Glycerol提取的蛋白序列与tigerfam数据库比对
```
# 将Glycerol提取的蛋白序列与tigerfam数据库比对
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 Glycerol/Glycerol_tigrfam.replace.tsv)  Glycerol/Glycerol.fa

E_VALUE=1e-10
NAME=Glycerol
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.txt --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/TIGERFAM/tigrfams.hmm  ${NAME}/${NAME}.fa

perl abstract1.pl Glycerol/Glycerol.tbl >Glycerol/Glycerol.abstract.tsv
tsv-filter --le 4:1e-50 Glycerol/Glycerol.abstract.tsv >Glycerol/Glycerol.cutoff.tsv
tsv-summarize  -g 8 --count Glycerol/Glycerol.cutoff.tsv
#查看以e值过滤后的结果
glycerol_kin:_glycerol_kinase   2393
XylB:_xylulokinase      3285
gntK_FGGY:_gluconate_kinase     229
L-ribulokinase:_ribulokinase    10
5C_CHO_kinase:_FGGY-family_pentulose_kinase     16
fuculo_kin_coli:_L-fuculokinase 14

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare1.pl Glycerol/Glycerol.cutoff.tsv > Glycerol/Glycerol_minevalue.tsv
tsv-summarize  -g 3 --count Glycerol/Glycerol_minevalue.tsv
#查看保留最显著的且符合description的结果
glycerol_kin:_glycerol_kinase   2318
XylB:_xylulokinase      1349
L-ribulokinase:_ribulokinase    1
gntK_FGGY:_gluconate_kinase     3
fuculo_kin_coli:_L-fuculokinase 4
5C_CHO_kinase:_FGGY-family_pentulose_kinase     13

#拼接属名等信息并统计拷贝数
cat Glycerol/Glycerol_minevalue.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"glycerol_kin:_glycerol_kinase" |
tsv-join -d 1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 3 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-summarize -g 3,4 --count |
keep-header -- tsv-sort -k3,3n >Glycerol/Glycerol_hmmscan_copy.tsv

#统计拷贝数的分布
tsv-summarize -g 3,2 --count  Glycerol/Glycerol_hmmscan_copy.tsv > Glycerol/Glycerol_hmmscan_GCF_copy.tsv
sed -i '1icopy\tgenus\tGCF'  Glycerol/Glycerol_hmmscan_GCF_copy.tsv

与预期的有比较大的区别
```

### 1.1.2将Glyoxalase提取的蛋白序列与tigerfam数据库比对
```
# 将Glyoxalase提取的蛋白序列与tigerfam数据库比对
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 Glyoxalase/Glyoxalase.replace.tsv)  Glyoxalase/Glyoxalase.fa

E_VALUE=1e-10
NAME=Glyoxalase
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.txt --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/TIGERFAM/tigrfams.hmm  ${NAME}/${NAME}.fa

perl abstract1.pl Glyoxalase/Glyoxalase.tbl >Glyoxalase/Glyoxalase.abstract.tsv
tsv-filter --le 4:1e-50 Glyoxalase/Glyoxalase.abstract.tsv >Glyoxalase/Glyoxalase.cutoff.tsv
tsv-summarize  -g 8 --count Glyoxalase/Glyoxalase.cutoff.tsv
#查看以e值过滤后的结果
glyox_I:_lactoylglutathione_lyase       1838
VI_Rhs_Vgr:_type_VI_secretion_system_Vgr_family_protein 1
vgr_GE:_Rhs_element_Vgr_protein 1

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare1.pl Glyoxalase/Glyoxalase.cutoff.tsv > Glyoxalase/Glyoxalase_minevalue.tsv
tsv-summarize  -g 3 --count Glyoxalase/Glyoxalase_minevalue.tsv
#查看保留最显著的且符合description的结果
glyox_I:_lactoylglutathione_lyase       1838
VI_Rhs_Vgr:_type_VI_secretion_system_Vgr_family_protein 1

#拼接属名等信息并统计拷贝数
cat Glyoxalase/Glyoxalase_minevalue.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"glyox_I:_lactoylglutathione_lyase" |
tsv-join -d 1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 3 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-summarize -g 3,4 --count |
keep-header -- tsv-sort -k3,3n >Glyoxalase/Glyoxalase_hmmscan_copy.tsv

#统计拷贝数的分布
tsv-summarize -g 3,2 --count  Glyoxalase/Glyoxalase_hmmscan_copy.tsv > Glyoxalase/Glyoxalase_hmmscan_GCF_copy.tsv
sed -i '1icopy\tgenus\tGCF'  Glyoxalase/Glyoxalase_hmmscan_GCF_copy.tsv

###与预期的有区别
```

### 1.1.3 将Guanine提取的蛋白序列与tigerfam数据库比对
```
# 将Guanine提取的蛋白序列与tigerfam数据库比对
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 Guanine/Guanine_tigrfam.replace.tsv)  Guanine/Guanine.fa

E_VALUE=1e-10
NAME=Guanine
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.txt --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/TIGERFAM/tigrfams.hmm  ${NAME}/${NAME}.fa


perl abstract1.pl Guanine/Guanine.tbl >Guanine/Guanine.abstract.tsv
tsv-filter --le 4:1e-50 Guanine/Guanine.abstract.tsv >Guanine/Guanine.cutoff.tsv
tsv-summarize  -g 8 --count Guanine/Guanine.cutoff.tsv
#查看以e值过滤后的结果
guan_deamin:_guanine_deaminase  1811
Se_ssnA:_putative_selenium_metabolism_protein_SsnA      52
hutF:_formiminoglutamate_deiminase

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare1.pl Guanine/Guanine.cutoff.tsv > Guanine/Guanine_minevalue.tsv
tsv-summarize  -g 3 --count Guanine/Guanine_minevalue.tsv
#查看保留最显著的且符合description的结果
guan_deamin:_guanine_deaminase  1808
Se_ssnA:_putative_selenium_metabolism_protein_SsnA      46
hutF:_formiminoglutamate_deiminase


#拼接属名等信息并统计拷贝数
cat Guanine/Guanine_minevalue.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"guan_deamin:_guanine_deaminase" |
tsv-join -d 1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 3 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-summarize -g 3,4 --count |
keep-header -- tsv-sort -k3,3n >Guanine/Guanine_hmmscan_copy.tsv


#统计拷贝数的分布
tsv-summarize -g 3,2 --count  Guanine/Guanine_hmmscan_copy.tsv > Guanine/Guanine_hmmscan_GCF_copy.tsv
sed -i '1icopy\tgenus\tGCF'  Guanine/Guanine_hmmscan_GCF_copy.tsv

###与预期的区别不大

```


### 1.1.4 将branched-chain提取的蛋白序列与tigerfam数据库比对
```
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
livcs:_branched-chain_amino_acid_transport_system_II_carrier_protein    1779

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare1.pl branched-chain/branched-chain.cutoff.tsv > branched-chain/branched-chain_minevalue.tsv

#拼接属名等信息并统计拷贝数
cat branched-chain/branched-chain_minevalue.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"branched-chain" |
tsv-join -d 1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 3 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-summarize -g 3,4 --count |
keep-header -- tsv-sort -k3,3n >branched-chain/branched-chain_hmmscan_copy.tsv

#统计拷贝数的分布
tsv-summarize -g 3,2 --count  branched-chain/branched-chain_hmmscan_copy.tsv > branched-chain/branched-chain_hmmscan_GCF_copy.tsv
sed -i '1icopy\tgenus\tGCF'  branched-chain/branched-chain_hmmscan_GCF_copy.tsv

###与预期结果区别不大

```

### 1.1.5 将MltB提取的蛋白序列与tigerfam数据库比对
```
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 MltB/MltB_tigrfam.replace.tsv) MltB/MltB.fa

E_VALUE=1e-10
NAME=MltB
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.txt --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/TIGERFAM/tigrfams.hmm  ${NAME}/${NAME}.fa

perl abstract1.pl MltB/MltB.tbl >MltB/MltB.abstract.tsv
tsv-filter --le 4:1e-50 MltB/MltB.abstract.tsv >MltB/MltB.cutoff.tsv
tsv-summarize  -g 8 --count MltB/MltB.cutoff.tsv
#查看以e值过滤后的结果
MltB:_lytic_murein_transglycosylase_B   3696
MltB_2:_lytic_murein_transglycosylase   3390

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare1.pl MltB/MltB.cutoff.tsv >MltB/MltB_minevalue.tsv
tsv-summarize  -g 3 --count MltB/MltB_minevalue.tsv
#查看保留最显著的且符合description的结果
MltB:_lytic_murein_transglycosylase_B   1790
MltB_2:_lytic_murein_transglycosylase   1989

#拼接属名等信息并统计拷贝数
cat MltB/MltB_minevalue.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"MltB:_lytic_murein_transglycosylase_B" |
tsv-join -d 1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 3 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-summarize -g 3,4 --count |
keep-header -- tsv-sort -k3,3n >MltB/MltB_hmmscan_copy.tsv

#统计拷贝数的分布
tsv-summarize -g 3,2 --count  MltB/MltB_hmmscan_copy.tsv > MltB/MltB_hmmscan_GCF_copy.tsv
sed -i '1icopy\tgenus\tGCF'  MltB/MltB_hmmscan_GCF_copy.tsv

###与预期结果区别不大

```


### 1.1.6 Tyrosine提取的蛋白序列与tigerfam数据库比对
```
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 Tyrosine/Tyrosine.replace.tsv) Tyrosine/Tyrosine.fa

E_VALUE=1e-10
NAME=Tyrosine
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.txt --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/TIGERFAM/tigrfams.hmm  ${NAME}/${NAME}.fa


perl abstract1.pl Tyrosine/Tyrosine.tbl >Tyrosine/Tyrosine.abstract.tsv
tsv-filter --le 4:1e-50 Tyrosine/Tyrosine.abstract.tsv >Tyrosine/Tyrosine.cutoff.tsv
tsv-summarize  -g 8 --count Tyrosine/Tyrosine.cutoff.tsv
#查看以e值过滤后的结果
tyrS:_tyrosine--tRNA_ligase     1912

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare1.pl Tyrosine/Tyrosine.cutoff.tsv >Tyrosine/Tyrosine_minevalue.tsv


#拼接属名等信息并统计拷贝数
cat Tyrosine/Tyrosine_minevalue.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"tyrS" |
tsv-join -d 1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 3 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-summarize -g 3,4 --count |
keep-header -- tsv-sort -k3,3n >Tyrosine/Tyrosine_hmmscan_copy.tsv

#统计拷贝数的分布
tsv-summarize -g 3,2 --count  Tyrosine/Tyrosine_hmmscan_copy.tsv > Tyrosine/Tyrosine_hmmscan_GCF_copy.tsv
sed -i '1icopy\tgenus\tGCF'  Tyrosine/Tyrosine_hmmscan_GCF_copy.tsv

###与预期结果区别不大

```

## 1.2使用pfam数据库
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

```





### 1.2.1 将heat_shock提取的蛋白序列与pfam数据库比对
```
# 将heat_shock提取的蛋白序列与pfam数据库比对
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 heat_shock/heat_shock_pfam.replace.tsv)  heat_shock/heat_shock.fa

E_VALUE=1e-10
NAME=heat_shock
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.txt --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/PFAM/Pfam-A.hmm  ${NAME}/${NAME}.fa


perl abstract1.pl heat_shock/heat_shock.tbl >heat_shock/heat_shock.abstract.tsv
tsv-filter --le 4:1e-50 heat_shock/heat_shock.abstract.tsv >heat_shock/heat_shock.cutoff.tsv
tsv-summarize  -g 8 --count heat_shock/heat_shock.cutoff.tsv
#查看以e值过滤后的结果
Hsp90_protein   1504
#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare1.pl heat_shock/heat_shock.cutoff.tsv >heat_shock/heat_shock_minevalue.tsv

#拼接属名等信息并统计拷贝数
cat heat_shock/heat_shock_minevalue.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"Hsp90_protein" |
tsv-join -d 1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 3 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-summarize -g 3,4 --count |
keep-header -- tsv-sort -k3,3n >heat_shock/heat_shock_hmmscan_copy.tsv

#统计拷贝数的分布
tsv-summarize -g 3,2 --count  heat_shock/heat_shock_hmmscan_copy.tsv > heat_shock/heat_shock_hmmscan_GCF_copy.tsv
sed -i '1icopy\tgenus\tGCF'  heat_shock/heat_shock_hmmscan_GCF_copy.tsv

###与预期结果有差别
```

### 1.2.2 将Proteasome提取的蛋白序列与pfam数据库比对
```
# 将proteasome提取的蛋白序列与pfam数据库比对
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 Proteasome/Proteasome_pfam.replace.tsv)  Proteasome/Proteasome.fa

E_VALUE=1e-10
NAME=Proteasome
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.txt --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/PFAM/Pfam-A.hmm  ${NAME}/${NAME}.fa

perl abstract1.pl Proteasome/Proteasome.tbl >Proteasome/Proteasome.abstract.tsv
#该family较为特别，e值全部在1e-20左右，不能进行过滤
 tsv-summarize  -g 8 --count Proteasome/Proteasome.abstract.tsv
#统计结果其对应的结构域
Proteasome_subunit      842

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare1.pl Proteasome/Proteasome.abstract.tsv >Proteasome/Proteasome_minevalue.tsv

#拼接属名等信息并统计拷贝数
cat Proteasome/Proteasome_minevalue.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"Proteasome_subunit" |
tsv-join -d 1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 3 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-summarize -g 3,4 --count |
keep-header -- tsv-sort -k3,3n >Proteasome/Proteasome_hmmscan_copy.tsv

#统计拷贝数的分布
tsv-summarize -g 3,2 --count  Proteasome/Proteasome_hmmscan_copy.tsv > Proteasome/Proteasome_hmmscan_GCF_copy.tsv
sed -i '1icopy\tgenus\tGCF'  Proteasome/Proteasome_hmmscan_GCF_copy.tsv

###与预期结果有差别



```
### 1.2.3 将cobn提取的蛋白序列与pfam数据库比对
```
# 将proteasome提取的蛋白序列与pfam数据库比对
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 cobn/cobn.pfam.replace.tsv)  cobn/cobn.fa

E_VALUE=1e-10
NAME=cobn
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.txt --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/PFAM/Pfam-A.hmm  ${NAME}/${NAME}.fa

perl abstract1.pl cobn/cobn.tbl >cobn/cobn.abstract.tsv
tsv-filter --le 4:1e-50 cobn/cobn.abstract.tsv >cobn/cobn.cutoff.tsv
tsv-summarize  -g 8 --count cobn/cobn.cutoff.tsv
#查看以e值过滤后的结果
CobN/Magnesium_Chelatase        1207
#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare1.pl cobn/cobn.cutoff.tsv >cobn/cobn_minevalue.tsv

#拼接属名等信息并统计拷贝数
cat cobn/cobn_minevalue.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"CobN/Magnesium_Chelatase" |
tsv-join -d 1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 3 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-summarize -g 3,4 --count |
keep-header -- tsv-sort -k3,3n >cobn/cobn_hmmscan_copy.tsv

#统计拷贝数的分布
tsv-summarize -g 3,2 --count cobn/cobn_hmmscan_copy.tsv > cobn/cobn_hmmscan_GCF_copy.tsv
sed -i '1icopy\tgenus\tGCF'  cobn/cobn_hmmscan_GCF_copy.tsv

###符合预期结果


```

## 1.3使用panther数据库
```
#下载panther数据库
mkdir -p ~/data/HMM/PANTHER
cd ~/data/HMM/PANTHER

wget -N ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-12.0.tar.gz
# wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-12.0.tar.gz.md5
# md5sum -c panther-data-12.0.tar.gz.md5
tar -xzvf panther-data-12.0.tar.gz
###注意：
1.panther.hmm文件中没有descrpiton，descptiron在names.tab文件中，并且与PTHR号一一对应
2.panther数据库已经提前格式化好了,不需要使用hmmpres
cd ~/data/Pseudomonas

```

### 1.3.1 将bac_Tyrosine提取的蛋白序列与pfam数据库比对
```
faops some PROTEINS/all.replace.fa <(tsv-select -f 2  bac_Tyrosine/bac_Tyrosine.replace.tsv) bac_Tyrosine/bac_Tyrosine.fa
#panther数据库太大，需要放在超算上跑
bsub -q mpi -n 24 -J "BAC" hmmscan -cpu 20 -E 1e-10 --domE 1e-10 \
--noali -o bac_Tyrosine_progress.txt --tblout bac_Tyrosine.tbl \
  ./panther/panther.hmm  ./bac_Tyrosine.fa