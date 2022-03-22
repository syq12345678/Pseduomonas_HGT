[TOC]


# 1.统计不同菌株中拷贝数的基因家族
注:
homologous superfamily:
1.来自cath-gene3D和superfamily数据库
2.使用underlying profile hidden Markov models而不是单一的模型来表示不同的结构家族
3.广泛使用的hmm数据库主要包括pfam pathern tigerfam和cdd
4.CATH和SCOP是目前最大的人工验证的蛋白质域结构层次分类数据库,这两个数据库都将蛋白质分为同源家族和超家族


| family | count | description | database | 显示family relationship |
| :--- | :--- |  :--- | :--- | :--- |
| IPR014311 | Guanine deaminase | Overlapping homologous superfamilies | panther tigerfam | 无 |
| IPR001404 | Heat shock protein Hsp90 family | Overlapping homologous superfamilies | pfam panther pirsf hampap | 无 |
| IPR004685 | Branched-chain amino acid transport system II carrier protein | none | pathern pfam tigerfam | 无|
| IPR004361 | Glyoxalase I | Overlapping homologous superfamilies | tigrfam | 无 | 
|           |              |             |
| IPR005999 | Glycerol kinase | Overlapping homologous superfamilies | hamap  tigerfam | 有 |
| IPR011757 | Lytic transglycosylase MltB | Overlapping homologous superfamilies | tigerfam panther | 有 |
| IPR001353 | Proteasome, subunit alpha/beta | Overlapping homologous superfamilies | pfam | 有 |
| IPR002307 | Tyrosine-tRNA ligase | Overlapping homologous superfamilies  |cdd tigrfam prints |  有 |
| IPR024088 | Tyrosine-tRNA ligase, bacterial-type | Overlapping homologous superfamilies | pathern | 有 |
| IPR003672 | CobN/magnesium chelatase | none |  pfam cdd pathern | 有 |
|           |            |              |
| IPR037532 | Peptidoglycan D,D-transpeptidase FtsI | Overlapping homologous superfamilies | hampvp |  |
| IPR024922 | Rubredoxin | none | pirsf |   |
| IPR000813 | 7Fe ferredoxin | none | prints |   |
| IPR007416 | YggL 50S ribosome-binding protein | none | pfam |  |




# 2.使用hmmsearch抓取相同domain的基因或者基因家族

## 2.1 IPR014311 鸟嘌呤脱氨酶 Guanine deaminase

```
#查看基因的hmm文件所在数据库
cd ~/data/Pseudomonas
cat STRAINS/Pseudom_aer_PAO1/*.tsv |
    grep "IPR014311"
#可以发现数据库中的登录号是
CDD     cd01303
PANTHER PTHR11271:SF6
TIGRFAM TIGR02967

#下载hmm文件
mkdir -p Guanine/HMM
curl -L https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR02967.1.HMM >Guanine/HMM/Guanine_tigrfam.hmm
curl -L http://www.pantherdb.org/panther/exportHmm.jsp?acc=PTHR11271:SF6 >Guanine/HMM/Guanine_panther.hmm

#在蛋白库里搜索该基因的结构域
E_VALUE=1e-20
for domain in Guanine_panther Guanine_tigrfam; do
    >&2 echo "==> domain [${domain}]"

    if [ -e Guanine/${domain}.replace.tsv ]; then
        continue;
    fi

    for GENUS in $(cat genus.lst); do
        >&2 echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw Guanine/HMM/${domain}.hmm - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '
                    />>\s+(\S+)/ or next;
                    $n = $1;
                    $s = $n;
                    $s =~ s/\.\d+//;
                    printf qq{%s\t%s_%s\n}, $n, $ENV{STRAIN}, $s;
                '
        done
    done \
        > Guanine/${domain}.replace.tsv

    >&2 echo
done



#统计不同数据库根据相同domain抓取出的基因或基因家族序列数目
wc -l Guanine/Guanine_tigrfam.replace.tsv
3583 Guanine/Guanine_tigrfam.replace.tsv
wc -l Guanine/Guanine_panther.replace.tsv
3634 Guanine/Guanine_panther.replace.tsv

tsv-join Guanine/Guanine_panther.replace.tsv \
    -f Guanine/Guanine_tigrfam.replace.tsv \
    > Guanine/Guanine.replace.tsv

wc -l Guanine/*.tsv
  3583 Guanine/Guanine.replace.tsv
  3634 Guanine/Guanine_panther.replace.tsv
  3583 Guanine/Guanine_tigrfam.replace.tsv

#根据ncbi的annotation查看有哪些结构域
cat Guanine/Guanine.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-summarize -g 3 --count
#总共有八种
guanine deaminase       1765
amidohydrolase family protein   175
amidohydrolase  20
cytosine deaminase      1
8-oxoguanine deaminase  830
formimidoylglutamate deiminase  1
TRZ/ATZ family hydrolase        789
chlorohydrolase/deaminase family protein        1
N-ethylammeline chlorohydrolase 1



#拼接结构域等信息
cat Guanine/Guanine.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-filter --str-in-fld 3:"guanine" |
tsv-filter --str-not-in-fld  3:"8-oxoguanine" | 
tsv-join -d  1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 4 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-select -f 1,5,3,2,4 >Guanine/Guanine_no_8-o_summary.tsv

#统计拷贝数
cat Guanine/Guanine_no_8-o_summary.tsv | 
tsv-summarize -g 5,2 --count |
keep-header -- tsv-sort -k3,3n >Guanine/Guanine_no_8-o_copy.tsv

```

## 2.2 IPR001404  Heat shock protein Hsp90 family热休克蛋白Hsp90家族

```
cd ~/data/Pseudomonas
cat STRAINS/Pseudom_aer_PAO1/*.tsv |
    grep "IPR001404"
#可以发现数据库中的登录号是
Pfam    PF00183
PANTHER PTHR11528
PIRSF   PIRSF002583 不能用
HAMAP   MF_00505 不能用

mkdir -p heat_shock/HMM
curl -L http://www.pantherdb.org/panther/exportHmm.jsp?acc=PTHR11528 >heat_shock/HMM/heat_shock_panther.hmm
curl -L https://pfam.xfam.org/family/PF00183/hmm > heat_shock/HMM/heat_shock_pfam.hmm

#在蛋白库里搜索该基因的结构域
E_VALUE=1e-20
for domain in heat_shock_panther heat_shock_pfam; do
    >&2 echo "==> domain [${domain}]"

    if [ -e heat_shock/${domain}.replace.tsv ]; then
        continue;
    fi

    for GENUS in $(cat genus.lst); do
        >&2 echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw heat_shock/HMM/${domain}.hmm - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '
                    />>\s+(\S+)/ or next;
                    $n = $1;
                    $s = $n;
                    $s =~ s/\.\d+//;
                    printf qq{%s\t%s_%s\n}, $n, $ENV{STRAIN}, $s;
                '
        done
    done \
        > heat_shock/${domain}.replace.tsv

    >&2 echo
done


wc -l heat_shock/heat_shock_panther.replace.tsv
1902 heat_shock/heat_shock_panther.replace.tsv
wc -l heat_shock/heat_shock_pfam.replace.tsv
1508 heat_shock/heat_shock_pfam.replace.tsv

tsv-join heat_shock/heat_shock_pfam.replace.tsv \
    -f heat_shock/heat_shock_panther.replace.tsv \
    > heat_shock/heat_shock.replace.tsv

wc -l heat_shock/*.tsv
1508 heat_shock/heat_shock.replace.tsv
  1902 heat_shock/heat_shock_panther.replace.tsv
  1508 heat_shock/heat_shock_pfam.replace.tsv
  
#查看有哪些结构域
cat heat_shock/heat_shock.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-summarize -g 3 --count
总共有六种
molecular chaperone HtpG        1498
heat shock protein 90   3
class III heat-shock protein (ATP-dependent molecular chaperone HSP90)  1
chaperone protein HtpG  4
chaperone protein       1
heat shock protein Hsp90        1


#拼接结构域等信息
cat heat_shock/heat_shock.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-filter --str-in-fld 3:"heat" |
tsv-join -d  1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 4 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-select -f 1,5,3,2,4 >heat_shock/heat_shock_summary.tsv

#统计拷贝数
cat heat_shock/heat_shock_summary.tsv | 
tsv-summarize -g 5,2 --count |
keep-header -- tsv-sort -k3,3n >heat_shock/heat_shock_copy.tsv
```

## 2.3 IPR004685 Branched-chain amino acid 支链氨基酸转运蛋白

```
#查看基因的hmm文件所在数据库
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

#在蛋白库里搜索该基因的结构域
E_VALUE=1e-20
for domain in branched.panther branched.tigerfam branched.pfam; do
    >&2 echo "==> domain [${domain}]"

    if [ -e branched-chain/${domain}.replace.tsv ]; then
        continue;
    fi

    for GENUS in $(cat genus.lst); do
        >&2 echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw branched-chain/HMM/${domain}.hmm - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '
                    />>\s+(\S+)/ or next;
                    $n = $1;
                    $s = $n;
                    $s =~ s/\.\d+//;
                    printf qq{%s\t%s_%s\n}, $n, $ENV{STRAIN}, $s;
                '
        done
    done \
        > branched-chain/${domain}.replace.tsv

    >&2 echo
done



#统计不同数据库根据相同domain抓取出的基因或基因家族序列数目
wc -l branched-chain/*.tsv
     1785 branched-chain/branched.panther.replace.tsv
     1785 branched-chain/branched.pfam.replace.tsv
     1785 branched-chain/branched.tigerfam.replace.tsv

tsv-join branched-chain/branched.tigerfam.replace.tsv \
    -f branched-chain/branched.pfam.replace.tsv \
    > branched-chain/branched.replace.tsv


#根据ncbi的annotation查看有哪些结构域
cat branched-chain/branched.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-summarize -g 3 --count
#总共有十种
branched-chain amino acid transport system II carrier protein   1775
branched chain amino acid transporter   1
low-affinity branched-chain amino acid transporter      1
branched-chain amino acid-Na+ symporter 1
branched-chain amino acid transport system carrier protein      1
branched chain amino acid transporter BrnQ      1
branched chain amino acid transport system II carrier protein   2
branched-chain amino acid transporter   1
branched-chain amino acid transport system 3 carrier protein    1
branched-chain amino acid ABC transporter       1


#拼接结构域等信息
cat branched-chain/branched.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-filter --str-in-fld 3:"branched-chain amino acid transport system II carrier protein" |
tsv-join -d  1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 4 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-select -f 1,5,3,2,4 >branched-chain/branched_summary.tsv

#统计拷贝数
cat branched-chain/branched_summary.tsv | 
tsv-summarize -g 5,2 --count |
keep-header -- tsv-sort -k3,3n >branched-chain/branched_copy.tsv

```

## 2.4.IPR004361  乙二醛酶 Glyoxalase I

```
#在interpro数据库中IPR004361对应的是Glyoxalase I
cd ~/data/Pseudomonas
cat STRAINS/Pseudom_aer_PAO1/*.tsv |
    grep "IPR004361"
#可以发现数据库中的登录号是TIGRFAM TIGR00068
mkdir -p Glyoxalase/HMM
curl -L https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR00068.2.HMM >Glyoxalase/HMM/Glyoxalase.hmm

#在蛋白库里搜索该基因的结构域
E_VALUE=1e-20
for domain in Glyoxalase; do
    >&2 echo "==> domain [${domain}]"

    if [ -e Glyoxalase/${domain}.replace.tsv ]; then
        continue;
    fi

    for GENUS in $(cat genus.lst); do
        >&2 echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw Glyoxalase/HMM/${domain}.hmm - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '
                    />>\s+(\S+)/ or next;
                    $n = $1;
                    $s = $n;
                    $s =~ s/\.\d+//;
                    printf qq{%s\t%s_%s\n}, $n, $ENV{STRAIN}, $s;
                '
        done
    done \
        > Glyoxalase/${domain}.replace.tsv

    >&2 echo
done

wc -l Glyoxalase/Glyoxalase.replace.tsv
2390 Glyoxalase/Glyoxalase.replace.tsv

#查看有哪些结构域
cat Glyoxalase/Glyoxalase.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-summarize -g 3 --count
总共有四种
lactoylglutathione lyase        2340
glyoxalase I    2
VOC family protein      46
glyoxalase      2


#拼接结构域等信息
cat Glyoxalase/Glyoxalase.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-filter --str-in-fld 3:"lactoylglutathione" |
tsv-join -d  1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 4 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-select -f 1,5,3,2,4 >Glyoxalase/Glyoxalase_lactoylglutathione_summary.tsv

#统计拷贝数
cat Glyoxalase/Glyoxalase_glyoxalase_summary.tsv | 
tsv-summarize -g 5,2 --count |
keep-header -- tsv-sort -k3,3n >Glyoxalase/Glyoxalase_glyoxalase_copy.tsv


```



## 2.5.IPR005999  Glycerol kinase甘油激酶

```
cd ~/data/Pseudomonas
cat STRAINS/Pseudom_aer_PAO1/*.tsv |
    grep "IPR005999"
#可以发现数据库中的登录号是
TIGRFAM TIGR01311 
mkdir -p Glycerol/HMM

curl -L https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR01311.1.HMM >Glycerol/HMM/Glycerol_tigrfam.hmm

#在蛋白库里搜索该基因的结构域

E_VALUE=1e-20
for domain in Glycerol_tigrfam; do
    >&2 echo "==> domain [${domain}]"

    if [ -e Glycerol/${domain}.replace.tsv ]; then
        continue;
    fi

    for GENUS in $(cat genus.lst); do
        >&2 echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E 1e-20 --domE 1e-20 --noali --notextw Glycerol/HMM/${domain}.hmm - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '
                    />>\s+(\S+)/ or next;
                    $n = $1;
                    $s = $n;
                    $s =~ s/\.\d+//;
                    printf qq{%s\t%s_%s\n}, $n, $ENV{STRAIN}, $s;
                '
        done
    done \
        > Glycerol/${domain}.replace.tsv

    >&2 echo
done


wc -l Glycerol/Glycerol_tigrfam.replace.tsv
3705 Glycerol/Glycerol_tigrfam.replace.tsv


#查看有哪些结构域
cat Glycerol/Glycerol_tigrfam.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-summarize -g 3 --count
总共有
glycerol kinase GlpK    1836
xylulokinase    770
glycerol kinase 336
glycerol kinase (sn-glycerol-3-phosphate generating)    1
D-gluconate kinase      1
xylulose kinase 3
hydroxylated metabolite kinase  1
L-ribulokinase  1
carbohydrate kinase     63
L-xylulose kinase       1
autoinducer-2 kinase    3
L-fuculokinase  5
putative sugar kinase YgcE      1
kinase  1
FGGY-family carbohydrate kinase 631
FGGY family pentulose kinase    6
Carbohydrate kinase, FGGY-like protein  1
putative pentose kinase 1
putative carbohydrate kinase    1
putative L-xylulose kinase      1
putative kinase 1
ribitol kinase  1
gluconate kinase        2
FGGY family carbohydrate kinase 34
sugar kinase    3

#拼接结构域等信息
cat Glycerol/Glycerol_tigrfam.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-filter --str-in-fld 3:"GlpK" |
tsv-join -d  1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 4 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-select -f 1,5,3,2,4 >Glycerol/Glycerol_summary.tsv

#统计拷贝数
cat Glycerol/Glycerol_summary.tsv | 
tsv-summarize -g 5,2 --count |
keep-header -- tsv-sort -k3,3n >Glycerol/Glycerol_copy.tsv
```

##2.6 IPR001353蛋白酶体Proteasome, subunit alpha/beta

```

cd ~/data/Pseudomonas
cat STRAINS/Pseudom_aer_PAO1/*.tsv |
    grep "IPR001353"
#可以发现数据库中的登录号是 Pfam    PF00227
mkdir -p Proteasome/HMM
curl -L http://pfam.xfam.org/family/PF00227/hmm > Proteasome/HMM/Proteasome_pfam.hmm 


#在蛋白库里搜索该基因的结构域
E_VALUE=1e-20
for domain in Proteasome_pfam; do
    >&2 echo "==> domain [${domain}]"

    if [ -e Proteasome/${domain}.replace.tsv ]; then
        continue;
    fi

    for GENUS in $(cat genus.lst); do
        >&2 echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw Proteasome/HMM/${domain}.hmm - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '
                    />>\s+(\S+)/ or next;
                    $n = $1;
                    $s = $n;
                    $s =~ s/\.\d+//;
                    printf qq{%s\t%s_%s\n}, $n, $ENV{STRAIN}, $s;
                '
        done
    done \
        > Proteasome/${domain}.replace.tsv

    >&2 echo
done

wc -l Proteasome/Proteasome_pfam.replace.tsv
842 Proteasome/Proteasome_pfam.replace.tsv

#查看有哪些结构域
cat Proteasome/Proteasome_pfam.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-summarize -g 3 --count
总共有
ATP-dependent protease subunit HslV     829
two-component ATP-dependent protease (N-terminal serine protease)       1
ATP-dependent protease peptidase subunit        5
ATP-dependent endopeptidase hsl proteolytic subunit HslV        1
ATP-dependent endopeptidase hsl proteolytic subunit     1
peptidase component of the HslVU protease       1
peptidase component of the HslUV protease       1
proteasome subunit alpha        1
proteasome subunit beta 1
HslU--HslV peptidase proteolytic subunit        1


#拼接结构域等信息
cat Proteasome/Proteasome_pfam.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-filter --str-in-fld 3:"ATP-dependent protease subunit HslV" |
tsv-join -d  1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 4 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-select -f 1,5,3,2,4 >Proteasome/Proteasome_summary.tsv

#统计拷贝数
cat Proteasome/Proteasome_summary.tsv | 
tsv-summarize -g 5,2 --count |
keep-header -- tsv-sort -k3,3n >Proteasome/Proteasome_copy.tsv
```


##2.7 IPR011757 Lytic transglycosylase MltB

```
cd ~/data/Pseudomonas
cat STRAINS/Pseudom_aer_PAO1/*.tsv |
    grep "IPR011757"
#可以发现数据库中的登录号是

mkdir -p MltB/HMM
curl -L https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR02282.1.HMM >MltB/HMM/MltB_tigrfam.hmm



#在蛋白库里搜索该基因的结构域
E_VALUE=1e-20
for domain in MltB_tigrfam; do
    >&2 echo "==> domain [${domain}]"

    if [ -e MltB/${domain}.replace.tsv ]; then
        continue;
    fi

    for GENUS in $(cat genus.lst); do
        >&2 echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw MltB/HMM/${domain}.hmm - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '
                    />>\s+(\S+)/ or next;
                    $n = $1;
                    $s = $n;
                    $s =~ s/\.\d+//;
                    printf qq{%s\t%s_%s\n}, $n, $ENV{STRAIN}, $s;
                '
        done
    done \
        > MltB/${domain}.replace.tsv

    >&2 echo
done

wc -l MltB/MltB_tigrfam.replace.tsv
3790 MltB/MltB_tigrfam.replace.tsv
  

#查看有哪些结构域
cat MltB/MltB_tigrfam.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-summarize -g 3 --count
#总共有八种
lytic murein transglycosylase B 1787
lytic murein transglycosylase   1990
membrane-bound lytic murein transglycosylase B  6
transglycosylase        2
soluble lytic transglycosylase B        1
murein hydrolase B      2
hypothetical protein PA3992     1
murein transglycosylase B       1

#拼接结构域等信息
cat MltB/MltB_tigrfam.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-filter --str-in-fld 3:"lytic murein transglycosylase B" |
tsv-filter --str-not-in-fld  3:"membrane-bound" | 
tsv-join -d  1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 4 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-select -f 1,5,3,2,4 >MltB/MltB_have_B_no_membrane_summary.tsv

#统计拷贝数
cat MltB/MltB_have_B_no_membrane_summary.tsv | 
tsv-summarize -g 5,2 --count |
keep-header -- tsv-sort -k3,3n >MltB/MltB_have_B_no_membrane_copy.tsv
```


## 2.8.IPR002307 络氨酸TRNA连接酶Tyrosine-tRNA ligase

```
cd ~/data/Pseudomonas
cat STRAINS/Pseudom_aer_PAO1/*.tsv |
    grep "IPR002307"
#可以发现数据库中的登录号是
TIGRFAM TIGR00234
CDD     cd00805
PRINTS  PR01040
mkdir -p Tyrosine/HMM
curl -L https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR00234.1.HMM >Tyrosine/HMM/Tyrosine.hmm
 
#在蛋白库里搜索该基因的结构域
E_VALUE=1e-20
for domain in Tyrosine; do
    >&2 echo "==> domain [${domain}]"

    if [ -e Tyrosine/${domain}.replace.tsv ]; then
        continue;
    fi

    for GENUS in $(cat genus.lst); do
        >&2 echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw Tyrosine/HMM/${domain}.hmm - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '
                    />>\s+(\S+)/ or next;
                    $n = $1;
                    $s = $n;
                    $s =~ s/\.\d+//;
                    printf qq{%s\t%s_%s\n}, $n, $ENV{STRAIN}, $s;
                '
        done
    done \
        > Tyrosine/${domain}.replace.tsv

    >&2 echo
done

wc -l Tyrosine/Tyrosine.replace.tsv
1916 Tyrosine/Tyrosine.replace.tsv

#查看有哪些结构域
cat Tyrosine/Tyrosine.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-summarize -g 3 --count
总共有
tyrosine--tRNA ligase   1907
tyrosyl-tRNA synthetase 9


#拼接结构域等信息
cat Tyrosine/Tyrosine.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-filter --str-in-fld 3:"tyrosine--tRNA ligase" |
tsv-join -d  1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 4 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-select -f 1,5,3,2,4 >Tyrosine/Tyrosine_summary.tsv

#统计拷贝数
cat Tyrosine/Tyrosine_summary.tsv | 
tsv-summarize -g 5,2 --count |
keep-header -- tsv-sort -k3,3n >Tyrosine/Tyrosine_copy.tsv
```

## 2.9 IPR024088  细菌中的络氨酸TRNA连接酶Tyrosine-tRNA ligase, bacterial-type

```
cd ~/data/Pseudomonas
cat STRAINS/Pseudom_aer_PAO1/*.tsv |
    grep "IPR024088"
#可以发现数据库中的登录号是
PANTHER PTHR11766
mkdir -p bac_Tyrosine/HMM
curl -L http://www.pantherdb.org/panther/exportHmm.jsp?acc=PTHR11766 >bac_Tyrosine/HMM/bac_Tyrosine.hmm
 
#在蛋白库里搜索该基因的结构域
E_VALUE=1e-20
for domain in bac_Tyrosine; do
    >&2 echo "==> domain [${domain}]"

    if [ -e bac_Tyrosine/${domain}.replace.tsv ]; then
        continue;
    fi

    for GENUS in $(cat genus.lst); do
        >&2 echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw bac_Tyrosine/HMM/${domain}.hmm - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '
                    />>\s+(\S+)/ or next;
                    $n = $1;
                    $s = $n;
                    $s =~ s/\.\d+//;
                    printf qq{%s\t%s_%s\n}, $n, $ENV{STRAIN}, $s;
                '
        done
    done \
        > bac_Tyrosine/${domain}.replace.tsv

    >&2 echo
done

wc -l bac_Tyrosine/bac_Tyrosine.replace.tsv
1916 bac_Tyrosine/bac_Tyrosine.replace.tsv

#查看有哪些结构域
cat bac_Tyrosine/bac_Tyrosine.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-summarize -g 3 --count
总共有
tyrosine--tRNA ligase   1907
tyrosyl-tRNA synthetase 9


#拼接结构域等信息
cat bac_Tyrosine/bac_Tyrosine.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-filter --str-in-fld 3:"tyrosine--tRNA ligase" |
tsv-join -d  1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 4 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-select -f 1,5,3,2,4 >bac_Tyrosine/bac_Tyrosine_summary.tsv

#统计拷贝数
cat bac_Tyrosine/bac_Tyrosine_summary.tsv | 
tsv-summarize -g 5,2 --count |
keep-header -- tsv-sort -k3,3n >bac_Tyrosine/bac_Tyrosine_copy.tsv
```

## 2.10  IPR003672 镁螯合酶 CobN/magnesium chelatase

```
cd ~/data/Pseudomonas
cat STRAINS/Pseudom_aer_PAO1/*.tsv |
    grep "IPR003672"
#可以发现数据库中的登录号是
Pfam    PF02514

mkdir -p cobn/HMM
curl -L http://pfam.xfam.org/family/PF02514/hmm >cobn/HMM/cobn.pfam.hmm

 
#在蛋白库里搜索该基因的结构域
E_VALUE=1e-20
for domain in cobn.pfam; do
    >&2 echo "==> domain [${domain}]"

    if [ -e cobn/${domain}.replace.tsv ]; then
        continue;
    fi

    for GENUS in $(cat genus.lst); do
        >&2 echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw cobn/HMM/${domain}.hmm - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '
                    />>\s+(\S+)/ or next;
                    $n = $1;
                    $s = $n;
                    $s =~ s/\.\d+//;
                    printf qq{%s\t%s_%s\n}, $n, $ENV{STRAIN}, $s;
                '
        done
    done \
        > cobn/${domain}.replace.tsv

    >&2 echo
done

wc -l cobn/cobn.pfam.replace.tsv
1211 cobn/cobn.pfam.replace.tsv

#查看有哪些结构域
cat cobn/cobn.pfam.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-summarize -g 3 --count
总共有
cobaltochelatase subunit CobN   1210
cobalamin biosynthesis protein CobN     1


#拼接结构域等信息
cat cobn/cobn.pfam.replace.tsv | tsv-select -f 2,1 |
tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 \
--append-fields 2 |
tsv-filter --str-in-fld 3:"cobaltochelatase subunit CobN" |
tsv-join -d  1 \
-f PROTEINS/strain.tsv -k 2 \
--append-fields 3 |
tsv-join -d 4 \
-f strains.taxon.tsv -k 1 \
--append-fields 4 | 
tsv-select -f 1,5,3,2,4 >cobn/cobn_summary.tsv

#统计拷贝数
cat cobn/cobn_summary.tsv | 
tsv-summarize -g 5,2 --count |
keep-header -- tsv-sort -k3,3n >cobn/cobn_copy.tsv
```