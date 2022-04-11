<!-- TOC -->

- [1.统计不同菌株中拷贝数的基因家族](#1统计不同菌株中拷贝数的基因家族)
- [2.使用hmmsearch抓取相同domain的基因或者基因家族](#2使用hmmsearch抓取相同domain的基因或者基因家族)
  - [2.1 IPR014311 鸟嘌呤脱氨酶 Guanine deaminase](#21-ipr014311-鸟嘌呤脱氨酶-guanine-deaminase)
  - [2.2 IPR001404  Heat shock protein Hsp90 family热休克蛋白Hsp90家族](#22-ipr001404--heat-shock-protein-hsp90-family热休克蛋白hsp90家族)
  - [2.3 IPR004685 Branched-chain amino acid 支链氨基酸转运蛋白](#23-ipr004685-branched-chain-amino-acid-支链氨基酸转运蛋白)
  - [2.4.IPR004361  乙二醛酶 Glyoxalase I](#24ipr004361--乙二醛酶-glyoxalase-i)
  - [2.5.IPR005999  Glycerol kinase甘油激酶](#25ipr005999--glycerol-kinase甘油激酶)
  - [2.6 IPR001353蛋白酶体Proteasome, subunit alpha/beta](#26-ipr001353蛋白酶体proteasome-subunit-alphabeta)
  - [2.7 IPR011757 Lytic transglycosylase MltB](#27-ipr011757-lytic-transglycosylase-mltb)
  - [2.8.IPR002307 络氨酸TRNA连接酶Tyrosine-tRNA ligase](#28ipr002307-络氨酸trna连接酶tyrosine-trna-ligase)
  - [2.9 IPR024088  细菌中的络氨酸TRNA连接酶Tyrosine-tRNA ligase, bacterial-type](#29-ipr024088--细菌中的络氨酸trna连接酶tyrosine-trna-ligase-bacterial-type)
  - [2.10  IPR003672 镁螯合酶 CobN/magnesium chelatase](#210--ipr003672-镁螯合酶-cobnmagnesium-chelatase)
- [3.使用hmmscan搜索结构域](#3使用hmmscan搜索结构域)
  - [3.1使用pfam数据库](#31使用pfam数据库)
    - [3.1.1 将heat_shock提取的蛋白序列与pfam数据库比对](#311-将heat_shock提取的蛋白序列与pfam数据库比对)
  - [3.2使用tigrfams数据库](#32使用tigrfams数据库)
    - [3.2.1将Glycerol提取的蛋白序列与tigerfam数据库比对](#321将glycerol提取的蛋白序列与tigerfam数据库比对)
    - [3.2.2将Glyoxalase提取的蛋白序列与tigerfam数据库比对](#322将glyoxalase提取的蛋白序列与tigerfam数据库比对)
    - [3.2.3 将Guanine提取的蛋白序列与tigerfam数据库比对](#323-将guanine提取的蛋白序列与tigerfam数据库比对)
    - [3.2.3 将branched-chain提取的蛋白序列与tigerfam数据库比对](#323-将branched-chain提取的蛋白序列与tigerfam数据库比对)
    - [3.2.4 将MltB提取的蛋白序列与tigerfam数据库比对](#324-将mltb提取的蛋白序列与tigerfam数据库比对)
    - [3.2.5 Tyrosine提取的蛋白序列与tigerfam数据库比对](#325-tyrosine提取的蛋白序列与tigerfam数据库比对)
  - [3.3使用panther数据库](#33使用panther数据库)

<!-- /TOC -->
# 1.统计不同菌株中拷贝数的基因家族
注:homologous superfamily:
* 1.来自cath-gene3D和superfamily数据库
* 2.使用underlying profile hidden Markov models而不是单一的模型来表示不同的结构家族
* 3.广泛使用的hmm数据库主要包括pfam pathern tigerfam和cdd
* 4.CATH和SCOP是目前最大的人工验证的蛋白质域结构层次分类数据库,这两个数据库都将蛋白质分为同源家族和超家族

| family | count | description | database | 显示family relationship |
| :--- | :--- |  :--- | :--- | :--- |
| IPR014311 | Guanine deaminase | Overlapping homologous superfamilies | panther tigerfam | 无 |
| IPR001404 | Heat shock protein Hsp90 family | Overlapping homologous superfamilies | pfam panther pirsf hamap | 无 |
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
| IPR037532 | Peptidoglycan D,D-transpeptidase FtsI | Overlapping homologous superfamilies | hamap| |
| IPR024922 | Rubredoxin | none | pirsf |   |
| IPR000813 | 7Fe ferredoxin | none | prints |   |
| IPR007416 | YggL 50S ribosome-binding protein | none | pfam |  |


# 2.使用hmmsearch抓取相同domain的基因或者基因家族
## 2.1 IPR014311 鸟嘌呤脱氨酶 Guanine deaminase
```BASH
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

#统计拷贝数的分布
tsv-summarize -g 3,2 --count Guanine_copy.tsv >Guanine_copy_GCF_num.tsv
sed -i '1icopy\tgenus\tGCF' Guanine_copy_GCF_num.tsv

```

## 2.2 IPR001404  Heat shock protein Hsp90 family热休克蛋白Hsp90家族

```BASH
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

#统计拷贝数的分布
tsv-summarize -g 3,2 --count heat_shock_copy.tsv >heat_shock_copy_GCF_num.tsv
sed -i '1icopy\tgenus\tGCF' heat_shock_copy_GCF_num.tsv

```

## 2.3 IPR004685 Branched-chain amino acid 支链氨基酸转运蛋白

```BASH
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

#统计拷贝数的分布
tsv-summarize -g 3,2 --count branched_copy.tsv >branched_copy_GCF_num.tsv
sed -i '1icopy\tgenus\tGCF' branched_copy_GCF_num.tsv

```

## 2.4.IPR004361  乙二醛酶 Glyoxalase I

```BASH
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

#统计拷贝数的分布
tsv-summarize -g 3,2 --count Glyoxalase_lactoylglutathione.copy.tsv >Glyoxalase_lactoylglutathione_copy_GCF_num.tsv
sed -i '1icopy\tgenus\tGCF' Glyoxalase_lactoylglutathione_copy_GCF_num.tsv

```



## 2.5.IPR005999  Glycerol kinase甘油激酶

```BASH
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

#统计拷贝数的分布
tsv-summarize -g 3,2 --count Glycerol_copy.tsv >Glycerol_copy_GCF_num.tsv
sed -i '1icopy\tgenus\tGCF' Glycerol_copy_GCF_num.tsv


```

##  2.6 IPR001353蛋白酶体Proteasome, subunit alpha/beta

```BASH

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


#统计拷贝数的分布
tsv-summarize -g 3,2 --count  Proteasome_copy.tsv > Proteasome_copy_GCF_num.tsv
sed -i '1icopy\tgenus\tGCF'  Proteasome_copy_GCF_num.tsv


```


##  2.7 IPR011757 Lytic transglycosylase MltB

```BASH
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

#统计拷贝数的分布
tsv-summarize -g 3,2 --count  MltB_no_B_copy.tsv > MltB_no_B_copy_GCF_num.tsv
sed -i '1icopy\tgenus\tGCF'  MltB_no_B_copy_GCF_num.tsv

```

## 2.8.IPR002307 络氨酸TRNA连接酶Tyrosine-tRNA ligase

```BASH
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

#统计拷贝数的分布
tsv-summarize -g 3,2 --count  Tyrosine_copy.tsv > Tyrosine_copy_GCF_num.tsv
sed -i '1icopy\tgenus\tGCF'  Tyrosine_copy_GCF_num.tsv


```

## 2.9 IPR024088  细菌中的络氨酸TRNA连接酶Tyrosine-tRNA ligase, bacterial-type

```BASH
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

#统计拷贝数的分布
tsv-summarize -g 3,2 --count  bac_Tyrosine_copy.tsv > bac_Tyrosine_copy_GCF_num.tsv
sed -i '1icopy\tgenus\tGCF'  bac_Tyrosine_copy_GCF_num.tsv

```

## 2.10  IPR003672 镁螯合酶 CobN/magnesium chelatase

```BASH
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

#统计拷贝数的分布
tsv-summarize -g 3,2 --count  cobn_copy.tsv > cobn_copy_GCF_num.tsv
sed -i '1icopy\tgenus\tGCF'  cobn_copy_GCF_num.tsv

```

# 3.使用hmmscan搜索结构域
## 3.1使用pfam数据库
```BASH
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
### 3.1.1 将heat_shock提取的蛋白序列与pfam数据库比对
```BASH
# 将heat_shock提取的蛋白序列与pfam数据库比对
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 heat_shock/heat_shock_pfam.replace.tsv)  heat_shock/heat_shock.fa

E_value=1e-10
NAME=heat_shock
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/PFAM/Pfam-A.hmm  heat_shock/heat_shock.fa

perl abstract.pl heat_shock/heat_shock.tbl >heat_shock/heat_shock.abstract.tsv
#查看heat_shock同家族蛋白以及涉及到的数据库登录号
cat  heat_shock/heat_shock.abstract.tsv | tsv-summarize -g 2,6  --count
PF00183.18      Hsp90_protein   1504
PF02518.26      Histidine_kinase-,_DNA_gyrase_B-,_and_HSP90-like_ATPase 1503
PF13589.6       Histidine_kinase-,_DNA_gyrase_B-,_and_HSP90-like_ATPase 1368
#查看domain以及domain的描述
cat  heat_shock/heat_shock.abstract.tsv | tsv-summarize -g 1,6 --count
HSP90   Hsp90_protein   1504
HATPase_c       Histidine_kinase-,_DNA_gyrase_B-,_and_HSP90-like_ATPase 1503
HATPase_c_3     Histidine_kinase-,_DNA_gyrase_B-,_and_HSP90-like_ATPase 1368

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare.pl heat_shock/heat_shock.abstract.tsv >heat_shock/heat_shock_minevalue.tsv
tsv-summarize -g -g 3 --count heat_shock/heat_shock_minevalue.tsv
#Hsp90_protein   1504
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
tsv-summarize -g 3,2 --count  heat_shock_hmmscan_copy.tsv > heat_shock_hmmscan_GCF_copy.tsv
sed -i '1icopy\tgenus\tGCF'  heat_shock_hmmscan_GCF_copy.tsv

###最终结果与annotion过滤的结果无较大区别

```

## 3.2使用tigrfams数据库
```BASH
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

### 3.2.1将Glycerol提取的蛋白序列与tigerfam数据库比对
```BASH
# 将Glycerol提取的蛋白序列与tigerfam数据库比对
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 Glycerol/Glycerol_tigrfam.replace.tsv)  Glycerol/Glycerol.fa

E_VALUE=1e-40
NAME=Glycerol
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.txt --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/TIGERFAM/tigrfams.hmm  ${NAME}/${NAME}.fa

perl abstract.pl Glycerol/Glycerol.tbl >Glycerol/Glycerol.abstract.tsv

#查看Glycerol同家族蛋白以及涉及到的数据库登录号
cat  Glycerol/Glycerol.abstract.tsv | tsv-summarize -g 2,6  --count
TIGR01311       glycerol_kin:_glycerol_kinase   3695
TIGR01312       XylB:_xylulokinase      3695
TIGR01314       gntK_FGGY:_gluconate_kinase     3695
TIGR02628       fuculo_kin_coli:_L-fuculokinase 3680
TIGR01315       5C_CHO_kinase:_FGGY-family_pentulose_kinase     3282
TIGR02627       rhamnulo_kin:_rhamnulokinase    2694
TIGR01234       L-ribulokinase:_ribulokinase    3218

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare.pl Glycerol/Glycerol.abstract.tsv Glycerol/Glycerol_minevalue.tsv
tsv-summarize  -g 3 --count Glycerol/Glycerol_minevalue.tsv
#
glycerol_kin:_glycerol_kinase   2318
XylB:_xylulokinase      1356
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

#####!!!!!注意：Glycerol kinase是Overlapping homologous superfamilies且有famliy relationship
family
IPR000577 Carbohydrate kinase, FGGY
IPR005999 Glycerol kinase
IPR006000 Xylulokinase

#统计拷贝数的分布
tsv-summarize -g 3,2 --count  Glycerol_hmmscan_copy.tsv > Glycerol_hmmscan_GCF_copy.tsv
sed -i '1icopy\tgenus\tGCF'  Glycerol_hmmscan_GCF_copy.tsv

###与annotaion注释结果区别较大，可能是由于glycerol_kin:_glycerol_kinase中包含FGGY-famliy未被过滤
###需要查看以上基因的结构域比较分析
```

### 3.2.2将Glyoxalase提取的蛋白序列与tigerfam数据库比对
```BASH
# 将Glyoxalase提取的蛋白序列与tigerfam数据库比对
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 Glyoxalase/Glyoxalase.replace.tsv)  Glyoxalase/Glyoxalase.fa

E_VALUE=1e-10
NAME=Glyoxalase
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.txt --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/TIGERFAM/tigrfams.hmm  ${NAME}/${NAME}.fa

perl abstract.pl Glyoxalase/Glyoxalase.tbl >Glyoxalase/Glyoxalase.abstract.tsv

#查看Glyoxalase同家族蛋白以及涉及到的数据库登录号
cat  Glyoxalase/Glyoxalase.abstract.tsv | tsv-summarize -g 2,6  --count
TIGR00068       glyox_I:_lactoylglutathione_lyase       2386
TIGR02295       HpaD:_3,4-dihydroxyphenylacetate_2,3-dioxygenase        10
TIGR03361       VI_Rhs_Vgr:_type_VI_secretion_system_Vgr_family_protein 1
TIGR01646       vgr_GE:_Rhs_element_Vgr_protein 1
TIGR03081       metmalonyl_epim:_methylmalonyl-CoA_epimerase    7

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare.pl Glyoxalase/Glyoxalase.abstract.tsv >Glyoxalase/Glyoxalase_minevalue.tsv
tsv-summarize  -g 3 --count Glyoxalase/Glyoxalase_minevalue.tsv
#查看需要过滤的基因
glyox_I:_lactoylglutathione_lyase       2385
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
tsv-summarize -g 3,2 --count  Glyoxalase_hmmscan_copy.tsv > Glyoxalase_hmmscan_GCF_copy.tsv
sed -i '1icopy\tgenus\tGCF'  Glyoxalase_hmmscan_GCF_copy.tsv

###与annotaion过滤的结果区别不大
```

### 3.2.3 将Guanine提取的蛋白序列与tigerfam数据库比对
```BASH
# 将Guanine提取的蛋白序列与tigerfam数据库比对
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 Guanine/Guanine_tigrfam.replace.tsv)  Guanine/Guanine.fa

E_VALUE=1e-10
NAME=Guanine
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.txt --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/TIGERFAM/tigrfams.hmm  ${NAME}/${NAME}.fa

perl abstract.pl Guanine/Guanine.tbl >Guanine/Guanine.abstract.tsv

#查看Guanine同家族蛋白以及涉及到的数据库登录号
cat Guanine/Guanine.abstract.tsv | tsv-summarize -g 2,6  --count
TIGR02967       guan_deamin:_guanine_deaminase  3566
TIGR02022       hutF:_formiminoglutamate_deiminase      3544
TIGR01224       hutI:_imidazolonepropionase     2323
TIGR03314       Se_ssnA:_putative_selenium_metabolism_protein_SsnA      3566
TIGR03178       allantoinase:_allantoinase      31
TIGR00221       nagA:_N-acetylglucosamine-6-phosphate_deacetylase       5
TIGR02033       D-hydantoinase:_dihydropyrimidinase     3
TIGR03583       EF_0837:_putative_amidohydrolase,_EF_0837/AHA_3915_family       12
TIGR02318       phosphono_phnM:_phosphonate_metabolism_protein_PhnM     13

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare.pl Guanine/Guanine.abstract.tsv >Guanine/Guanine_minevalue.tsv
tsv-summarize  -g 3 --count Guanine/Guanine_minevalue.tsv
#查看需要过滤的基因
guan_deamin:_guanine_deaminase  3168
Se_ssnA:_putative_selenium_metabolism_protein_SsnA      161
hutF:_formiminoglutamate_deiminase      237
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

###与annotaion过滤的结果差别较大
可能是tigrfam数据库的原因
有可能未过滤TRZ_ATZ famliy hydrolase和9-oxoguanine deaminase

###最后解决方法：比对各个结果文件，发现需要筛选e值小于1e-90和score大于300的


```


### 3.2.3 将branched-chain提取的蛋白序列与tigerfam数据库比对
```BASH
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 branched-chain/branched.tigerfam.replace.tsv) branched-chain/branched-chain.fa

E_VALUE=1e-10
NAME=branched-chain
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.txt --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/TIGERFAM/tigrfams.hmm  ${NAME}/${NAME}.fa

perl abstract.pl branched-chain/branched-chain.tbl >branched-chain/branched-chain.abstract.tsv

#查看branched同家族蛋白以及涉及到的数据库登录号
cat branched-chain/branched-chain.abstract.tsv | tsv-summarize -g 2,6  --count
TIGR00796       livcs:_branched-chain_amino_acid_transport_system_II_carrier_protein    1779

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare.pl branched-chain/branched-chain.abstract.tsv >branched-chain/branched-chain_minevalue.tsv
tsv-summarize  -g 3 --count branched-chain/branched-chain_minevalue.tsv

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

###与annotaion过滤的结果差别不大

```

### 3.2.4 将MltB提取的蛋白序列与tigerfam数据库比对
```BASH
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 MltB/MltB_tigrfam.replace.tsv) MltB/MltB.fa

E_VALUE=1e-10
NAME=MltB
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.txt --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/TIGERFAM/tigrfams.hmm  ${NAME}/${NAME}.fa

perl abstract.pl MltB/MltB.tbl >MltB/MltB.abstract.tsv

#查看MltB同家族蛋白以及涉及到的数据库登录号
cat MltB/MltB.abstract.tsv | tsv-summarize -g 2,6  --count
TIGR02282       MltB:_lytic_murein_transglycosylase_B   3779
TIGR02283       MltB_2:_lytic_murein_transglycosylase   3779
TIGR02869       spore_SleB:_spore_cortex-lytic_enzyme   2

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare.pl MltB/MltB.abstract.tsv >MltB/MltB_minevalue.tsv
tsv-summarize  -g 3 --count MltB/MltB_minevalue.tsv
#结果
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

###与annotaion过滤的结果差别不大

```


### 3.2.5 Tyrosine提取的蛋白序列与tigerfam数据库比对
```BASH
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 Tyrosine/Tyrosine.replace.tsv) Tyrosine/Tyrosine.fa

E_VALUE=1e-10
NAME=Tyrosine
hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali \
-o ${NAME}/${NAME}_progress.txt --tblout ${NAME}/${NAME}.tbl  \
   ~/data/HMM/TIGERFAM/tigrfams.hmm  ${NAME}/${NAME}.fa

perl abstract.pl Tyrosine/Tyrosine.tbl >Tyrosine/Tyrosine.abstract.tsv

#查看MltB同家族蛋白以及涉及到的数据库登录号
cat Tyrosine/Tyrosine.abstract.tsv | tsv-summarize -g 2,6  --count
TIGR00234       tyrS:_tyrosine--tRNA_ligase     1912
TIGR00233       trpS:_tryptophan--tRNA_ligase   651

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl compare.pl Tyrosine/Tyrosine.abstract.tsv >Tyrosine/Tyrosine_minevalue.tsv
tsv-summarize  -g 3 --count Tyrosine/Tyrosine_minevalue.tsv
#结果
tyrS:_tyrosine--tRNA_ligase     1912

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

###与annotaion过滤的结果差别不大


```

## 3.3使用panther数据库
```BASH
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

# 将bac_Tyrosine提取的蛋白序列与pfam数据库比对
faops some PROTEINS/all.replace.fa <(tsv-select -f 2  bac_Tyrosine/bac_Tyrosine.replace.tsv) bac_Tyrosine/bac_Tyrosine.fa

#panther数据库太大，需要放在超算上跑

```





