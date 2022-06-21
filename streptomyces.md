<!-- TOC -->

- [1.下载链霉菌属streptomyces的基因组](#1下载链霉菌属streptomyces的基因组)
- [2.处理下载后的文件](#2处理下载后的文件)
- [3.下载所有hmm文件(ncbi的文件及信息)](#3下载所有hmm文件ncbi的文件及信息)
- [4.使用链霉菌属的基因组蛋白序列建库,使用上述hmm文件去比对](#4使用链霉菌属的基因组蛋白序列建库使用上述hmm文件去比对)
- [4.使用antismash预测代谢基因簇](#4使用antismash预测代谢基因簇)
- [安装big-scape](#安装big-scape)
- [5.ppanggolin分析泛基因组](#5ppanggolin分析泛基因组)

<!-- /TOC -->
# 1.下载链霉菌属streptomyces的基因组
```bash
cd /mnt/d
mkdir -p pan/strep
cd /mnt/d/pan/strep
#下载refer数据库所有的链霉菌的assembly的详细信息
GENUS=$(
    nwr member Streptomyces -r genus |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        grep -v " sp."|
        sed '1d' |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN ($GENUS)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > assembly_one.tsv
    wc -l  assembly_one.tsv #237
#下载refer数据库中所有的链霉菌的assembly
GENUS=$(
    nwr member Streptomyces -r genus |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        grep -v " sp." |
        sed '1d' |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND genus_id IN ($GENUS)
         AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite  >assembly_two.tsv #236
#给菌株添加GCF名字
cat assembly_two.tsv |
    grep -v '^#' |
    tsv-uniq |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "1,2,3" -s '\t' -min 5 --shortsub |
    (echo -e '#name\tftp_path\torganism\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[3]}++; # ftp_path
        $seen{$F[3]} > 1 and next;
        $seen{$F[5]}++; # abbr_name
        $seen{$F[5]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\n}, $F[5], $F[3], $F[1], $F[4];
        ' |
    keep-header -- sort -k3,3 -k1,1 \
    >assembly_three.tsv #237
sed -i 's/^S_/Str_/g' assembly_three.tsv
#准备下载材料
perl ~/Scripts/withncbi/taxon/assembly_prep.pl  -f assembly_three.tsv -o ASSEMBLY
bash ASSEMBLY/assembly_three.rsync.sh
#检查下载文件数
find ASSEMBLY -maxdepth 2  -name "*_genomic.fna.gz" | wc -l  #708
find ASSEMBLY -maxdepth 2  -name "*_genomic.gbff.gz" | wc -l   #236
find ASSEMBLY -maxdepth 2  -name "*_protein.faa.gz" | wc -l  #236
#或者md5检查
bash ASSEMBLY/assembly_three.collect.sh
cat ASSEMBLY/rsync.tsv |
    tsv-select -f 1 |
    parallel -j 4 --keep-order '
        echo "==> {}"
        cd ASSEMBLY/{}
        md5sum --check md5checksums.txt
    ' |
    grep -v ": OK"

```
# 2.处理下载后的文件
```bash
find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    sort |
    grep 'ASSEMBLY/' |
    wc -l
# 236

find ASSEMBLY -type f -name "*_protein.faa.gz" |
    wc -l
# 236

cat ASSEMBLY/rsync.tsv |cut -f 1 > strains.lst 
cat strains.lst  |  wc -l
# 236

mkdir -p PROTEINS
for STRAIN in $(cat strains.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz
done \
    > PROTEINS/all.pro.fa

cat PROTEINS/all.pro.fa |
    perl -nl -e '
        BEGIN { our %seen; our $h; }

        if (/^>/) {
            $h = (split(" ", $_))[0];
            $seen{$h}++;
            $_ = $h;
        }
        print if $seen{$h} == 1;
    ' \
    > PROTEINS/all.uniq.fa

# counting proteins
cat PROTEINS/all.pro.fa |
    grep "^>" |
    wc -l
#1658462

# annotations may be different
cat PROTEINS/all.uniq.fa |
    grep "^>" |
    wc -l
#1226527

for STRAIN in $(cat strains.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
        grep "^>" |
        cut -d" " -f 1 |
        sed "s/^>//" |
        STRAIN=${STRAIN} perl -nl -e '
            $n = $_;
            $s = $n;
            $s =~ s/\.\d+//;
            printf qq{%s\t%s_%s\t%s\n}, $n, $ENV{STRAIN}, $s, $ENV{STRAIN};
        ' \
    > PROTEINS/${STRAIN}.replace.tsv

    cut -f 2,3 PROTEINS/${STRAIN}.replace.tsv >> PROTEINS/all.strain.tsv

    faops replace -s ASSEMBLY/${STRAIN}/*_protein.faa.gz <(cut -f 1,2 PROTEINS/${STRAIN}.replace.tsv) stdout

    rm PROTEINS/${STRAIN}.replace.tsv
done \
    > PROTEINS/all.replace.fa

cat PROTEINS/all.replace.fa |
    grep "^>" |
    wc -l
#1658462
```

# 3.下载所有hmm文件(ncbi的文件及信息)
```bash
#256 aa    #AadA family aminoglycoside 3''-O-nucleotidyltransferase   #Domain cutoff 375.0
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF012157.0.HMM -O aadA.hmm
#397aa    #multidrug efflux RND transporter periplasmic adaptor subunit AcrA    #Domain cutoff 765.0   
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF011604.1.HMM -O acrA.hmm 
#1049aa  #efflux RND transporter permease AcrB       #Identity cutoff (%) 94  Target coverage cutoff (%) 90  
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF011706.2.HMM -O acrB.hmm 
#1036aa  #multidrug efflux RND transporter permease AcrD       #Domain cutoff 1854.0
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF007842.1.HMM -O acrD.hmm
#383aa   #ADC family extended-spectrum class C beta-lactamase   #Domain cutoff 800.0
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000425.2.HMM -O ADC.hmm
#465aa   #multidrug efflux RND transporter AdeABC outer membrane channel subunit AdeC   #domain cutoff 975
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF033142.1.HMM -O adeC.hmm
#480aa  #multidrug efflux RND transporter AdeIJK outer membrane channel subunit AdeK   #domain cutoff 900
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF033143.1.HMM -O adeK.hmm
#239aa #efflux system response regulator transcription factor AdeR   #domain cutoff 340
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF012227.0.HMM -O adeR.hmm
#353aa #two-component sensor histidine kinase AdeS  #domain cutoff 500
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF012226.0.HMM -O adeS.hmm
#256aa #aminoglycoside O-nucleotidyltransferase ANT(4')-Ib  #domain cutoff 500
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000079.1.HMM -O ANT(4)-Ib.hmm
#259aa    APH(3')-VI family aminoglycoside O-phosphotransferase  #domain cutoff 529
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF033062.0.HMM -O APH(3)-VI.hmm
#272aa aminoglycoside O-phosphotransferase APH(3'')-Ia  #domain cutoff 600
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF032894.1.HMM -O APH(3)-Ia.hmm
#264aa aminoglycoside O-phosphotransferase APH(3')-IIIa #domain cutoff 530
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF033064.0.HMM -O APH(3)-IIIa.hmm
#321aa bile acid:sodium symporter family transporter #domain cutoff 326
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR00841.1.HMM -O basS.hmm
#281aa penicillin-hydrolyzing class A beta-lactamase BlaZ #domain cutoff 580 
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF033139.1.HMM -O blaZ.hmm
#291aa CTX-M family extended-spectrum class A beta-lactamase #domain cutoff 570
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF033089.1.HMM -O CTX-M.hmm
#287aa GES family class A beta-lactamase #domain cutoff 635
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF012103.0.HMM -O GES.hmm
#292aa IMI family carbapenem-hydrolyzing class A beta-lactamase     #domain cutoff  625
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000400.2.HMM -O IMI.hmm
#245aa IMP family subclass B1 metallo-beta-lactamase     #domain cutoff   460
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF012147.1.HMM -O IMP.hmm
#293aa KPC family carbapenem-hydrolyzing class A beta-lactamase   #domain cutoff   625
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF012141.0.HMM -O KPC.hmm
#541aa MCR-1 family phosphoethanolamine--lipid A transferase    #domain cutoff    1100
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000465.1.HMM -O MCR-1.hmm
#540aa  MCR-3 family phosphoethanolamine--lipid A transferase #domain cutoff 1225
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF033408.5.HMM -O MCR-3.hmm
#539aa MCR-9 family phosphoethanolamine--lipid A transferase   #domain cutoff 1240
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF033836.1.HMM -O MCR-9.hmm
#668aa PBP2a family beta-lactam-resistant peptidoglycan transpeptidase MecA  #domain cutoff
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000409.1.HMM -O mecA.hmm
#123aa mecA-type methicillin resistance repressor MecI  #domain cutoff 250
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000243.1.HMM  -O mceI.hmm
#451aa multidrug efflux MATE transporter MepA   850
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000131.1.HMM -O mepA.hmm
#153aa peptide-methionine (S)-S-oxide reductase MsrA 108
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR00401.1.HMM -O msrA.hmm
#270aa  NDM family subclass B1 metallo-beta-lactamase    #domain cutoff 580
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000259.2.HMM -O NDM.hmm
#391aa multidrug efflux RND transporter periplasmic adaptor subunit OqxA    domain cutoff 875
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000272.1.HMM -O oqxA.hmm
#1050aa multidrug efflux RND transporter permease subunit OqxB   domain cutoff 2200
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000037.1.HMM -O oqxB.hmm
#276aa OXA-1 family class D beta-lactamase   domain cutoff 580
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000388.2.HMM -O  OXA-1.hmm
#275aa OXA-2 family class D beta-lactamase domain cutoff 590
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000267.2.HMM -O   OXA-2.hmm
#275aa OXA-10 family class D beta-lactamase   domain cutoff 580
wget -c  https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000386.2.HMM  -O OXA-10.hmm
#273aa OXA-23 family carbapenem-hydrolyzing class D beta-lactamase   domain cutoff 600
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000266.2.HMM -O OXA-23.hmm
#265aa OXA-48 family class D beta-lactamase   domain cutoff 580
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000387.2.HMM -O OXA-48.hmm
#262aa OXA-50 family oxacillin-hydrolyzing class D beta-lactamase domain cutoff 580
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000432.2.HMM -O OXA-50.hmm
#274aa OXA-51 family carbapenem-hydrolyzing class D beta-lactamase domain cutoff 600 
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000268.2.HMM -O OXA-51.hmm
#397aa PDC family class C beta-lactamase  #domain cutoff 500 
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000422.6.HMM -O PDC.hmm
# 308aa PER family extended-spectrum class A beta-lactamase domain cutoff 650
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000389.2.HMM -O PER.hmm
# 216aa two-component system response regulator PmrA domain cutoff 330
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF007928.2.HMM -O pmrA.hmm
# 514aa QacA/B family quaternary ammonium compound efflux MFS transporter  doamin cutoff 900
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000089.1.HMM -O qacA.hmm
#218aa QnrA family quinolone resistance pentapeptide repeat protein domain cutoff 450
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000071.4.HMM -O qnrA.hmm
#214aa QnrB family quinolone resistance pentapeptide repeat protein domain cutoff 460
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000420.1.HMM -O qnrB.hmm
#214AA QnrD family quinolone resistance pentapeptide repeat protein  doamin cutoff 425
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000139.2.HMM -O qnrD.hmm
#218aa QnrS family quinolone resistance pentapeptide repeat protein    domain cutoff 460
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000056.3.HMM -O qnrS.hmm
#285aa SHV family class A beta-lactamase domain cutoff 640
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000285.3.HMM -O SHV.hmm
#278aa  sulfonamide-resistant dihydropteroate synthase Sul1   doamin cutoff 600
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000294.1.HMM -O sul1.hmm
#271aa sulfonamide-resistant dihydropteroate synthase Sul2  domain cutoff 575
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000295.1.HMM -O sul2.hmm
#263aa sulfonamide-resistant dihydropteroate synthase Sul3 domain cutoff 475
wget -c  https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000296.1.HMM -O sul3.hmm
#286aa TEM family class A beta-lactamase domain cutoff 625
wget n-c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000531.2.HMM -O TEM.hmm
# 420aa tetracycline efflux MFS transporter TetA(P) domain cutoff 750
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000043.1.HMM -O tetA(P).hmm
#513aa tetracycline efflux ABC transporter TetAB subunit A domain cut off 10000
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000506.1.HMM -O tetA.hmm
#197aa tetracyline resistance-associated transcriptional repressor TetC domain cutoff 425
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000040.3.HMM -O tetC.hmm
#342aa D-alanine--(R)-lactate ligase domain cutoff 550
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000206.1.HMM -O vanA.hmm
# 322aa VanH-A/VanH-Pt family D-lactate dehydrogenase  domain cutoff 700
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000371.1.HMM -O vanHA.hmm
#231aa VanA-type vancomycin resistance DNA-binding response regulator VanR   domain cutoff 500
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF000401.1.HMM -O vanA.hmm
#293aa vancomycin resistance histidine kinase VanS   domain cutoff 400 
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF033091.0.HMM -O vanS.hmm
#199aa D-Ala-D-Ala dipeptidase VanX domain cutoff 320
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF033115.1.HMM -O vanX.hmm
#164aa   glycopeptide resistance protein VanZ-A   domain cutoff 250
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF033125.2.HMM -O vanZ.hmm
#266AA VIM family subclass B1 metallo-beta-lactamase    domain cutoff 500 
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/NF012100.0.HMM -O VIM.hmm
```

# 4.使用链霉菌属的基因组蛋白序列建库,使用上述hmm文件去比对
```bash
cd /mnt/d/pan/strep
#使用hmmsearch查找hmm文件对应的蛋白序列
for filename in hmm/*.hmm
do
base=$(basename $filename .hmm)
echo $base
hmmsearch --noali --notextw -E 1e-5 --domE 1e-5  hmm/$base.hmm   PROTEINS/all.replace.fa  >hmm_result/$base.hmmer.result.tsv
done

for filename in hmm_result/*.tsv
do
base=$(basename $filename .tsv)
echo $base
cat hmm_result/$base.tsv | grep '>>' | perl -nl -e '/>>\s+(\S+)/ and print $1' >hmm_result/$base.txt
done

#提取hmmsearch对应的序列
for filename in hmm_result/*.txt
do
base=$(basename $filename .txt)
echo $base
faops some  PROTEINS/all.replace.fa hmm_result/$base.txt  hmm_result/$base.fa
done

#查看上述查找序列对应的domain 
for filename in hmm_result/*.fa
do
base=$(basename $filename .fa)
echo $base
hmmscan --cpu 4 -E 1e-5  --domE 1e-5  --noali -o hmm_result/$base.pfam.txt  --tblout hmm_result/$base.pfam.tbl   ~/data/HMM/PFAM/Pfam-A.hmm  hmm_result/$base.fa
done 

for filename in hmm_result/*.fa
do
base=$(basename $filename .fa)
echo $base
perl  ~/data/Pseudomonas/abstract1.pl hmm_result/$base.pfam.tbl   >hmm_result/$base.pfam.json
done 
```

# 4.使用antismash预测代谢基因簇
* 次级代谢产物包括聚酮PK(type ⅠⅡⅢ),非核糖体肽，萜烯，羊毛硫抗生素，细菌素，β-内酰胺，氨基糖苷类，铁载体，丁内酯，吲哚类，黑素类，磷酸糖脂，四氢嘧啶等
* 其中聚酮PKS和非核糖体肽NRPS是目前研究最多的两类次级代谢产物，PKS和NRPS复合酶的基因主要是由连续的模块构成，每个模块具有它各自的功能，含有对应的结构域，一个基因可能会由多个模块（module）。一种次级代谢产物的合成由多个基因共同控制
* PKS的主要结构是酰基转移酶（AT）、酮基合酶（KS）、酰基载体蛋白（ACP），另外可能还会有一些非必需的结构域，1~3个修饰酮基的结构域，酮基还原酶（KR）、脱水酶（DH）、烯酰基还原酶（ER）。NPRS的主要结构域是由腺苷酰化结构域（A）、缩合结构域（C）、肽酰载体蛋白结构域（PCP/T），还会有些非必需的结构域：差向异构化（E）、N甲基化、氧化等修饰结构域。
* antiSMASH的流程是（gbk，embl，fasta）文件输入---基因预测（glimmer3）----基因簇确定(标志基因的hmm)----得到输出文件（domain分析，化学结构预测，smCOG分析,clusterblasst)
* COG和eggnog都是蛋白同源簇数据库，但是eggnog包含COG,KOG,NOG,KEGG,GO,SMART,PFAM等信息.
```bash
#安装antismash所需要的依赖
conda create -n antismash
conda activate antismash
conda install hmmer2 hmmer diamond fasttree prodigal blast muscle glimmerhmm
conda install meme==4.11.2
conda activate antismash
#安装antismash
conda activate antismash
wget https://dl.secondarymetabolites.org/releases/6.0.0/antismash-6.0.0.tar.gz
tar -zxf antismash-6.0.0.tar.gz
pip3 install ./antismash-6.0.0
#报错 ImportError: cannot import name 'Markup' from 'jinja2'
#解决  pip install jinjia2   vi /home/syq/miniconda3/envs/antismash/lib/python3.8/site-packages/jinja2/__init__.py     from markupsafe import Markup
#
#查看是否安装成功
antismash -h
#下载数据库(包括pfam  resfams tigrfam   clusterblast   mibig 共五个数据库)
download-antismash-databases
```
* antismash支持.fasta .gbk .embl作为格式文件输入，推荐以genbank文件输入，该文件包含编码蛋白基因的注释信息，不需要调用glimmer3和glimmerhmm来进行基因预测后再进行次级代谢基因簇的鉴定
* --fullhmmer 运行全基因组hmmer分析,--cassis使用cassis算法预测基因簇边界,--cf-borders-only仅注释现有簇,--cf-create-clusters寻找额外簇,--clusterhmmer运行簇的hmmer分析,--smcog-trees寻找簇的直系同源群,--tta-threshold运行tta密码子检测模块,--cb-genrnal：将预测的簇与antiSMASH现有的簇进行比较,--cb-subclusters将已经识别的基因簇和已知负责合成前体物质的子簇进行比对,--cb-knowclusters将已经识别的基因簇和mibig数据库中的已知基因簇进行比对，--asf运行活性中心检测模块--pfam2go运行pfam模块
```bash
cd /mnt/d/pan/strep
mkdir -p gbff
#移动文件
for name in $(cat strains.lst)
do
echo $name
cp ASSEMBLY/$name/*_genomic.gbff.gz ASSEMBLY/$name/$name.gbff.gz
done

for name in $(cat strains.lst)
do
echo $name
mv ASSEMBLY/$name/$name.gbff.gz ./gbff
done
#find ASSEMBLY -maxdepth 2  -name "*_genomic.gbff.gz"  | parallel --no-run-if-empty --linebuffer -k -j 4 'cp {}  ./gbff'
#批量解压gz
gunzip  gbff/*.gz
#antismash注释(信息最全面)
#cat strains.lst | parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 'antismash --cb-general   gbff/{1}.gbff -c 6 '

for name in $(cat strains.lst)
do
echo $name
antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go  gbff/$name.gbff -c 8
done
```
* 结果解析
* NRPS非核糖体多肽，betalactone β-丁内酯,lanthipeptide羊毛硫肽,terpene萜烯,butyrolactone γ-丁内酯,siderophore铁载体

* 结果可视化软件  CORASON和BIG-SCAPE

# 安装big-scape
```bash
conda install numpy scipy scikit-learn hmmer biopython fasttree networkx
cd  ~
wget https://github.com/medema-group/BiG-SCAPE/archive/refs/heads/master.zip
unzip BiG-SCAPE-master.zip
mv BiG-SCAPE-master bigscape


```




# 5.ppanggolin分析泛基因组
```bash
#提取genbank文件
cat strains.taxon.tsv | grep "Pseudom_aeru" >Pse_aeru_strain_GCF.tsv
wc -l Pse_aeru_strain_GCF.tsv

for f in $(cat Pse_aeru_strain_GCF.tsv)
do
cp  ASSEMBLY/$f/*.gbff.gz ASSEMBLY/$f/$f.gbff.gz
done

for f in $(cat Pse_aeru_strain_GCF.tsv)
do
mv  ASSEMBLY/$f/$f.gbff.gz /mnt/d/pan
done

#将文件名和路径合并为一个列表，方便泛基因组分析
for filename in *.gbff
do
base=$(basename $filename .gbff)
echo -e "$base\t$base.gbff" >>pan_name_dir.tsv
done 

#进行pangeomone分析，默认覆盖度为0.8，默认相似性为0.8
ppanggolin workflow --anno ./pan_name_dir.tsv --cpu 8 -o ./result


```
