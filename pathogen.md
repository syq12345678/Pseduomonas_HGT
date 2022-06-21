
# 1.下载常见病原菌的基因组

# 1.1下载铜绿假单胞菌的基因组
```bash
cd /mnt/d
mkdir -p pan/pathogen/pseudomonas
cd /mnt/d/pan/pathogen/pseudomonas
#下载refer数据库所有的链霉菌的assembly的详细信息
GENUS=$(
    nwr member Pseudomonas -r genus |
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
    wc -l  assembly_one.tsv #1000
#下载refer数据库中所有的链霉菌的assembly
GENUS=$(
    nwr member Pseudomonas -r genus |
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
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite  >assembly_two.tsv #999
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
    >assembly_three.tsv #1000
sed -i 's/^P_/Pse_/g' assembly_three.tsv
#准备下载材料
perl ~/Scripts/withncbi/taxon/assembly_prep.pl  -f assembly_three.tsv -o ASSEMBLY
bash ASSEMBLY/assembly_three.rsync.sh
```

# 2.处理下载后的文件
```bash
find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    sort |
    grep 'ASSEMBLY/' |
    wc -l
# 391

find ASSEMBLY -type f -name "*_protein.faa.gz" |
    wc -l
# 391

cat ASSEMBLY/rsync.tsv |cut -f 1 > strains.lst 
cat strains.lst  |  wc -l
# 391

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
#2371874

# annotations may be different
cat PROTEINS/all.uniq.fa |
    grep "^>" |
    wc -l
#225191

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
#2371874
```

# 3.使用铜绿假单胞菌的基因组蛋白序列建库,使用上述hmm文件去比对
```bash
cd /mnt/d/pan/pathogen/pseudomonas
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
faops some  PROTEINS/all.replace.fa hmm_result/$base.txt hmm_result/$base.fa
done

#使用hmmscan查看序列对应的domain
for filename in hmm_result/*.fa
do
base=$(basename $filename .fa)
echo $base
hmmscan --cpu 4 -E 1e-5 --domE 1e-5 --noali -o hmm_result/$base.pfam.txt --tblout hmm_result/$base.pfam.tbl ~/data/HMM/PFAM/Pfam-A.hmm  hmm_result/$base.fa
done



```