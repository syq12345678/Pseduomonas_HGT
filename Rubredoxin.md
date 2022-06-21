* [hmm网址](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/NF012522/)
# 1.使用hmmsearch
```bash
#查看branched-chain的hmm文件所在数据库
cd ~/data/Pseudomonas
cat STRAINS/Pseudom_aer_PAO1/*.tsv |
    grep "IPR024922"
#可以发现数据库中的登录号是
PIRSF  PIRSF000071
#在PIRSF网站中输入PIRSF000071，发现对应的Pfam id PF00301
#下载hmm文件 (47aa)
mkdir -p Rubr/HMM
curl -L https://pfam.xfam.org/family/PF00301/hmm >Rubr/HMM/Rubr.pfam.hmm
#使用Hmm文件抓取对应的序列
hmmsearch -E 1e-5 --domE 1e-5 --noali --notextw Rubr/HMM/Rubr.pfam.hmm PROTEINS/all.replace.fa  | grep ">>" | perl -nl -e 'm{>>\s+(\S+)} or next;$n = $1;print"$n";' >Rubr/Rubr.pfam.replace.tsv
#查看序列数量 2235
wc -l Rubr/Rubr.pfam.replace.tsv
#查看抓取的序列长度分布范围
faops some  PROTEINS/all.replace.fa Rubr/Rubr.pfam.replace.tsv Rubr/Rubr.pfam.replace.fa
faops size Rubr/Rubr.pfam.replace.fa | tsv-summarize -g 2 --count
faops size Rubr/Rubr.pfam.replace.fa | cut -f 2 | sort -n | uniq
#最短48，最长858
```

# 2.使用hmmscan
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

#查看抓取出来的序列对应的domain
hmmscan --cpu 4 -E 1e-5 --domE 1e-5 --noali -o Rubr/Rubr_progress.pfam.txt --tblout Rubr/Rubr.pfam.tbl  ~/data/HMM/PFAM/Pfam-A.hmm  Rubr/Rubr.pfam.replace.fa
#处理结果
perl abstract1.pl Rubr/Rubr.pfam.tbl >Rubr/Rubr.abstract.pfam.tsv #2927
tsv-summarize  -g 8 --count Rubr/Rubr.abstract.pfam.tsv 
#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
cut -f 3 Rubr/Rubr.abstract.pfam.tsv | sort -n | uniq | wc -l  #2235
#仅保留e值最小的匹配项
perl compare1.pl  Rubr/Rubr.abstract.pfam.tsv  >  Rubr/Rubr.minevalue.pfam.tsv  #2235
tsv-summarize  -g 3 --count Rubr/Rubr.minevalue.pfam.tsv
#保留description符合项
cat Rubr/Rubr.minevalue.pfam.tsv |tsv-filter --str-in-fld 3:"Rubredoxin"    #2186
cat Rubr/Rubr.abstract.pfam.tsv | tsv-select -f 2,3,5,8 | grep -f <(cat Rubr/Rubr.minevalue.pfam.tsv |tsv-filter --str-in-fld 3:"Rubredoxin" ) |tsv-summarize -g 1 --count
#PF00301.20      2186
```