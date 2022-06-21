
# 1.处理结果文件
```bash
#合并结果
cd /mnt/d/WGCNA/rma_WGCNA/condition_thhree_combine/result
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
echo $base
braB=$(cat select/$base.select.txt| grep "PA3046"| cut -f 3)
braZ=$(cat select/$base.select.txt | grep "PA1842"| cut -f 3)
cat select/$base.select.txt| grep "$braB" >> PA3046.select.txt
cat select/$base.select.txt | grep "$braZ" >> PA1842.select.txt
done

#准备频次密度统计图的绘图文件
cd /mnt/d/WGCNA/rma_WGCNA/condition_three_combine/result
cat PA3046.select.txt | grep -v "grey" | tsv-summarize  -g 1 --count >yggl1_whole_number.tsv
cat PA1842.select.txt | grep -v "grey" | tsv-summarize  -g 1 --count >yggl2_whole_number.tsv
cut -f 2 yggl1_whole_number.tsv | sed  '1iratio' >yggl1_density_ggplot.tsv
cut -f 2 yggl2_whole_number.tsv | sed '1iratio' >yggl2_density_ggplot.tsv
```


# 2.绘制正态密度分布直方图
```bash
#绘制yggl1的正态密度分布直方图
setwd("D:/WGCNA/rma_WGCNA/condition_three_combine/result")
library(ggplot2)
library(showtext)#用showtext包输入需要的字体
font_add('Arial','/Library/Fonts/Arial.ttf')
showtext_auto()
library(patchwork)
ratio<-read.table("yggl1_density_ggplot.tsv",header=T)
set.seed(1000)
df <- data.frame(ratio)
data_mean=mean(df$ratio) 
data_mean  #24980.79
data_sd=sd(df$ratio) 
data_sd   # 2267.087
two_sd_plus=mean(df$ratio)+data_sd*2
two_sd_plus    #29514.96

p1<-ggplot(df, aes(x = ratio)) +geom_histogram(aes(y =..density..),breaks = seq(min(df$ratio), max(df$ratio), length =100),colour ="#C6DBEF", fill ="#08519C", size = 0.1,alpha=0.8) + #使用density代替y轴
stat_function(fun = dnorm, args = list(mean = data_mean, sd = data_sd,color ="black", size = 1))+
geom_density(color = "black", size = 1)+  #直方图上加密度曲线，也可以用geom_density()
geom_vline(aes(xintercept=data_mean), color="red", linetype="dashed", size=1)+ #添加均值线和两个标准差线
geom_vline(aes(xintercept=two_sd_plus), color="red", linetype="dashed", size=1)+
annotate(geom="text",fontface="bold",family="Arial",color="black",x=32000,y=0.00015,label="29514.96",size=9)+   
 #通过annotate函数添加注释，将geom变量设置为"text"和"rect"分别代表文字与矩形，并且可以调整位置，大小颜色
annotate(geom="text",fontface="bold",family="Arial",color="black",x=10000,y=0.000200,label="μ+3σ=24980.79+4534.17",size=9)+
labs(x="frequency",y = "density")+
theme_bw()+
theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
theme(axis.text=element_text(size=24,family="Arial",face = "bold"),axis.title=element_text(size=24,family="Arial",face="bold"))
ggsave(file="yggl1_density_picture.pdf",plot=p1,width=10,height=8)


#绘制yggl2的正态密度分布直方图
setwd("D:/WGCNA/rma_WGCNA/condition_three_combine/result")
library(ggplot2)
library(showtext)#用showtext包输入需要的字体
options(scipen = 200)
font_add('Arial','/Library/Fonts/Arial.ttf')
showtext_auto()
library(patchwork)
ratio1<-read.table("yggl2_density_ggplot.tsv",header=T)
set.seed(1000)
df1 <- data.frame(ratio1)
data_mean1=mean(df1$ratio) 
data_mean1   #24118.65
data_sd1=sd(df1$ratio) 
data_sd1    #2019.258
two_sd_plus1=mean(df1$ratio)+data_sd1*2
two_sd_plus1  # 28157.16

p2<-ggplot(df1, aes(x = ratio)) +geom_histogram(aes(y =..density..),breaks = seq(min(df1$ratio), max(df1$ratio), length =100),colour ="#C6DBEF", fill ="#08519C", size = 0.1,alpha=0.8) + #使用density代替y轴
stat_function(fun = dnorm, args = list(mean = data_mean1, sd = data_sd1,color ="black", size = 1))+  #直方图上加密度曲线，也可以用geom_density()
geom_vline(aes(xintercept=data_mean1), color="red", linetype="dashed", size=1)+ #添加均值线和两个标准差线
geom_density(color = "black", size = 1)+  #直方图上加密度曲线，也可以用geom_density()
geom_vline(aes(xintercept=two_sd_plus1), color="red", linetype="dashed", size=1)+
annotate(geom="text",fontface="bold",family="Arial",color="black",x=32000,y=0.00015,label="28157.16",size=9)+   
#通过annotate函数添加注释，将geom变量设置为"text"和"rect"分别代表文字与矩形，并且可以调整位置，大小颜色
annotate(geom="text",fontface="bold",family="Arial",color="black",x=10000,y=0.000200,label="μ+3σ=24118.65+4038.52",size=9)+
labs(x="frequency",y = "density")+
theme_bw()+
theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
theme(axis.text=element_text(size=24,family="Arial",face = "bold"),axis.title=element_text(size=24,family="Arial",face="bold"))
ggsave(file="yggl2_density_picture.pdf",plot=p2,width=10,height=8)
```

# 3.和yggl1，yggl2共表达的基因
```bash
# three_sigma yggl1 31782.05  yggl2 30176.42   yggl1是4个  yggl2是8个
# two_sigma  yggl1 29514.96  yggl2 28157.16  yggl1是123个  yggl2是135个
cd /mnt/d/WGCNA/rma_WGCNA/condition_three_combine/result

cat PA3046.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | tsv-filter --ge 2:29515 | cut -f 1 >yggl1_two_sigma_neighbor_high_fre.tsv
cat PA1842.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | tsv-filter --ge 2:28158 | cut -f 1 >yggl2_two_sigma_neighbor_high_fre.tsv

cat PA3046.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | tsv-filter --ge 2:31783 | cut -f 1 >yggl1_three_sigma_neighbor_high_fre.tsv
cat PA1842.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | tsv-filter --ge 2:30177 | cut -f 1 >yggl2_three_sigma_neighbor_high_fre.tsv

mv yggl1_three_sigma_neighbor_high_fre.tsv yggl2_three_sigma_neighbor_high_fre.tsv  yggl1_two_sigma_neighbor_high_fre.tsv  yggl2_two_sigma_neighbor_high_fre.tsv /mnt/d/WGCNA/venn

#将yggl1_neighber_high_fre.tsv的gene id转换locus_tag id,使用david
cd /mnt/d/WGCNA/venn
wget https://david.ncifcrf.gov/data/download/conv_28132519CBFB1655275871960.txt -O yggl1_two_sigma_PA.tsv #120,有3个id未知
wget https://david.ncifcrf.gov/data/download/conv_0CB67D61B0231655276772607.txt -O yggl2_two_sigma_PA.tsv #134，有1个id未知
cut -f 1 yggl1_two_sigma_PA.tsv | grep "PA" | perl -alne '$_=~/PA(.*?)\_/;$name=$1;print"PA$name"' >yggl1_two_sigma_PA_id.tsv
cut -f 1 yggl2_two_sigma_PA.tsv | grep "PA" | perl -alne '$_=~/PA(.*?)\_/;$name=$1;print"PA$name"' >yggl2_two_sigma_PA_id.tsv
cp yggl1_two_sigma_PA_id.tsv yggl1.tsv
cp yggl2_two_sigma_PA_id.tsv yggl2.tsv
#绘制mean+2个sigma的韦恩图
cd /mnt/d/WGCNA/venn
plotr venn yggl1.tsv yggl2.tsv
```

# 1.取PF00466 Ribosomal_L10(rplJ)(PA4272)和PF00410 Ribosomal_S8(rpsH)（PA4249)共表达基因的结果
```bash
#合并结果
cd /mnt/d/WGCNA/rma_WGCNA/condition_three_combine/result
for filename in select/*.select.txt
do
base=$(basename $filename .select.txt)
echo $base
braB=$(cat select/$base.select.txt| grep "PA4272)"| cut -f 3)
braZ=$(cat select/$base.select.txt | grep "PA4249"| cut -f 3)
cat select/$base.select.txt| grep "$braB" >> PA4272.select.txt
cat select/$base.select.txt | grep "$braZ" >> PA4249.select.txt
done

#准备频次密度统计图的绘图文件
cd /mnt/d/WGCNA/rma_WGCNA/condition_three_combine/result
cat PA4272.select.txt | grep -v "grey" | tsv-summarize  -g 1 --count >L10_whole_number.tsv
cat PA4249.select.txt | grep -v "grey" | tsv-summarize  -g 1 --count >S8_whole_number.tsv
cut -f 2 L10_whole_number.tsv | sed  '1iratio' >L10_density_ggplot.tsv
cut -f 2 S8_whole_number.tsv | sed '1iratio' >S8_density_ggplot.tsv
```


# 2.绘制正态密度分布直方图
```bash
#绘制L10的正态密度分布直方图
setwd("D:/WGCNA/rma_WGCNA/condition_three_combine/result")
library(ggplot2)
library(showtext)#用showtext包输入需要的字体
font_add('Arial','/Library/Fonts/Arial.ttf')
showtext_auto()
library(patchwork)
ratio<-read.table("L10_density_ggplot.tsv",header=T)
set.seed(1000)
df <- data.frame(ratio)
data_mean=mean(df$ratio) 
data_mean  #27903.24
data_sd=sd(df$ratio) 
data_sd   #2228.299
two_sd_plus=mean(df$ratio)+data_sd*2
two_sd_plus    # 32359.84

p1<-ggplot(df, aes(x = ratio)) +geom_histogram(aes(y =..density..),breaks = seq(min(df$ratio), max(df$ratio), length =100),colour ="#C6DBEF", fill ="#08519C", size = 0.1,alpha=0.8) + #使用density代替y轴
stat_function(fun = dnorm, args = list(mean = data_mean, sd = data_sd,color ="black", size = 1))+
geom_density(color = "black", size = 1)+  #直方图上加密度曲线，也可以用geom_density()
geom_vline(aes(xintercept=data_mean), color="red", linetype="dashed", size=1)+ #添加均值线和两个标准差线
geom_vline(aes(xintercept=two_sd_plus), color="red", linetype="dashed", size=1)+
annotate(geom="text",fontface="bold",family="Arial",color="black",x=35000,y=0.00015,label="32359.84",size=9)+   
 #通过annotate函数添加注释，将geom变量设置为"text"和"rect"分别代表文字与矩形，并且可以调整位置，大小颜色
annotate(geom="text",fontface="bold",family="Arial",color="black",x=12000,y=0.000200,label="μ+3σ=27903.24+4456.60",size=9)+
labs(x="frequency",y = "density")+
theme_bw()+
theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
theme(axis.text=element_text(size=24,family="Arial",face = "bold"),axis.title=element_text(size=24,family="Arial",face="bold"))
ggsave(file="L10_two_sigma.pdf",plot=p1,width=10,height=8)


#绘制S8的正态密度分布直方图
setwd("D:/WGCNA/rma_WGCNA/condition_three_combine/result")
library(ggplot2)
library(showtext)#用showtext包输入需要的字体
options(scipen = 200)
font_add('Arial','/Library/Fonts/Arial.ttf')
showtext_auto()
library(patchwork)
ratio1<-read.table("S8_density_ggplot.tsv",header=T)
set.seed(1000)
df1 <- data.frame(ratio1)
data_mean1=mean(df1$ratio) 
data_mean1   #25162.64
data_sd1=sd(df1$ratio) 
data_sd1    #2247.653
two_sd_plus1=mean(df1$ratio)+data_sd1*2
two_sd_plus1  #29657.95

p2<-ggplot(df1, aes(x = ratio)) +geom_histogram(aes(y =..density..),breaks = seq(min(df1$ratio), max(df1$ratio), length =100),colour ="#C6DBEF", fill ="#08519C", size = 0.1,alpha=0.8) + #使用density代替y轴
stat_function(fun = dnorm, args = list(mean = data_mean1, sd = data_sd1,color ="black", size = 1))+  #直方图上加密度曲线，也可以用geom_density()
geom_vline(aes(xintercept=data_mean1), color="red", linetype="dashed", size=1)+ #添加均值线和两个标准差线
geom_density(color = "black", size = 1)+  #直方图上加密度曲线，也可以用geom_density()
geom_vline(aes(xintercept=two_sd_plus1), color="red", linetype="dashed", size=1)+
annotate(geom="text",fontface="bold",family="Arial",color="black",x=32000,y=0.00015,label="29657.95",size=9)+   
#通过annotate函数添加注释，将geom变量设置为"text"和"rect"分别代表文字与矩形，并且可以调整位置，大小颜色
annotate(geom="text",fontface="bold",family="Arial",color="black",x=10000,y=0.000200,label="μ+3σ=25162.64+4495.31",size=9)+
labs(x="frequency",y = "density")+
theme_bw()+
theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
theme(axis.text=element_text(size=24,family="Arial",face = "bold"),axis.title=element_text(size=24,family="Arial",face="bold"))
ggsave(file="S8_two_sigma.pdf",plot=p2,width=10,height=8)
```

# 3.和yggl1，yggl2共表达的基因
```bash
# two_sigma  L10 32359.84  L10是156个   S8 29657.95  S8是123个
cd /mnt/d/WGCNA/rma_WGCNA/condition_two_combine/result
cat PA4272.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | tsv-filter --ge 2:32360 | cut -f 1 >L10_two_sigma_neighbor_high_fre.tsv
cat PA4249.select.txt | grep -v "grey" | cut -f 1 | tsv-summarize -g 1 --count | tsv-filter --ge 2:29658 | cut -f 1 >S8_two_sigma_neighbor_high_fre.tsv
mv L10_two_sigma_neighbor_high_fre.tsv  S8_two_sigma_neighbor_high_fre.tsv /mnt/d/WGCNA/venn

#将yggl1_neighber_high_fre.tsv的gene id转换locus_tag id,使用david
cd /mnt/d/WGCNA/venn
wget https://david.ncifcrf.gov/data/download/conv_472BDA4A97A91655374131218.txt -O S8_two_sigma_PA.tsv #120  有3个id未知
wget https://david.ncifcrf.gov/data/download/conv_472BDA4A97A91655374152039.txt -O L10_two_sigma_PA.tsv  #153 #有3个id未知

cut -f 1 S8_two_sigma_PA.tsv | grep "PA" | perl -alne '$_=~/PA(.*?)\_/;$name=$1;print"PA$name"' >S8_two_sigma_PA_id.tsv
cut -f 1 L10_two_sigma_PA.tsv | grep "PA" | perl -alne '$_=~/PA(.*?)\_/;$name=$1;print"PA$name"' >L10_two_sigma_PA_id.tsv
cp S8_two_sigma_PA_id.tsv  S8.tsv
cp L10_two_sigma_PA_id.tsv   L10.tsv
#绘制mean+2个sigma的韦恩图
cd /mnt/d/WGCNA/venn
plotr venn S8.tsv  L10.tsv
```