#packages----
.libPaths("/home/user02/R/x86_64-pc-linux-gnu-library/4.3")
library(ggplot2)
library(reshape2)
library(ggpubr)
library(patchwork)
library(IOBR)
library(ggVennDiagram)
library(vegan)
library(gridExtra)
library(ALDEx2)
library(phyloseq)
library(pheatmap)
library(tidyverse)
library(edgeR)
library(DESeq2)
library(Hmisc)
library(clusterProfiler)
library(enrichplot)
library(ggrepel)


wjscptlwkr<-wjscptlwkrdmfqkdwlgkqslek #실수로 전체 시작 눌렀을때 오류로 중지시키기

library(ttr) #for geom_split_violin()

used.data<-'/data2/msi/Rdata/EMT_bulk.RData'

#data processing----
patient.info<-read.table('/home/user02/srr_file/SraRunTable_GSE179640.txt',sep = ',',header = T)

buf<-readxl::read_xlsx('/data2/msi/endometriosis/endometriosis/41556_2022_961_MOESM2_ESM.xlsx')
buf1<-patient.info[c('Run','source_name','condition','sample_location','tissue','Method')]
buf2<-merge(buf1,buf,by.x='source_name',by.y='PID',all.x=T)

patient.info.bulk<-buf2[buf2$Method=='RNA-seq',]
patient.info.sc<-buf2[buf2$Method!='RNA-seq',]

total.data<-patient.info.bulk[c('Run','source_name','condition','sample_location','tissue','Age','Ethnicity','Pregnancy','rASRM Stage','Oral contraceptive treatment','Endometrium histology')]
total.data$tissue<-ifelse(total.data$condition=='Control','Control',total.data$tissue)
s.data<-data.frame()

for (i in rownames(total.data)){
  buf<-paste0('/data2/msi/kraken2/pluspf/bracken_',total.data[i,'Run'],'.report.bracken')
  buf<-read.table(buf,sep = '\t',header=T)
  buf<-buf[c('name','new_est_reads')]
  colnames(buf)<-c('name',total.data[i,'Run'])
  if (nrow(s.data)==0){
    s.data<-buf
  }else{
    s.data<-merge(s.data,buf,by = 'name',all = T)
  }
} #merge bracken data
s.data<-s.data[s.data$name!='Homo sapiens',]

s.data$name<-unlist(lapply(s.data$name,FUN = function(x){gsub(' ','_',x)}))
rownames(s.data)<-s.data$name
s.data<-s.data[,-1]
s.data<-data.frame(t(s.data))

s.data[startsWith(colnames(s.data),'Bacteria_s_Fuso')]

g.data$name<-unlist(lapply(g.data$name,FUN = function(x){gsub(' ','_',x)}))
rownames(g.data)<-g.data$name
g.data<-g.data[,-1]
g.data<-data.frame(t(g.data))


g.data[is.na(g.data)]<-0
s.data[is.na(s.data)]<-0


#make taxanomy reference
mpa<-data.frame()
for (i in rownames(total.data)){
  buf<-paste0('/data2/msi/endometriosis/endometriosis/kraken/',total.data[i,'Run'],'.trim.fastq.report.tsv')
  buf<-read.table(buf,sep = '\t',header=F)
  for (j in buf$V1){
    k<-unlist(strsplit(j,'[|]'))
    k<-unlist(strsplit(k[length(k)],'__'))
    j<-c(j,k)
    if (nrow(mpa)==0){mpa<-data.frame(j)}
    else{mpa<-cbind(mpa,j)}
  }
}
mpa<-data.frame(t(mpa))
mpa$X3<-gsub(' ','_',mpa$X3)
mpa$X3<-gsub('-','.',mpa$X3)
mpa$X3<-gsub('[[]','X.',mpa$X3)
mpa$X3<-gsub('[]]','.',mpa$X3)
mpa$X3<-gsub('[(]','.',mpa$X3)
mpa$X3<-gsub('[)]','.',mpa$X3)
mpa<-mpa[!duplicated(mpa),]

#labeling
buf<-data.frame(colnames(s.data))

buf1<-merge(buf,mpa,all.x=T,by.x='colnames.s.data.',by.y='X3')
table(is.na(buf1$X1))
rownames(buf1)<-buf1$colnames.s.data.
colnames(s.data)<-buf1[colnames(s.data),'X1']

lefse.sdata<-data.frame(t(s.data))

ccol<-c()
for (j in colnames(s.data)){
  i<-unlist(strsplit(j,'[|]'))
  i<-lapply(i,FUN = function(x){unlist(strsplit(x,'__'))})
  if (length(i)==1){
    ccol<-c(ccol,j)
  } else if (i[1][[1]][2]=='Viruses'){
    ccol<-c(ccol,paste0('Virus_',i[length(i)][[1]][1],'_',i[length(i)][[1]][2]))
  } else if(i[1][[1]][2]=='Bacteria'){
    ccol<-c(ccol,paste0('Bacteria_',i[length(i)][[1]][1],'_',i[length(i)][[1]][2]))
  } else if(i[2][[1]][2]=='Fungi'){
    ccol<-c(ccol,paste0('Fungi_',i[length(i)][[1]][1],'_',i[length(i)][[1]][2]))
  } else if(i[1][[1]][2]=='Eukaryota'){
    ccol<-c(ccol,paste0('Eukaryota_',i[length(i)][[1]][1],'_',i[length(i)][[1]][2]))
  } else if(i[1][[1]][2]=='Archaea'){
    ccol<-c(ccol,paste0('Archaea_',i[length(i)][[1]][1],'_',i[length(i)][[1]][2]))
  } else(ccol<-c(ccol,j))
}

colnames(s.data)<-ccol


s.data['Run']<-rownames(s.data)

total.data.s<-merge(total.data,s.data,by='Run',all=T)


#classified by fungi, bacteria, virus
Fun<-colnames(total.data.s)[startsWith(colnames(total.data.s),'F')]
Bac<-colnames(total.data.s)[startsWith(colnames(total.data.s),'B')]
Vir<-colnames(total.data.s)[startsWith(colnames(total.data.s),'V')]

total.data.sf<-total.data.s[,c(colnames(total.data.s)[1:11],Fun)]
total.data.sb<-total.data.s[,c(colnames(total.data.s)[1:11],Bac)]
total.data.sv<-total.data.s[,c(colnames(total.data.s)[1:11],Vir)]

#########################################여기까지###########################################
##richness and evenness check----
alpha.s<-total.data.s[,c('source_name','condition','tissue')]
alpha.s$Chao1<-estimateR(total.data.s[,12:(ncol(total.data.s))])[1,]
alpha.s$shannon<-vegan::diversity(total.data.s[,12:(ncol(total.data.s))],index="shannon")

alpha.s$Chao1_f<-estimateR(total.data.sf[,12:(ncol(total.data.sf))])[1,]
alpha.s$shannon_f<-vegan::diversity(total.data.sf[,12:(ncol(total.data.sf))],index="shannon")

alpha.s$Chao1_b<-estimateR(total.data.sb[,12:(ncol(total.data.sb))])[1,]
alpha.s$shannon_b<-vegan::diversity(total.data.sb[,12:(ncol(total.data.sb))],index="shannon")


###alpha-----
M.alpha.s<-melt(alpha.s,id.vars = c("tissue",'source_name','condition'))
M.alpha.s$variable<-factor(M.alpha.s$variable,levels=c('Chao1','Chao1_b','Chao1_f','shannon','shannon_b','shannon_f'))


p.s.diversity.chao<-ggplot(M.alpha.s[M.alpha.s$variable%in%c('Chao1','Chao1_b'),],aes(x=condition,y=value,fill=condition))+geom_boxplot()+theme_pubr()+
  facet_wrap(~variable,nrow=1)+
  stat_compare_means(method='t.test',label = 'p',label.x.npc = 0,label.y.npc = 0.05, family='serif',comparisons=list(c('Control','Endometriosis')))+
  theme(text=element_text(family="serif",face = 'bold'), axis.text = element_text(face='plain'),
        axis.title.x = element_blank(),axis.text.x=element_blank(), axis.title.y=element_text(face='plain'),
        legend.title=element_blank(),legend.text=element_text(face='plain'))+
  labs(y='score')
p.s.diversity.shannon<-ggplot(M.alpha.s[M.alpha.s$variable%in%c('shannon','shannon_b'),],aes(x=condition,y=value,fill=condition))+geom_boxplot()+theme_pubr()+
  facet_wrap(~variable,nrow=1)+
  stat_compare_means(method='t.test',label = 'p',label.x.npc = 0,label.y.npc = 0.05, family='serif',comparisons=list(c('Control','Endometriosis')))+
  theme(text=element_text(family="serif",face = 'bold'), axis.text = element_text(face='plain'),
        axis.title.x = element_blank(),axis.text.x=element_blank(), axis.title.y=element_text(face='plain'),
        legend.title=element_blank(),legend.text=element_text(face='plain'))+
  labs(y='score')

p.s.diversity.chao+p.s.diversity.shannon

p.s.diversity.tb<-ggplot(M.alpha.s,aes(x=tissue,y=value,fill=tissue))+geom_boxplot()+theme_pubr()+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(method='t.test',label = 'p',label.x.npc = 0,label.y.npc = 0.05, family='serif',comparisons=list(c("Control", "Endometrium"),c("Endometrium", "Ovary"),c("Endometrium", "Peritoneum"),c('Ovary','Peritoneum')))+
  theme(text=element_text(family="serif",face = 'bold'), axis.text = element_text(face='plain'),
        axis.title.x = element_blank(),axis.text.x=element_blank(), axis.title.y=element_text(face='plain'),
        legend.title=element_blank(),legend.text=element_text(face='plain'))+
  labs(y='score',tag='')

p.s.diversity.bb+p.s.diversity.tb



#visualization
p.s.diversity.bi<-ggplot(M.alpha.s,aes(x=condition,y=value,fill=condition))+geom_boxplot()+theme_pubr()+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(method='t.test',label = 'p',label.x.npc = 0,label.y.npc = 0.05, family='serif',comparisons=list(c('Control','Endometriosis')))+
  theme(text=element_text(family="serif",face = 'bold'), axis.text = element_text(face='plain'),
        axis.title.x = element_blank(),axis.text.x=element_blank(), axis.title.y=element_text(face='plain'),
        legend.title=element_blank(),legend.text=element_text(face='plain'))+
  labs(y='score',tag='Species')
p.s.diversity.ti<-ggplot(M.alpha.s,aes(x=tissue,y=value,fill=tissue))+geom_boxplot()+theme_pubr()+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(method='t.test',label = 'p',label.x.npc = 0,label.y.npc = 0.05, family='serif',comparisons=list(c("Control", "Endometrium"),c("Endometrium", "Ovary"),c("Endometrium", "Peritoneum"),c('Ovary','Peritoneum')))+
  theme(text=element_text(family="serif",face = 'bold'), axis.text = element_text(face='plain'),
        axis.title.x = element_blank(),axis.text.x=element_blank(), axis.title.y=element_text(face='plain'),
        legend.title=element_blank(),legend.text=element_text(face='plain'))+
  labs(y='score',tag='')

p.g.diversity.bi<-ggplot(M.alpha.g,aes(x=condition,y=value,fill=condition))+geom_boxplot()+theme_pubr()+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(method='t.test',label = 'p',label.x.npc = 0,label.y.npc = 0.05, family='serif',comparisons=list(c('Control','Endometriosis')))+
  theme(text=element_text(family="serif",face = 'bold'), axis.text = element_text(face='plain'),
        axis.title.x = element_blank(),axis.text.x=element_blank(), axis.title.y=element_text(face='plain'),
        legend.title=element_blank(),legend.text=element_text(face='plain'))+
  labs(y='score',tag='Genus')
p.g.diversity.ti<-ggplot(M.alpha.g,aes(x=tissue,y=value,fill=tissue))+geom_boxplot()+theme_pubr()+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(method='t.test',label = 'p',label.x.npc = 0,label.y.npc = 0.05, family='serif',comparisons=list(c("Control", "Endometrium"),c("Endometrium", "Ovary"),c("Endometrium", "Peritoneum"),c('Ovary','Peritoneum')))+
  theme(text=element_text(family="serif",face = 'bold'), axis.text = element_text(face='plain'),
        axis.title.x = element_blank(),axis.text.x=element_blank(), axis.title.y=element_text(face='plain'),
        legend.title=element_blank(),legend.text=element_text(face='plain'))+
  labs(y='score',tag='')




###beta====
beta.dist.g<-vegdist(total.data.g[,c(6:ncol(total.data.g))],method='bray')
beta.dist.s<-vegdist(total.data.s[,c(6:ncol(total.data.s))],method='bray')
beta.dist.gb<-vegdist(total.data.gb[,c(6:ncol(total.data.gb))],method='bray')
beta.dist.sb<-vegdist(total.data.sb[,c(6:ncol(total.data.sb))],method='bray')

mds.g<-metaMDS(beta.dist.g)
mds.s<-metaMDS(beta.dist.s)
mds.gb<-metaMDS(beta.dist.gb)
mds.sb<-metaMDS(beta.dist.sb)

mds.data.g<-as.data.frame(mds.g$points)
mds.data.s<-as.data.frame(mds.s$points)
mds.data.gb<-as.data.frame(mds.gb$points)
mds.data.sb<-as.data.frame(mds.sb$points)

mds.data.g$type<-total.data.g$tissue
mds.data.s$type<-total.data.s$tissue
mds.data.gb$type<-total.data.gb$tissue
mds.data.sb$type<-total.data.sb$tissue
mds.data.gb$bi<-total.data.gb$condition
mds.data.sb$bi<-total.data.sb$condition


mds.data.g$condition<-ifelse(mds.data.g$type=='Control','Control','Endometriosis')
mds.data.s$condition<-ifelse(mds.data.s$type=='Control','Control','Endometriosis')

pv.g<-adonis2(total.data.g[,c(6:ncol(total.data.g))]~tissue,total.data.g)$`Pr(>F)`[1] #p value calculate (Pr(>F))
pv.s<-adonis2(total.data.s[,c(6:ncol(total.data.s))]~tissue,total.data.s)$`Pr(>F)`[1]

pv.gb<-adonis2(total.data.gb[,c(6:ncol(total.data.gb))]~tissue,total.data.gb)$`Pr(>F)`[1] #p value calculate (Pr(>F))
pv.sb<-adonis2(total.data.sb[,c(6:ncol(total.data.sb))]~tissue,total.data.sb)$`Pr(>F)`[1]


p.g.beta<-ggplot(mds.data.g,aes(x=MDS1,y=MDS2,color=type))+geom_point()+theme_pubr()+
  theme(text=element_text(family="serif"),legend.title=element_blank())+
  annotate("text", x=min(mds.data.g$MDS1)+0.1, y=max(mds.data.g$MDS2), label=paste0("p = ",pv.g),family='serif')
p.s.beta<-ggplot(mds.data.s,aes(x=MDS1,y=MDS2,color=type))+geom_point()+theme_pubr()+
  theme(text=element_text(family="serif"),legend.title=element_blank())+
  annotate("text", x=min(mds.data.s$MDS1)+0.1, y=max(mds.data.s$MDS2), label=paste0("p = ",pv.s),family='serif')

p.gb.beta<-ggplot(mds.data.gb,aes(x=MDS1,y=MDS2,color=type))+geom_point()+theme_pubr()+
  theme(text=element_text(family="serif"),legend.title=element_blank())+
  annotate("text", x=min(mds.data.g$MDS1)+0.1, y=max(mds.data.g$MDS2), label=paste0("p = ",pv.gb),family='serif')
p.sb.beta<-ggplot(mds.data.sb,aes(x=MDS1,y=MDS2,color=type))+geom_point()+theme_pubr()+
  theme(text=element_text(family="serif"),legend.title=element_blank())+
  annotate("text", x=min(mds.data.s$MDS1)+0.1, y=max(mds.data.s$MDS2), label=paste0("p = ",pv.sb),family='serif')

grid.arrange(grid.arrange(p.s.diversity.chao,p.s.diversity.shannon)+plot_layout(ncol=1),p.gb.beta,p.sb.beta)+plot_layout(widths=c(1,1))

###visualization----
p.s.diversity.chao+p.s.diversity.shannon+p.sb.beta+plot_layout(widths=c(1,1,1.4))
#p.s.diversity.bi+p.s.diversity.ti+p.s.beta+p.g.diversity.bi+p.g.diversity.ti+p.g.beta+plot_layout(widths=c(1,1.4,0.8),heights = c(1,1))



####다른 정보에 따른 microbiome 차이====
total.data.sb$`Bacteria_s_Bacillus cereus`
buf<-total.data.s[c('source_name','condition','tissue','Age','Ethnicity','Pregnancy','rASRM Stage','Oral contraceptive treatment','Endometrium histology','Bacteria_s_Bacillus cereus')]

buf1<-ggplot(buf,aes(x=Age,y=`Bacteria_s_Bacillus cereus`,col=condition))+geom_point()
buf2<-ggplot(buf,aes(x=Ethnicity,y=`Bacteria_s_Bacillus cereus`,col=condition))+geom_point()
buf3<-ggplot(buf,aes(x=Pregnancy,y=`Bacteria_s_Bacillus cereus`,col=condition))+geom_point()
buf5<-ggplot(buf,aes(x=tissue,y=`Bacteria_s_Bacillus cereus`,col=`Endometrium histology`))+geom_point()
buf7<-ggplot(buf,aes(x=`Oral contraceptive treatment`,y=`Bacteria_s_Bacillus cereus`,col=condition))+geom_point()
buf1+buf2+buf3+buf5+buf7 #전체적 확인

buf<-total.data.s[total.data.s$condition=='Control',c('source_name','condition','tissue','Age','Ethnicity','Pregnancy','rASRM Stage','Oral contraceptive treatment','Endometrium histology','Bacteria_s_Bacillus cereus')]
buf4<-ggplot(buf,aes(x=source_name,y=`Bacteria_s_Bacillus cereus`,col=`Endometrium histology`))+geom_point()

buf<-total.data.s[total.data.s$`Endometrium histology`=='exogenous hormonal effect',c('source_name','condition','tissue','Age','Ethnicity','Pregnancy','rASRM Stage','Endometrium histology','Bacteria_s_Bacillus cereus')]

buf6<-ggplot(buf,aes(x=source_name,y=`Bacteria_s_Bacillus cereus`,col=condition))+geom_point()

buf4+buf6

#different abundance analysis----
##Aldex----

set.seed(4232)
aldex.gf<-aldex(t(total.data.gf[,6:ncol(total.data.gf)]), total.data.gf$condition)
aldex.gb<-aldex(t(total.data.gb[,6:ncol(total.data.gb)]), total.data.gb$condition)
aldex.s<-aldex(t(total.data.s[,6:ncol(total.data.s)]), total.data.s$condition)

aldex.sb<-aldex(t(total.data.sb[,6:ncol(total.data.sb)]), total.data.sb$condition) #결과 x


biomelist[order(biomelist)]==rownames(res.g)[order(rownames(res.g))]
res.g<-res.gastric %>%  filter(wi.ep<0.05) %>% filter(abs(diff.btw)>1.5) #123
res.g<-res.g[order(res.g$diff.btw,decreasing = T),]
res.g$diff.btw<-res.g$diff.btw*(-1)
res.g$group<-ifelse(res.g$diff.btw>0,'Cancer','Normal')
res.g$group<-factor(res.g$group,levels=c('Normal','Cancer'))

res.g$row<-rownames(res.g)

ggbarplot(res.g,x='row',y='diff.btw',rotate=TRUE,fill='group',size = 0.3)+scale_fill_manual(values=Color)+theme(text=element_text(family='serif'),legend.title=element_blank(),axis.title.y = element_blank(),axis.text.y=element_text(face = 'bold'))

##Lefse====
library("phyloseq")
library(microbiomeMarker)


###genus====
taxamet<-data.frame(Kingdom="", Phylum='', Class='', Order='', Family='', Genus='')
taxamet<-taxamet[-1,]
for (i in rownames(lefse.gdata)){
  for (j in unlist(strsplit(i,'[|]'))){
    slp<-unlist(strsplit(j,'__'))
    if (slp[1]=='d'){taxamet[i,'Kingdom']<-paste0('D',slp[2])}
    else if (slp[1]=='k'){taxamet[i,'Kingdom']<-slp[2]}
    else if (slp[1]=='p'){taxamet[i,'Phylum']<-slp[2]}
    else if (slp[1]=='c'){taxamet[i,'Class']<-slp[2]}
    else if (slp[1]=='o'){taxamet[i,'Order']<-slp[2]}
    else if (slp[1]=='f'){taxamet[i,'Family']<-slp[2]}
    else if (slp[1]=='g'){taxamet[i,'Genus']<-slp[2]}
  }
}

lefse.gdata.b<-lefse.gdata[startsWith(rownames(lefse.gdata),'d__Bacteria'),]
lefse.gdata.f<-lefse.gdata[startsWith(rownames(lefse.gdata),'d__Eukaryota|k__Fungi'),]
lefse.gdata.v<-lefse.gdata[startsWith(rownames(lefse.gdata),'d__Viruses'),]

lefse.sdata.b<-lefse.sdata[startsWith(rownames(lefse.sdata),'d__Bacteria'),]
lefse.sdata.f<-lefse.sdata[startsWith(rownames(lefse.sdata),'d__Eukaryota|k__Fungi'),]
lefse.sdata.v<-lefse.sdata[startsWith(rownames(lefse.sdata),'d__Viruses'),]

taxamet.b<-taxamet[startsWith(rownames(taxamet),'d__Bacteria'),]
taxamet.f<-taxamet[startsWith(rownames(taxamet),'d__Eukaryota|k__Fungi'),]
taxamet.v<-taxamet[startsWith(rownames(taxamet),'d__Viruses'),]

#total
OTU = otu_table(as.matrix(lefse.gdata), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxamet))
SAM<-sample_data(data.frame(total.data,row.names = 'Run'))
physeq.g = phyloseq(OTU, TAX, SAM)

rownames(physeq.g@otu_table)<-1:nrow(physeq.g@otu_table)
rownames(physeq.g@tax_table)<-1:nrow(physeq.g@tax_table)

#bacteria
OTU = otu_table(as.matrix(lefse.gdata.b), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxamet.b))
physeq.gb <- phyloseq(OTU, TAX, SAM)

rownames(physeq.gb@otu_table)<-1:nrow(physeq.gb@otu_table)
rownames(physeq.gb@tax_table)<-1:nrow(physeq.gb@tax_table)

#fungi
OTU = otu_table(as.matrix(lefse.gdata.f), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxamet.f))
physeq.gf <- phyloseq(OTU, TAX, SAM)

rownames(physeq.gf@otu_table)<-1:nrow(physeq.gf@otu_table)
rownames(physeq.gf@tax_table)<-1:nrow(physeq.gf@tax_table)

#virus
OTU = otu_table(as.matrix(lefse.gdata.v), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxamet.v))
physeq.gv <- phyloseq(OTU, TAX, SAM)

rownames(physeq.gv@otu_table)<-1:nrow(physeq.gv@otu_table)
rownames(physeq.gv@tax_table)<-1:nrow(physeq.gv@tax_table)


###Species====
taxamet<-data.frame(Kingdom="", Phylum='', Class='', Order='', Family='', Genus='',Species='')
taxamet<-taxamet[-1,]
for (i in rownames(lefse.sdata)){
  for (j in unlist(strsplit(i,'[|]'))){
    slp<-unlist(strsplit(j,'__'))
    if (slp[1]=='d'){taxamet[i,'Kingdom']<-paste0('D',slp[2])}
    else if (slp[1]=='k'){taxamet[i,'Kingdom']<-slp[2]}
    else if (slp[1]=='p'){taxamet[i,'Phylum']<-slp[2]}
    else if (slp[1]=='c'){taxamet[i,'Class']<-slp[2]}
    else if (slp[1]=='o'){taxamet[i,'Order']<-slp[2]}
    else if (slp[1]=='f'){taxamet[i,'Family']<-slp[2]}
    else if (slp[1]=='g'){taxamet[i,'Genus']<-slp[2]}
    else if (slp[1]=='s'){taxamet[i,'Species']<-slp[2]}
  }
}

taxamet.b<-taxamet[startsWith(rownames(taxamet),'d__Bacteria'),]
taxamet.f<-taxamet[startsWith(rownames(taxamet),'d__Eukaryota|k__Fungi'),]
taxamet.v<-taxamet[startsWith(rownames(taxamet),'d__Viruses'),]

#total
OTU = otu_table(as.matrix(lefse.sdata), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxamet))
SAM<-sample_data(data.frame(total.data,row.names = 'Run'))
physeq.s = phyloseq(OTU, TAX, SAM)

rownames(physeq.s@otu_table)<-1:nrow(physeq.s@otu_table)
rownames(physeq.s@tax_table)<-1:nrow(physeq.s@tax_table)

#bacteria
OTU = otu_table(as.matrix(lefse.sdata.b), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxamet.b))
physeq.sb <- phyloseq(OTU, TAX, SAM)

rownames(physeq.sb@otu_table)<-1:nrow(physeq.sb@otu_table)
rownames(physeq.sb@tax_table)<-1:nrow(physeq.sb@tax_table)

#fungi
OTU = otu_table(as.matrix(lefse.sdata.f), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxamet.f))
physeq.sf <- phyloseq(OTU, TAX, SAM)

rownames(physeq.sf@otu_table)<-1:nrow(physeq.sf@otu_table)
rownames(physeq.sf@tax_table)<-1:nrow(physeq.sf@tax_table)

#virus
OTU = otu_table(as.matrix(lefse.sdata.v), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxamet.v))
physeq.sv <- phyloseq(OTU, TAX, SAM)

rownames(physeq.sv@otu_table)<-1:nrow(physeq.sv@otu_table)
rownames(physeq.sv@tax_table)<-1:nrow(physeq.sv@tax_table)


#visualization
mm_test.gt <- run_lefse( physeq.g,
                         wilcoxon_cutoff = 0.05,
                         norm = "CPM", 
                         group = "tissue",
                         kw_cutoff = 0.01,
                         taxa_rank = 'all',
                         multigrp_strat = TRUE
)
mm_test.gb <- run_lefse( physeq.g,
                         wilcoxon_cutoff = 0.05,
                         norm = "CPM", 
                         group = "condition",
                         kw_cutoff = 0.01,
                         taxa_rank = 'all',
                         multigrp_strat = TRUE
)
mm_test.sb <- run_lefse( physeq.s,
                        wilcoxon_cutoff = 0.05,
                        norm = "CPM", 
                        group = "condition",
                        kw_cutoff = 0.01,
                        taxa_rank = 'all',
                        multigrp_strat = TRUE
)
mm_test.st <- run_lefse( physeq.s,
                         wilcoxon_cutoff = 0.05,
                         norm = "CPM", 
                         group = "tissue",
                         kw_cutoff = 0.01,
                         taxa_rank = 'all',
                         multigrp_strat = TRUE
)


mm_test.sbb<- run_lefse( physeq.sb,
                         wilcoxon_cutoff = 0.05,
                         norm = "CPM", 
                         group = "condition",
                         kw_cutoff = 0.05,
                         taxa_rank = 'all',
                         multigrp_strat = TRUE
)
mm_test.stb<- run_lefse( physeq.sb,
                         wilcoxon_cutoff = 0.05,
                         norm = "CPM", 
                         group = "tissue",
                         kw_cutoff = 0.05,
                         taxa_rank = 'all',
                         multigrp_strat = TRUE
)

tof<-c()
for (i in mm_test.sbb@marker_table$feature){
  if(length(unlist(strsplit(i,'[|]')))==7){tof<-c(tof,T)}
  else{tof<-c(tof,F)}
}
mm_test.sbb@marker_table$feature[!tof]


  
plot_ef_bar(mm_test.sbb)+plot_ef_bar(mm_test.sbb@marker_table[tof])


buf<-total.data.s[,c('Bacteria_s_Bacillus cereus','condition')]
ggplot(buf,aes(x=condition,y=`Bacteria_s_Bacillus cereus`))+geom_point()

((plot_ef_bar(mm_test.sbf)+labs(tag='Fungi'))+(plot_ef_bar(mm_test.sb)+labs(tag='All')))




#limma with Species====
##limma----
library(limma)
#limma.eset<-total.data.sf[,6:ncol(total.data.sf)]
limma.eset<-total.data.sb[,6:ncol(total.data.sb)]
design<-data.frame(ifelse(total.data.sf$condition=='Control',0,1),ifelse(total.data.sf$condition=='Control',1,0))
colnames(design)<-c('Endometriosis','Control')
fit <- lmFit(t(limma.eset),design)
cont <- makeContrasts(Endometriosis-Control,levels=design)
fit.cont <- contrasts.fit(fit,cont)
fit.cont <- eBayes(fit.cont)
res <- topTable(fit.cont,number=Inf)
#res.fliter<-res[res$adj.P.Val<0.05,]
res<-res[abs(res$logFC)>1,]
res.filter<-res[res$P.Value<0.05,]

res$row<-rownames(res)
res$group<-ifelse(res$logFC>0,'Endometriosis','Control')
res<-sort_df(res,vars='logFC')
res<-res[abs(res$logFC)>1,]

res.filter$row<-rownames(res.filter)
res.filter$group<-ifelse(res.filter$logFC>0,'Endometriosis','Control')
res.filter<-sort_df(res.filter,vars='logFC')

p.limma<-ggbarplot(res,x='row',y='logFC',fill='group',rotate=TRUE,size = 0.3)+theme(text=element_text(family='serif'),legend.title=element_blank(),axis.title.y = element_blank(),axis.text.y=element_text(size=10))+labs(tag='limma')#+scale_y_continuous(limits =c(-6,2))

p.limma.filter<-ggbarplot(res.filter,x='row',y='logFC',fill='group',rotate=TRUE,size = 0.3)+theme(text=element_text(family='serif'),legend.title=element_blank(),axis.title.y = element_blank(),axis.text.y=element_text(size=10))+labs(tag='limma')+scale_y_continuous(limits =c(-6,0))

p.limma+p.limma.filter

##edgeR----
library(edgeR)
phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  #ref: https://github.com/nearinj/Comparison_of_DA_microbiome_methods/blob/main/Pipeline_scripts/Tool_scripts/Run_edgeR.R
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}

test <- phyloseq_to_edgeR(physeq = physeq.s, group=physeq.s@sam_data$condition)
et = exactTest(test)
tt = topTags(et, n=nrow(test$table), adjust.method="fdr", sort.by="PValue")
edger.res <- tt@.Data[[1]]

test <- phyloseq_to_edgeR(physeq = physeq.sb, group=physeq.sf@sam_data$condition)
et = exactTest(test)
tt = topTags(et, n=nrow(test$table), adjust.method="fdr", sort.by="PValue")
edger.res.b <- tt@.Data[[1]]

edger.res.b$row<-edger.res.b$Species
edger.res.b$group<-ifelse(edger.res.b$logFC>0,'Endometriosis','Control')
edger.res.b<-sort_df(edger.res.b,vars='logFC')
edger.res.b<-edger.res.b[edger.res.b$PValue<0.05,]

p.edger<-ggbarplot(edger.res.b,x='row',y='logFC',fill='group',rotate=TRUE,size = 0.3)+theme(text=element_text(family='serif'),legend.title=element_blank(),axis.title.y = element_blank(),axis.text.y=element_text(size=10))+labs(tag='edgeR')

p.limma.filter+p.edger

#limma with Genus====
##limma----
library(limma)
limma.eset<-total.data.gf[,6:ncol(total.data.gf)]
design<-data.frame(ifelse(total.data.gf$condition=='Control',0,1),ifelse(total.data.gf$condition=='Control',1,0))
colnames(design)<-c('Endometriosis','Control')
fit <- lmFit(t(limma.eset),design)

##edgeR----
library(edgeR)
phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  #ref: https://github.com/nearinj/Comparison_of_DA_microbiome_methods/blob/main/Pipeline_scripts/Tool_scripts/Run_edgeR.R
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}

test <- phyloseq_to_edgeR(physeq = physeq.sf, group=physeq.sf@sam_data$condition)
et = exactTest(test)
tt = topTags(et, n=nrow(test$table), adjust.method="fdr", sort.by="PValue")
edger.res.sf <- tt@.Data[[1]]

edger.res.sf$row<-edger.res.sf$Species
edger.res.sf$group<-ifelse(edger.res.sf$logFC>0,'Endometriosis','Control')
edger.res.sf<-sort_df(edger.res.sf,vars='logFC')

buf<-edger.res.sf[edger.res.sf$PValue<0.05,]

p.edger<-ggbarplot(edger.res.sf,x='row',y='logFC',fill='group',rotate=TRUE,size = 0.3)+theme(text=element_text(family='serif'),legend.title=element_blank(),axis.title.y = element_blank(),axis.text.y=element_text(size=10))+labs(tag='edgeR')

p.limma+p.edger



#after upside result====
##preprocessing----
total.data.s.filter<-total.data.s[!total.data.s$`Endometrium histology`%in%c('exogenous hormonal effect','inactive endometrium','n.a'),]

#total.data.s.filter<-total.data.s.filter[,c(rep(T,11),colSums(total.data.s.filter[12:ncol(total.data.s.filter)])!=0)]

total.data.s.filter<-total.data.s.filter[total.data.s.filter$tissue!='Endometrium',]


##richness evenness----
alpha.filter<-total.data.s.filter[,c('source_name','condition','tissue')]
alpha.filter$Chao1<-estimateR(total.data.s.filter[,(ncol(total.data)+1):(ncol(total.data.s.filter))])[2,]
alpha.filter$shannon<-vegan::diversity(total.data.s.filter[,(ncol(total.data)+1):(ncol(total.data.s.filter))],index="shannon")

buf<-colnames(total.data.s.filter)[ncol(total.data)+1:ncol(total.data.s.filter)]
Fun<-startsWith(colnames(total.data.s.filter),'F')
Vir<-startsWith(colnames(total.data.s.filter),'V')
Bac<-startsWith(colnames(total.data.s.filter),'B')
Euc<-startsWith(colnames(total.data.s.filter),'d__')

alpha.filter$Chao1_f<-estimateR(total.data.s.filter[,Fun])[2,]
alpha.filter$Chao1_b<-estimateR(total.data.s.filter[,Bac])[2,]

alpha.filter$shannon_f<-vegan::diversity(total.data.s.filter[,Fun],index="shannon")
alpha.filter$shannon_b<-vegan::diversity(total.data.s.filter[,Bac],index="shannon")

##C O L O R====
Color=c('Control'='#01D11A','Endometriosis'='#CB1106')
Color1=c('Control'='#01D11A99','Endometriosis'='#CB110699')
Color3<-c(Color,'Patient_Endometrium'='#FFD700')
abb.Color<-Color
abb.Color1<-Color1
names(abb.Color)<-c('HC','EM')
names(abb.Color1)<-c('HC','EM')

###alpha-----


M.alpha.filter<-melt(alpha.filter,id.vars = c("tissue",'source_name','condition'))

M.alpha.filter$variable<-factor(M.alpha.filter$variable,levels=c('Chao1','Chao1_b','Chao1_f','shannon','shannon_b','shannon_f'))
M.alpha.filter$condition<-factor(ifelse(M.alpha.filter$condition=='Control','HC','EM'),levels=c('HC','EM'))

p.filter.diversity<-ggplot(M.alpha.filter[M.alpha.filter$variable%in%c('Chao1','shannon'),],aes(x=condition,y=value,fill=condition,color=condition))+geom_boxplot()+theme_pubr()+scale_fill_manual(values = abb.Color1)+
  facet_wrap(~variable,nrow=1,scales = 'free')+
  stat_compare_means(method='wilcox.test',label = 'p',label.x.npc = 0.6,label.y.npc = 0.9, family='arial',)+
  theme(text=element_text(family="arial",face = 'bold',size=10), axis.text = element_text(face='plain'),
        axis.title.x = element_blank(),axis.text.x=element_blank(), axis.title.y=element_text(face='plain'),
        legend.title=element_blank(),legend.text=element_text(face='plain'))+
  labs(y='score')+scale_color_manual(values=abb.Color)

shapiro.test(M.alpha.filter[M.alpha.filter$variable%in%c('Chao1'),'value'])
shapiro.test(M.alpha.filter[M.alpha.filter$variable%in%c('shannon'),'value'])

###beta====

beta.dist.filter<-vegdist(total.data.s.filter[Bac+Vir+Fun+Euc==1],method='bray')
mds.filter<-metaMDS(beta.dist.filter)

mds.data.filter<-as.data.frame(mds.filter$points)

mds.data.filter$type<-total.data.s.filter$condition
mds.data.filter$type<-factor(ifelse(mds.data.filter$type=='Control','HC','EM'),levels=c('HC','EM'))

#set.seed(6212)
set.seed(4230)
pv.filter<-adonis2(total.data.s.filter[Bac+Vir+Fun+Euc==1]~condition,total.data.s.filter, method='bray')$`Pr(>F)`[1] #p value calculate (Pr(>F))
pv.filter3<-adonis2(total.data.s.filter[Bac+Vir+Fun+Euc==1]~condition3,total.data.s.filter, method='bray')$`Pr(>F)`[1] #p value calculate (Pr(>F))


p.filter.beta<-ggplot(mds.data.filter, aes(x=MDS1,y=MDS2,color=type))+geom_point()+theme_pubr()+
  theme(text=element_text(family="arial",size=10),legend.title=element_blank())+scale_color_manual(values = c(abb.Color))+
  annotate("text", x=min(mds.data.filter$MDS1)+0.4, y=max(mds.data.filter$MDS2), label=paste0("p = ",pv.filter),family='arial')


###visualization----

p.filter.diversity+p.filter.beta+plot_layout(widths=c(1,1.8))



##different abundance analysis====

###lefse----
library("phyloseq")
library(microbiomeMarker)

lefse.filter<-lefse.sdata[total.data.s.filter$Run]
lefse.filter.b<-lefse.filter[startsWith(rownames(lefse.filter),'d__Bacteria'),]
lefse.filter.f<-lefse.filter[startsWith(rownames(lefse.filter),'d__Eukaryota|k__Fungi'),]
lefse.filter.v<-lefse.filter[startsWith(rownames(lefse.filter),'d__Viruses'),]
lefse.filter.o<-lefse.filter[!(rownames(lefse.filter)%in%rownames(lefse.filter.b)),]

taxamet<-data.frame(Kingdom="", Phylum='', Class='', Order='', Family='', Genus='',Species='')
taxamet<-taxamet[-1,]
for (i in rownames(lefse.filter)){
  for (j in unlist(strsplit(i,'[|]'))){
    slp<-unlist(strsplit(j,'__'))
    if (slp[1]=='d'){taxamet[i,'Kingdom']<-paste0('D',slp[2])}
    else if (slp[1]=='k'){taxamet[i,'Kingdom']<-slp[2]}
    else if (slp[1]=='p'){taxamet[i,'Phylum']<-slp[2]}
    else if (slp[1]=='c'){taxamet[i,'Class']<-slp[2]}
    else if (slp[1]=='o'){taxamet[i,'Order']<-slp[2]}
    else if (slp[1]=='f'){taxamet[i,'Family']<-slp[2]}
    else if (slp[1]=='g'){taxamet[i,'Genus']<-slp[2]}
    else if (slp[1]=='s'){taxamet[i,'Species']<-slp[2]}
  }
}

taxamet.b<-taxamet[startsWith(rownames(taxamet),'d__Bacteria'),]
taxamet.f<-taxamet[startsWith(rownames(taxamet),'d__Eukaryota|k__Fungi'),]
taxamet.v<-taxamet[startsWith(rownames(taxamet),'d__Viruses'),]
taxamet.o<-taxamet[!(rownames(taxamet)%in%rownames(taxamet.b)),]

#total

OTU = otu_table(as.matrix(lefse.filter), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxamet))
SAM<-sample_data(data.frame(total.data.s.filter,row.names = 'Run'))
physeq.filter = phyloseq(OTU, TAX, SAM)

rownames(physeq.filter@otu_table)<-1:nrow(physeq.filter@otu_table)
rownames(physeq.filter@tax_table)<-1:nrow(physeq.filter@tax_table)

#bacteria
OTU = otu_table(as.matrix(lefse.filter.b), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxamet.b))

physeq.filter.b <- phyloseq(OTU, TAX, SAM)

rownames(physeq.filter.b@otu_table)<-1:nrow(physeq.filter.b@otu_table)
rownames(physeq.filter.b@tax_table)<-1:nrow(physeq.filter.b@tax_table)

#fungi
OTU = otu_table(as.matrix(lefse.filter.f), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxamet.f))
physeq.filter.f <- phyloseq(OTU, TAX, SAM)

rownames(physeq.filter.f@otu_table)<-1:nrow(physeq.filter.f@otu_table)
rownames(physeq.filter.f@tax_table)<-1:nrow(physeq.filter.f@tax_table)

#other
OTU = otu_table(as.matrix(lefse.filter.o), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxamet.o))
SAM = sample_data(data.frame(total.data.s.filter[1:11],row.names = 'Run'))
physeq.filter.o <- phyloseq(OTU, TAX, SAM)

rownames(physeq.filter.o@otu_table)<-1:nrow(physeq.filter.o@otu_table)
rownames(physeq.filter.o@tax_table)<-1:nrow(physeq.filter.o@tax_table)

#visualization
mm_test.filter <- run_lefse( physeq.filter,
                         norm = "CPM", 
                         kw_cutoff=0.1,
                         group = "condition",
                         taxa_rank = 'all',
                         multigrp_strat = TRUE
)
plot_ef_bar(mm_test.filter)

mm_test.filter.b <- run_lefse( physeq.filter.b,
                         kw_cutoff=0.1,
                         lda_cutoff = 2.5,
                         norm = "CPM", 
                         group = "condition",
                         taxa_rank = 'all',
                         multigrp_strat = TRUE
)
plot_ef_bar(mm_test.filter.b)

total.data.s.filter[endsWith(colnames(total.data.s.filter),'tum')]

tof<-c()
for (i in mm_test.filter.b@marker_table$feature){
  if(length(unlist(strsplit(i,'[|]')))==7){tof<-c(tof,T)}
  else{tof<-c(tof,F)}
}

plot_ef_bar(mm_test.filter.b)+plot_ef_bar(mm_test.filter.b@marker_table[tof])


buf<-mm_test.filter.b@marker_table[tof]
buf$feature<-gsub('s__','',buf$feature)
buf$enrich_group<-factor(ifelse(buf$enrich_group=='Control','HC','EM'),levels=c('HC','EM'))

plot_ef_bar(buf)+theme_pubr()+scale_fill_manual(values = Color1)+theme(text=element_text(family="serif",face = 'plain'), axis.text = element_text(face='plain'), axis.title.y=element_text(face='plain'),legend.title=element_blank(),legend.text=element_text(face='plain'))+labs(y='score')+scale_color_manual(values=Color)
p.lefse<-plot_ef_bar(buf)+theme_pubr()+scale_fill_manual(values = abb.Color1)+theme(text=element_text(family="arial",face = 'plain',size=10), axis.text = element_text(face='plain'), axis.title.y=element_text(face='plain'),legend.title=element_blank(),legend.text=element_text(face='plain'))+labs(y='')

#visualization diversity + lefse====
p.filter.diversity+p.filter.beta+p.lefse+plot_layout(widths=c(2.4,1.6,1.7))
#p.filter.diversity+p.filter.beta+plot_layout(heights = c(1,2))

####b.cereus count check----
buf<-mm_test.filter.b@marker_table[tof]
buf2<-buf$feature
buf2<-gsub('k__D','d__',buf2)

buf<-lefse.filter
buf1<-colSums(buf)
for (i in colnames(buf)){
  buf[i]<-buf[i]/buf1[i]
}


cereus.count<-data.frame(t(buf[buf2,]))
buf<-mm_test.filter.b@marker_table[tof]
for (i in 1:length(buf$feature)){
  buf$feature[i]<-strsplit(buf$feature[i],'s__')[[1]][2]
}
colnames(cereus.count)<-buf$feature
cereus.count$type<-patient.info[patient.info$Run %in% rownames(cereus.count),'condition']
cereus.count$type<-factor(ifelse(cereus.count$type=='Control','HC','EM'),levels=c('HC','EM'))
cereus.count.melt<-melt(cereus.count,id.vars = 'type')
cereus.count.melt<-cereus.count.melt[cereus.count.melt$variable %in%]

p.lefse+ggplot(cereus.count.melt, aes(x=type,y=value,fill=type,color=type))+geom_boxplot()+
  facet_wrap(~variable,scales = 'free')+theme_pubr()+scale_fill_manual(values=abb.Color1)+scale_color_manual(values=abb.Color)+
  stat_compare_means(method='kruskal.test',label='p',family='arial',label.y.npc = 0.8)+
  theme(text = element_text(family='arial',size=12),axis.title.x = element_blank(),axis.text.x = element_blank())+
  labs(color='',fill='',y='relative abundance')+plot_layout(ncol=1)


#immune cell data processing----
make.timer_sample.csv<-function(){
  buf1<-read.csv('/data2/msi/endometriosis/endometriosis/gene_count_matrix.csv') # by stringtie
  buf1<-buf1[buf1$gene_id!="<class 'str'>",]
  buf1['g']<-sapply(buf1$gene_id,FUN = function(x){
    y<-unlist(strsplit(x,'[|]'))
    if (is.na(y[2])){return(y[1])}
    else{return(y[2])}
  })
  
  buf3<-aggregate(buf1[,2:25], by=list(g=buf1$g), FUN=sum)
  rownames(buf3)<-buf3$g
  buf3<-buf3[,-1]
  buf3<-buf3[!grepl("^MT-", rownames(buf3)), ]
  write.csv(buf3,'/data2/msi/endometriosis/endometriosis/timer_sample.csv')
}

##all preprosess====
rna.expression<-read.csv('/data2/msi/endometriosis/endometriosis/timer_sample.csv')
rownames(rna.expression)<-rna.expression$X
rna.expression<-rna.expression[,-1]


#normalization
buf<-rna.expression
buf<-DESeqDataSetFromMatrix(buf,total.data,~condition)
buf<-estimateSizeFactors(buf)
rna.expression<-counts(buf,normalized=T)

write.csv(rna.expression,'/data2/msi/endometriosis/endometriosis/timer_sample_normalization.csv')

##filter preprosess====
rna.expression<-read.csv('/data2/msi/endometriosis/endometriosis/timer_sample.csv')
rownames(rna.expression)<-rna.expression$X
rna.expression<-rna.expression[,-1]

rna.expression<-rna.expression[,total.data.s.filter$Run]
write.csv(rna.expression,'/data2/msi/endometriosis/endometriosis/timer_sample.filter.csv')

#normalization
buf<-rna.expression
buf<-DESeqDataSetFromMatrix(buf,total.data.s.filter,~condition)
buf<-estimateSizeFactors(buf)
rna.expression<-counts(buf,normalized=T)

write.csv(rna.expression,'/data2/msi/endometriosis/endometriosis/timer_sample.filter_normalization.csv')

#immune cell data====
rna.expression.filter<-read.csv('/data2/msi/endometriosis/endometriosis/timer_sample.filter.csv',row.names = 'X')
rna.expression.filter.normalization<-read.csv('/data2/msi/endometriosis/endometriosis/timer_sample.filter_normalization.csv',row.names = 'X')

##cibersort(교수님 픽!)====
###timer2.0====
immune.ovary.based<-read.csv('/data2/msi/endometriosis/endometriosis/estimation_ovary_all_normalization.csv')
immune.ovary.filter<-read.csv('/data2/msi/endometriosis/endometriosis/estimation_ovary.filter_normalization.csv')

table(unlist(lapply(immune.ovary.filter$cell_type,FUN = function(x){unlist(strsplit(x, split = "_"))[2]})))
cibersort<-unlist(lapply(immune.ovary.filter$cell_type,FUN = function(x){unlist(strsplit(x, split = "_"))[2]=='CIBERSORT'}))
cibersort.abs<-unlist(lapply(immune.ovary.filter$cell_type,FUN = function(x){unlist(strsplit(x, split = "_"))[2]=='CIBERSORT-ABS'}))
epic<-unlist(lapply(immune.ovary.filter$cell_type,FUN = function(x){unlist(strsplit(x, split = "_"))[2]=='EPIC'}))
mcpcounter<-unlist(lapply(immune.ovary.filter$cell_type,FUN = function(x){unlist(strsplit(x, split = "_"))[2]=='MCPCOUNTER'}))
quantiseq<-unlist(lapply(immune.ovary.filter$cell_type,FUN = function(x){unlist(strsplit(x, split = "_"))[2]=='QUANTISEQ'}))
timer<-unlist(lapply(immune.ovary.filter$cell_type,FUN = function(x){unlist(strsplit(x, split = "_"))[2]=='TIMER'}))
xcell<-unlist(lapply(immune.ovary.filter$cell_type,FUN = function(x){unlist(strsplit(x, split = "_"))[2]=='XCELL'}))
cell_type<-unlist(lapply(immune.ovary.filter$cell_type,FUN = function(x){unlist(strsplit(x, split = "_"))[1]}))
cell_type<-gsub(' ','_',cell_type) #rowname mapping 할때 쓰임 (###1,2) (중복된 index가 존재해서 바로 rowname에 넣지 못함)
immune.ovary.filter<-immune.ovary.filter[,-1]

###corr with b.cereus----
Fun<-startsWith(colnames(total.data.s.filter),'Fun')
Vir<-startsWith(colnames(total.data.s.filter),'Vi')
Bac<-startsWith(colnames(total.data.s.filter),'Bac')
Arc<-startsWith(colnames(total.data.s.filter),'Arc')
Euk<-startsWith(colnames(total.data.s.filter),'Euk')
buf<-total.data.s.filter
buf<-buf[Fun | Vir |Bac | Euk | Arc]
buf<-buf/rowSums(buf)
buf$Run<-total.data.s.filter$Run
buf<-buf[c('Run',"Bacteria_s_Bacillus cereus")]

get('immune.ovary.filter.cibersort.abs')<-c()
method=c('cibersort','cibersort.abs','epic','mcpcounter','quantiseq','timer','xcell')
i=1
buf1<-immune.ovary.filter[get(method[i]),]
rownames(buf1)<-cell_type[get(method[i])]
buf1<-data.frame(t(buf1))
buf1$Run<-rownames(buf1)
buf1<-merge(total.data,buf1,by='Run',all.y=T)

immune.filter.cibersort<-merge(buf1,buf,by='Run',all.y=T)#1
immune.filter.cibersort.abs<-merge(buf1,buf,by='Run',all.y=T)#2
immune.filter.epic<-merge(buf1,buf,by='Run',all.y=T)#3
immune.filter.mcpcounter<-merge(buf1,buf,by='Run',all.y=T)#4
immune.filter.quantiseq<-merge(buf1,buf,by='Run',all.y=T)#5
immune.filter.timer<-merge(buf1,buf,by='Run',all.y=T)#6
immune.filter.xcell<-merge(buf1,buf,by='Run',all.y=T)#7


buf$type<-total.data.s.filter$tissue
buf$type<-ifelse(buf$type=='Control','Control',ifelse(buf$type=='Endometrium','Patient_Endometrium','Endometriosis'))

#b.cereus의 정상과 환자사이의 endometrium에서의 비교
ggplot(buf[buf$type!='Control',],aes(x=type,y=`Bacteria_s_Bacillus cereus`,fill=type))+geom_boxplot()+theme_pubr()+
  theme(text=element_text(family="arial",size=12),legend.position = 'none')+scale_fill_manual(values=Color3)+stat_compare_means(label = 'p.format')+
  labs(x='',y='relative abundance (%)')

library(Hmisc)
library(ggplot2)
library(ggpubr)
library(reshape2)

f.cor.data<-function(df){
  corr.raw.data<-df[4:10,12:ncol(df)]
  corr.matrix<-rcorr(as.matrix(corr.raw.data),type='spearman')
  cor.df.cibersort.filter<-data.frame(x='',y='',r=0,p=0)[-1,]
  
  for (i in colnames(corr.matrix$r)){
    for (j in rownames(corr.matrix$r))
      cor.df.cibersort.filter[nrow(cor.df.cibersort.filter)+1,]<-c(i,j,corr.matrix$r[i,j],corr.matrix$P[i,j])
  }
  cor.df.cibersort.filter$r<-as.numeric(cor.df.cibersort.filter$r)
  cor.df.cibersort.filter$p<-as.numeric(cor.df.cibersort.filter$p)
  
  cor.df.cibersort.filter<-cor.df.cibersort.filter[order(cor.df.cibersort.filter$r),]
  
  buf<-cor.df.cibersort.filter
  buf<-buf[buf$x=='Bacteria_s_Bacillus cereus',]
  buf<-buf[buf$y!='Bacteria_s_Bacillus cereus',]
  #buf$p.sym<-ifelse(buf$p<0.1,round(buf$p,3),'ns')
  buf$p.sym<-ifelse(buf$p<0.1,round(buf$p,3),round(buf$p,3))
  buf$p.sym<-ifelse(buf$p.sym==0,'<0.001',buf$p.sym)
  buf$type<-ifelse(buf$r>0,'po','ne')
  buf$Y<-gsub('_',' ',buf$y)
  buf<-buf[order(buf$r,decreasing = F),]
  #buf$Y<-factor(buf$Y,level=buf$Y)
  return(buf)
}
corr.df<-data.frame()
for (i in 1:length(method)){
  buf<-f.cor.data(get(paste0('immune.filter.',method[i])))
  buf$method<-method[i]
  if (nrow(corr.df)==0){corr.df<-buf}
  else{corr.df<-rbind(corr.df,buf)}
}

corr.df<-corr.df[!is.na(corr.df$r),]
corr.df<-corr.df[order(corr.df$r,decreasing = T),]
corr.df$sig<-ifelse(corr.df$p<0.1,corr.df$y,'')

  

ggplot(corr.df, aes(y=-log10(p),x=r,color=method,label=sig))+geom_point()+theme_pubr()+
  geom_text_repel(nudge_y = 1) 

buf<-corr.df[corr.df$method=='epic',]
buf<-buf[order(buf$r),]
buf$Y<-factor(buf$Y,levels = buf$Y)
ggplot(buf, aes(y=Y,x=r,fill=type))+geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values=c('red','blue'))+
  theme_pubr()+theme(text=element_text(family="arial"),legend.position = 'none')+geom_text(aes(label=p.sym,x=-r/abs(r)/12,y=Y))+
  labs(y='')


buf1<-immune.filter.cibersort[c(3,5,12:ncol(immune.filter.cibersort))]
immune.melt<-melt(buf1,id.vars = c('condition','tissue'))
immune.melt$variable<-gsub('_',' ',immune.melt$variable)
immune.melt<-immune.melt[immune.melt$variable!='Bacteria s Bacillus cereus',]
immune.melt$variable<-factor(immune.melt$variable,levels=buf$Y)
ggplot(immune.melt,aes(x=condition,y=value,fill=condition))+geom_boxplot()+facet_wrap(~variable,scales = 'free',ncol=6)+stat_compare_means(label='p',label.y.npc = 0.9)+theme_pubr()+scale_fill_manual(values=Color)


#corr with gene====
library(DESeq2)

{
  ensg<-c()
  gene<-c()
  buf1<-read.csv('/data2/msi/endometriosis/endometriosis/gene_count_matrix.csv') # by stringtie
  buf1<-buf1[buf1$gene_id!="<class 'str'>",]
  for (x in buf1$gene_id){
    y<-unlist(strsplit(x,'[|]'))
    if (is.na(y[2])){
      ensg<-c(ensg,y[1])
      gene<-c(gene,y[1])}
    else{
      ensg<-c(ensg,y[1])
      gene<-c(gene,y[2])}
  }
  buf1$ensg<-ensg
  buf1$gene<-gene
  buf1<-buf1[c(26,27,2:25)]
  buf1<-buf1[!grepl("^MT-", buf1$gene), ]
  table(rowSums(buf1[3:26]))
  rna.expression.ensg<-buf1
} #make.timer_sample_ensg


rna.expression.ensg.filter<-rna.expression.ensg[c('ensg','gene',total.data.s.filter$Run)]

##by ensg code====


#reference가 factor레벨에서 높아야함 (여기서는 control)
total.data.s.filter$condition<-as.factor(total.data.s.filter$condition)
buf<-rna.expression.ensg.filter
buf<-buf[-1]
buf<-aggregate(.~gene,data=buf,sum)

rownames(buf)<-buf$gene
buf<-buf[-1]
buf<-buf[rowSums(buf==0)<8,]

buf1<-data.frame(rownames(buf))
buf1<-merge(buf1,rna.expression.ensg.filter[c('ensg','gene')],all.x=T,by.x='rownames.buf.',by.y='ensg')

dds.ensg<-DESeqDataSetFromMatrix(countData = buf, colData = total.data.s.filter, design = ~condition)

dds.ensg<-DESeq(dds.ensg)
dds.ensg.res<-results(dds.ensg, alpha=0.05)

plotMA(dds.ensg)

EnhancedVolcano(dds.ensg.res,
                lab = rownames(dds.ensg.res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 2,
                title='Endometriosis vs Control')

buf<-dds.ensg.res[abs(dds.ensg.res$log2FoldChange)>2 & !is.na(dds.ensg.res$padj),]
buf<-buf[buf$padj<0.05,]
for (i in rownames(buf)){
  if (buf[i,'log2FoldChange']>0){show(i)}
}

##by gene name----
###remove sncRNA version(dds.cereus)
buf<-rna.expression.ensg.filter
buf$gene[startsWith(buf$gene,'RPS')]

buf<-buf[!startsWith(buf$gene,'RNU'),]
buf<-buf[!startsWith(buf$gene,'SNOR'),]
buf<-buf[!startsWith(buf$gene,'Y_'),]
buf<-buf[!startsWith(buf$gene,'SCARNA'),]
buf<-buf[!startsWith(buf$gene,'RPS'),]
buf<-buf[!startsWith(buf$gene,'RPL'),]


for (i in rownames(buf)){
  buf[i,'gene']<-strsplit(buf[i,'gene'],'-')[[1]][1]
}

buf<-aggregate(buf[-c(1,2)], by=list(gene=buf$gene),FUN=sum)
rownames(buf)<-buf$gene
buf<-buf[-1]
buf<-buf[rowSums(buf)!=0,]


dds.cereus<-DESeqDataSetFromMatrix(countData = buf, colData = total.data.s.filter, design = ~condition)
dds.cereus<-DESeq(dds.cereus)
dds.res.cereus<-results(dds.cereus,alpha=0.05)

EnhancedVolcano(dds.res.cereus,
                lab = rownames(dds.res.cereus),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 2,
                title='Endometriosis vs Control', subtitle = 'Volcano Plot')+theme(text=element_text(family = 'serif'))

gene_list<-dds.res.cereus$log2FoldChange
names(gene_list)<-rownames(dds.res.cereus)
gene_list<-gene_list[order(gene_list,decreasing = T)]


##
{vegf<-'VEGFA
KDR
SH2D2A
PLCG1
PLCG2
PRKCA
PRKCB
PRKCG
SPHK1
SPHK2
HRAS
KRAS
NRAS
RAF1
MAP2K1
MAP2K2
MAPK1
MAPK3
PLA2G4E
PLA2G4A
JMJD7-PLA2G4B
PLA2G4B
PLA2G4C
PLA2G4D
PLA2G4F
PPP3CA
PPP3CB
PPP3CC
PPP3R1
PPP3R2
NFATC2
PTGS2
PTK2
SHC2
PXN
CDC42
MAPK11
MAPK12
MAPK13
MAPK14
MAPKAPK2
MAPKAPK3
HSPB1
SRC
PIK3CA
PIK3CD
PIK3CB
PIK3R1
PIK3R2
PIK3R3
RAC1
RAC2
RAC3
AKT1
AKT2
AKT3
NOS3
CASP9
BAD'
  vegf<-unlist(strsplit(vegf,'\n'))
}#vegf관련 gene


dds.res.df<-as.data.frame(dds.res)
dds.res.df[is.na(dds.res)]<-1
dds.res.df<-dds.res.df[dds.res.df$log2FoldChange>2 & dds.res.df$padj<0.05,]

for (i in rownames(dds.res.df)){
  show(i)
}
dds.res.df<-as.data.frame(dds.res)
dds.res.df[is.na(dds.res)]<-1
dds.res.df<-dds.res.df[dds.res.df$log2FoldChange>2 & dds.res.df$padj<0.05,]

#corr with bacteria====
buf<-rna.expression.ensg.filter
buf<-buf[rowSums(buf[3:12])!=0,]
buf[nrow(buf)+1,]<-c(0,0,total.data.s.filter$`Bacteria_s_Bacillus cereus`)
buf[nrow(buf),1]<-'bacteria'
buf[nrow(buf),2]<-'Bacillus cereus'

rownames(buf)<-buf$ensg
buf<-buf[c(-1,-2)]
correlation.gene.filter<-rcorr(t(buf),type = 'spearman')

cor.df.gene.filter<-data.frame(x='',y='',r=0,p=0)
cor.df.gene.filter<-cor.df.gene.filter[-1,]
for (i in colnames(correlation.gene.filter$r)){
    cor.df.gene.filter[nrow(cor.df.gene.filter)+1,]<-c(i,'Bacillus cereus',correlation.gene.filter$r[i,'bacteria'],correlation.gene.filter$P[i,'bacteria'])
}
cor.df.gene.filter$r<-as.numeric(cor.df.gene.filter$r)
cor.df.gene.filter$p<-as.numeric(cor.df.gene.filter$p)

##vlookup start##
buf<-rna.expression.ensg.filter[c('ensg','gene')]
rownames(buf)<-rna.expression.ensg.filter$ensg
buf<-buf[cor.df.gene.filter$x,]
cor.df.gene.filter$x<-buf$gene
##vlookup end##
cor.df.gene.filter<-cor.df.gene.filter[order(cor.df.gene.filter$r),]
buf<-cor.df.gene.filter

buf<-buf[buf$p<0.05,]

ggplot(buf,aes(x=x,y=r,col=r))+geom_point()+scale_color_gradient2(low = "red", high = "blue",mid='gray',limits=c(-1,1))+theme(axis.text.x=element_blank())

buf1<-buf[buf$r>0,]
buf2<-buf[buf$r<0,]
write.csv(buf1,'~/positive.csv')
write.csv(buf2,'~/negative.csv')

#GO with B.cereus====
buf<-as.data.frame(dds.res)
buf<-buf[buf1$x,]

dds.res.df<-dds.res.df[order(dds.res.df$log2FoldChange,decreasing = T),]
gene_list<-dds.res.df$log2FoldChange
names(gene_list)<-rownames(dds.res.df)

gene_list[names(gene_list) %in% buf1$x]

library(org.Hs.eg.db)
gse.cereus<-gseGO(geneList=gene_list[names(gene_list) %in% buf1$x], 
                 ont ="ALL", 
                 keyType = "ALIAS", 
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Hs.eg.db, 
                 pAdjustMethod = "none")



require(DOSE)
p1<-dotplot(gse.cereus, showCategory=7,title='GO',split='ONTOLOGY') + facet_grid(ONTOLOGY~.sign, scale="free")+scale_fill_gradient(low='red',high='#ffafaf') + theme(axis.text.y=element_text(size=10),text = element_text(family = 'serif'))

ids<-bitr(names(gene_list), fromType = "ALIAS", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
df <- dds.res.cereus[rownames(dds.res.cereus)%in%ids$ALIAS,]
ids<-ids[ids$ALIAS %in% rownames(df),]
table(rownames(df)==ids$ALIAS)
df$Y = ids[which(rownames(df),ids$ALIAS),'ENTREZID']


df<-dds.res.cereus[rownames(dds.res.cereus)%in%buf1$x,]
kegg_gene_list <- df$log2FoldChange
names(kegg_gene_list) <- ids$ENTREZID
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
p2<-dotplot(kk2, showCategory=7,title='KEGG') + facet_grid(.~.sign, scale="free")+scale_fill_gradient(low='red',high='#ffafaf') + theme(axis.text.y=element_text(size=10),text = element_text(family = 'serif'))

#Gene enrichment analysis====
##Gene Ontology====
library('org.Hs.eg.db', character.only = TRUE)


dds.ensg.res.df<-dds.ensg.res.df[order(dds.ensg.res.df$log2FoldChange,decreasing = T),]
gene_list<-dds.ensg.res.df$log2FoldChange
names(gene_list)<-rownames(dds.ensg.res.df)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

dds.res.df<-dds.res.df[order(dds.res.df$log2FoldChange,decreasing = T),]

gene_list<-dds.res.df[rownames(dds.res.df)%in%buf1$x,]$log2FoldChange
names(gene_list)<-rownames(dds.res.df)

gse.alias<-gseGO(geneList=gene_list, 
                 ont ="ALL", 
                 keyType = "ALIAS", 
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Hs.eg.db, 
                 pAdjustMethod = "none")

require(DOSE)
p1<-dotplot(gse, showCategory=7,title='GO',split='ONTOLOGY') + facet_grid(ONTOLOGY~.sign, scale="free")+scale_fill_gradient(low='red',high='#ffafaf') + theme(axis.text.y=element_text(size=10),text = element_text(family = 'serif'))


##kegg pathway====
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df <- dds.res.df[ids$ALIAS,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df$Y = ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

p2<-dotplot(kk2, showCategory=7,title='KEGG') + facet_grid(.~.sign, scale="free")+scale_fill_gradient(low='red',high='#ffafaf') + theme(axis.text.y=element_text(size=10),text = element_text(family = 'serif'))


p1+p2+plot_layout(widths=c(2.5,1.5)) #horizenal
p1/p2+plot_layout(heights =c(5,1)) #vertical

kk2
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)

ridgeplot(kk2)

save.image(used.data)
#single cell----
#seurat previous ver.4 install (not ver.5)
#install.packages('Seurat', repos = c('https://satijalab.r-universe.dev', 'https://cloud.r-project.org'))
#remotes::install_version(package = 'Seurat', version = package_version('4.4.0'))
#remotes::install_version(package = 'SeuratObject', version = package_version('4.1.4'))
#seurat
library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
#doublet
library(DoubletFinder)
library(ROCR)
#singleR
library(SingleR)
library(celldex)
#visualization
library(ggplot2)
library(ggpubr)
library(reshape2)
library(stringr)


##data preprocessing 1(base paper)====
used.data<-'/data2/msi/Rdata/EMT.scRNA.RData' #based paper data
load(used.data)

path<-'/data2/msi/endometriosis/endometriosis/raw_data/matrix/'
file.list<-list.files(path)
file.list<-file.list[endsWith(file.list,'h5')]

unlist(strsplit(file.list[1],'_'))
mr<-function(x){
  cell.count<-ncol(x)
  if(cell.count<1000){
    return(0.008*cell.count/1000)
  }else if(cell.count<2000){
    return(0.007*(cell.count-1000)/1000+0.008)
  }else if(cell.count<3000){
    return(0.008*(cell.count-2000)/1000+0.015)
  }else if(cell.count<4000){
    return(0.007*(cell.count-3000)/1000+0.023)
  }else if(cell.count<5000){
    return(0.008*(cell.count-4000)/1000+0.03)
  }else if(cell.count<6000){
    return(0.008*(cell.count-5000)/1000+0.038)
  }else if(cell.count<7000){
    return(0.007*(cell.count-6000)/1000+0.046)
  }else if(cell.count<8000){
    return(0.008*(cell.count-7000)/1000+0.053)
  }else if(cell.count<9000){
    return(0.007*(cell.count-8000)/1000+0.061)
  }else if(cell.count<10000){
    return(0.012*(cell.count-9000)/1000+0.068)
  }else{
    return(0.013*(cell.count-10000)/1000+0.08)
  }
}
paramSweep_v3_Seurat5<-function (seu, PCs = 1:10, sct = FALSE, num.cores = 1) {
  require(Seurat)
  require(fields)
  pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.3, by = 0.05)
  min.cells <- round(nrow(seu@meta.data)/(1 - 0.05) - nrow(seu@meta.data))
  pK.test <- round(pK * min.cells)
  pK <- pK[which(pK.test >= 1)]
  orig.commands <- seu@commands
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 
                                                 10000, replace = FALSE)]
    data <- seu@assays$RNA$counts[, real.cells]
    n.real.cells <- ncol(data)
  }
  if (nrow(seu@meta.data) <= 10000) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA$counts
    n.real.cells <- ncol(data)
  }
  if (num.cores > 1) {
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)), FUN = parallel_paramSweep_v3, 
                        n.real.cells, real.cells, pK, pN, data, orig.commands, 
                        PCs, sct, mc.cores = num.cores)
    stopCluster(cl)
  }
  else {
    output2 <- lapply(as.list(1:length(pN)), FUN = parallel_paramSweep_v3, 
                      n.real.cells, real.cells, pK, pN, data, orig.commands, 
                      PCs, sct)
  }
  sweep.res.list <- list()
  list.ind <- 0
  for (i in 1:length(output2)) {
    for (j in 1:length(output2[[i]])) {
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, 
                                  sep = "_"))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}
doubletFinder_v3_SeuratV5 <- function (seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, 
                                       sct = FALSE, annotations = NULL){
  require(Seurat)
  require(fields)
  require(KernSmooth)
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    return(seu)
  }
  if (reuse.pANN == FALSE) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA$counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating", n_doublets, "artificial doublets...", 
                sep = " "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    if (!is.null(annotations)) {
      stopifnot(typeof(annotations) == "character")
      stopifnot(length(annotations) == length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    orig.commands <- seu@commands
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets, normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method, 
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor, 
                                     margin = orig.commands$NormalizeData.RNA@params$margin)
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets, 
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method, 
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span, 
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max, 
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function, 
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function, 
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin, 
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method, 
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures, 
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff, 
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                                 model.use = orig.commands$ScaleData.RNA$model.use, 
                                 do.scale = orig.commands$ScaleData.RNA$do.scale, 
                                 do.center = orig.commands$ScaleData.RNA$do.center, 
                                 scale.max = orig.commands$ScaleData.RNA$scale.max, 
                                 block.size = orig.commands$ScaleData.RNA$block.size, 
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                              npcs = length(PCs), rev.pca = orig.commands$RunPCA.RNA$rev.pca, 
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, 
                              verbose = FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                 ncol = 1))
    if (!is.null(annotations)) {
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                             ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if (!is.null(annotations)) {
        for (ct in unique(annotations)) {
          neighbor_types[i, ] <- table(doublet_types1[neighbors - 
                                                        n_real.cells]) + table(doublet_types2[neighbors - 
                                                                                                n_real.cells])
          neighbor_types[i, ] <- neighbor_types[i, ]/sum(neighbor_types[i, 
          ])
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data), 
                                                                    1]
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    if (!is.null(annotations)) {
      colnames(neighbor_types) = levels(doublet_types1)
      for (ct in levels(doublet_types1)) {
        seu@meta.data[, paste("DF.doublet.contributors", 
                              pN, pK, nExp, ct, sep = "_")] <- neighbor_types[, 
                                                                              ct]
      }
    }
    return(seu)
  }
}
singlecell<-function(x){
  pbmc.data<-Read10X_h5(paste0('/home/user02/raw_data1/seongik/GSE179640/matrix/',x))
  buf<-unlist(strsplit(x,'_'))
  name<-paste0(buf[2],buf[3])
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = name)
  #mitochondria and hemoglobin gene expressions more than 10% of total genes
  #pbmc<-PercentageFeatureSet(pbmc,pattern='^HB[^(P)]',col.name = 'percent.hb') #없는 gene 제외
  pbmc<-PercentageFeatureSet(pbmc,pattern='^MT-',col.name='percent.mt')
  #pbmc[['percent.mthb']]<-pbmc[['percent.mt']]+pbmc[['percent.hb']]
  # + 100<RNAcount<20000
  pbmc <- subset(pbmc, subset= nCount_RNA > 1000 & nCount_RNA<100000 & percent.mt < 25)
  
  # standard deviation lower than 1 per cell
  #pbmc[['sd']]<-apply(pbmc@assays$RNA@counts,2,function(x){sd(x)})
  #pbmc<-subset(pbmc,subset=sd<1)
  
  # 1st normalization
  pbmc<-NormalizeData(pbmc)
  pbmc<-FindVariableFeatures(pbmc,selection.method = 'vst')
  pbmc<-ScaleData(pbmc)
  #pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
  pbmc<-RunPCA(pbmc)
  
  #remove doublet
  #sweep.res.list_myData<-paramSweep_v3_Seurat5(pbmc,PCs=1:20,sct=F)
  sweep.res.list_myData<-paramSweep_v3(pbmc,PCs=1:20,sct=F)
  sweep.stats_myData<-summarizeSweep(sweep.res.list_myData,GT=F)
  bcmvn_myData<-find.pK(sweep.stats_myData)
  bcmaxcol<-which.max(bcmvn_myData$BCmetric)
  pk<-as.numeric(attr(bcmvn_myData[bcmaxcol,'pK'],'levels')[bcmaxcol])
  m<-mr(pbmc)
  
  #doublet
  nExp_poi<-round(m*nrow(pbmc@meta.data))
  #pbmc<-doubletFinder_v3_SeuratV5(pbmc,PCs = 1:20,pN = 0.25,pK = pk,nExp = nExp_poi,reuse.pANN = F,sct=F)
  pbmc<-doubletFinder_v3(pbmc,PCs = 1:20,pN = 0.25,pK = pk,nExp = nExp_poi,reuse.pANN = F,sct=F)
  colnames(pbmc@meta.data)[ncol(pbmc@meta.data)]<-'doublet'
  pbmc<-subset(pbmc, doublet=='Singlet')
  
  #saveRDS(pbmc,file=paste0('/data2/msi/rds/',name,'.rds'))
  
  #return(c(pbmc,m))
  return(pbmc)
}

{
c1<-singlecell(file.list[1])
c2<-singlecell(file.list[2])
c3<-singlecell(file.list[3])
e1e<-singlecell(file.list[4])
e2e<-singlecell(file.list[5])
e3e<-singlecell(file.list[6])
e4e<-singlecell(file.list[7])
e5e<-singlecell(file.list[8])
e6e<-singlecell(file.list[10])
e7e<-singlecell(file.list[13])
e8e<-singlecell(file.list[14])
e9e<-singlecell(file.list[16])
e1p<-singlecell(file.list[20])
e2p<-singlecell(file.list[22])
e3p<-singlecell(file.list[24])
e4p<-singlecell(file.list[26])
e5p<-singlecell(file.list[28])
e6p<-singlecell(file.list[9])
e7p<-singlecell(file.list[12])
e9p<-singlecell(file.list[30])
e1pa<-singlecell(file.list[21])
e2pa<-singlecell(file.list[23])
e3pa<-singlecell(file.list[25])
e4pa<-singlecell(file.list[27])
e5pa<-singlecell(file.list[29])
e9pa1<-singlecell(file.list[31])
e9pa2<-singlecell(file.list[32])
e7o<-singlecell(file.list[11])
e9o<-singlecell(file.list[15])
e10o<-singlecell(file.list[33])
e11o<-singlecell(file.list[17])
eor1<-singlecell(file.list[18])
eor3<-singlecell(file.list[19])
} #each sample preprocessing



mergelist<-c(c2,c3,
             e1e,e2e,e3e,e4e,e5e,e6e,e7e,e8e,e9e,
             e1p,e2p,e3p,e4p,e5p,e6p,e7p,e9p,
             e1pa,e2pa,e3pa,e4pa,e5pa,e9pa1,e9pa2,
             e7o,e9o,e10o,e11o)

emt.merge<-merge(c1, y=mergelist,add.cell.id=as.character(1:31)) # raw 111178 -> preprocessing 106910

rm(c1,c2,c3,e1e,e2e,e3e,e4e,e5e,e6e,e7e,e8e,e9e,e1p,e2p,e3p,e4p,e5p,e6p,e7p,e9p,e1pa,e2pa,e3pa,e4pa,e5pa,e9pa1,e9pa2,e7o,e9o,e10o,e11o,mergelist)

{
##data preprocessing 2(another paper)====

used.data<-'/data2/msi/Rdata/EMT.scRNA1.RData' #other paper val data
load(used.data)

path<-'/data2/msi/endometriosis/singlecell1_GSE203191/h5_matrix/'
file.list<-list.files(path)
file.list<-file.list[endsWith(file.list,'h5')]

mr<-function(x){
  cell.count<-ncol(x)
  if(cell.count<1000){
    return(0.008*cell.count/1000)
  }else if(cell.count<2000){
    return(0.007*(cell.count-1000)/1000+0.008)
  }else if(cell.count<3000){
    return(0.008*(cell.count-2000)/1000+0.015)
  }else if(cell.count<4000){
    return(0.007*(cell.count-3000)/1000+0.023)
  }else if(cell.count<5000){
    return(0.008*(cell.count-4000)/1000+0.03)
  }else if(cell.count<6000){
    return(0.008*(cell.count-5000)/1000+0.038)
  }else if(cell.count<7000){
    return(0.007*(cell.count-6000)/1000+0.046)
  }else if(cell.count<8000){
    return(0.008*(cell.count-7000)/1000+0.053)
  }else if(cell.count<9000){
    return(0.007*(cell.count-8000)/1000+0.061)
  }else if(cell.count<10000){
    return(0.012*(cell.count-9000)/1000+0.068)
  }else{
    return(0.013*(cell.count-10000)/1000+0.08)
  }
}

singlecell<-function(p,x){
  pbmc.data<-Read10X_h5(paste0(p,x))
  buf<-unlist(strsplit(x,'_'))
  name<-paste0(buf[2],buf[3])
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = name)
  #mitochondria and hemoglobin gene expressions more than 10% of total genes
  #pbmc<-PercentageFeatureSet(pbmc,pattern='^HB[^(P)]',col.name = 'percent.hb') #없는 gene 제외
  pbmc<-PercentageFeatureSet(pbmc,pattern='^MT-',col.name='percent.mt')
  #pbmc[['percent.mthb']]<-pbmc[['percent.mt']]+pbmc[['percent.hb']]
  pbmc$log10GenesPerUMI<-log10(pbmc$nFeature_RNA)/log10(pbmc$nCount_RNA)
  pbmc$GenesPerUMI<-pbmc$nFeature_RNA/pbmc$nCount_RNA
  
  # + 100<RNAcount<20000
  pbmc <- subset(pbmc, subset= nCount_RNA > 500 & nCount_RNA<50000 & percent.mt < 10 & nFeature_RNA < 6000)
  
  # standard deviation lower than 1 per cell
  #pbmc[['sd']]<-apply(pbmc@assays$RNA@counts,2,function(x){sd(x)})
  #pbmc<-subset(pbmc,subset=sd<1)
  
  # 1st normalization
  pbmc<-NormalizeData(pbmc)
  pbmc<-FindVariableFeatures(pbmc,selection.method = 'vst')
  pbmc<-ScaleData(pbmc)
  #pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
  pbmc<-RunPCA(pbmc)
  
  #remove doublet
  #sweep.res.list_myData<-paramSweep_v3_Seurat5(pbmc,PCs=1:20,sct=F)
  sweep.res.list_myData<-paramSweep_v3(pbmc,PCs=1:20,sct=F)
  sweep.stats_myData<-summarizeSweep(sweep.res.list_myData,GT=F)
  bcmvn_myData<-find.pK(sweep.stats_myData)
  bcmaxcol<-which.max(bcmvn_myData$BCmetric)
  pk<-as.numeric(attr(bcmvn_myData[bcmaxcol,'pK'],'levels')[bcmaxcol])
  m<-mr(pbmc)
  
  #doublet
  nExp_poi<-round(m*nrow(pbmc@meta.data))
  #pbmc<-doubletFinder_v3_SeuratV5(pbmc,PCs = 1:20,pN = 0.25,pK = pk,nExp = nExp_poi,reuse.pANN = F,sct=F)
  pbmc<-doubletFinder_v3(pbmc,PCs = 1:20,pN = 0.25,pK = pk,nExp = nExp_poi,reuse.pANN = F,sct=F)
  colnames(pbmc@meta.data)[ncol(pbmc@meta.data)]<-'doublet'
  pbmc<-subset(pbmc, doublet=='Singlet')
  
  #saveRDS(pbmc,file=paste0('/data2/msi/rds/',name,'.rds'))
  
  #return(c(pbmc,m))
  return(pbmc)
}

singlecell<-function(p,x){
  pbmc.data<-Read10X_h5(paste0(p,x))
  buf<-unlist(strsplit(x,'_'))
  name<-paste0(buf[2],buf[3])
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = name)
  #mitochondria and hemoglobin gene expressions more than 10% of total genes
  #pbmc<-PercentageFeatureSet(pbmc,pattern='^HB[^(P)]',col.name = 'percent.hb') #없는 gene 제외
  pbmc<-PercentageFeatureSet(pbmc,pattern='^MT-',col.name='percent.mt')
  return(pbmc)
} #not preprocessing
#seurat v5 experiment

{
  wME<-singlecell(path,file.list[1])
  ME1<-singlecell(path,file.list[2])
  ME2<-singlecell(path,file.list[3])
  ME3<-singlecell(path,file.list[4])
  ME4<-singlecell(path,file.list[5])
  ME5<-singlecell(path,file.list[6])
  ME6<-singlecell(path,file.list[7])
  ME7<-singlecell(path,file.list[8])
  w2<-singlecell(path,file.list[9])
  w3<-singlecell(path,file.list[10])
  w1<-singlecell(path,file.list[11])
} #data preprocessing

pbmclist<-c(ME1,ME2,ME3,ME4,ME5,ME6,ME7,w1,w2,w3)
me.list<-c('wME','ME1','ME2','ME3','ME4','ME5','ME6','ME7','w1','w2','w3')
emt.merge<-merge(wME, y=pbmclist, add.cell.ids = me.list)

buf<-GetAssayData(emt.merge, slot="counts", assay="RNA")
buf<-buf[rowSums(buf)>2,]
emt.merge<-CreateSeuratObject(buf)

cell.info<-read.csv('/data2/msi/endometriosis/singlecell1_GSE203191/h5_matrix/gse203191_endo_meta_frame.csv')
table(cell.info$run)
cell.info$temp<-gsub('tissue','',cell.info$run)
cell.info$temp<-gsub('holeME','',cell.info$temp)
cell.info$temp<-gsub('_','',cell.info$temp)
cell.info$cellname<-paste0(cell.info$temp,'_',cell.info$barcode)

emt.merge$cellname<-colnames(emt.merge)
emt.merge<-subset(emt.merge, subset= cellname %in% cell.info$cellname)

rownames(cell.info)<-cell.info$cellname
cell.info<-cell.info[colnames(emt.merge),]
table(rownames(cell.info)==colnames(emt.merge))

emt.merge<-PercentageFeatureSet(emt.merge,pattern='^MT-',col.name = 'percent.mt')
emt.merge$type<-cell.info$pheno
emt.merge$condition<-cell.info$prep
emt.merge$clusterID<-cell.info$clusterID

emt.merge<-subset(emt.merge,subset= percent.mt<10)

VlnPlot(emt.merge,features = c('nCount_RNA','nFeature_RNA','percent.mt'),group.by = 'orig.ident')


##data structure====

p.sc.feature1 <- FeatureScatter(emt.merge, feature1 = "nCount_RNA", feature2 = "percent.mt",raster = F)
p.sc.feature2 <- FeatureScatter(emt.merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster=F)
p.sc.feature1+p.sc.feature2

emt.merge<-NormalizeData(emt.merge)
emt.merge<-FindVariableFeatures(emt.merge, selection.method = "vst")
all.genes <- rownames(emt.merge)
emt.merge <- ScaleData(emt.merge, features = all.genes)

BoxPlot<-function(obj,mk,group.by='seurat_cluster',x.label='type',comparisons=list()){
  data<-data.frame(obj@assays$RNA@counts)
  show('not contain gene: ')
  show(mk[!(mk %in% rownames(data))])
  mk<-mk[mk %in% rownames(data)]
  if(length(mk)==1){buf<-data[mk,]}
  else{buf<-data.frame(t(data[mk,]))}
  if (group.by=='seurat_cluster'){buf[x.label]<-obj@active.ident }
  else{buf[x.label]<-obj[[group.by]]}
  rm(data)
  buf<-melt(buf,id.vars=x.label)
  return(ggplot(buf,aes(x=get(x.label),y=value,fill=get(x.label)))+geom_boxplot()+facet_wrap(~variable,scales='free')+stat_compare_means(label='p.format',method = 'wilcox.test',comparisons = comparisons)+theme_pubr()+labs(x='',y='count',fill='')+theme(legend.position="none"))
}
} #버려진 데이터터

#emt.merge[['type']]<-ifelse(startsWith(emt.merge@meta.data$orig.ident,'C'),'Control','Endometriosis')

emt.merge <- RunPCA(emt.merge, features = VariableFeatures(emt.merge),pcs.compute=30,seed.use=12345)
VizDimLoadings(emt.merge, dims = 1:2, reduction = "pca")
DimPlot(emt.merge, reduction = "pca",raster = F)

emt.merge <- FindNeighbors(emt.merge, dims = 1:20)
emt.merge <- FindClusters(emt.merge,resolution=3)
emt.merge <- RunTSNE(emt.merge,dims=1:30,seed.use=12345)
emt.merge <- RunUMAP(emt.merge,dims=1:30)
DimPlot(emt.merge, reduction = "tsne",raster = F,label=T)+DimPlot(emt.merge, reduction = "umap",raster = F,label=T)

Idents(emt.merge)<-'clusterID'
DimPlot(emt.merge, reduction = "umap",raster = F,label=T)+DimPlot(emt.merge, reduction = "umap",raster = F,label=T)

FeaturePlot(emt.merge,features = c('CD3D','CD3E','CD4','CD8A'),raster=F)+DimPlot(emt.merge,reduction = 'umap',label=T,raster=F)


emt.marker<-FindAllMarkers(emt.merge1,logfc.threshold = 0.25)

Idents(emt.merge) = "RNA_snn_res.3"





#data1(base paper)----

emt.merge[['type1']]<-emt.merge$type
emt.merge[['type1']]<-ifelse(endsWith(emt.merge$orig.ident,'sis'),'Endometrium',emt.merge$type1)
emt.merge[['type1']]<-ifelse(endsWith(emt.merge$orig.ident,'P'),'Peritoneum',emt.merge$type1)
emt.merge[['type1']]<-ifelse(endsWith(emt.merge$orig.ident,'A'),'Peritoneum_adg',emt.merge$type1)
emt.merge[['type1']]<-ifelse(endsWith(emt.merge$orig.ident,'O'),'Ovary',emt.merge$type1)

DimPlot(emt.merge, reduction = "umap",raster = F,label=F,group.by = 'type1')+DimPlot(emt.merge, reduction = "umap",raster = F,label=F,group.by = 'type')+DimPlot(emt.merge, reduction = "umap",raster = F,label=T)

#https://www.nature.com/articles/s41467-022-33982-7
markers<-c('CD4','CD8A','CD8B') #T cell
markers<-c('EPCAM','CLDN3','SMIM22',
           'MUC1','C2orf88','CDH1') #epithelial
markers<-c('VWF','CLEC14A', 'EMCN',
           'ADGRL4','ECSCR','SOX18') #endothelial
markers<-c('CD3E', 'CD7',
           'IKZF3','TBC1D10C') #lymphocyte
markers<-c('AIF1','CYBB','CLEC7A',
           'SPI1','MARCH1') #myeloid (30,70 cluster를 여기 넣어야 하나??)
markers<-c('PDGFRA','PDGFRB','ISLR',
           'MXRA8','EMILIN1','PRRX1') #stromal
markers<-c('CD3E', 'CD7', 'IKZF3',
           'PDGFRA','PDGFRB','ISLR',
           'EPCAM','CLDN3','SMIM22',
           'AIF1','CYBB','CLEC7A',
           'VWF','CLEC14A', 'EMCN'
           )
markers<-c('COL1A1','COL3A1','COL1A2')

DotPlot(emt.merge,features = markers)+theme_pubr()+theme(axis.text.x=element_text(angle=90, hjust=1),text = element_text(family='serif'),legend.position = 'right')+labs(x='Gene',y='Cell Type')
FeaturePlot(emt.merge,markers[c(3,5,7,10,13)],raster=F,ncol=5)
DimPlot(emt.merge,raster=F,label=T)
VlnPlot(emt.merge,markers,raster=F)

cluster<- c('lymphocyte','stromal','stromal','lymphocyte','stromal', #0~4
            'stromal','epithelial','myeloid','endothelial','lymphocyte', #5~9
            'stromal','lymphocyte','lymphocyte','stromal','stromal', #10~14
            'stromal','endothelial','endothelial','myeloid','epithelial', #15~19
            'epithelial','stromal','lymphocyte','stromal','myeloid', #20~24
            'epithelial','epithelial','epithelial','stromal','stromal', #25~29
            'lymphocyte','stromal','stromal','endothelial','stromal', #30~34
            'stromal','stromal','myeloid','stromal','endothelial', #35~39
            'epithelial','myeloid','stromal','stromal','lymphocyte', #40~44
            'stromal','epithelial','myeloid','','endothelial', #45~49
            'epithelial','epithelial','stromal','myeloid','lymphocyte', #50~54
            'epithelial','stromal','stromal','lymphocyte','lymphocyte', #55~59
            'myeloid','lymphocyte','myeloid','','epithelial', #60~64
            'lymphocyte','epithelial','myeloid','epithelial','', #65~69
            'lymphocyte','','stromal','epithelial','', #70~74
            'epithelial','lymphocyte' #75~76
)
cluster<-ifelse(cluster=='','unknown',cluster)
names(cluster) <- levels(emt.merge)
emt.merge <- RenameIdents(emt.merge, cluster)

DimPlot(emt.merge,raster = F,label=T)
Idents(emt.merge) = "RNA_snn_res.3"

emt.T<-emt.merge[,emt.merge@active.ident=='lymphocyte']

emt.T<-FindVariableFeatures(emt.T, selection.method = "vst")
all.genes <- rownames(emt.T)
emt.T <- ScaleData(emt.T, features = all.genes)

emt.T <- RunPCA(emt.T, features = VariableFeatures(emt.T),pcs.compute=30,seed.use=12345)

emt.T <- FindNeighbors(emt.T, dims = 1:20)
emt.T <- FindClusters(emt.T,resolution=6)
emt.T <- RunTSNE(emt.T,dims=1:30,seed.use=12345)
emt.T <- RunUMAP(emt.T,dims=1:30,seed.use=12345)

DimPlot(emt.T,label=T,raster=F)
DimPlot(emt.T,label=T,raster=F,reduction='tsne')

FeaturePlot(emt.T,features=c('temp'))
emt.T[['temp']]<-ifelse(emt.T$RNA_snn_res.6%in%c('0'),1,0)

FeaturePlot(emt.T,markers,ncol=2)

cluster<- c('CD8 T','NK','CD4 naive','CTL','CD8 T', #0~4
            'NK','CD4 memory','B cell','CTL','CD8 T', #5~9
            'CD8 T','CTL','CD4 memory','CD8 T','CD4 naive', #10~14
            'B cell','B cell','T reg','T reg','CD4 memory', #15~19
            'CD8 T','NK','CD4 memory','CD8 T','NK', #20~24
            'plasma','NK','NK','CD8 T','CD8 T', #25~29
            'NK','CD4 naive','CTL','NK','CD8 T', #30~34
            'CD8 T','CD4 memory','CD8 T','CD8 T','NK', #35~39
            'ILC','NK','','CD8 T','CD8 T', #40~44
            'CD8 T','CD8 naive','T reg','NK','CD8 T', #45~49
            '','NK','CD4 memory','','NK', #50~54
            'CD8 T','NK','NK','CD8 T','CD4 naive', #55~59
            'CD8 naive','NK','NK','','CD8 T', #60~64
            'NK','plasma','B cell','B cell','CD8 T', #65~69
            'CD8 T','','CD8 T','B cell','CD8 T' #70~74
)
cluster<-ifelse(cluster=='','unknown',cluster)
names(cluster) <- levels(emt.T)
emt.T<-RenameIdents(emt.T,cluster)
DimPlot(emt.T,label=T)
Idents(emt.T)='RNA_snn_res.6'
table(emt.T$RNA_snn_res.6)

markers<-list('NK'=c('ITGAX','NCAM1'), #nk
           'T reg'=c('FOXP3','RTKN2'), #t reg
           'T cell'=c('CD4','CD8B'), #cd4/cd8
           'memory'='CD40LG', #memory
           'naive'=c('CCR7','LEF1'), #naive cell
           'CTL'=c('FGFBP2','SPON2'), #ctl
           'ILC'=c('IL7R','KIT'), #ilc
           'B cell'=c('CD40','BANK1'), #B cell
           'plasma'=c('IGKC','FCRL5') #plasma
)

{
markers<-c("IL7R",
    "CD40LG",
    "DUSP1",
    "FOS",
    "JUN",
    "TSC22D3",
    "TNFAIP3",
    "SPOCK2",
    "RORA",
    "CD69",
    "RGS1",
    "ERN1",
    "NFKBIA",
    "CXCR4",
    "GPR183",
    "GPR171",
    "AC058791.1",
    "CD4",
    "S100A4",
    "BTG1",
    "AQP3",
    "CD2",
    "KLRB1",
    "TNFRSF25",
    "FKBP5",
    "PBXIP1",
    "FOSB",
    "TOB1",
    "DPP4",
    "PRDM1",
    "JAML",
    "TRAC",
    "KLF6",
    "TNFSF13B",
    "CRYBG1",
    "PDE4D",
    "CD6",
    "CD3D",
    "CCR6",
    "LTB",
    "PIK3IP1",
    "SERINC5",
    "PTGER4",
    "CD3G",
    "LGALS3",
    "TXNIP",
    "CCL20",
    "ANXA1",
    "CDC14A",
    "SARAF",
    "ZFP36L2",
    "LEPROTL1",
    "PDE3B",
    "IFI44",
    "CAMK4",
    "ODF2L",
    "CXCR6",
    "JUNB",
    "ANKRD28",
    "RPLP1",
    "MAF",
    "TPT1",
    "INPP4B",
    "S100A6",
    "BCL2",
    "PAG1",
    "VIM",
    "NFKBIZ",
    "MYBL1",
    "NABP1",
    "AC016831.5",
    "AAK1",
    "ELK3",
    "RCAN3",
    "TMEM173",
    "SLAMF1",
    "CITED2",
    "OXNAD1",
    "BCL11B",
    "TRAT1",
    "S100A11",
    "EML4",
    "AHNAK",
    "TTC39C",
    "B3GALT2",
    "PPP2R5C",
    "SYNE2",
    "IFNGR1",
    "TC2N",
    "ADAM19",
    "SLC4A7",
    "RUNX2",
    "RNF19A",
    "DGKA",
    "PPP1R15A",
    "IL2",
    "EMB",
    "FKBP11",
    "TNF",
    "MGAT4A",
    "TNFAIP8",
    "IL18R1",
    "CD5",
    "ITM2B",
    "CD226",
    "S100A10",
    "PERP",
    "RPLP0",
    "PDCD4",
    "IKZF1",
    "ANKRD12",
    "LPIN2",
    "CD28",
    "IL4I1",
    "VMP1",
    "CLEC2D",
    "CCNH",
    "CYTIP",
    "SCML4",
    "RGS2",
    "CTSB",
    "KIAA1551",
    "GOLGA8A",
    "FLT3LG",
    "CTLA4",
    "GOLGA8B",
    "ANK3",
    "LINC00513",
    "TIMP1",
    "ARHGAP15",
    "ABRACL",
    "IFI44L",
    "THEM4",
    "KDSR",
    "AL357060.1",
    "KIF5C",
    "CREM",
    "TLE4",
    "TRBC2",
    "FXYD5",
    "ITK",
    "MALAT1",
    "CERK",
    "AC114760.2",
    "LINC01871",
    "MCL1",
    "VSIR",
    "TRADD",
    "PRKCA",
    "TRBC1",
    "PLEC",
    "GSTK1",
    "FTH1",
    "CDKN1B",
    "GALM",
    "PTGER2",
    "GPR155",
    "OCIAD2",
    "EEF1A1",
    "MAL",
    "SORL1",
    "DYNLT3",
    "RETREG1",
    "TCF7",
    "MCAM",
    "CD44",
    "SLC4A10",
    "NCK2",
    "CHRM3-AS2",
    "KLF9",
    "PTPN13",
    "ATF7IP2",
    "LTK",
    "LINC00649",
    "TNFSF8",
    "PBX4",
    "DUSP16",
    "SYTL2",
    "RPS12",
    "RNF157",
    "RGS10",
    "PPP1R2",
    "FRMD4B",
    "RNASET2",
    "AC020916.1",
    "JUND",
    "CIB1",
    "NMRK1",
    "BORCS5",
    "GABPB1-AS1",
    "TESPA1",
    "SENP7",
    "ACP6",
    "SESN1",
    "RNF149",
    "CCR2",
    "SIT1",
    "ZC3HAV1",
    "CCDC107",
    "CORO1B",
    "PFKFB3",
    "GPRIN3",
    "HIST1H1D",
    "PHTF2",
    "SUSD3",
    "RGS14",
    "SATB1",
    "CHD7",
    "CTSH",
    "PRNP",
    "FBLN7",
    "ANTXR2",
    "SRSF7",
    "PCNX1",
    "LAT",
    "CTSL",
    "SEPT6",
    "CD3E",
    "MPRIP",
    "AC020571.1",
    "TGFBR2",
    "STK17B",
    "SNHG8",
    "THEMIS",
    "RPL10",
    "FYB1",
    "RPS4X",
    "RPL12",
    "N4BP2L2",
    "ZC3H8",
    "CD46",
    "TMEM63A",
    "CCNL1",
    "RFLNB",
    "AL096865.1",
    "NOSIP",
    "AP3M2",
    "LDHB",
    "ACSL6",
    "TMEM123",
    "PDE4B",
    "ZFP36L1",
    "GALNT10",
    "RPL34",
    "CYB561",
    "RPL32",
    "RPL39",
    "MYC",
    "MZT2A")
} #cd4 memory
{
  markers<-c("CD8A",
             "GZMK",
             "CD8B",
             "CD3G",
             "CD3D",
             "TRGC2",
             "CD2",
             "GZMH",
             "KLRG1",
             "CCL4L2",
             "TRAC",
             "IL32",
             "CCL5",
             "CXCR6",
             "MT2A",
             "TRBC2",
             "CD27",
             "LAG3",
             "HMGB2",
             "CCL4",
             "THEMIS",
             "S100A6",
             "SIT1",
             "STMN1",
             "CD84",
             "LYAR",
             "VIM",
             "ITM2A",
             "PTMS",
             "LCK",
             "S100A10",
             "SH2D1A",
             "CD3E",
             "GZMM",
             "KIAA1551",
             "SH3BGRL3",
             "DUSP2",
             "PPP2R5C",
             "RARRES3",
             "LINC01871",
             "CD99",
             "GIMAP4",
             "ITM2C",
             "CYTOR",
             "TUBA4A",
             "PTTG1",
             "GRAP2",
             "CLEC2D",
             "CCR5",
             "F2R",
             "GPR174",
             "CISH",
             "PTPRC",
             "AHNAK",
             "SLF1",
             "PSMB9",
             "LAT",
             "CRTAM",
             "TMSB4X",
             "TC2N",
             "CXCR4",
             "IKZF3",
             "MYL12B",
             "HLA-C",
             "APOBEC3G",
             "CD6",
             "PPDPF",
             "CD69",
             "LGALS3",
             "ANXA2",
             "C12orf75",
             "WNK1",
             "CXCR3",
             "S100A4",
             "HCST",
             "BCL11B",
             "ARL4C",
             "RGS1",
             "CKLF",
             "MIAT",
             "LINC00649",
             "ANXA1",
             "RNF213",
             "HLA-F",
             "SYNE2",
             "CD52",
             "RGS10",
             "PDCD1",
             "MAF",
             "DTHD1",
             "TRBC1",
             "TUBB",
             "GUK1",
             "MT1X",
             "GIMAP1",
             "SLAMF7",
             "GIMAP7",
             "HLA-A",
             "FYB1",
             "NIN",
             "KLF6",
             "ANXA5",
             "ATP5F1E",
             "CST7",
             "TNIP3",
             "PTP4A2",
             "JPT1",
             "TRAT1",
             "RAB27A",
             "PLA2G16",
             "ACTG1",
             "CD48",
             "ATXN1",
             "CCL3L1",
             "ISG15",
             "APBB1IP",
             "LYST",
             "TIGIT",
             "JAML",
             "B2M",
             "ITGA4",
             "EMB",
             "EMP3",
             "STK17A",
             "PSMA7",
             "A2M-AS1",
             "APOBEC3C",
             "CD44",
             "RGL4",
             "CRIP1",
             "TRG-AS1",
             "ITGAL",
             "KIF21A",
             "MT1F",
             "GPRIN3",
             "ACTB",
             "TMA7",
             "DEK",
             "IL10RA",
             "CLEC2B",
             "TMEM50A",
             "SH3KBP1",
             "LINC01934",
             "OASL",
             "GPR171",
             "HLA-B",
             "LEPROTL1",
             "H2AFZ",
             "MYL12A",
             "UBL5",
             "MFSD10",
             "TPRG1",
             "IDH2",
             "SLC2A4RG",
             "MYH9",
             "HMGB1",
             "UBE2S",
             "FABP5",
             "STK17B",
             "CIB1",
             "OPTN",
             "ITGB2",
             "DUT",
             "CRYBG1",
             "BLOC1S1",
             "CAPN2",
             "CD5",
             "CYTIP",
             "H2AFV",
             "MGST3",
             "PRDM1",
             "PFN1",
             "PTGER4",
             "ARHGAP26",
             "DNAJC1",
             "CDC42EP3",
             "GNG2",
             "INPP4B",
             "NDUFS5",
             "SEPT1",
             "RBL2",
             "BICDL1",
             "CDC42SE2",
             "CEMIP2",
             "CALM1",
             "TOB1",
             "GBP1",
             "OST4",
             "CFL1",
             "COX8A",
             "NLRC5",
             "LCP1",
             "PDLIM2",
             "SRSF7",
             "ZNF683",
             "KMT2A",
             "RASAL3",
             "ARHGAP18",
             "GSTK1",
             "SUB1",
             "TERF2IP",
             "HIST1H1D",
             "ST8SIA1",
             "TAP1",
             "ATP1B1",
             "RHOH",
             "TPM3",
             "ARHGAP25",
             "CYBA",
             "RNF19A",
             "HMGN2",
             "CSTB",
             "ATP5MD",
             "CTSC",
             "CASP1",
             "PPP2R2B",
             "FLNA",
             "CAPZB",
             "MRPL41",
             "CALM2",
             "MYO1G",
             "SYNE1",
             "PITPNC1",
             "TUBB4B",
             "OCIAD2",
             "MCUB",
             "PTGES3",
             "HLA-E",
             "NELFCD",
             "IFNGR1",
             "DRAP1",
             "COPS9",
             "SSBP4",
             "ZNF101",
             "CORO1A",
             "PLEC",
             "PYHIN1",
             "CCDC167",
             "ABRACL",
             "REEP5",
             "ANKRD12",
             "TXN",
             "EZR",
             "TMPO",
             "CLIC1",
             "COX7A2",
             "MXD4",
             "DYNLT1",
             "C9orf78",
             "ARPC2",
             "S100A11",
             "ANXA2R",
             "ZBTB38",
             "FGFR1OP2",
             "NFATC2",
             "HSP90AA1",
             "BTN3A2",
             "PAXX",
             "LINC00861",
             "RUNX3",
             "PTPN7",
             "ATP6V0E2"
  )
} #cd8 t
{
  markers<-c(
    "FGFBP2",
    "SPON2",
    "FCGR3A",
    "KLRF1",
    "PLAC8",
    "CLIC3",
    "CX3CR1",
    "PRF1",
    "EFHD2",
    "KLRD1",
    "KLF2",
    "NKG7",
    "AKR1C3",
    "PTGDS",
    "PLEK",
    "ADGRG1",
    "GZMB",
    "CST7",
    "ZEB2",
    "ARL4C",
    "ITGB2",
    "CD247",
    "CEP78",
    "S1PR5",
    "FCER1G",
    "TTC38",
    "C1orf162",
    "PTPN12",
    "ABHD17A",
    "CHST2",
    "MYOM2",
    "RAP1B",
    "CEBPD",
    "GPR65",
    "TYROBP",
    "GZMH",
    "CD160",
    "C1orf21",
    "GZMM",
    "GSAP",
    "BIN2",
    "PRSS23",
    "PYHIN1",
    "LITAF",
    "TBX21",
    "SYNE1",
    "SH3BP5",
    "MYBL1",
    "IFITM2",
    "RIPOR2",
    "TXK",
    "CYBA",
    "MYO1F",
    "KLRB1",
    "MBP",
    "ITGAM",
    "JAK1",
    "GK5",
    "XBP1",
    "GNG2",
    "TFDP2",
    "PTGDR",
    "HSH2D",
    "CD300A",
    "TPST2",
    "CCND3",
    "APMAP",
    "CTBP2",
    "LINC00861",
    "CDC42SE1",
    "TGFBR3",
    "MATK",
    "HOPX",
    "IL18RAP",
    "SH2D1B",
    "TXNIP",
    "CMC1",
    "FLNA",
    "MAPK1",
    "IGFBP7",
    "LAIR2",
    "METRNL",
    "SPN",
    "LPCAT1",
    "C1orf56",
    "KLF3",
    "FCRL6",
    "CCL3",
    "UCP2",
    "SAMD3",
    "PRMT2",
    "NCR1",
    "MTSS1",
    "ADRB2",
    "HIPK2",
    "ATM",
    "SYTL3",
    "RHOBTB3",
    "UTRN",
    "HAVCR2",
    "AES",
    "ID2",
    "PLEKHG3",
    "IGF2R",
    "IRF1",
    "CD226",
    "THEMIS2",
    "OSBPL5",
    "PXN",
    "PTPN4",
    "ICAM2",
    "GLRX",
    "AOAH",
    "DDIT4",
    "MSN",
    "CYTH1",
    "CALM1",
    "SERPINB6",
    "YPEL1",
    "FGR",
    "PTPRE",
    "GFOD1",
    "GIMAP7",
    "ZBTB16",
    "ADD3",
    "STK38",
    "LYAR",
    "EIF3G",
    "MLC1",
    "HLA-E",
    "DBI",
    "VAV3",
    "CD47",
    "EMP3",
    "PDGFD",
    "SYNE2",
    "CAST",
    "FCRL3",
    "PCSK5",
    "PLEKHA1",
    "DIP2A",
    "AGTPBP1",
    "STK10",
    "MYL12A",
    "LINC02384",
    "CXXC5",
    "USP28",
    "LYN",
    "RASA3",
    "C12orf75",
    "DOK2",
    "ITGAL",
    "FAM49B",
    "SYTL1",
    "TES",
    "NDUFB7",
    "GPATCH8",
    "TCF25",
    "DENND2D",
    "IQGAP2",
    "NEIL1",
    "HLA-F",
    "DAB2",
    "PFN1",
    "RAB29",
    "SAMHD1",
    "SIGLEC7",
    "AKNA",
    "HCST",
    "SSBP3",
    "PIP4K2A",
    "ABI3",
    "RAP1GAP2",
    "RAP2B",
    "FEZ1",
    "HLA-C",
    "NCR3",
    "TMSB10",
    "ADGRE5",
    "AREG",
    "ASCL2",
    "TLE4",
    "CFLAR",
    "CCL4",
    "RASGRP2",
    "APOL6",
    "MCTP2",
    "PPP2R5C",
    "RNF125",
    "ZAP70",
    "SORL1",
    "CCDC28B",
    "LAIR1",
    "UBB",
    "BTN3A2",
    "RAB11FIP1",
    "SLC9A3R1",
    "CARD16",
    "CTSW",
    "CDKN2D",
    "RNF213",
    "UBE2F",
    "DHRS7",
    "RUNX3",
    "PRDM1",
    "PRR5L",
    "PLEKHF1",
    "IFITM1",
    "VCL",
    "SH2D2A",
    "GIMAP4",
    "NEAT1",
    "MMP23B",
    "RASSF4",
    "HDDC2",
    "LIMD2",
    "ARHGEF3",
    "ARPC2",
    "TGFB1",
    "G6PD",
    "CEMIP2",
    "SRPK2",
    "OSTF1",
    "AKAP13",
    "MIDN",
    "CELF2",
    "CXCR2",
    "ZNF600",
    "RAB9A",
    "NFATC2",
    "PRKCH",
    "CDC42",
    "FAM107B",
    "S1PR1",
    "SERPINB1",
    "B4GALT4",
    "SUN2",
    "CYB5B",
    "SRGN",
    "MYO1G",
    "TSPAN32",
    "NDUFB2",
    "KLF13",
    "PTPRA",
    "ZNF276",
    "FUT11",
    "DSTN",
    "FYN",
    "ARPC5L",
    "HMGN3",
    "SELPLG",
    "LLGL2",
    "MTPN",
    "CNOT6L",
    "RNF166",
    "SLAMF7",
    "SYTL2",
    "PRR5",
    "RNF165",
    "SSH2",
    "CD244",
    "GAB3",
    "B2M",
    "YPEL3",
    "KLHDC4",
    "IDI1",
    "IQGAP1",
    "RARRES3",
    "LBH",
    "RNF115",
    "SLFN5",
    "PITPNC1",
    "GTF3C1",
    "RAC2",
    "PPM1L",
    "DGKZ",
    "FGL2",
    "CD53",
    "SCLT1",
    "TBXAS1",
    "MAPRE2",
    "ZCCHC17",
    "GBP5",
    "PILRB",
    "SLC15A4",
    "DIAPH1",
    "F2R",
    "PDIA3",
    "ZNF595",
    "MEX3C",
    "IFITM3",
    "NFATC3",
    "TRGC1",
    "PSMB9",
    "ARHGAP26",
    "DGKD",
    "TAOK3",
    "PPIL2",
    "PSMB8",
    "SLC44A2",
    "NFE2L2",
    "SLCO4C1",
    "RASSF1",
    "HES6",
    "NMUR1",
    "CD7",
    "GBP4",
    "CMKLR1",
    "ADCY7",
    "TLN1",
    "MVD",
    "ZBTB44",
    "SETD3",
    "PPIA",
    "AK5",
    "CISD3",
    "PPP1R18",
    "TNIK",
    "SIPA1",
    "MGAT4A",
    "PRKCB",
    "RORA",
    "MIEN1",
    "SH3BP2",
    "KLRG1",
    "NPC1",
    "SCP2",
    "CFL1",
    "LAT2",
    "APLP2",
    "VASP",
    "ENPP4",
    "SUGP2",
    "CMTR2",
    "UCHL5",
    "LINC00623",
    "HLA-A",
    "STARD3NL",
    "SERF2",
    "VPS37B",
    "HELZ",
    "AGTRAP",
    "BTN3A1",
    "SBK1",
    "ARPC5",
    "RHOC",
    "ANXA4",
    "POLR3GL",
    "CASP10",
    "CAP1",
    "IER2",
    "QKI",
    "TNFRSF1B",
    "SAMD9",
    "CD164",
    "ADAMTS1",
    "DTNBP1",
    "RBM38",
    "CCDC88C",
    "CLINT1",
    "PTGER2",
    "RAB37",
    "RPA2",
    "HMOX2",
    "SIRT2",
    "ADD1",
    "PCNX4",
    "NCOA1",
    "STX8",
    "CTDSP1",
    "CTBP1",
    "INPP4A",
    "RAMP1",
    "CCDC85B",
    "CD320",
    "PSME1",
    "SH3KBP1",
    "ITCH",
    "MIB2",
    "PTPN2",
    "ADGRG5",
    "ESYT2",
    "P2RY8",
    "MAP3K8",
    "FMNL1",
    "CAPNS1",
    "NCALD",
    "ANXA6",
    "RPS6KA1",
    "BRD7",
    "KCTD10",
    "RNF126",
    "FHL3",
    "ARHGAP30",
    "PIK3R5",
    "STK26",
    "SIGIRR",
    "TRAPPC10",
    "RRBP1",
    "YWHAZ",
    "GNPTAB",
    "GOLGA4",
    "PSMB10",
    "CEBPB",
    "RAB8A",
    "NDUFA12",
    "BZW1",
    "RPE",
    "S1PR4",
    "DTHD1",
    "TRDC",
    "LGALS9B",
    "SYAP1",
    "PKN1",
    "ARHGAP25",
    "PPP2R5A",
    "PTPN18",
    "LYST",
    "C1orf174",
    "RASA2",
    "PLCB1",
    "REL",
    "WSB1",
    "SLCO3A1",
    "CLIC1",
    "CD55",
    "ADAM8",
    "ERBB2",
    "CASP8",
    "ANXA1",
    "TSPOAP1",
    "GALC",
    "PAXX",
    "CPT1A",
    "LPGAT1",
    "EBP",
    "CHST12",
    "BST2",
    "ORAI1",
    "NPRL2",
    "ERBIN",
    "FOSL2",
    "MLLT6",
    "CEP126",
    "KCNAB2",
    "NFIL3",
    "RHBDF2",
    "PTPN6",
    "KIR2DL1",
    "PRKACB",
    "TCEAL8",
    "SLAMF6",
    "RAB8B",
    "HLA-B",
    "KIFAP3",
    "CDK2AP2",
    "AC245297.3",
    "TMEM71",
    "VPS4B",
    "PDZD4",
    "TNFRSF1A",
    "SEPT7",
    "ITGA6",
    "RIN3",
    "NOTCH1",
    "CTSC",
    "H3F3A",
    "AP005482.1",
    "DECR1",
    "YES1",
    "RHOG",
    "FAM43A",
    "NFAT5",
    "TLE1",
    "ADAM10",
    "C11orf21",
    "H3F3B",
    "TRAF3IP3",
    "RGS3",
    "RBM39",
    "ATOX1",
    "DOCK11",
    "RRAS2",
    "PAFAH2",
    "FNDC3B",
    "GIMAP1",
    "TET2",
    "FAM8A1",
    "KIR3DL1",
    "SLC27A3",
    "GLTP",
    "TROVE2",
    "GSN",
    "SELL",
    "DNAJC8",
    "RFTN1",
    "PPP1CA",
    "FASLG",
    "FKBP11",
    "PPP1R12A",
    "CRBN",
    "GATAD2A",
    "NPRL3",
    "IFI16",
    "CCDC69",
    "TENT5A",
    "C9orf139",
    "AP2A1",
    "DNAJB6",
    "ETS1",
    "SLFN13",
    "RGS19",
    "A2M-AS1",
    "AXIN1",
    "TBC1D10C",
    "CCDC82",
    "PSME2",
    "RASSF3",
    "CA5B",
    "XCL2",
    "FBXO6",
    "DYNLT1",
    "ARPC4",
    "SELENOT",
    "SRSF9",
    "PRDX5",
    "PHLDB2",
    "NUB1",
    "LFNG",
    "TERF1",
    "TFEB",
    "FBXW5",
    "VAMP5",
    "JADE1",
    "TSEN54",
    "TNFAIP8L2",
    "WDR45",
    "FLI1",
    "DUSP2",
    "TSPAN2",
    "LBR",
    "RSBN1L",
    "BCL9L",
    "MAT2B",
    "PPP2R2B",
    "SERBP1",
    "ZFPM1",
    "CDC42EP3",
    "WHAMM",
    "RTN4",
    "HENMT1",
    "VPS16",
    "ACTR3",
    "PCBP1",
    "EPB41",
    "MYLIP",
    "SASH3",
    "PAIP2",
    "TAGAP",
    "CPQ",
    "FOXN3",
    "GALNT3",
    "SLC20A1",
    "ARSG",
    "POLR2G",
    "DNAJC1",
    "TMEM181",
    "RAB27B",
    "ZNF22",
    "IFNGR1",
    "MRPL10",
    "CASP4",
    "CAPN15",
    "SPART",
    "ENPP5",
    "USP34",
    "PIK3AP1",
    "BANF1",
    "GNAI2",
    "PRELID1",
    "PRR13",
    "WDR37",
    "TNFRSF14",
    "ISG20",
    "PREX1",
    "HIPK1",
    "UHMK1",
    "RAB7A",
    "ACAA2",
    "ATXN7L3B",
    "SLC25A20",
    "ABCF1",
    "CCM2",
    "TINF2",
    "CYB561D2",
    "CSK",
    "FRG1",
    "CALCOCO1",
    "ADA",
    "ARNTL",
    "SUPT4H1",
    "NCL",
    "CALR",
    "RALBP1",
    "CD2BP2",
    "COLQ",
    "FES",
    "ILKAP",
    "TAP1",
    "STAT6",
    "LCP1",
    "STAT4",
    "SSX2IP",
    "MORC3",
    "CRAMP1"
  )
} #ctl
{markers<-c("AREG",
            "LST1",
            "IL4I1",
            "TNFSF13B",
            "PCDH9",
            "KIT",
            "IL7R",
            "FOS",
            "LTB",
            "ALDOC",
            "KRT86",
            "KRT81",
            "TNFRSF25",
            "TNFRSF4",
            "KLRB1",
            "NFKBIA",
            "LINC00299",
            "DUSP1",
            "ENPP1",
            "LIF",
            "GSN",
            "ZFP36L1",
            "H2AFY",
            "MDFIC",
            "BST2",
            "AHR",
            "IL1R1",
            "RORC",
            "PLAT",
            "MPG",
            "SPRY1",
            "CCL20",
            "PRR5",
            "AQP3",
            "TMIGD2",
            "AFF3",
            "CAT",
            "CXCR4",
            "IL23R",
            "B4GALNT1",
            "CD69",
            "NRIP1",
            "MFGE8",
            "CD81",
            "CTSH",
            "TCF7",
            "TLE1",
            "AC007952.4",
            "YWHAH",
            "FOSB",
            "TTC39C-AS1",
            "TOX2",
            "CD300LF",
            "JUN",
            "CITED2",
            "NCOA7",
            "TIMP1",
            "PHACTR2",
            "TNFAIP3",
            "C20orf204",
            "HERPUD1",
            "SLA",
            "SRI",
            "TNFRSF18",
            "DLL1",
            "IL18R1",
            "IGFBP4",
            "CA2",
            "PEBP1",
            "CXXC5",
            "TRBC1",
            "KIAA1324",
            "TSC22D3",
            "DDIT4",
            "CCR6",
            "IFI44",
            "MIS18BP1",
            "SCN1B",
            "SVIL",
            "SLC4A10",
            "ZG16B",
            "SKIL",
            "TEX30",
            "TGM2",
            "AC245014.3",
            "JAML",
            "APOL4",
            "TIAM1",
            "LTC4S",
            "TPT1",
            "TNFSF11",
            "SMAP1",
            "CD164",
            "CASP3",
            "SKAP2",
            "ANP32E",
            "FXYD5",
            "PPP1R9A",
            "BLVRA",
            "LTA4H",
            "ITM2B",
            "GOLGA8A",
            "TMEM123",
            "NFKB1",
            "IER2",
            "RGS1",
            "GNAS",
            "PLXNC1",
            "PTPN22",
            "ABHD15",
            "TNFSF12",
            "PNP",
            "JMY",
            "TP53I13",
            "FTH1",
            "PRPF6",
            "IKZF2",
            "ELK3",
            "ID2",
            "TNF",
            "NEO1",
            "SEC11A",
            "RORA",
            "HMGN3",
            "ERN1",
            "SLC16A3",
            "DCAF11",
            "DYNLT3",
            "DHRS3",
            "MPV17L",
            "CPNE7",
            "SPOCK2",
            "GPR183",
            "PTGER3",
            "COL4A4",
            "MGAT4A",
            "PFKFB3",
            "LMO4",
            "RPLP1",
            "RRBP1",
            "NCF4",
            "FURIN",
            "UNC93B1",
            "KIAA1211L",
            "METTL9",
            "CIDEB",
            "VWA5A",
            "TGIF1",
            "BTG2",
            "COL9A2",
            "RBPJ",
            "MAN1A1",
            "GRAMD2B",
            "CITED4",
            "DST",
            "STN1",
            "ARHGAP10",
            "RUNX2",
            "BEX2",
            "TMEM243",
            "REEP3",
            "RETREG1",
            "SPINK2",
            "RPS6KA3",
            "REPIN1",
            "RTKN2",
            "GOLGA8B",
            "NAA38",
            "B3GALT5",
            "CMTM6",
            "PLSCR1",
            "ATP10D",
            "BZW1",
            "MAF",
            "BEX4",
            "FES",
            "NCR3",
            "SOCS1",
            "NAP1L1",
            "YWHAQ",
            "OSTC",
            "DPP4",
            "RAB37",
            "SELENOS",
            "PLXDC2",
            "ENG",
            "SOST",
            "IFI44L",
            "CCDC88A",
            "FKBP11",
            "SCML4",
            "ZFAND5",
            "MATN2",
            "S100A13",
            "PLEKHA5",
            "CEP83",
            "LAPTM4B",
            "MBOAT7",
            "MGAT5",
            "CTNNB1",
            "ANP32A",
            "AL034397.3",
            "NECTIN2",
            "SOX18",
            "CFH",
            "MAML2",
            "CD2AP",
            "CENPV",
            "HDAC9",
            "AC103591.3",
            "PRAM1",
            "ADAM10",
            "PREX1",
            "ODF2L",
            "CERS6",
            "SEPT11",
            "FLOT1",
            "TMED8",
            "ADI1",
            "CCNG2",
            "FARS2",
            "TRDC",
            "RTN4",
            "SLC43A2",
            "HNRNPA0",
            "SUSD6",
            "EGR1",
            "HMGN1",
            "LTB4R",
            "PTGDR",
            "SCX",
            "ITGB1",
            "SPART",
            "BEX5",
            "NSMCE1",
            "OTUD5",
            "LTBP3",
            "FUCA2",
            "APP",
            "SKP1",
            "RALBP1",
            "CTR9",
            "GBP4",
            "ELOVL6",
            "LTA",
            "CIRBP",
            "TMEM14C",
            "ABRACL",
            "EEF1A1",
            "BCL6",
            "TBC1D31",
            "VAMP8",
            "EML4",
            "CHD7",
            "KMT5A",
            "TSPO",
            "SIGLEC7",
            "PGLS",
            "RIPK2",
            "SRSF5",
            "SH2D1B",
            "SSBP2",
            "IFNGR2",
            "MID1IP1",
            "SPTLC2",
            "ALAS1",
            "ANXA7",
            "DUSP6",
            "TNFSF4",
            "TMEM167A",
            "BIRC3",
            "RNASET2",
            "HSPB1",
            "PCID2",
            "CFAP36",
            "IFNGR1",
            "CSF2",
            "NUDT14",
            "SCRN1",
            "MT-ND6",
            "PBXIP1",
            "PER1",
            "EIPR1",
            "TOX",
            "HIC1",
            "HYI",
            "RHEB",
            "ABCB1",
            "DDX3X",
            "IL17RE",
            "IL2",
            "MTRNR2L12",
            "SNHG8",
            "PLPP1",
            "BHLHE40",
            "BCAS1",
            "CERK",
            "IL4R",
            "DOCK5",
            "CEP170",
            "MAP2K6",
            "HEXA",
            "CD52",
            "CHRM3-AS2",
            "SCARF1",
            "ST3GAL5",
            "ZNF667-AS1",
            "FAM189A2",
            "ATP10A",
            "TNS3",
            "FNBP1",
            "TMEM156",
            "SEM1",
            "AHI1",
            "MMP25",
            "GYG1",
            "PCMTD1",
            "CTDSP2",
            "PERP",
            "RAB11FIP1",
            "CCNI",
            "MAP3K8",
            "REV3L",
            "SYPL1",
            "AC245060.5",
            "PPP3CA",
            "SLC1A4",
            "LDLRAD4",
            "GRAMD1A",
            "MAP3K7",
            "MIF",
            "LAIR1",
            "NMT2",
            "HSP90AB1",
            "FUCA1",
            "HNRNPA1",
            "AC020916.1",
            "PGGHG",
            "CYSLTR1",
            "GPR171",
            "RAP1B",
            "ZBTB16",
            "ANKRD28",
            "EIF1AX",
            "PKDCC",
            "SLC25A39",
            "TPMT",
            "SPATS2L",
            "PRNP",
            "VDAC1",
            "RERE",
            "SLC40A1",
            "LPXN",
            "DNMT3A",
            "SEPT2",
            "ZNRF1",
            "IFT57",
            "AGPAT5",
            "TPI1",
            "ATP8B4",
            "MTHFS",
            "PPFIBP1",
            "TBL1XR1",
            "IRF2BP2",
            "KIFC3",
            "CBL",
            "HNRNPH3",
            "SLC25A36",
            "HOTAIRM1",
            "SYDE2",
            "ITGB7",
            "UBB",
            "RAPGEF6",
            "MICU3",
            "VPS28",
            "RPL8",
            "LAYN",
            "RASSF8",
            "GDE1",
            "SEPHS2",
            "KLHL5",
            "PTPN18",
            "SERP1",
            "SRGN",
            "FOXC1",
            "SPINT2",
            "LDLRAD3",
            "MAP2K1",
            "BAZ2B",
            "CKLF",
            "NOTCH2",
            "SMIM3",
            "HSBP1L1",
            "INSIG2",
            "DMTF1",
            "PABPC1",
            "SERTAD2",
            "ELOVL5",
            "KDM6B",
            "DEK",
            "N4BP2",
            "PLIN3",
            "PLEKHO1",
            "CD59",
            "NFATC1",
            "HIPK2",
            "KCTD9",
            "TCEAL8",
            "RGS9",
            "FYTTD1",
            "LZTS1",
            "SORCS1",
            "SH3TC1",
            "LEPROT",
            "TSPAN4",
            "ESPN",
            "GPX4",
            "VEGFA",
            "PPP1R15A",
            "ANXA1",
            "CHD9",
            "LPAR1",
            "NUDT7",
            "COMMD3",
            "ILK",
            "RPN2",
            "CYP26A1",
            "SPECC1",
            "CD40LG",
            "ATP6AP2",
            "RNF207",
            "AL356272.1",
            "PRMT9",
            "NRP1",
            "TMEM230",
            "WDR83OS",
            "ZNF217",
            "SESN1",
            "TMEM106B",
            "SDC4",
            "ADAM19",
            "ZCCHC2",
            "INPP4A",
            "SPSB3",
            "ESD",
            "TRGC1",
            "ETV6",
            "FSD1",
            "HEBP2",
            "LZTFL1",
            "CLEC11A",
            "H3F3A",
            "MGST2",
            "USP53",
            "SLC31A2",
            "SELENOT",
            "ABCC1",
            "ATP2B1-AS1",
            "H3F3B",
            "ISOC1",
            "SOX4",
            "AC022182.2",
            "PAPSS1",
            "CREG1",
            "APIP",
            "PAK2",
            "STK26",
            "INPP1",
            "TRAM2",
            "ARL3",
            "MRPL18",
            "STX17",
            "DENND5A",
            "RPL15",
            "MVB12B",
            "IRAK3",
            "CLN5",
            "NUTF2",
            "LSP1",
            "F2R",
            "GPR82",
            "PAWR",
            "RPS6KA5",
            "LINGO4",
            "PDZD2",
            "AP1S1",
            "EIF3E",
            "TOR1A",
            "TEC",
            "PHPT1",
            "CAMTA1",
            "OXNAD1",
            "GGT1",
            "CLNS1A",
            "TBXAS1",
            "SSH1",
            "WDFY2",
            "PLAUR",
            "KMT5B",
            "RNF146",
            "TSPYL1",
            "ST3GAL1",
            "NCR2",
            "HOXA7",
            "CDC14A",
            "MPST",
            "BIVM",
            "RGS16",
            "TMEM63A",
            "ECH1",
            "ANAPC16",
            "RASA2",
            "CLEC7A",
            "NME3",
            "SCCPDH",
            "SLC49A3",
            "RHOB",
            "HTATSF1",
            "M6PR",
            "SPN",
            "TANC1",
            "ACER3",
            "TALDO1",
            "PLEKHN1",
            "CBX5",
            "GGNBP2",
            "FAM241A",
            "LINC02256",
            "RTN3",
            "RAC1",
            "CHKA",
            "IFNAR1",
            "C1GALT1",
            "BTBD7",
            "CD96",
            "BRD1",
            "NFIL3",
            "PTPN9",
            "HOXA9",
            "RAB11FIP5",
            "AF117829.1",
            "IL12RB1",
            "MUM1",
            "NDUFS8",
            "SAMD12",
            "NUDT4",
            "GRINA",
            "DYNLL1",
            "TAGLN2",
            "MED30",
            "RCAN1",
            "CAMK1",
            "B3GAT2",
            "TBC1D17",
            "GPR65",
            "EIF1",
            "RPS24",
            "SEC14L1",
            "METTL21A",
            "GOLGA8R",
            "SGIP1",
            "SRSF2",
            "PIM2",
            "MAST4",
            "CNDP2",
            "LAMP1",
            "RALB",
            "GALNT1",
            "NFE2L3",
            "AL096865.1",
            "G3BP1",
            "MCAM",
            "DNAJB12",
            "TTC39C",
            "ZBED4",
            "HNRNPUL1",
            "PNRC1",
            "RANBP17",
            "GYPC",
            "DUSP22",
            "TXK",
            "PSMA1",
            "OCIAD1",
            "ARMCX2",
            "TMBIM6",
            "CLNK",
            "C1QTNF1",
            "COX7A2L",
            "TLN1",
            "MARCH3",
            "TGFBI",
            "ARHGAP5",
            "IL2RA",
            "TCIRG1",
            "CD200R1",
            "TESPA1",
            "NPM1",
            "TMED2",
            "ERO1A",
            "HSD17B4",
            "ARL14EP",
            "IFRD1",
            "USP47",
            "ABCA1",
            "NR1D1",
            "GTF2B",
            "HOXA10-AS",
            "COQ8A",
            "OLA1",
            "MAPKAPK3",
            "RPS4X",
            "PITPNM1",
            "LACTB",
            "ARF3",
            "MIB1",
            "RPL37A",
            "SEC62",
            "PRICKLE1",
            "ZNF609",
            "HSPD1",
            "RNF165",
            "CST3",
            "RNF152",
            "THOC7",
            "AAMDC",
            "DYNLRB1",
            "DAD1",
            "HOXB2",
            "MYO7A",
            "A4GALT",
            "TANC2",
            "AAK1",
            "PHTF2",
            "SRSF6",
            "RHOBTB3",
            "TFAM",
            "PYURF",
            "GAPDH",
            "RPL5",
            "FGD6",
            "MINDY2",
            "PRKG1",
            "CD33",
            "RPL22",
            "KDM5B",
            "ZNF652",
            "ARF4",
            "NR1D2",
            "TNFAIP8",
            "IMMT",
            "PSTPIP2",
            "MPZL3",
            "GPRC5C",
            "PHB",
            "RASSF8-AS1",
            "TCEAL4",
            "AKT1S1",
            "CD44",
            "MCF2L2",
            "FBXL3",
            "CDC42",
            "SSBP3",
            "UBAC2",
            "RAB32",
            "RACK1",
            "DYNC2H1",
            "SIK3",
            "OGFRL1",
            "MLXIP",
            "EPCAM",
            "HNRNPA2B1",
            "FAM84B",
            "EPS8L2",
            "ATF7IP2",
            "PPP2CA",
            "NBN",
            "DLD",
            "CRTAP",
            "LACC1",
            "CD63",
            "PRKAB1",
            "RIN3",
            "MTCL1",
            "ATP6V1G1",
            "WHRN",
            "FIBP",
            "ZNF462",
            "AC013264.1",
            "RNFT1",
            "IRF4",
            "MT-ND3",
            "NCL",
            "EPB41L2",
            "PPP4R1",
            "MKNK1",
            "TMEM220",
            "HNRNPK",
            "CAPN12",
            "GOLGA8N",
            "DACH1",
            "NDUFC2",
            "BTBD11",
            "HDLBP",
            "SMAD2",
            "EDF1",
            "EPDR1",
            "ECHS1",
            "NBL1",
            "INPP5B",
            "AC090152.1",
            "RHOC",
            "ZFYVE9",
            "AC109446.3",
            "ELL2",
            "PIK3CD",
            "LAT2",
            "YBX1",
            "ATP2B1",
            "ZYG11B",
            "Z93241.1",
            "CYFIP1",
            "ARL4A",
            "HOXA5",
            "SERPINB6",
            "RPA1",
            "KLF11",
            "TST",
            "PLCB1",
            "FAM169A",
            "CCNC",
            "DENND3",
            "TMEM9",
            "KDSR",
            "EXOSC6",
            "ATF5",
            "POC1B",
            "CAMK4",
            "TMEM59")} #ilc
{markers<-c("IGKC",
            "IGLC2",
            "IGLC3",
            "IGHG1",
            "IGHA1",
            "IGHG3",
            "IGHG2",
            "JCHAIN",
            "IGHG4",
            "IGLL5",
            "MZB1",
            "IGHM",
            "IGHGP",
            "IGHA2",
            "IGLV2-11",
            "SSR4",
            "XBP1",
            "DERL3",
            "PRDX4",
            "SEC11C",
            "HERPUD1",
            "FKBP11",
            "HSP90B1",
            "IGHV3-30",
            "IGKV4-1",
            "FCRL5",
            "FKBP2",
            "SSR3",
            "TNFRSF17",
            "IGLV6-57",
            "IGLV3-1",
            "SPCS3",
            "JSRP1",
            "ERLEC1",
            "BTG2",
            "CD79A",
            "SEL1L",
            "UBE2J1",
            "TENT5C",
            "RRBP1",
            "DNAJB9",
            "SPCS2",
            "SELENOK",
            "LMAN1",
            "MANF",
            "TP53INP1",
            "PIM2",
            "SDC1",
            "IGKV3-20",
            "POU2AF1",
            "SDF2L1",
            "PLPP5",
            "HSPA5",
            "MYDGF",
            "PTP4A3",
            "SELENOM",
            "ITM2C",
            "TXNDC11",
            "ELL2",
            "SELENOS",
            "TMEM59",
            "CRELD2",
            "LY96",
            "SEC61B",
            "SPCS1",
            "PDK1",
            "GAS6",
            "TXNDC15",
            "H1FX",
            "CLPTM1L",
            "SPAG4",
            "CREB3L2",
            "RPN2",
            "ST6GAL1",
            "CPEB4",
            "PDIA6",
            "PDIA4",
            "EAF2",
            "AC012236.1",
            "SEL1L3",
            "KDELR2",
            "TMEM258",
            "SIL1",
            "FAM30A",
            "ISG20",
            "TXNDC5",
            "SPATS2",
            "CPNE5",
            "CTHRC1",
            "PECAM1",
            "FTL",
            "HDLBP",
            "SERP1",
            "PPIB",
            "DNAAF1",
            "MEF2C",
            "KDELR1",
            "DNAJC3",
            "IGLC7",
            "MAN1A1",
            "LINC02362",
            "CKAP4",
            "MEI1",
            "KRTCAP2",
            "IGHV4-34",
            "IGHV1-69D",
            "CCPG1",
            "CTSH",
            "SEC14L1",
            "CYBA",
            "NPC2",
            "ALG5",
            "IRF4",
            "FAM92B",
            "CAV1",
            "HM13",
            "GNG7",
            "ARF4",
            "CADM1",
            "NANS",
            "RASSF6",
            "CD38",
            "RABAC1",
            "TPD52",
            "FBXW7",
            "RALGPS2",
            "ITGA6",
            "TRAM1",
            "TBCEL",
            "MARCKS",
            "PAIP2B",
            "SRPRB",
            "CHPF",
            "TRIB1",
            "GADD45A",
            "DDOST",
            "METTL7A",
            "U62317.5",
            "ARSA",
            "FNDC3B",
            "COBLL1",
            "ATF5",
            "CHID1",
            "SLAMF7",
            "LINC01480",
            "FNDC3A",
            "BRSK1",
            "HSPA13",
            "HIST1H1C",
            "DERL1",
            "IGF1",
            "SRM",
            "TMED9",
            "CLIC4",
            "QPRT",
            "TRAM2",
            "ERN1",
            "RPN1",
            "DENND5B",
            "RAB30",
            "SEC61G",
            "PSAP",
            "ZNF706",
            "CALU",
            "CCDC88A",
            "SOCS3",
            "BHLHE41",
            "ST6GALNAC4",
            "QPCT",
            "PNOC",
            "TOR3A",
            "AC078883.1",
            "FCGR2B",
            "VOPP1",
            "UAP1",
            "MANEA",
            "TMED10",
            "LMAN2",
            "SYVN1",
            "GMPPB",
            "SERPINI1",
            "TMEM208",
            "Z93930.2",
            "LRRC59",
            "SEC24A",
            "ANKRD28",
            "SAR1B",
            "PGM3",
            "NT5DC2",
            "CDK14",
            "TMEM205",
            "ERGIC3",
            "B4GALT7",
            "MDK",
            "KCNN3",
            "CEP128",
            "HIST1H2BG",
            "CCR10",
            "CHST15",
            "FICD")} #plasma


markers<-list('a'=c('KLRB1','CD81'),'b'=c('KIT','CCR6','CD69'),'c'=c('IL7R','IL8R1'))
emt.T$temp<-emt.T@active.ident
levels(emt.T)<-c('NK','T reg','CD4 memory','CD4 naive','CD8 naive','CD8 T','CTL','ILC','B cell','plasma','unknown')

DotPlot(emt.T,features = c(markers))+theme_pubr()+theme(text = element_text(family='serif'),legend.position = 'right')+labs(x='Gene',y='Cell Type')+RotatedAxis()

DotPlot(emt.T,features = c(markers))

buf<-FindMarkers(emt.T,ident.1 = 'plasma')
buf<-buf[buf$p_val_adj<0.05,]
buf<-buf[buf$p_val<0.05,]
buf1<-buf[buf$avg_log2FC>(0),]
buf1<-buf1[order(buf1$avg_log2FC,decreasing = T),]
for (i in rownames(buf1)){show(i)}

##visualization----
buf<-data.frame(table(emt.merge@active.ident[emt.merge$type=='Control']))
buf3<-data.frame(table(emt.merge@active.ident[emt.merge$type!='Control']))
buf['count']<-buf$Freq/sum(buf$Freq)
buf3['count']<-buf3$Freq/sum(buf3$Freq)
buf['type']<-'Control'
buf3['type']<-'Endometriosis'
count.all<-rbind(buf,buf3)
buf<-data.frame(table(emt.T@active.ident[emt.T$type=='Control']))
buf3<-data.frame(table(emt.T@active.ident[emt.T$type!='Control']))
buf['count']<-buf$Freq/sum(buf$Freq)
buf3['count']<-buf3$Freq/sum(buf3$Freq)
buf['type']<-'Control'
buf3['type']<-'Endometriosis'
count.tcell<-rbind(buf,buf3)
p.abundance.all<-ggplot(count.all,aes(x=type,y=count,fill=type))+geom_bar(stat = 'identity')+facet_wrap(~Var1,scales='free')+scale_fill_manual(values = Color)+theme_pubr()+labs(fill='',y='relative abundance')+theme(text=element_text(family='serif'))
p.abundance.tcell<-ggplot(count.tcell,aes(x=type,y=count,fill=type))+geom_bar(stat = 'identity')+facet_wrap(~Var1,scales='free')+scale_fill_manual(values = Color)+theme_pubr()+labs(fill='',y='relative abundance')+theme(text=element_text(family='serif'), axis.text.x=element_blank())

DimPlot(emt.merge,label=T)+DimPlot(emt.T,label=T)
p.abundance.all+p.abundance.tcell+plot_layout(widths=c(4,5))

buf<-data.frame(patient=emt.merge$orig.ident)
buf$type<-emt.merge@active.ident
buf$count<-1
buf1<-aggregate(buf$count, by=list(patient=buf$patient,type=buf$type), FUN=sum)
buf1$group<-ifelse(startsWith(buf1$patient,'C'),'Control','Endometriosis')
count.groupby.all<-buf1
buf<-aggregate(buf$count, by=list(patient=buf$patient), FUN=sum)
for (i in rownames(count.groupby.all)){
  count.groupby.all[i,'x']<-count.groupby.all[i,'x']/buf[buf$patient==count.groupby.all[i,'patient'],'x']
}
ggplot(count.groupby.all,aes(x=group,y=x,color=group))+geom_point()+facet_wrap(~type,scales='free')+scale_color_manual(values=Color)+stat_compare_means(family='serif',label = 'p.format',label.y.npc =0.9)+theme_pubr()+theme(text=element_text(family='serif'),axis.text.x=element_blank())+labs(color='',y='relative abundance')

buf<-data.frame(patient=emt.T$orig.ident)
buf$type<-emt.T@active.ident
buf$count<-1
buf1<-aggregate(buf$count, by=list(patient=buf$patient,type=buf$type), FUN=sum)
buf1$group<-ifelse(startsWith(buf1$patient,'C'),'Control','Endometriosis')
buf1<-buf1[buf1$type!='unknown',]
count.groupby.tcell<-buf1
buf<-aggregate(buf$count, by=list(patient=buf$patient), FUN=sum)
for (i in rownames(count.groupby.tcell)){
  count.groupby.tcell[i,'x']<-count.groupby.tcell[i,'x']/buf[buf$patient==count.groupby.tcell[i,'patient'],'x']
}


ggplot(count.groupby.tcell,aes(x=group,y=x,fill=group,color=group))+geom_boxplot()+facet_wrap(~type,scales='free',ncol=5)+scale_color_manual(values=Color1)+stat_compare_means(family='serif',label = 'p.format',label.y.npc =0.9)+theme_pubr()+theme(text=element_text(family="serif",face = 'bold'), axis.text = element_text(face='plain'),axis.title.x = element_blank(),axis.text.x=element_blank(), axis.title.y=element_text(face='plain'),legend.text=element_text(face='plain'))+labs(color='',fill='',y='relative abundance')+scale_fill_manual(values=Color)


DimPlot(emt.merge,label=T,raster=F)+DimPlot(emt.T,label=T,raster=F)

count.groupby.all$count['lymphocyte']

{
#data2(another paper)====
Idents(emt.merge)<-'clusterID'
DimPlot(emt.merge, reduction = "umap",raster = F,label=T)+DimPlot(emt.merge, reduction = "umap",raster = F,label=T)


emt.T<-subset(emt.merge, subset= clusterID %in% c('B','CD8T1','uNK1','CD4T','CD8T2','uNK2','CD8T3','CD8T4'))

emt.T<-NormalizeData(emt.T)
emt.T<-FindVariableFeatures(emt.T, selection.method = "vst")
all.genes <- rownames(emt.T)
emt.T <- ScaleData(emt.T, features = all.genes)

emt.T <- RunPCA(emt.T, features = VariableFeatures(emt.T),pcs.compute=30,seed.use=12345)

emt.T <- FindNeighbors(emt.T, dims = 1:20)
emt.T <- FindClusters(emt.T,resolution=6)
emt.T <- RunTSNE(emt.T,dims=1:30,seed.use=12345)
emt.T <- RunUMAP(emt.T,dims=1:30,seed.use=12345)



DotPlot(emt.T,features = markers)+theme_pubr()+theme(text = element_text(family='serif'),legend.position = 'right')+labs(x='Gene',y='Cell Type')+RotatedAxis()

FeaturePlot(emt.T,c('FGFBP2','SPON2'), cols = c("#DDDDDD","red"))+DimPlot(emt.T,label=T,group.by = 'RNA_snn_res.5')
table(emt.T$seurat_clusters)



emt.T$temp1<-emt.T$RNA_snn_res.5 %in% c('12')

FeaturePlot(emt.T,features=c('temp','temp1'))

markers<-c('CD4','CD8A','FOXP3','RTKN2')

DotPlot(emt.T,features = markers)+theme_pubr()+theme(text = element_text(family='serif'),legend.position = 'right')+labs(x='Gene',y='Cell Type')+RotatedAxis()


cluster<-c('CTL','','T','NK','T', #00~04
           'B cell','','NK','T','NK', #05~09
           '','T','NK','NK','NK', #10~14
           'NK','B cell','CD8','T','T', #15~19
           'NK','CD8','NK','B cell','T', #20~24
           'CD8','NK','NK','T','T', #25~29
           'CD8','CD8','NK','T','T', #30~34
           'NK','NK','T','B cell','NK', #35~39
           'CD8','','B cell','NK','CD8', #40~44
           'B cell','NK','T','NK','NK', #45~49
           'CD8','NK','NK','NK','T', #50~54
           '','T','B cell','NK','NK', #55~59
           'NK','NK','T','T','CD8', #60~64
           'NK','ILC','CD8','','T', #65~69
           '' #70
           )

cluster<-ifelse(cluster=='','unknown',cluster)
#cluster<-ifelse(cluster=='',0:70,cluster)
cluster<-ifelse(cluster=='T',0:70,cluster)
Idents(emt.T)='RNA_snn_res.6'
names(cluster) <- levels(emt.T)
emt.T<-RenameIdents(emt.T,cluster)
DimPlot(emt.T,label=T)


DimPlot(emt.T,label=T,group.by = 'RNA_snn_res.3')+DimPlot(emt.T,label=T,group.by = 'RNA_snn_res.5')
DimPlot(emt.merge,label=T)/DimPlot(emt.T,label=T)

markers<-list('NK'=c('ITGAX','NCAM1'), #
              'CTL'=c('FGFBP2','SPON2'),
              'T reg'=c('FOXP3','RTKN2'),
              'T cell'=c('CD3D','CD4','CD8B'),
              'memory'='CD40LG',
              'naive'=c('CCR7','LEF1'),
              'ILC'=c('IL7R','KIT'),
              'B cell'=c('CD19','MS4A1') #
)
buf<-FindMarkers(emt.T,ident.1 = '2')
buf<-buf[buf$p_val<0.05,]
buf<-buf[buf$p_val_adj<0.05,]
buf1<-buf[buf$avg_log2FC>0.5,]
buf1<-buf1[order(buf1$avg_log2FC,decreasing = T),]
for (i in rownames(buf1)){show(i)}
{markers<-c("IL7R",
            "TNFRSF4",
            "LTB",
            "LMNA",
            "GPR183",
            "SPOCK2",
            "ZFP36L2",
            "ZFP36",
            "ARL4C",
            "CREM",
            "TNFAIP3",
            "JUNB",
            "FTH1",
            "TRAC",
            "AC058791.1",
            "ICOS",
            "CALM1",
            "JUN",
            "ZC3H12D",
            "MAF",
            "KLF6",
            "GADD45B",
            "FOS",
            "CD6",
            "SOCS3",
            "ARID5B",
            "VIM",
            "RORA",
            "KLF2",
            "FAM107B",
            "CRYBG1",
            "TRBC1",
            "DDIT4",
            "EML4",
            "S100A4",
            "ANXA1",
            "PPP2R5C",
            "KDM6B",
            "BIRC3",
            "CYTOR",
            "TNFRSF25",
            "SYNE2",
            "PTGER4",
            "TUBB4B",
            "BICDL1",
            "ADAM19",
            "FNBP1",
            "IL32",
            "PABPC1",
            "PGAP1",
            "AHNAK",
            "SESN3",
            "AC016831.5",
            "CD2",
            "BCL11B",
            "BTG2",
            "TOB1",
            "CXCR4",
            "ZC3HAV1",
            "ICAM2",
            "IRF1",
            "FLNA",
            "NFKBIA",
            "TMEM173",
            "ZNF331",
            "CHMP1B",
            "AQP3",
            "VMP1",
            "CD48",
            "ANKRD12",
            "PMEPA1",
            "RPS16",
            "AC016831.1",
            "PLK3",
            "AL121944.1",
            "ANK3",
            "SARAF",
            "SOCS1",
            "BTG1",
            "CD3G",
            "ATP1B3",
            "PDE4D",
            "RNF125",
            "CRIP1",
            "STK24",
            "CD5",
            "EMB",
            "PDE4B",
            "EZR",
            "CD3D",
            "RNASET2",
            "TIPARP",
            "CCR6",
            "EEF1A1",
            "NOP53",
            "CDKN1A",
            "BAZ1A",
            "AP3M2",
            "S100A10",
            "RPL3",
            "YPEL5",
            "RPL4",
            "RGCC",
            "PAG1",
            "HIVEP2",
            "PIM3",
            "BATF",
            "CD28",
            "GSTK1",
            "MFHAS1",
            "DUSP16",
            "LINC00513",
            "IDS",
            "AC016831.7",
            "IFITM1",
            "GPR171",
            "TRBC2",
            "KLF3",
            "GLIPR1"
)}

##visualization----
cluster<-c('myeloid','lymphocyte','lymphocyte','lymphocyte','lymphocyte',
           'myeloid','lymphocyte','stromal','lymphocyte','lymphocyte',
           'epithelial','myeloid','myeloid','epithelial','myeloid',
           'epithelial','lymphocyte','endothelial-like','unknown')
names(cluster)<-levels(emt.merge)
emt.merge<-RenameIdents(emt.merge,cluster)
DimPlot(emt.merge,label=T)+DimPlot(emt.T,label=T)

table(emt.merge@active.ident) %>% group_by(type)
unique(emt.merge$type)
buf<-data.frame(table(emt.merge@active.ident[emt.merge$type=='Control']))
buf2<-data.frame(table(emt.merge@active.ident[emt.merge$type=='Symptomatic']))
buf3<-data.frame(table(emt.merge@active.ident[emt.merge$type=='Diagnosed']))
buf['count']<-buf$Freq/sum(buf$Freq)
buf2['count']<-buf2$Freq/sum(buf2$Freq)
buf3['count']<-buf3$Freq/sum(buf3$Freq)
buf['type']<-'Control'
buf2['type']<-'Symptomatic'
buf3['type']<-'Endometriosis'
count.all<-rbind(buf,buf2,buf3)
buf<-data.frame(table(emt.T@active.ident[emt.T$type=='Control']))
buf2<-data.frame(table(emt.T@active.ident[emt.T$type=='Symptomatic']))
buf3<-data.frame(table(emt.T@active.ident[emt.T$type=='Diagnosed']))
buf['count']<-buf$Freq/sum(buf$Freq)
buf2['count']<-buf2$Freq/sum(buf2$Freq)
buf3['count']<-buf3$Freq/sum(buf3$Freq)
buf['type']<-'Control'
buf2['type']<-'Symptomatic'
buf3['type']<-'Endometriosis'
count.tcell<-rbind(buf,buf2,buf3)

count.all$type<-factor(count.all$type,levels=c('Control','Symptomatic','Endometriosis'))
count.tcell$type<-factor(count.tcell$type,levels=c('Control','Symptomatic','Endometriosis'))
Color=c('Control'='#007B3D99','Symptomatic'='#e2e20099','Endometriosis'='#D5131799')
p.abundance.all<-ggplot(count.all,aes(x=type,y=count,fill=type))+geom_bar(stat = 'identity')+facet_wrap(~Var1,scales='free')+scale_fill_manual(values = Color)+theme_pubr()+labs(fill='',y='relative abundance',x='')+theme(text=element_text(family='serif'),axis.text.x = element_blank())
p.abundance.tcell<-ggplot(count.tcell,aes(x=type,y=count,fill=type))+geom_bar(stat = 'identity')+facet_wrap(~Var1,scales='free')+scale_fill_manual(values = Color)+theme_pubr()+labs(fill='',y='relative abundance',x='')+theme(text=element_text(family='serif'),axis.text.x = element_blank())

p.abundance.all+p.abundance.tcell+plot_layout(widths=c(4,5))

unique(cell.info[cell.info$pheno=='Diagnosed',]$subjectID)

buf<-cell.info[colnames(emt.merge),]
buf$pheno<-gsub('Diagnosed','D',(gsub('Symptomatic','S',gsub('Control','C',buf$pheno))))
buf$subjectID<-paste0(buf$subjectID,'_',buf$pheno)
table(rownames(buf)==colnames(emt.merge))
emt.merge$patient<-buf$subjectID

buf<-data.frame(patient=emt.merge$patient)
buf$type<-emt.merge@active.ident
buf$count<-1
buf1<-aggregate(buf$count, by=list(patient=buf$patient,type=buf$type), FUN=sum)
buf1$group<-ifelse(endsWith(buf1$patient,'C'),'Control',buf1$patient)
buf1$group<-ifelse(endsWith(buf1$patient,'S'),'Symptomatic',buf1$group)
buf1$group<-ifelse(endsWith(buf1$patient,'D'),'Endometriosis',buf1$group)
count.groupby.all<-buf1
buf<-aggregate(buf$count, by=list(patient=buf$patient), FUN=sum)
for (i in rownames(count.groupby.all)){
  count.groupby.all[i,'x']<-count.groupby.all[i,'x']/buf[buf$patient==count.groupby.all[i,'patient'],'x']
}
ggplot(count.groupby.all,aes(x=group,y=x,color=group))+geom_point()+facet_wrap(~type,scales='free')+scale_color_manual(values=Color)+stat_compare_means(family='serif',label = 'p.format')+theme_pubr()+theme(text=element_text(family='serif'))+
  labs(color='',y='relative abundance')

buf<-cell.info[colnames(emt.T),]
buf$pheno<-gsub('Diagnosed','D',(gsub('Symptomatic','S',gsub('Control','C',buf$pheno))))
buf$subjectID<-paste0(buf$subjectID,'_',buf$pheno)
table(rownames(buf)==colnames(emt.T))
emt.T$patient<-buf$subjectID

buf<-data.frame(patient=emt.T$patient)
buf$type<-emt.T@active.ident
buf$count<-1
buf1<-aggregate(buf$count, by=list(patient=buf$patient,type=buf$type), FUN=sum)
buf1$group<-ifelse(endsWith(buf1$patient,'C'),'Control',buf1$patient)
buf1$group<-ifelse(endsWith(buf1$patient,'S'),'Symptomatic',buf1$group)
buf1$group<-ifelse(endsWith(buf1$patient,'D'),'Endometriosis',buf1$group)
count.groupby.tcell<-buf1
buf<-aggregate(buf$count, by=list(patient=buf$patient), FUN=sum)
for (i in rownames(count.groupby.tcell)){
  count.groupby.tcell[i,'x']<-count.groupby.tcell[i,'x']/buf[buf$patient==count.groupby.tcell[i,'patient'],'x']
}

count.groupby.tcell$group<-factor(count.groupby.tcell$group,c('Control','Symptomatic','Endometriosis'))

ggplot(count.groupby.tcell,aes(x=group,y=x,color=group))+geom_point()+facet_wrap(~type,scales='free',ncol=4)+scale_color_manual(values=Color)+stat_compare_means(family='serif',label = 'p.format',comparisons = list(c('Control','Symptomatic'),c('Symptomatic','Endometriosis'),c('Control','Endometriosis')))+theme_pubr()+theme(text=element_text(family='serif'),axis.text.x=element_blank())+labs(color='',y='relative abundance',x='')


buf<-names(table(count.groupby.tcell$x)[table(count.groupby.tcell$x)>1])
count.groupby.tcell[count.groupby.tcell$x %in% buf,]

rownames(emt.T)[startsWith(rownames(emt.T),'FCGR')]

markers<-c('LINC00996','ITGAX','NCAM1', #NK
           'FGFBP2','FCGR3A', #CTL
           'CD3D','CD8A','CD8B', #T cell
           'TNFRSF4', 'LTB', #CD4 memory
           'b', #CD4 naive
           'c', #CD8 memory
           'd', #CD8 naive
           'HLA-DQB1','HLA-DQA1','BANK1' #B cell
           )
{
markers<-c("TNFRSF4",
           "IL7R",
           "LTB",
           "CCL20",
           "FTH1",
           "GPR183",
           "SPOCK2",
           "TRAC",
           "CD2",
           "CD3D",
           "TRBC1",
           "TRBC2",
           "ICOS",
           "TNFRSF25",
           "BATF",
           "HBA2",
           "NINJ1",
           "IL32",
           "TMEM173",
           "KLF2",
           "HBB",
           "ZFP36",
           "RPS12",
           "CD48",
           "CREM",
           "PASK",
           "CD6",
           "RORA",
           "VIM",
           "LMNA",
           "S100A10",
           "EEF1A1",
           "KLF6",
           "SARAF",
           "RPS16",
           "DDIT4",
           "PPP2R5C",
           "HBA1",
           "LGALS3",
           "SOD1",
           "JUNB",
           "LEPROTL1",
           "TNFAIP3",
           "RPL13",
           "ADAM19",
           "MAF",
           "FAM107B",
           "NFKBIA",
           "RPS6",
           "RPL3",
           "ARL4C",
           "S100A4",
           "ANXA1",
           "SESN3",
           "RPL11",
           "LDHB",
           "RPL10",
           "ARID5B",
           "ZC3H12D",
           "RPL36A",
           "TPT1",
           "RPS25",
           "CALM1",
           "CLEC2D",
           "CD5",
           "RPL4",
           "ICAM2",
           "ISG20",
           "FNBP1",
           "CCR7",
           "JUN",
           "BIRC3",
           "NOP53",
           "GSTK1",
           "AQP3",
           "RPL8",
           "CRYBG1",
           "CD3G",
           "CD28",
           "FAM177A1",
           "ZFP36L2",
           "RPS3",
           "RPL34",
           "GPR171",
           "MT2A",
           "RILPL2",
           "ANK3",
           "CIB1",
           "RPL9",
           "RPL21",
           "SOCS3",
           "RPS18",
           "CXCR4",
           "BICDL1",
           "TNFAIP8",
           "RPL32",
           "FLT3LG",
           "TC2N",
           "LCK",
           "RPL38",
           "DUSP16",
           "AC058791.1",
           "CCR6",
           "RPL36",
           "ANKRD12",
           "RPS4X",
           "TRAT1",
           "PRDX1",
           "RPL30",
           "SELENOK",
           "RPS10",
           "RPS28",
           "RPS2",
           "C12orf57",
           "PABPC1",
           "RPL12",
           "CYTOR",
           "FXYD5",
           "RPL39",
           "RPS8",
           "RPS21",
           "RPS20",
           "AL118516.1",
           "PIK3IP1",
           "RPS3A",
           "GLIPR1",
           "RPL19",
           "RPS14",
           "CORO1B",
           "RPL18A",
           "YIPF5",
           "HINT1",
           "RPL14",
           "RPSA",
           "CYCS",
           "RPL27A",
           "SELL",
           "RPS29",
           "PPA1",
           "AP3M2",
           "RNASET2",
           "PTGER4",
           "RPS27A",
           "YWHAB",
           "EEF1B2",
           "PLK3",
           "IL2RA",
           "ATP1B3",
           "RPL10A",
           "BTG1",
           "RPL22",
           "PIM3",
           "RPS27",
           "ARHGAP15",
           "EML4",
           "TOB1",
           "EIF1",
           "RPS13",
           "DUSP1",
           "RPL28",
           "RPS9",
           "RCAN3",
           "SERP1",
           "RPLP0",
           "IL4I1",
           "CERK",
           "RPL18",
           "RPL13A",
           "SLAMF1",
           "RPL7",
           "FRMD4B",
           "RPL5",
           "RPL23A",
           "EEF1D",
           "NOSIP",
           "GADD45B",
           "BCL11B",
           "RPLP2",
           "FURIN",
           "TRADD",
           "AC026979.2",
           "CCDC167",
           "RPL41",
           "PRDX2",
           "CDC14A",
           "AC245014.3",
           "EMP3",
           "RPL6",
           "KDSR",
           "TSHZ2",
           "MAL",
           "EIF3E",
           "RPL35A",
           "RPS15A",
           "GNAS",
           "IL23A",
           "RPS23",
           "SNHG15",
           "RPL29",
           "BAZ1A",
           "RPS5",
           "RPS26",
           "CDC37",
           "CAMK4",
           "H2AFV",
           "SVIP",
           "CDC42SE2",
           "PGAP1",
           "NAMPT",
           "AC044849.1",
           "CFAP36",
           "TNIP2",
           "NAP1L1",
           "TOMM7",
           "PBXIP1",
           "TMEM238",
           "RPS15",
           "NDUFS5",
           "RPS17",
           "AES",
           "EIF3K",
           "NAF1",
           "MGAT4A",
           "EMB",
           "TTC39C",
           "RPL26",
           "RGCC",
           "SNHG25",
           "PRDM1",
           "RPL27",
           "BEX2",
           "PEBP1",
           "NCK2",
           "FAM118A",
           "PDE3B",
           "APRT",
           "PAG1",
           "PFDN5",
           "LTA",
           "YBX1",
           "STK24",
           "SGK1",
           "NDUFV2",
           "TUBB2A",
           "HNRNPA1",
           "FYB1",
           "P2RY8",
           "S1PR1",
           "KLF3",
           "UBA52",
           "NPM1",
           "MCUB",
           "VMP1",
           "CMSS1",
           "RPL36AL",
           "FOS",
           "AREG",
           "SLC39A8",
           "ANP32B",
           "UQCRB",
           "SNHG9",
           "DNPH1",
           "AL121944.1",
           "SYNE2",
           "PMEPA1",
           "NSA2",
           "SC5D",
           "IDS",
           "NACA",
           "RPL17",
           "C21orf91",
           "RPL24",
           "RGS10",
           "IL6ST"
) #cd8 naive
} # cd4 memory
{
markers<-c("SELL",
           "SARAF",
           "CCR7",
           "HBA1",
           "LTB",
           "HBA2",
           "BATF",
           "ADTRP",
           "GIMAP7",
           "RGS10",
           "HBB",
           "MAL",
           "RPS12",
           "RPS13",
           "IL7R",
           "LDHB",
           "CD48",
           "FHIT",
           "LEF1",
           "NOSIP",
           "RPL32",
           "RPL11",
           "RPL5",
           "RPL30",
           "RPS3A",
           "EEF1B2",
           "RPL22",
           "ICOS",
           "NCK2",
           "CD3D",
           "ARID5A",
           "LEPROTL1",
           "RPS10",
           "RPS25",
           "KLF2",
           "PIK3IP1",
           "PRMT2",
           "RPS6",
           "C12orf57",
           "PDCD4",
           "RPL13",
           "IKZF1",
           "RPL19",
           "TRAC",
           "RPL34",
           "CMTM8",
           "RPL9",
           "RPL21",
           "RPS14",
           "RPS5",
           "RPS28",
           "CD2",
           "CLEC2D",
           "PNRC1",
           "GYPC",
           "NPM1",
           "RPL38",
           "RPS8",
           "RPL14",
           "TAF4B",
           "SOD1",
           "EIF3E",
           "EEF1A1",
           "RPL36A",
           "RPL4",
           "RPS23",
           "TSHZ2",
           "RPL10A",
           "NACA",
           "TRBC1",
           "NDFIP1",
           "RPL10",
           "RPS27A",
           "SBDS",
           "MSL3",
           "RPL35A",
           "RPL3",
           "RPS21",
           "BEX2",
           "RPL18A",
           "CD55",
           "NOP53",
           "OCIAD2",
           "CAMK4",
           "PASK",
           "AIF1",
           "RPL39",
           "TMEM243",
           "CD6",
           "RPS3",
           "IL6ST",
           "LYPD3",
           "RPL29",
           "RPS9",
           "AES",
           "NELL2",
           "RPS20",
           "RPL18",
           "RPSA",
           "SNHG8",
           "CYTIP",
           "RPL12",
           "RPL36",
           "RPL6",
           "RPL7",
           "TNFRSF25",
           "RPS4X",
           "GIMAP1",
           "SOD2",
           "RPS2",
           "RPS15A",
           "RPS27",
           "MCUB",
           "RPS29",
           "RPS18",
           "RNASET2",
           "RPLP2",
           "SOCS3",
           "PIM2",
           "CYCS",
           "RPL8",
           "BICDL1",
           "KLF3",
           "FLT3LG",
           "BAZ1A",
           "TOMM7",
           "GBP1",
           "TMEM123",
           "TMSB10",
           "UBE2B",
           "PPA1",
           "H2AFZ",
           "GPR183",
           "HAPLN3",
           "STAT1",
           "SVIP",
           "EIF3H",
           "FYB1",
           "AC026979.2",
           "SF1",
           "TMEM173",
           "TRBC2",
           "RPS15",
           "CRYBG1",
           "CFL2",
           "SCML1",
           "NSA2",
           "TPT1",
           "PDE3B",
           "IRF1",
           "ITK",
           "APRT",
           "BEX3",
           "RPL17",
           "DDX21",
           "RPL24",
           "SEPT1",
           "CRTC3",
           "RPL37",
           "CIB1",
           "FAM177A1",
           "SERINC5",
           "SPOCK2",
           "CD3G",
           "ANP32B",
           "IFIT3",
           "SC5D",
           "SNHG25",
           "RPL28",
           "RPL31",
           "AP3M2",
           "ANXA2R",
           "UQCRB",
           "LINC00861",
           "FOXP1",
           "BTG1",
           "HINT1",
           "ATP5PO",
           "EIF1",
           "ITGA6",
           "EEF1D",
           "DEAF1",
           "RSL1D1",
           "RPS7",
           "BCL11B",
           "NSMCE3",
           "RPS17",
           "CD28",
           "GTF2B",
           "CCND3",
           "TGFBR2",
           "EIF3K",
           "UBALD2",
           "ANK3",
           "LCK",
           "KRAS",
           "C1QBP",
           "RIPOR2",
           "PFDN5",
           "FAU",
           "FKBP5",
           "RPL27A",
           "DNPH1",
           "MGAT4A",
           "RPL7A",
           "SESN3",
           "RETREG1",
           "SNHG9",
           "SATB1",
           "RPL27",
           "ARMH1",
           "SH3YL1",
           "RACK1",
           "S1PR1",
           "CSGALNACT1",
           "ALKBH7",
           "COX4I1",
           "CHRM3-AS2",
           "PRDX2",
           "RPL23",
           "SELENOM",
           "RPL26",
           "TCF7",
           "CCDC57",
           "MARCKSL1",
           "TOMM20",
           "CD69",
           "GSTK1",
           "RPL23A",
           "DDX5",
           "PPP2R5C",
           "KDSR",
           "RCAN3",
           "LIMD2",
           "SMDT1",
           "FAM118A",
           "HNRNPA1",
           "RPS16",
           "FAM184A",
           "EIF2S3",
           "CD27",
           "ANKRD12",
           "SUSD3",
           "RPL41",
           "RPL13A",
           "TRABD2A",
           "CHD7",
           "EIF3F",
           "UXT",
           "C1orf162",
           "IFITM1",
           "GLRX",
           "JUNB",
           "RASA3",
           "RPS24",
           "CALM1",
           "LTA",
           "PITPNC1",
           "FXYD5",
           "RASGRP2",
           "COX7C",
           "LPAR6",
           "ASXL1",
           "OXNAD1",
           "WDR74",
           "STMN1",
           "BTF3",
           "SOCS2",
           "NSMCE1",
           "YBX1",
           "EIF3L",
           "TESPA1",
           "MXI1",
           "TRAT1",
           "UBA52",
           "BCAS2",
           "SUB1",
           "NDUFV2",
           "USP3",
           "REXO2",
           "PLAC8",
           "CDKN2D",
           "NDUFAF4",
           "BIRC2",
           "PABPC1",
           "LRRC75A",
           "MYC",
           "SNRPG",
           "AREG",
           "ERP29",
           "SMAP2",
           "PDK1",
           "SNHG7",
           "DUT",
           "CISH",
           "RAPGEF6",
           "AK6",
           "RHOH",
           "MX1",
           "CALM2",
           "BCL7B",
           "NXT1",
           "ISG20",
           "ARL4C",
           "GIMAP4",
           "TRIM22",
           "HELB",
           "PCED1B",
           "RSL24D1",
           "MARCH9",
           "TNFAIP8",
           "DDX18",
           "CITED4",
           "PEBP1",
           "SSH2",
           "ALG13",
           "COMMD6",
           "MX2",
           "LINS1",
           "YBX3",
           "TC2N",
           "C6orf48",
           "ZNF101",
           "FAM241A",
           "RPS26",
           "RIC3",
           "HNRNPDL",
           "CIRBP",
           "CD5",
           "TRAPPC6A",
           "RPLP0",
           "CDC14A",
           "PDE4B",
           "IFI44L",
           "IRF9",
           "TAF9",
           "TXNIP",
           "LAT",
           "CNBP",
           "DGKA",
           "NAP1L1",
           "HSBP1L1",
           "BCL3",
           "ATM",
           "SERP1",
           "PRKCA",
           "TRADD",
           "UBXN1",
           "SLC30A7",
           "SLC25A6",
           "SBNO2",
           "EIF4B",
           "CMSS1",
           "NDUFC1",
           "ZC3HAV1"
)
} # cd4 naive
{
markers<-c(
  "CCL5",
  "CD8A",
  "CD3D",
  "GZMH",
  "TUBA4A",
  "CD8B",
  "TRGC2",
  "GZMK",
  "IFNG",
  "GZMM",
  "LYAR",
  "FTH1",
  "CYBA",
  "TRBC2",
  "PPP2R5C",
  "CD3G",
  "TRAC",
  "LINC01871",
  "PPP1R14B",
  "LAG3",
  "CCL20",
  "ANXA1",
  "CXCR4",
  "S100A10",
  "CST7",
  "KLRG1",
  "TRBC1",
  "FGFBP2",
  "SPOCK2",
  "CLIC3",
  "TNF",
  "DUSP2",
  "PAXX",
  "C12orf75",
  "IFITM2",
  "GNG2",
  "KLRB1",
  "RARRES3",
  "DNAJB1",
  "IL7R",
  "HLA-B",
  "TC2N",
  "RGCC",
  "VIM",
  "CREM",
  "GCHFR",
  "HSPA8",
  "TGFB1",
  "HLA-C",
  "NDUFS5",
  "ALOX5AP",
  "LMNA",
  "SRSF7",
  "B2M",
  "CD2",
  "RGS2",
  "CCL4L2",
  "CALM1",
  "LCK",
  "DUSP1",
  "GADD45B",
  "ABHD17A",
  "UBE2S",
  "IL32",
  "RORA",
  "APOBEC3G",
  "PMAIP1",
  "MIAT",
  "ZEB2",
  "FAM177A1",
  "NINJ1",
  "GUK1",
  "CCNH",
  "OAZ1",
  "PSMA7",
  "C1orf21",
  "HLA-A",
  "S100A4",
  "OASL",
  "GABARAPL1",
  "MT2A",
  "NCR3",
  "IFITM1",
  "CD99",
  "AC058791.1",
  "TIGIT",
  "SRSF2",
  "UBB",
  "GPR65",
  "TNFSF9",
  "SERTAD1",
  "ITGB2",
  "FLNA",
  "GNAS",
  "GTF3C1",
  "CD48",
  "TRMT112",
  "CLEC2D",
  "MAF",
  "RPS3",
  "H2AFV",
  "ITM2C",
  "CXCR3",
  "HSPE1",
  "EFHD2",
  "AL118516.1",
  "IER5L",
  "ERN1",
  "ODC1",
  "SELENOK",
  "ANXA5",
  "ORAI1",
  "TNFAIP3",
  "ADRB2",
  "LUZP1",
  "ZFP36",
  "H2AFZ",
  "ISG20",
  "ARHGDIB",
  "TMEM50A",
  "H1FX",
  "SH2D1A",
  "SKP1",
  "LEPROTL1",
  "PLIN2",
  "CIB1",
  "SYTL2",
  "LDHA",
  "LGALS3",
  "PITPNC1",
  "YWHAB",
  "CAMK2N1",
  "ARL4C",
  "MYL12B",
  "HCST",
  "CALM2",
  "RNF19A",
  "KLF2",
  "LITAF",
  "TRAT1",
  "YWHAQ",
  "PRNP",
  "CD6",
  "FKBP2",
  "IDS",
  "PTMA",
  "EML4",
  "CRIP1",
  "LAT",
  "JUN",
  "YIPF5",
  "NDUFB8",
  "TUBB4B",
  "SEC11C",
  "ZFP36L2",
  "METRNL",
  "RPS12",
  "DYNLL1",
  "ISCU",
  "TUBB2A",
  "TOB1",
  "MAP1LC3B",
  "RAB5IF",
  "RNF166",
  "RPL36AL",
  "SAP18",
  "CDC42SE2",
  "PRDX1",
  "H3F3B",
  "YPEL5",
  "PREX1",
  "ITGA4",
  "UBC",
  "CD5",
  "FKBP11",
  "JUNB",
  "GPR183",
  "JUND",
  "CNBP",
  "TTC39C",
  "SNRPN",
  "PTMS",
  "CD3E",
  "RPS27",
  "CERK",
  "ARF4",
  "CDC42EP3",
  "RABAC1",
  "SH3BGRL3",
  "DNAJA1",
  "RPS26",
  "HLA-F",
  "PSME1",
  "JMJD6",
  "POLR2L",
  "SLC9A3R1",
  "AC044849.1",
  "OTULIN",
  "F2R",
  "CHMP1B",
  "GLIPR1",
  "TUBA1B",
  "DUSP4",
  "FCGR3A",
  "JPT1",
  "RPL10",
  "SERP1",
  "MAGOH",
  "PTPRC",
  "SSR4",
  "MCUB",
  "NFATC2",
  "UBL5",
  "RBM8A",
  "EMB",
  "MXRA7",
  "GYPC",
  "RPL23A",
  "RPS16",
  "FNBP1",
  "AKIRIN1",
  "MYL6",
  "PTTG1",
  "SLC3A2",
  "FAM107B",
  "PSMB9",
  "CCSER2",
  "SRRT",
  "SYNE2",
  "GCLM",
  "EIF1",
  "CLTB",
  "POLR2K",
  "UBE2A",
  "LYST",
  "SPCS1",
  "C9orf78",
  "RPL28",
  "PNP",
  "APMAP",
  "LBH",
  "STX11",
  "UQCRB",
  "RPL3",
  "FDX1",
  "PSME2",
  "TGFBR3",
  "RPS14",
  "HNRNPUL1",
  "SVIP",
  "FYN",
  "VAMP5",
  "CYCS",
  "NDUFA6",
  "ANXA2R",
  "C21orf91",
  "EEF1D",
  "SLC7A5",
  "TNIP2",
  "DYNLL2",
  "ZC3H12A",
  "RPL41",
  "TECR",
  "SOCS1",
  "AC243960.1",
  "SH3KBP1",
  "LCP1",
  "MYO1F",
  "CCL4",
  "KRTCAP2",
  "DNAJB9",
  "TAF10",
  "NAA50",
  "RPS15A",
  "TSC22D3",
  "CACYBP",
  "RAB27A",
  "FTL",
  "EIF3K",
  "DSTN",
  "ARHGAP26",
  "NDUFA4",
  "SEC61B",
  "MOB2",
  "REEP5",
  "NENF",
  "FABP5",
  "C16orf54",
  "PPDPF",
  "RPL21",
  "PLK3",
  "RPL26",
  "ABI3",
  "MGAT4A",
  "TSPYL2",
  "CDK2AP2",
  "EMD",
  "AHNAK",
  "PPIB",
  "BTG1",
  "CCDC167",
  "SOD1",
  "P2RY8",
  "ARPC2",
  "FAM129A",
  "HIST2H2AA4",
  "AC016831.5",
  "EEF1A1",
  "ATP6V0E1",
  "DUSP16",
  "HLA-E",
  "MRPL41",
  "HNRNPLL",
  "XBP1",
  "RARA",
  "GPR171",
  "TGIF1"
)
} #cd8 memory


markers<-c('SELL','CCR7','CD4','CD8A')
markers<-toupper(markers)

DotPlot(emt.T,features=c(markers),group.by = 'RNA_snn_res.3')+theme(axis.text.x=element_text(angle=90, hjust=1))
DotPlot(emt.T,features=c(markers))+theme(axis.text.x=element_text(angle=90, hjust=1))

DimPlot(emt.T,label=T,group.by = 'RNA_snn_res.3')

buf<-FindMarkers(emt.T,ident.1 = 'CD8 naive')
buf<-buf[buf$p_val_adj<0.05,]
buf<-buf[buf$p_val<0.05,]
buf1<-buf[buf$avg_log2FC>(0),]
buf1<-buf1[order(buf1$avg_log2FC,decreasing = T),]
for (i in rownames(buf1)){show(i)}
} #버려진 데이터

#single cell KEGG, GO ----

library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(IOBR)
library(ggplot2)
library(reshape)
library(ggpubr)
library(vegan)
library(gridExtra)
library(ALDEx2)
library(microbial)
library(microbiome)
library(microbiomeMarker)
library(phyloseq)
library(pheatmap)
library(tidyverse)
library(edgeR)
library(DESeq2)
library(Hmisc)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library('org.Hs.eg.db', character.only = TRUE)
used.data<-"/data2/msi/Rproject/EMT.scRNA.RData"
load(used.data)

##compare with other lymphocyte----
###GO(lymphocyte)----
DimPlot(emt.T,label=T)
buf<-FindMarkers(emt.T,ident.1 = 'CD4 memory')
buf<-buf[buf$p_val<0.05,]
buf<-buf[buf$p_val_adj<0.05,]
buf<-buf[order(buf$avg_log2FC,decreasing = T),]
buf1<-buf[buf$avg_log2FC>0.5,]
buf2<-buf[buf$avg_log2FC<(-0.5),]


gene_list<-buf1$avg_log2FC
gene_list.neg<-buf2$avg_log2FC
names(gene_list)<-rownames(buf1)
names(gene_list.neg)<-rownames(buf2)

gse.alias <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ALIAS", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")
gse.alias.neg <- gseGO(geneList=gene_list.neg, 
                   ont ="ALL", 
                   keyType = "ALIAS", 
                   nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = org.Hs.eg.db, 
                   pAdjustMethod = "none")


p1<-clusterProfiler::dotplot(gse.alias, showCategory=10,title='GO',split='ONTOLOGY') + facet_grid(ONTOLOGY~., scale="free")+scale_fill_gradient(low='red',high='#ffafaf') + theme(axis.text.y=element_text(size=10),text = element_text(family = 'serif'))
p1.neg<-clusterProfiler::dotplot(gse.alias.neg, showCategory=10,title='GO',split='ONTOLOGY') + facet_grid(ONTOLOGY~., scale="free")+scale_fill_gradient(low='red',high='#ffafaf') + theme(axis.text.y=element_text(size=10),text = element_text(family = 'serif'))

p1

###kegg pathway(lymphocyte)====
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted

ids<-bitr(names(gene_list), fromType = "ALIAS", toType = "ENTREZID", OrgDb=org.Hs.eg.db)

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df <- buf[ids$ALIAS,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df$Y = ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df$avg_log2FC

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

ids<-bitr(names(gene_list.neg), fromType = "ALIAS", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
df <- buf[ids$ALIAS,]
df$Y = ids$ENTREZID
kegg_gene_list.neg <- df$avg_log2FC
names(kegg_gene_list.neg) <- df$Y
kegg_gene_list.neg<-na.omit(kegg_gene_list.neg)
kegg_gene_list.neg = sort(kegg_gene_list.neg, decreasing = TRUE)

kegg <- gseKEGG(geneList     = kegg_gene_list,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
kegg.neg <- gseKEGG(geneList     = kegg_gene_list.neg,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")


p2<-dotplot(kk1, showCategory=10,title='KEGG') +scale_fill_gradient(low='red',high='#ffafaf') + theme(axis.text.y=element_text(size=10),text = element_text(family = 'serif'))
p2.neg<-dotplot(kk1.neg, showCategory=10,title='KEGG') +scale_fill_gradient(low='red',high='#ffafaf') + theme(axis.text.y=element_text(size=10),text = element_text(family = 'serif'))

p2


p1+p2+plot_layout(widths=c(2.5,1.5)) #horizenal
p1/p2+plot_layout(heights =c(5,2)) #vertical

p1.neg/p2.neg+plot_layout(heights =c(5,2)) #vertical

ridgeplot(kk2)

buf<-names(emt.T@active.ident[emt.T@active.ident=='CD4 memory'])

table(buf %in% names(emt.merge@active.ident))

##compare with other cell====
###GO(all cell)====
emt.merge$RNA_snn_res.0.8<-emt.merge@active.ident
levels(emt.merge$RNA_snn_res.0.8)<-c(levels(emt.merge$RNA_snn_res.0.8),'CD4 memory')
emt.merge$RNA_snn_res.0.8[buf]<-factor('CD4 memory')
Idents(emt.merge)<-'RNA_snn_res.0.8'

buf<-FindMarkers(emt.merge,ident.1 = 'CD4 memory')
buf<-buf[order(buf$avg_log2FC,decreasing = T),]
buf<-buf[buf$p_val<0.05,]
buf<-buf[buf$p_val_adj<0.05,]
buf1<-buf[buf$avg_log2FC>0.5,]
buf2<-buf[buf$avg_log2FC<(-0.5),]
gene_list<-buf1$avg_log2FC
gene_list.neg<-buf2$avg_log2FC
names(gene_list)<-rownames(buf1)
names(gene_list.neg)<-rownames(buf2)

gse.allcell <- gseGO(geneList=gene_list, 
                   ont ="ALL", 
                   keyType = "ALIAS", 
                   nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = org.Hs.eg.db, 
                   pAdjustMethod = "none")

p1<-clusterProfiler::dotplot(gse.allcell, showCategory=10,title='GO',split='ONTOLOGY') + facet_grid(ONTOLOGY~., scale="free")+scale_fill_gradient(low='red',high='#ffafaf') + theme(axis.text.y=element_text(size=10),text = element_text(family = 'serif'))
p1.neg<-clusterProfiler::dotplot(gse.allcell, showCategory=10,title='GO',split='ONTOLOGY') + facet_grid(ONTOLOGY~., scale="free")+scale_fill_gradient(low='red',high='#ffafaf') + theme(axis.text.y=element_text(size=10),text = element_text(family = 'serif'))


###kegg(all cell)----
ids<-bitr(names(gene_list), fromType = "ALIAS", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
df <- buf[ids$ALIAS,]
df$Y = ids$ENTREZID
kegg_gene_list <- df$avg_log2FC
names(kegg_gene_list) <- df$Y
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
kk2.neg <- gseKEGG(geneList     = kegg_gene_list,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

p2.neg<-dotplot(kk2.neg, showCategory=10,title='KEGG') + facet_grid(~., scale='free')+scale_fill_gradient(low='red',high='#ffafaf') + theme(axis.text.y=element_text(size=10),text = element_text(family = 'serif'))


p1+p2+plot_layout(widths=c(2.5,1.5)) #horizenal
p1/p2+plot_layout(heights =c(5,2)) #vertical

p1.neg/p2.neg+plot_layout(heights =c(5,2)) #vertical


#brief report----
library(ggplot2)
library(ggpubr)
library(reshape2)
stain=c('-','±','+','++','+++')
stain=c('0','1','2','3','4')
buf<-data.frame(
  staining=stain,
  B.cereus_control=c(10,0,0,0,0),
  B.cereus_endometriosis=c(5,3,19,21,20)
)
buf1<-colSums(buf[2:ncol(buf)])
for (i in colnames(buf)[2:ncol(buf)]){
  buf[i]<-buf[i]/buf1[i]
}
buf$staining<-factor(buf$staining,levels = stain)
buf<-melt(buf,id.vars = 'staining')
#buf<-buf[(buf$variable %in% c('B.cereus_endometriosis','B.cereus_control')),]
buf$variable<-ifelse(buf$variable=='B.cereus_endometriosis','Endometriosis','Control')
buf$variable<-factor(buf$variable,levels = c('Endometriosis','Control'))
ggplot(buf,aes(x='',y=value,fill=staining))+
  geom_bar(position = 'fill',stat='identity',color="black")+
  coord_polar('y')+
  theme_void()+
  scale_fill_brewer(palette = 'Blues')+
  theme(text=element_text(family='arial',size=16),legend.title=element_text(size=12))+
  geom_text(aes(x=1.65,label=paste0(round(value*100,1), '%')),size=4,position = position_stack(vjust=0.5))+
  facet_wrap(~variable,ncol=2)+labs(fill='Histopathologic\nscore') #2배 안하고 저장

#Staining analysis----
library(readxl)
library(ggplot2)
library(ggpubr)
patient.stain<-read_xlsx('/home/user02/raw_data1/seongik/Supplementary Tables.xlsx',sheet = 'Supplement Table 2')
patient.stain<-patient.stain[patient.stain$Sample!='A67-B',]
patient.stain$Sample[patient.stain$Sample=='A67-A']<-'A67'

image.j<-read_excel('ihc-imagej.xlsx',sheet = 'image J data')
colnames(image.j)<-image.j[3,]
image.j<-image.j[-c(1,2,3,4,15,83,85),-c(5,6)]
#                   번호  bac, tissue에 대한 rid값    b rid/t rid   area에 대해 rid와 같이 수행
colnames(image.j)<-c('No','Bacillus_rid','tissue_rid','b.t_rid','Bacillus_area','tissue_area','b.t_area')
#rid=raw integrated density, area=단순 면적
image.j[,-1]<-t(apply(FUN = as.numeric,X=image.j[,-1],MARGIN = 1))
patient.stain<-cbind(patient.stain,image.j[,-1])


buf<-patient.stain[patient.stain$type=='Endometriosis',]
buf$B.cereus<-factor(buf$`Staining level`,levels = c('-','±','+','++','+++'))
buf$stain<-factor(buf$`Histopathologic score`)

buf$group<-ifelse(buf$B.cereus %in% c('-','±'),'low','high')
buf$group<-factor(buf$group,levels = c('low','high'))

buf$infertility<-ifelse(buf$`Infertility Y/N`=='Y',1,0)
#buf<-buf[buf$B.cereus!='-',]
colnames(buf)


Color4<-c('#FFa5a5','#CB1106')

value<-'AMH' #size AMH
p1<-ggplot(buf[buf$type=='Endometriosis',],aes(x=B.cereus, y=get(value),fill=B.cereus,color=B.cereus))+geom_boxplot()+theme_pubr(legend='none')+labs(x='Histopathologic score',y=value)+
  scale_fill_manual(values = paste0(Color2,'99'))+scale_color_manual(values = Color2)+
  theme(text=element_text(family='arial',size=12))
p2<-ggplot(buf[buf$type=='Endometriosis',],aes(x=group, y=get(value), fill=group, color=group))+geom_boxplot()+theme_pubr(legend='none')+labs(y=value)+
  scale_fill_manual(values = paste0(Color4,'99'))+scale_color_manual(values = Color4)+
  theme(text=element_text(family='arial',size=12))+labs(x='')#+stat_compare_means(label='p') #0.28

p1+p2

value<-'Gravidity' #Gravidity Parity
p1<-ggplot(buf,aes(x=B.cereus, y=get(value),fill=B.cereus,color=B.cereus))+geom_violin()+theme_pubr(legend='none')+labs(x='Histopathologic score',y=value)+
  scale_fill_manual(values = paste0(Color2,'99'))+scale_color_manual(values = Color2)+
  theme(text=element_text(family='arial',size=12))
p2<-ggplot(buf[buf$type=='Endometriosis',],aes(x=group, y=get(value), fill=group, color=group))+geom_violin()+theme_pubr(legend='none')+labs(y=value)+
  scale_fill_manual(values = paste0(Color4,'99'))+scale_color_manual(values = Color4)+
  theme(text=element_text(family='arial',size=12))+stat_compare_means(label='p')

p1+p2

value<-'infertility'
p1<-ggplot(buf,aes(fill=`Infertility Y/N`,x=B.cereus))+geom_bar(position='dodge')+labs(y=value)+
  scale_fill_manual(values = paste0(Color,'99'))+theme_pubr()+
  theme(text=element_text(family='arial',size=12))
p2<-ggplot(buf,aes(fill=`Infertility Y/N`,x=group))+geom_bar(position='fill')+labs(y=value)+
  scale_fill_manual(values = paste0(Color,'99'))+theme_pubr()+
  theme(text=element_text(family='arial',size=12))
p1<-ggplot(buf,aes(fill=`Infertility Y/N`,x=group))+geom_bar(position='fill')+labs(y=value)+
  scale_fill_manual(values = paste0(Color,'99'))+theme_pubr()+
  theme(text=element_text(family='arial',size=12))
p1+p2

ggplot(buf,aes(x=`Infertility Y/N`,fill=group))+geom_bar()+labs(y=value)+
  theme_pubr()+
  theme(text=element_text(family='arial',size=12))

ggplot(buf,aes(fill=`Infertility Y/N`,x=group))+geom_bar(position='fill')+
  coord_polar('y')+labs(y=value)+
  scale_fill_manual(values = paste0(Color,'99'))+theme_pubr()+facet_wrap(~group)+
  theme(text=element_text(family='arial',size=12))

ggplot(buf,aes(x='',y=value,fill=staining))+
  geom_bar(position = 'fill',stat='identity',color="black")+
  coord_polar('y')+
  theme_void()+
  scale_fill_brewer(palette = 'Blues')+
  theme(text=element_text(family='arial',size=16),legend.title=element_text(size=12))+
  geom_text(aes(x=1.65,label=paste0(round(value*100,1), '%')),size=4,position = position_stack(vjust=0.5))+
  facet_wrap(~variable,ncol=2)+labs(fill='Histopathologic\nscore') #2배 안하고 저장

ggplot(buf,aes(x='',y=stain,fill=stain))+
  geom_bar(position = 'fill',stat='identity',color="black")+
  coord_polar('y')+
  theme_void()+
  scale_fill_brewer(palette = 'Blues')+
  theme(text=element_text(family='arial',size=16),legend.title=element_text(size=12))+
  facet_wrap(~type,ncol=2)+labs(fill='B. cereus\nstaining intensity') #label text 수정정
##t-test----
library(car)

leveneTest(age ~ type,data=patient.stain)
t.test(age ~ type, data=patient.stain, var.equal=T, conf.level = 0.95) 
sd(patient.stain[patient.stain$type == "Normal",]$age)
sd(patient.stain[patient.stain$type != "Normal",]$age)

leveneTest(height ~ type,data=patient.stain)
t.test(height ~ type, data=patient.stain, var.equal=T, conf.level = 0.95) 
sd(patient.stain[patient.stain$type == "Normal",]$height)
sd(patient.stain[patient.stain$type != "Normal",]$height)

leveneTest(weight ~ type,data=patient.stain)
t.test(weight ~ type, data=patient.stain, var.equal=T, conf.level = 0.95) 
sd(patient.stain[patient.stain$type == "Normal",]$weight)
sd(patient.stain[patient.stain$type != "Normal",]$weight)

leveneTest(BMI ~ type,data=patient.stain)
t.test(BMI ~ type, data=patient.stain, var.equal=T, conf.level = 0.95)
sd(patient.stain[patient.stain$type == "Normal",]$BMI)
sd(patient.stain[patient.stain$type != "Normal",]$BMI)

ggplot(patient.stain,aes(x=type,y=age))+geom_boxplot()+stat_compare_means()

# correlation staining----
coef <- round(coef(model), 2)
text(2, 70,  paste("Y = ", coef[1], "+", coef[2], "x"))

Color2<-c('#01D11A','#FFa5a5','#FF6D6D','#FF1D1D','#CB1106')
(ggplot(buf[buf$`Histopathologic score`!=0,],aes(x=stain, y=size,fill=B.cereus,color=B.cereus))+geom_boxplot()+theme_pubr(legend='none')+labs(x='Histopathologic score',y='Size'))+
  scale_fill_manual(values = paste0(Color2,'99'))+scale_color_manual(values = Color2)+
  theme(text=element_text(family='arial',size=12))

Color2<-c('#FFa5a5','#CB1106')


ggplot(buf,aes(x=stain, y=b.t_rid,fill=stain))+geom_boxplot()+theme_pubr()+
  labs(x='group')+scale_fill_manual(values = paste0(Color2,'aa'))+scale_color_manual(values = Color2)+
  theme(text=element_text(family='arial',size=12))

Color2<-c('#01D11A','#fff700','#ff9500','#ff1100','#CB1106')

value<-'AMH'
ymax=20
bb<-summary(lm(get(value)~b.t_rid,buf))
p.corr.withsize.btrid<-ggplot(buf,aes(x=b.t_rid, y=get(value),color=stain))+geom_point()+theme_pubr()+
  labs(x='b/t RawIntDen',y=value)+scale_fill_manual(values = paste0(Color2,'aa'))+scale_color_manual(values = Color2)+
  theme(text=element_text(family='arial',size=12))+
  geom_abline(intercept=bb$coefficients[1], slope=bb$coefficients[2],linewidth=1)+ #회귀선
  #geom_smooth(method='lm',se=F)+
  annotate("text", x=max(buf['b.t_rid'])*0.7, y=ymax, label = paste0('R^2: ',round(bb$r.squared,3), ' / P_val: ',round(bb$coefficients[2,4],3)), size=6)


bb<-summary(lm(get(value)~b.t_area,buf))
p.corr.withsize.btarea<-ggplot(buf,aes(x=b.t_area, y=get(value),color=stain))+geom_point()+theme_pubr()+
  labs(x='b/t area',y=value)+scale_fill_manual(values = paste0(Color2,'aa'))+scale_color_manual(values = Color2)+
  theme(text=element_text(family='arial',size=12))+
  geom_abline(intercept=bb$coefficients[1], slope=bb$coefficients[2],linewidth=1)+ #회귀선
  #geom_smooth(method='lm',se=F)+
  annotate("text", x=max(buf['b.t_area'])*0.7, y=ymax, label = paste0('R^2: ',round(bb$r.squared,3), ' / P_val: ',round(bb$coefficients[2,4],3)), size=6)

bb<-summary(lm(get(value)~Bacillus_rid,buf))
p.corr.withsize.brid<-ggplot(buf,aes(x=Bacillus_rid, y=get(value),color=stain))+geom_point()+theme_pubr()+
  labs(x='Bacillus RawIntDen',y=value)+scale_fill_manual(values = paste0(Color2,'aa'))+scale_color_manual(values = Color2)+
  theme(text=element_text(family='arial',size=12))+
  geom_abline(intercept=bb$coefficients[1], slope=bb$coefficients[2],linewidth=1)+ #회귀선
  #geom_smooth(method='lm',se=F)+
  annotate("text", x=max(buf['Bacillus_rid'])*0.7, y=ymax, label = paste0('R^2: ',round(bb$r.squared,3), ' / P_val: ',round(bb$coefficients[2,4],3)), size=6)

bb<-summary(lm(get(value)~Bacillus_area,buf))
p.corr.withsize.barea<-ggplot(buf,aes(x=Bacillus_area, y=get(value),color=stain))+geom_point()+theme_pubr()+
  labs(x='Bacillus area',y=value)+scale_fill_manual(values = paste0(Color2,'aa'))+scale_color_manual(values = Color2)+
  theme(text=element_text(family='arial',size=12))+
  geom_abline(intercept=bb$coefficients[1], slope=bb$coefficients[2],linewidth=1)+ #회귀선
  #geom_smooth(method='lm',se=F)+
  annotate("text", x=max(buf['Bacillus_area'])*0.7, y=ymax, label = paste0('R^2: ',round(bb$r.squared,3), ' / P_val: ',round(bb$coefficients[2,4],3)), size=6)

p.corr.withsize.btrid+p.corr.withsize.brid+p.corr.withsize.btarea+p.corr.withsize.barea+plot_layout(ncol=2)



buf1<-melt(buf,measure.vars = c('size','b.t_rid','b.t_area'))

buf1<-buf1[!(buf1$stain==0 & buf1$variable=='size'),]

ggplot(buf1,aes(x=stain, y=value,fill=stain,color=stain))+geom_boxplot()+theme_pubr(legend='none')+
  labs(x='group')+stat_compare_means(comparisons = list(c('0','1'),c('0','2'),c('0','3'),c('0','4')),label = 'p.signif' )+facet_wrap(~variable,scales = 'free')+labs(y='Proportion')+theme(text=element_text(family='arial',size=12))+
  scale_fill_manual(values = paste0(Color2,'aa'))+scale_color_manual(values = Color2)


###microbiome pca====
total.data.microbiome<-data.frame(t(total.data.s.filter[12:233]))
colnames(total.data.microbiome)<-total.data.s.filter$Run

total.data.microbiome<-total.data.microbiome[rowSums(total.data.microbiome!=0)>1,]
total.data.s.filter$condition3<-ifelse(total.data.s.filter$tissue=='Endometrium','Patient_Endometrium',total.data.s.filter$condition)
dds.bulk.microbiome<-DESeqDataSetFromMatrix(total.data.microbiome,total.data.s.filter,~condition3)

table(rowSums(counts(dds.bulk.microbiome)>10))
Filter<-rowSums(counts(dds.bulk.microbiome)>10)>0
dds.bulk.microbiome<-dds.bulk.microbiome[Filter,]
dds.bulk.microbiome.norm<-varianceStabilizingTransformation(dds.bulk.microbiome)
p.microbiome.pca.filter<-plotPCA(dds.bulk.microbiome.norm,intgroup='condition3')+theme_pubr()+
  theme(text=element_text(family="arial",size=12),legend.title=element_blank())+
  geom_text_repel(aes(label= ifelse(dds.bulk.gene.norm$condition=='Control','Endometrium',dds.bulk.gene.norm$tissue)),color='black')
p.microbiome.pca.filter
#deg and pca plot====
total.data.gene<-data.frame()
for (i in total.data.s.filter$Run){
  file=paste0('/data2/msi/endometriosis/endometriosis/stringtie/',i,'.txt')
  buf<-read.table(file,sep = '\t',header=T)
  buf$Gene.Name<-ifelse(buf$Gene.Name=='.','-',buf$Gene.Name)
  buf$Gene<-paste0(buf$Gene.ID,'/',buf$Gene.Name)
  buf<-buf[c('Gene','FPKM')]
  colnames(buf)<-c('Gene',i)
  if (nrow(total.data.gene)==0){total.data.gene<-buf}
  else{total.data.gene<-merge(total.data.gene,buf,by='Gene',all = T)}
}

total.data.gene<-read.csv('/data2/msi/endometriosis/endometriosis/gene_count_matrix.csv',row.names = 'gene_id')
total.data.gene<-total.data.gene[total.data.s.filter$Run]
total.data.gene<-total.data.gene[-nrow(total.data.gene),]
buf<-rownames(total.data.gene)
buf<-strsplit(buf,'[|]')
for (i in 1:length(buf)){
  buf[[i]]<-buf[[i]][length(buf[[i]])]
}
total.data.gene$gene<-unlist(buf)

total.data.gene<-aggregate(.~gene,total.data.gene,sum)
rownames(total.data.gene)<-total.data.gene$gene
total.data.gene<-total.data.gene[,-1]
total.data.s.filter$condition3<-ifelse(total.data.s.filter$tissue=='Endometrium','Patient_Endometrium',total.data.s.filter$condition)
dds.bulk.gene<-DESeqDataSetFromMatrix(total.data.gene,total.data.s.filter,~condition3)
dds.bulk.tissue<-DESeqDataSetFromMatrix(total.data.gene,total.data.s.filter,~tissue)

Filter<-rowSums(counts(dds.bulk.gene)>10)>2
dds.bulk.gene<-dds.bulk.gene[Filter,]
dds.bulk.gene.norm<-vst(dds.bulk.gene)


Filter<-rowSums(counts(dds.bulk.tissue)>10)>2
dds.bulk.tissue<-dds.bulk.tissue[Filter,]
dds.bulk.tissue.norm<-vst(dds.bulk.tissue)



p.gene.pca.filter<-plotPCA(dds.bulk.tissue.norm,intgroup='condition3')+theme_pubr()+
  theme(text=element_text(family="arial",size=12),legend.title=element_blank())+
  geom_text_repel(aes(label= ifelse(dds.bulk.gene.norm$condition=='Control','Endometrium',dds.bulk.gene.norm$tissue)),color='black')
p.gene.pca.filter

dds.bulk.gene.norm$`Endometrium histology`

dds.bulk.gene<-DESeq(dds.bulk.gene)
dds.bulk.tissue<-DESeq(dds.bulk.tissue)

res.bulk.ConVSPat<-results(dds.bulk.gene,contrast = c('condition3','Patient_Endometrium','Control'))
res.bulk.umVSosis<-results(dds.bulk.gene,contrast = c('condition3','Endometriosis','Patient_Endometrium'))
res.bulk.gene<-results(dds.bulk.gene,contrast = c('condition3','Endometriosis','Control'))

res.tissue.gene<-results(dds.bulk.tissue,contrast = c('tissue','Endometrium','Peritoneum'))

res.tissue.gene$diffexpressed <- "NO"
res.tissue.gene$diffexpressed[res.tissue.gene$log2FoldChange > cut_lfc & res.tissue.gene$padj < cut_pvalue] <- "UP"
res.tissue.gene$diffexpressed[res.tissue.gene$log2FoldChange < -cut_lfc & res.tissue.gene$padj < cut_pvalue] <- "DOWN"
table(res.tissue.gene$diffexpressed)
deg.EnPe<-rownames(res.tissue.gene)[res.tissue.gene$diffexpressed!='NO']

ggVennDiagram(list(ConEn=deg.ConEn,ConOv=deg.ConOv,ConPe=deg.ConPe,EnOv=deg.EnOv,EnPe=deg.EnPe,OvPe=deg.OvPe))
##volcano plot----
cut_lfc <- 2
cut_pvalue <- 0.05
res.bulk.gene<-res.bulk.ConVSPat

res.bulk.gene$diffexpressed <- "None"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res.bulk.gene$diffexpressed[res.bulk.gene$log2FoldChange > cut_lfc & res.bulk.gene$padj < cut_pvalue] <- "Endometrium"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res.bulk.gene$diffexpressed[res.bulk.gene$log2FoldChange < -cut_lfc & res.bulk.gene$padj < cut_pvalue] <- "Control"
ggplot(data=data.frame(res.bulk.gene), aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) +
  geom_point() + theme_minimal() +
  theme(text=element_text(family="arial",size=12),legend.title=element_blank())+
  scale_color_manual(values=c("blue", "red", "black")) +
  geom_vline(xintercept=c(-cut_lfc, cut_lfc), col="red") +
  geom_hline(yintercept=-log10(cut_pvalue), col="red")#+
  geom_text_repel(aes(label= ifelse(res.bulk.gene$diffexpressed!='None',rownames(res.bulk.gene),'')),color='black')

library(ggplot2)
library(ggpubr)
library(ggrepel)
  
table(res.bulk.gene$diffexpressed)

buf<-total.data.s.filter$`Bacteria_s_Bacillus cereus`
buf<-buf/rowSums(total.data.s.filter[,Bac+Vir+Euc+Fun==1],)
ggplot(total.data.s.filter,aes(y=buf, x=condition3,fill=condition3,color=condition3))+geom_boxplot()+stat_compare_means(comparisons = list(c('Control','Endometriosis'),c('Endometriosis','Patient_Endometrium'),c('Control','Patient_Endometrium')))+theme_pubr()+
  labs(x='',y='relative abundance',fill='',color='')+scale_fill_manual(values = paste0(Color3,'99'))+scale_color_manual(values=Color3)

gn.selected<-res.bulk.gene$log2FoldChange>cut_lfc & res.bulk.gene$padj<cut_pvalue
buf3<-rownames(res.bulk.gene)[!is.na(gn.selected) & gn.selected]
for (i in buf4[buf4 %in% buf1]){show(i)}

library(ggVennDiagram)
ggVennDiagram(list(control_p=buf3,control_n=buf4, endometriosis_p=buf1, endometriosis_n=buf2))





#Visualization set----
#color
Color=c('Control'='#01D11A','Endometriosis'='#CB1106')
Color1=c('Control'='#01D11A99','Endometriosis'='#CB110699')
#diversity

p.filter.diversity<-ggplot(M.alpha.filter[M.alpha.filter$variable%in%c('Chao1','shannon'),],aes(x=condition,y=value,fill=condition,color=condition))+geom_boxplot()+theme_pubr()+scale_fill_manual(values = Color1)+
  facet_wrap(~variable,nrow=1,scales = 'free')+
  stat_compare_means(method='wilcox.test',size=3.5, family='arial',comparisons=list(c('Control','Endometriosis')))+
  theme(text=element_text(family="arial",face = 'bold',size=12), axis.text = element_text(face='plain'),
        axis.title.x = element_blank(),axis.text.x=element_blank(), axis.title.y=element_text(face='plain'),
        legend.title=element_blank(),legend.text=element_text(face='plain'))+
  labs(y='Score')+scale_color_manual(values=Color)

M.alpha.filter$condition3<-ifelse(M.alpha.filter$condition=='Control','Control',
                                  ifelse(M.alpha.filter$tissue=='Endometrium','Patient_Endometrium','Endometriosis'))
p.filter.diversity3<-ggplot(M.alpha.filter[M.alpha.filter$variable%in%c('Chao1','shannon'),],aes(x=condition3,y=value,fill=condition3,color=condition3))+geom_boxplot()+theme_pubr()+scale_fill_manual(values = paste0(Color3,'aa'))+
  facet_wrap(~variable,nrow=1,scales = 'free')+
  stat_compare_means(method='wilcox.test',size=3.5, family='arial',comparisons=list(c('Control','Endometriosis'),c('Control','Patient_Endometrium'),c('Patient_Endometrium','Endometriosis')))+
  theme(text=element_text(family="arial",face = 'bold',size=12), axis.text = element_text(face='plain'),
        axis.title.x = element_blank(),axis.text.x=element_blank(), axis.title.y=element_text(face='plain'),
        legend.title=element_blank(),legend.text=element_text(face='plain'))+
  labs(y='Score')+scale_color_manual(values=Color3)

p.filter.beta<-ggplot(mds.data.filter, aes(x=MDS1,y=MDS2,color=type))+geom_point()+theme_pubr()+
  theme(text=element_text(family="arial",size=12),legend.title=element_blank())+scale_color_manual(values = Color)+
  annotate("text", x=min(mds.data.filter$MDS1)+0.2, y=max(mds.data.filter$MDS2), label=paste0("p = ",pv.filter),family='arial',size=3.5)

mds.data.filter$type3<-total.data.s.filter$condition3
p.filter.beta3<-ggplot(mds.data.filter, aes(x=MDS1,y=MDS2,color=type3))+geom_point()+theme_pubr()+
  theme(text=element_text(family="arial",size=12),legend.title=element_blank())+scale_color_manual(values = Color3)+
  annotate("text", x=min(mds.data.filter$MDS1)+0.2, y=max(mds.data.filter$MDS2), label=paste0("p = ",pv.filter3),family='arial',size=3.5)

p.filter.diversity3+p.filter.beta3

#lefse
buf<-mm_test.filter.b@marker_table[tof]
buf$feature<-gsub('s__','',buf$feature)
p.lefse<-plot_ef_bar(buf)+theme_pubr()+scale_fill_manual(values = Color1)+theme(text=element_text(family="arial",face = 'plain',size=12), axis.text = element_text(face='plain'), axis.title.y=element_text(face='plain'),legend.title=element_blank(),legend.text=element_text(face='plain'))+labs(y='')+scale_color_manual(values=Color)

p.filter.diversity+p.filter.beta+p.lefse+plot_layout(widths=c(5,6,8))


p.lefse

cereus.count$type3<-factor(cereus.count$type3,levels=c('Control','Patient_Endometrium','Endometriosis'))
cereus.count$history<-dds.bulk.gene$`Endometrium histology`
cereus.count$patient<-dds.bulk.gene$source_name
cereus.count$tissue<-dds.bulk.gene$tissue
cereus.count.melt.history<-melt(cereus.count, id.vars = c('type','history','type3','patient','tissue'))
#bacillus cereus counts
ggplot(cereus.count.melt.history[cereus.count.melt.history$variable%in%c('Bacillus cereus'),], aes(x=history,y=value,color=tissue,shape=type3))+geom_point(size=10)+facet_wrap(~variable,scales = 'free')+theme_pubr()+theme(text = element_text(family='arial'))+labs(color='',fill='',y='Relative abundance (%)')

str(cereus.count.melt.history)

ggplot(cereus.count.melt.history[cereus.count.melt.history$variable%in%c('Bacillus cereus'),], aes(x=type3,y=value,fill=type3,color=type3))+geom_boxplot()+facet_wrap(~variable,scales = 'free')+theme_pubr()+scale_fill_manual(values=paste0(c(Color3[1],Color3[3],Color3[2]), 'aa'))+scale_color_manual(values=Color3)+theme(text = element_text(family='arial'),axis.title.x = element_blank())+labs(color='',fill='',y='Relative abundance (%)')
stat_compare_means(method='wilcox.test',label='p',family='arial',comparisons = list(c('Control','Patient_Endometrium'),c('Patient_Endometrium','Endometriosis'),c('Control','Endometriosis')))

ggplot(cereus.count.melt[cereus.count.melt$variable%in%c('Bacillus cereus','Tenacibaculum sp. GRR-S3-23'),], aes(x=type,y=value,fill=type,color=type))+geom_boxplot()+facet_wrap(~variable,scales = 'free')+theme_pubr()+scale_fill_manual(values=Color1)+scale_color_manual(values=Color)+stat_compare_means(method='wilcox.test',label='p',family='arial')+theme(text = element_text(family='arial'),axis.title.x = element_blank())+labs(color='',fill='',y='Relative abundance (%)')

hif1a<-cereus.count
hif1a$hif1a<-unlist(c(total.data.gene['HIF1A',]))
hif1a$hif1a<-log2(hif1a$hif1a+1)
hif1a.melt<-melt(hif1a, id.vars = c('type','history','type3','patient','tissue'))

ggplot(hif1a,aes(x=hif1a,y=`Bacillus cereus`,color=history,shape=type))+geom_point(size=3)+theme_pubr()


#cell cycle gene----
Color=c('Control'='#33B211','Endometriosis'='#CB1106')
library(readxl)
files<-list.files('/data2/msi/endometriosis/endometriosis/stringtie')
files<-files[startsWith(files,'SRR') & endsWith(files,'txt')]
#files<-paste0('/data2/msi/endometriosis/endometriosis/stringtie/',files)

total.data.gene<-data.frame()
for (i in files){ #with ensg
  buf<-read_tsv(paste0('/data2/msi/endometriosis/endometriosis/stringtie/',i))
  buf<-data.frame(buf)
  buf<-buf[c('Gene.ID','FPKM')]
  buf<-aggregate(buf[2],list(buf$`Gene.ID`),sum)
  colnames(buf)<-c('Gene',gsub('.txt','',i))
  if (sum(dim(total.data.gene))==0){ total.data.gene<- buf} else{total.data.gene<-merge(total.data.gene,buf,by='Gene',all=T)}
}

rownames(total.data.gene)<-total.data.gene$Gene
total.data.gene<-total.data.gene[-1,-1]
total.data.gene<-total.data.gene[,total.data.s.filter$Run]
table(rowSums(total.data.gene==0)<8)

cellcycle.table<-read_xlsx('/data2/msi/Rproject/endometriosis/cell_cycle_gene.xlsx')
cellcycle.table<-cellcycle.table[(cellcycle.table$ENSG %in% rownames(total.data.gene)),]
cellcycle.G1S<-cellcycle.table[cellcycle.table$Phase=='G1-S',]
cellcycle.G2M<-cellcycle.table[cellcycle.table$Phase=='G2-M',]


gene.expression<-total.data.gene
gene.expression$Phase<-'no'
for (i in 1:nrow(gene.expression)){
  gene.expression[i,'Phase']<-ifelse(rownames(gene.expression)[i] %in% cellcycle.table$ENSG,cellcycle.table[cellcycle.table$ENSG==rownames(gene.expression)[i],'Phase'],'no')
}

gene.expression$ensg<-rownames(gene.expression)
gene.expression<-gene.expression[gene.expression$Phase%in%c('G1-S','G2-M'),]
gene.expression$gene<-rownames(gene.expression)
for (i in 1:nrow(gene.expression)){
  gene.expression[i,'gene']<-ifelse(rownames(gene.expression)[i] %in% cellcycle.table$ENSG,cellcycle.table[cellcycle.table$ENSG==rownames(gene.expression)[i],'Symbol'],'no')
}
gene.expression<-gene.expression[rowSums(gene.expression[1:10]==0)<8,]

gene.expression.melt<-melt(gene.expression,id.vars = c('Phase','ensg','gene'))
gene.expression.melt$type<-ifelse(gene.expression.melt$variable %in% c('SRR19075268','SRR19075269','SRR19075270'),'Control','Endometriosis')
gene.expression.melt$log<-log2(gene.expression.melt$value+1)
gene.expression.melt<-gene.expression.melt[gene.expression.melt$Phase%in%c('G1-S','G2-M'),]

gene.expression1<-data.frame(t(total.data.gene))
gene.expression1<-gene.expression1[colnames(gene.expression1)%in%rownames(gene.expression)]
gene.expression1$type<-total.data.s.filter$condition

cellcycle.gene<-c()
for (i in colnames(gene.expression1)[-ncol(gene.expression1)]){
  kruskal.result<-kruskal.test(get(i)~type,gene.expression1)
  if (kruskal.result$p.value<0.06){cellcycle.gene<-c(cellcycle.gene,i)}
  }
  
ggplot(gene.expression.melt[gene.expression.melt$ensg%in%cellcycle.gene,],aes(x=gene,y=log,fill=type))+geom_split_violin()


a<-ggplot(gene.expression.melt[gene.expression.melt$ensg%in%cellcycle.gene[1:16],],aes(x=gene,y=log,fill=type,color=type))+geom_boxplot()+theme_pubr()+theme(axis.text.x=element_text(angle=90, hjust=1))+scale_fill_manual(values=paste0(Color,'aa'))+scale_color_manual(values=Color)+labs(x='',y='log2(FPKM)',fill='',color='')
b<-ggplot(gene.expression.melt[gene.expression.melt$ensg%in%cellcycle.gene[17:31],],aes(x=gene,y=log,fill=type,color=type))+geom_boxplot()+theme_pubr()+theme(axis.text.x=element_text(angle=90, hjust=1),legend.position = 'none')+labs(x='',fill='',y='log2(FPKM)')+scale_fill_manual(values=paste0(Color,'aa'))+scale_color_manual(values=Color)

ggplot(gene.expression.melt[gene.expression.melt$ensg%in%cellcycle.gene,],aes(x=type,y=log,fill=type,color=type))+geom_boxplot()+theme_pubr()+theme(axis.text.x=element_blank(),legend.position = 'right')+labs(x='',y='log2(FPKM)',fill='',color='')+scale_fill_manual(values=paste0(Color,'aa'))+scale_color_manual(values=Color)+facet_wrap(~gene,ncol=10,scales='free')+stat_compare_means()

#p.cellcycle.vln<-a/b/c
p.cellcycle.box<-a/b
p.cellcycle.box



#deg, go, kegg----
vol_plot<-function(res.bulk.gene,cut_lfc,cut_pvalue,upgroup,downgroup){
  res.bulk.gene$diffexpressed <- "None"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  res.bulk.gene$diffexpressed[res.bulk.gene$log2FoldChange > cut_lfc & res.bulk.gene$padj < cut_pvalue] <- upgroup
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  res.bulk.gene$diffexpressed[res.bulk.gene$log2FoldChange < -cut_lfc & res.bulk.gene$padj < cut_pvalue] <- downgroup
  #res.bulk.gene$diffexpressed[rownames(res.bulk.gene) %in% rownames(gene.expression)] <- "Cell Cycle"
  ggplot(data=data.frame(res.bulk.gene), aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) +
    geom_point() + theme_minimal() +
    theme(text=element_text(family="arial",size=12),legend.title=element_blank())+
    scale_color_manual(values=c('Control'="blue", 'Endometriosis'="red", 'None'="black", 'Cell Cycle'='green')) +
    geom_vline(xintercept=c(-cut_lfc, cut_lfc), col="red") +
    geom_hline(yintercept=-log10(cut_pvalue), col="red")+
    geom_text_repel(aes(label= ifelse(res.bulk.gene$diffexpressed=='Cell Cycle',gene.expression[rownames(gene.expression),'gene'],'')),color='black')
}
##count matrix----
library(DESeq2)

###without normalization=====
{
  ensg<-c()
  gene<-c()
  buf1<-read.csv('/data2/msi/endometriosis/endometriosis/gene_count_matrix.csv') # by stringtie
  buf1<-buf1[buf1$gene_id!="<class 'str'>",]
  for (x in buf1$gene_id){
    y<-unlist(strsplit(x,'[|]'))
    if (is.na(y[2])){
      ensg<-c(ensg,y[1])
      gene<-c(gene,y[1])}
    else{
      ensg<-c(ensg,y[1])
      gene<-c(gene,y[2])}
  }
  buf1$ensg<-ensg
  buf1$gene<-gene
  buf1<-buf1[c(26,27,2:25)]
  buf1<-buf1[!grepl("^MT-", buf1$gene), ]
  table(rowSums(buf1[3:26]))
  rna.expression.ensg<-buf1
} #make.timer_sample_ensg

rna.expression.ensg.filter<-rna.expression.ensg[c('ensg','gene',total.data.s.filter$Run)]

files<-list.files('/data2/msi/endometriosis/endometriosis/stringtie')
files<-files[startsWith(files,'SRR') & endsWith(files,'txt')]
#files<-paste0('/data2/msi/endometriosis/endometriosis/stringtie/',files)

###with normalization----
total.data.gene<-data.frame()
for (i in files){ #with ensg
  buf<-read_tsv(paste0('/data2/msi/endometriosis/endometriosis/stringtie/',i))
  buf<-data.frame(buf)
  buf<-buf[c('Gene.ID','FPKM')]
  buf<-aggregate(buf[2],list(buf$`Gene.ID`),sum)
  colnames(buf)<-c('Gene',gsub('.txt','',i))
  if (sum(dim(total.data.gene))==0){ total.data.gene<- buf} else{total.data.gene<-merge(total.data.gene,buf,by='Gene',all=T)}
}
total.data.gene<-total.data.gene[c('Gene',total.data.s.filter$Run)]


##deg====
ensg.to.gene<-data.frame(rna.expression.ensg[c('ensg','gene')],row.names = 'ensg')
#reference가 factor레벨에서 높아야함 (여기서는 control)
total.data.s.filter$condition<-as.factor(total.data.s.filter$condition)

###deg with gene----
buf<-rna.expression.ensg.filter
buf<-buf[-1]
buf<-aggregate(.~gene,data=buf,sum)

rownames(buf)<-buf$gene
buf<-buf[-1]
buf<-buf[rowSums(buf==0)<8,]

dds.gene<-DESeqDataSetFromMatrix(countData = buf, colData = total.data.s.filter, design = ~condition)

dds.gene<-DESeq(dds.gene)
dds.gene.res<-results(dds.gene, alpha=0.05)

plotMA(dds.gene)

vol_plot(dds.gene.res,1.5,0.05,'Endometriosis','Control')

buf<-data.frame(dds.gene.res)
buf<-buf[!is.na(buf$padj),]
buf<-buf[abs(buf$log2FoldChange)>2&buf$padj<0.05,]
buf<-data.frame(rownames(buf))
buf<-merge(buf,rna.expression.ensg.filter[c('gene','ensg')],by.x='rownames.buf.',by.y='gene',all.x=T)
buf<-merge(buf,total.data.gene,by.x='ensg',by.y='Gene',all.x=T)
buf1<-total.data.s.filter[Arc|Bac|Euk|Fun|Vir]
buf1<-buf1/rowSums(buf1)
buf1<-buf1['Bacteria_s_Bacillus cereus']
buf[nrow(buf)+1,]<-c('Bacillus','Bacillus cereus',buf1$`Bacteria_s_Bacillus cereus`)


###endometriosis gene corr with b.cereus----
rownames(buf)<-buf$ensg
corr.data<-rcorr(t(buf[3:12]),type = 'spearman')



cor.df.gene.filter<-data.frame(x='',y='',r=0,p=0)
cor.df.gene.filter<-cor.df.gene.filter[-1,]
for (i in colnames(corr.data$r)){
  for (j in rownames(corr.data$r)){
  cor.df.gene.filter[nrow(cor.df.gene.filter)+1,]<-c(i,j,corr.data$r[i,j],corr.data$P[i,j])
  }
}
cor.df.gene.filter$r<-as.numeric(cor.df.gene.filter$r)
cor.df.gene.filter$p<-as.numeric(cor.df.gene.filter$p)

cor.df.gene.filter<-cor.df.gene.filter[order(cor.df.gene.filter$r),]

buf<-cor.df.gene.filter
buf<-buf[buf$y=='Bacillus',]
buf<-buf[buf$x!='Bacillus',]
buf<-buf[!is.na(buf$p),]
buf$type<-ifelse(buf$r>0,'po','ne')
buf$p.sym<-ifelse(buf$p<0.05,round(buf$p,3),'ns')
buf1<-data.frame(rna.expression.ensg.filter[c('ensg','gene')],row.names = 'ensg')
buf$X<-buf1[buf$x,]

bacillus.corr.gene<-buf
bacillus.corr.gene[bacillus.corr.gene$p.sym!='ns','X']


ggplot(bacillus.corr.gene[bacillus.corr.gene$p.sym!='ns',], aes(y=x,x=r,fill=type))+ geom_bar(position="dodge", stat="identity")+scale_fill_manual(values=c('red','blue'))+
  theme_pubr()+theme(text=element_text(family="arial",size=12),legend.position = 'none')+geom_text(aes(label=p.sym,x=-r/abs(r)/12))+
  labs(y='')


#
###deg with ensg----
buf<-rna.expression.ensg.filter
rownames(buf)<-buf$ensg
buf<-buf[-c(2,1)]
buf<-buf[rowSums(buf==0)<8,]

dds.ensg<-DESeqDataSetFromMatrix(countData = buf, colData = total.data.s.filter, design = ~condition)

dds.ensg<-DESeq(dds.ensg)
dds.ensg.res<-results(dds.ensg, alpha=0.05)

plotMA(dds.ensg)

EnhancedVolcano(dds.ensg.res,
                lab = rownames(dds.ensg.res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 2,
                title='Endometriosis vs Control')

for (i in rownames(dds.ensg.res)){show(i)}

###Gene Ontology====
library('org.Hs.eg.db', character.only = TRUE)



buf<-dds.ensg.res[bacillus.corr.gene$x,]
gene_list<-buf$log2FoldChange
names(gene_list)<-rownames(buf)
gse.ensg <- gseGO(geneList=gene_list[order(gene_list,decreasing=T)], 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")



require(DOSE)
p.go.cereus<-dotplot(gse.ensg, showCategory=7,title='GO',split='ONTOLOGY') + facet_grid(ONTOLOGY~.sign, scale="free")+scale_fill_gradient(low='red',high='#ffafaf') + theme(axis.text.y=element_text(size=10),text = element_text(family = 'serif'))

p.go.cereus
###kegg pathway====
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df <- dds.ensg.res[ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df$Y = ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg <- gseKEGG(geneList     = kegg_gene_list,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

p.kegg.cereus<-dotplot(kegg, showCategory=7,title='KEGG') + facet_grid(.~.sign, scale="free")+scale_fill_gradient(low='red',high='#ffafaf') + theme(axis.text.y=element_text(size=10),text = element_text(family = 'serif'))


p.go.cereus+p.kegg.cereus

##b.cereus corr with specific gene----
rna.expression.ensg.filter
total.data.gene
###endometriosis gene====
endometriosis.gene<-c('IL8','ENA18','HNPA3','TGFB','IL10','IL12')
endometriosis.gene<-c("CXCL9","TENT5C","GADD45G","GH1","IGHG1","JCHAIN","IL23A","PTPRC","SERPINA3")
cxcr.gene<-ensg.to.gene$gene[startsWith(ensg.to.gene$gene,'CXCR')]
ccl.gene<-ensg.to.gene$gene[startsWith(ensg.to.gene$gene,'CCL')]
#https://www.researchgate.net/figure/An-overview-of-the-EMT-process-and-classic-gene-markers-The-EMT-is-a-dynamic-and_fig1_362300465
Emt.gene<-c('CDH1','CLDN1','CLDN3','COL1A1','CTNNB1','DSC','EPCAM','KRT7','KRT14','TJP1')
eMt.gene<-c('ACTA2','CDH2','FIB1','FN1','ITGB1','ITGB3','MMP2','MMP3','MMP9','PXN','S100A11','SNAI1','SNAI2','TBRIII','TWIST1','TGFB1','TGFB2','VIM','ZEB1','ZEB2')
#https://www.gsea-msigdb.org/gsea
emt.gene<-c("ABI3BP","ACTA2","ADAM12","ANPEP","APLP1","AREG","BASP1","BDNF","BGN","BMP1","CADM1","CALD1","CALU","CAP2","CAPG","CD44","CD59","CDH11","CDH2","CDH6","COL11A1","COL12A1","COL16A1","COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL5A3","COL6A2","COL6A3","COL7A1","COL8A2","COMP","COPA","CRLF1","CCN2","CTHRC1","CXCL1","CXCL12","CXCL6","CCN1","DAB2","DCN","DKK1","DPYSL3","DST","ECM1","ECM2","EDIL3","EFEMP2","ELN","EMP3","ENO2","FAP","FAS","FBLN1","FBLN2","FBLN5","FBN1","FBN2","FERMT2","FGF2","FLNA","FMOD","FN1","FOXC2","FSTL1","FSTL3","FUCA1","FZD8","GADD45A","GADD45B","GAS1","GEM","GJA1","GLIPR1","COLGALT1","GPC1","GPX7","GREM1","HTRA1","ID2","IGFBP2","IGFBP3","IGFBP4","IL15","IL32","IL6","CXCL8","INHBA","ITGA2","ITGA5","ITGAV","ITGB1","ITGB3","ITGB5","JUN","LAMA1","LAMA2","LAMA3","LAMC1","LAMC2","P3H1","LGALS1","LOX","LOXL1","LOXL2","LRP1","LRRC15","LUM","MAGEE1","MATN2","MATN3","MCM7","MEST","MFAP5","MGP","MMP1","MMP14","MMP2","MMP3","MSX1","MXRA5","MYL9","MYLK","NID2","NNMT","NOTCH2","NT5E","NTM","OXTR","PCOLCE","PCOLCE2","PDGFRB","PDLIM4","PFN2","PLAUR","PLOD1","PLOD2","PLOD3","PMEPA1","PMP22","POSTN","PPIB","PRRX1","PRSS2","PTHLH","PTX3","PVR","QSOX1","RGS4","RHOB","SAT1","SCG2","SDC1","SDC4","SERPINE1","SERPINE2","SERPINH1","SFRP1","SFRP4","SGCB","SGCD","SGCG","SLC6A8","SLIT2","SLIT3","SNAI2","SNTB1","SPARC","SPOCK1","SPP1","TAGLN","TFPI2","TGFB1","TGFBI","TGFBR3","TGM2","THBS1","THBS2","THY1","TIMP1","TIMP3","TNC","TNFAIP3","TNFRSF11B","TNFRSF12A","TPM1","TPM2","TPM4","VCAM1","VCAN","VEGFA","VEGFC","VIM","WIPF1","WNT5A")
migration.gene<-c("MSN","EPAS1","TXNIP","JUP","TGM2","THBS1","THBS1","THBS1","THBS1","CDH1","CDH1","TIMP3","TIMP3","TIMP3","IGFBP7","IGFBP7","CSTB","ARHGDIB","CCN1","VIM","CLDN4","TGFBI","PODXL","KRT18","BST2","KRT19","GJA1","EPB41L2","GALNT1","EPCAM","EFEMP1","CAPG","SRGN","SRGN","CD55","ST14","MEST","WWTR1","NEDD9","DKK3","ARL4C","ARL4C","NNMT","NNMT","TACSTD2","CTSH","FABP5","FXYD3","TRIM29","TNFAIP2","GSTM3","EPS8","SERPINE1","EFNB2","EFNB2","AXL","PALM2AKAP2","CLDN7","SPINT1","SIRPA","SOX9","SOX9","RHOBTB3","FBLN1","LOXL2","SLPI","CHST15","ANXA8","DAPK1","S100A4","CDH3","LAD1","DYNLT3","GALNT3","PPL","SCNN1A","ADIRF","LAMA3","EMP3","IGFBP6","ARHGAP29","CXADR","SH3YL1","LPAR1","TPM2","GNG11","AKR1C1","GLIPR1","S100A2","S100P","FOSL1","FGF2","CXCL1","UGT1A10","EEF1A2","VCAN","KCNK1","KCNK1","MT1G","PTGS2","","UGCG","PTGER4","ITGB4","ITGB4","KRT17","GJB3","SMURF2","CDA","SLC27A2","SLC27A2","WNT5A","NMU","UGT1A10","CST6","UGT1A10","PIR","HMGA2","H2BC12L","MT1X","UGT1A10","CLU","LGALS3","KRT7","COL18A1","NR2F2","PLIN2","COL6A2","AKR1C3","AGR2","TUBB6","KLF5","TSPAN4","TFPI2","TFPI2","PLAGL1","TM4SF1","TM4SF1","RBPMS","GATA3","GATA3","AKR1C1","CD24","","PKP3","MAGEA3","OLR1","MAPK13","IGFBP3","FYN","IL1A","RGS20","PTGES","FN1","AKAP12","KYNU","CCN1","TRIM29","TGM2","MT1H","TGM2","AKR1C1","FN1","COL4A2","COL4A2","COL4A1","MTUS1","IGFBP3","MT2A","MAP1B","KRT17","TUBA4A","LIMCH1","LIMCH1","LIMCH1","FN1","CRYBG1","DTX4","ARID5B","MAOA","MT1E","AHNAK2","SNAI2","SACS","TMEM30B","BICC1","F2RL1","RAC2","ASAH1","DCBLD2","GLIPR1","DKK3","TNFRSF21","MAGEA3","TM4SF1","UGT1A10","GJB3","VGLL1","RHOBTB3","JAG1","CD24","FN1","AKR1C1","ST14","MT1F","GLUL","KYNU","S100A6","PERP","RAB31","RAB31","LIMA1","TSPAN13","SQOR","ZNF22","RAB25","S100A14","GRAMD2B","PDGFC","IGF2BP2","DSE","TNFRSF21","PLAC8","FKBP11","FKBP11","ESRP1","RNF128","GALNT14","ERAP2","EHF","GPR87","VCAN","UGCG","ADAMTS1","SFN","AQP3"
)

{adhension.gene<-c("	CD58	",
            "	CD2	",
            "	VSIR	",
            "	IGSF11	",
            "	CD80	",
            "	CD274	",
            "	CD28	",
            "	CD86	",
            "	CTLA4	",
            "	ICOSLG	",
            "	ICOS ligand [KO:K06710]	",
            "	ICOS	",
            "	HLA-DMA	",
            "	HLA-DMB	",
            "	HLA-DOA	",
            "	HLA-DOB	",
            "	HLA-DPA1	",
            "	HLA-DPB1	",
            "	HLA-DQA1	",
            "	HLA-DQA2	",
            "	HLA-DQB1	",
            "	HLA-DRA	",
            "	HLA-DRB1	",
            "	HLA-DRB3	",
            "	HLA-DRB4	",
            "	HLA-DRB5	",
            "	CD4	",
            "	HLA-A	",
            "	HLA-B	",
            "	HLA-C	",
            "	HLA-F	",
            "	HLA-G	",
            "	HLA-E	",
            "	CD8A	",
            "	CD8B	",
            "	CD8B2	",
            "	PDCD1LG2	",
            "	CD276	",
            "	VTCN1	",
            "	PDCD1	",
            "	CD40	",
            "	CD40LG	",
            "	ALCAM	",
            "	CD6	",
            "	PVR	",
            "	CD226	",
            "	NECTIN2	",
            "	TIGIT	",
            "	ITGAL	",
            "	ITGB2	",
            "	ICAM1	",
            "	ICAM2	",
            "	ICAM3	",
            "	CD22	",
            "	PTPRC	",
            "	SIGLEC1	",
            "	SPN	",
            "	NECTIN3	",
            "	CLDN4	",
            "	CLDN3	",
            "	CLDN7	",
            "	CLDN19	",
            "	CLDN16	",
            "	CLDN14	",
            "	CLDN15	",
            "	CLDN17	",
            "	CLDN20	",
            "	CLDN11	",
            "	CLDN18	",
            "	CLDN22	",
            "	CLDN5	",
            "	CLDN10	",
            "	CLDN8	",
            "	CLDN6	",
            "	CLDN2	",
            "	CLDN1	",
            "	CLDN9	",
            "	CLDN23	",
            "	CLDN34	",
            "	CLDN25	",
            "	CLDN24	",
            "	CLDN12	",
            "	OCLN	",
            "	F11R	",
            "	JAM2	",
            "	JAM3	",
            "	ESAM	",
            "	CDH5	",
            "	PECAM1	",
            "	CD99	",
            "	CD99L2	",
            "	ITGAM	",
            "	SELPLG	",
            "	SELP	",
            "	ITGA4	",
            "	ITGB1	",
            "	ITGA9	",
            "	VCAM1	",
            "	ITGB7	",
            "	MADCAM1	",
            "	SELL	",
            "	CD34	",
            "	GLG1	",
            "	SELE	",
            "	NECTIN1	",
            "	CDH2	",
            "	NCAM1	",
            "	NCAM2	",
            "	L1CAM	",
            "	CADM1	",
            "	NEGR1	",
            "	NTNG1	",
            "	LRRC4C	",
            "	NTNG2	",
            "	LRRC4	",
            "	PTPRF	",
            "	LRRC4B	",
            "	SDC1	",
            "	SDC2	",
            "	SDC3	",
            "	SDC4	",
            "	ITGAV	",
            "	ITGB8	",
            "	ITGA8	",
            "	NRXN1	",
            "	NRXN2	",
            "	NRXN3	",
            "	NLGN1	",
            "	NLGN2	",
            "	NLGN3	",
            "	NLGN4X	",
            "	NLGN4Y	",
            "	PTPRS	",
            "	SLITRK1	",
            "	SLITRK2	",
            "	SLITRK4	",
            "	SLITRK5	",
            "	SLITRK6	",
            "	PTPRD	",
            "	SLITRK3	",
            "	CADM3	",
            "	NRCAM	",
            "	CNTN1	",
            "	PTPRM	",
            "	CNTN2	",
            "	NFASC	",
            "	CNTNAP1	",
            "	CNTNAP2	",
            "	MPZ	",
            "	MAG	",
            "	CDH1	",
            "	VCAN	",
            "	ITGA6	",
            "	CDH3	",
            "	CDH4	",
            "	CDH15	",
            "	NEO1	"
)
  adhension.gene<-gsub('\t','',adhension.gene)} #adhension.gene by kegg


#####################필요없어$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
cxcr.gene%in%ensg.to.gene$gene
Emt.gene%in%ensg.to.gene$gene
eMt.gene%in%ensg.to.gene$gene

dds.ensg.res[rownames(dds.ensg.res)%in%rownames(ensg.to.gene)[ensg.to.gene$gene%in%endometriosis.gene],] #0
dds.ensg.res[rownames(dds.ensg.res)%in%rownames(ensg.to.gene)[ensg.to.gene$gene%in%cxcr.gene],] #0
buf<-dds.ensg.res[rownames(dds.ensg.res)%in%rownames(ensg.to.gene)[ensg.to.gene$gene%in%emt.gene],] #56
buf[!is.na(buf$padj) & buf$padj<0.1,] #1
buf1<-ggplot(data.frame(buf),aes(x=log2FoldChange))+geom_histogram()+annotate('text',x=20,y=15,label=paste0('mean: ',round(mean(buf$log2FoldChange),3)))+labs(title='epithelial mesenchymal transition')
buf<-dds.ensg.res[rownames(dds.ensg.res)%in%rownames(ensg.to.gene)[ensg.to.gene$gene%in%migration.gene],] #58
buf[!is.na(buf$padj) & buf$padj<0.1,] #1
buf2<-ggplot(data.frame(buf),aes(x=log2FoldChange))+geom_histogram()+annotate('text',x=3,y=8,label=paste0('mean: ',round(mean(buf$log2FoldChange),3)))+labs(title='migration')
buf<-dds.ensg.res[rownames(dds.ensg.res)%in%rownames(ensg.to.gene)[ensg.to.gene$gene%in%adhension.gene],] #22
buf[!is.na(buf$padj) & buf$padj<0.1,] #1
buf3<-ggplot(data.frame(buf),aes(x=log2FoldChange))+geom_histogram()+annotate('text',x=23,y=7,label=paste0('mean: ',round(mean(buf$log2FoldChange),3)))+labs(title='adhension')


buf1/buf2/buf3
###########################필요없어 끝 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#corr gene with b.cereus----
dds.ensg.res[rownames(dds.ensg.res)%in%rownames(ensg.to.gene)[ensg.to.gene$gene%in%endometriosis.gene],]
bacillus.corr.gene$gene<-'none'
bacillus.corr.gene$gene<-ifelse(bacillus.corr.gene$X %in% emt.gene, 'emt', bacillus.corr.gene$gene)
bacillus.corr.gene$gene<-ifelse(bacillus.corr.gene$X %in% migration.gene, 'migration', bacillus.corr.gene$gene)
bacillus.corr.gene$gene<-ifelse(bacillus.corr.gene$X %in% adhension.gene, 'adhension', bacillus.corr.gene$gene)


buf<-bacillus.corr.gene[bacillus.corr.gene$gene!='none',]

buf<-merge(buf,total.data.gene[total.data.gene$Gene %in% buf$x,],by.x='x',by.y='Gene')
buf1<-melt(buf,id.vars = c('gene','X','p.sym'))
buf1<-buf1[-(1:20),]
buf1$value<-as.numeric(buf1$value)
buf1$type<-ifelse(buf1$variable %in% c('SRR19075268','SRR19075269','SRR19075270'),'Control','Endometriosis')
ggplot(buf1,aes(x=type,y=value,fill=type))+geom_boxplot()+facet_wrap(~gene+X,scales='free')+theme_pubr()+stat_compare_means(label.y.npc = 0.9)
ggplot(buf,aes(x=y,y=X,color=r))+geom_point(size=20)+theme_pubr()+
  scale_color_gradient2(high='blue',low='red',mid='gray',limits=c(-1,1))+geom_text(label=buf$p.sym,color='black')


buf<-data.frame(dds.gene.res)
buf<-buf[!is.na(buf$padj),]
buf<-buf[buf$padj<0.1,]
buf<-data.frame(rownames(buf))
buf<-merge(buf,rna.expression.ensg.filter[c('gene','ensg')],by.x='rownames.buf.',by.y='gene',all.x=T)
buf<-merge(buf,total.data.gene,by.x='ensg',by.y='Gene',all.x=T)
buf1<-total.data.s.filter[Arc|Bac|Euk|Fun|Vir]
buf1<-buf1/rowSums(buf1)
buf1<-buf1['Bacteria_s_Bacillus cereus']
buf[nrow(buf)+1,]<-c('Bacillus','Bacillus cereus',buf1$`Bacteria_s_Bacillus cereus`)


###endometriosis gene corr with b.cereus----
rownames(buf)<-buf$ensg
corr.data1<-rcorr(t(buf[3:12]),type = 'spearman')
dim(corr.data1$r)
table(rownames(corr.data1$r) %in% rownames(ensg.to.gene))
cor.df.gene.filter<-data.frame(x='',y='',r=0,p=0)
cor.df.gene.filter<-cor.df.gene.filter[-1,]
corr.data1$n<-corr.data1$n[c('Bacillus','ENSG00000115541'),]



for (i in colnames(corr.data1$n)){
  for (j in rownames(corr.data1$n)){
    cor.df.gene.filter[nrow(cor.df.gene.filter)+1,]<-c(i,j,corr.data1$r[i,j],corr.data1$P[i,j])
  }
}
buf<-cor.df.gene.filter
buf<-buf[!is.na(buf$p) & buf$y=='Bacillus' & buf$x!='Bacillus',]
buf$X<-ensg.to.gene[buf$x,'gene']

bacillus.corr.gene<-buf
bacillus.corr.gene$gene<-'none'
bacillus.corr.gene$gene<-ifelse(bacillus.corr.gene$X %in% emt.gene, 'emt', bacillus.corr.gene$gene)
bacillus.corr.gene$gene<-ifelse(bacillus.corr.gene$X %in% migration.gene, 'migration', bacillus.corr.gene$gene)
bacillus.corr.gene$gene<-ifelse(bacillus.corr.gene$X %in% adhension.gene, 'adhension', bacillus.corr.gene$gene)


buf<-bacillus.corr.gene[bacillus.corr.gene$gene!='none',]
buf$r<-as.numeric(buf$r)
buf$p<-as.numeric(buf$p)

buf<-merge(buf,total.data.gene[total.data.gene$Gene %in% buf$x,],by.x='x',by.y='Gene')
buf1<-melt(buf,id.vars = c('gene','X','p.sym'))
buf1<-buf1[-(1:20),]
buf1$value<-as.numeric(buf1$value)
buf1$type<-ifelse(buf1$variable %in% c('SRR19075268','SRR19075269','SRR19075270'),'Control','Endometriosis')
ggplot(buf1,aes(x=type,y=value,fill=type))+geom_boxplot()+facet_wrap(~gene+X,scales='free')+theme_pubr()+stat_compare_means(label.y.npc = 0.9)
ggplot(buf,aes(x=y,y=X,color=r))+geom_point(size=3)+theme_pubr()+facet_wrap(~gene)+
  scale_color_gradient2(high='blue',low='red',mid='gray',limits=c(-1,1))#+geom_text(label=round(buf$p,3),color='black')

buf1<-ggplot(buf[buf$gene=='emt',],aes(x=y,y=X,color=r))+geom_point(size=3)+theme_pubr()+labs(title='epithelial mesenchymal transition')+
  scale_color_gradient2(high='blue',low='red',mid='gray',limits=c(-1,1))#+geom_text(label=round(buf$p,3),color='black')
buf2<-ggplot(buf[buf$gene=='adhension',],aes(x=y,y=X,color=r))+geom_point(size=3)+theme_pubr()+labs(title='adhension')+
  scale_color_gradient2(high='blue',low='red',mid='gray',limits=c(-1,1))#+geom_text(label=round(buf$p,3),color='black')
buf3<-ggplot(buf[buf$gene=='migration',],aes(x=y,y=X,color=r))+geom_point(size=3)+theme_pubr()+labs(title='migration')+
  scale_color_gradient2(high='blue',low='red',mid='gray',limits=c(-1,1))#+geom_text(label=round(buf$p,3),color='black')
  
buf1+buf2+buf3
ensg.to.gene[startsWith(ensg.to.gene$gene,'IL1'),]
rna.expression.ensg.filter[rna.expression.ensg.filter$gene=='IL1B',]
buf$p<-round(buf$p,3)
write.csv(buf,'main_talbe.csv')


#infertility amh====
buf<-patient.stain#[patient.stain$type=='Endometriosis',]
buf$B.cereus<-buf$`Staining level`
buf$B.cereus<-ifelse(buf$type=='Normal','Control',buf$B.cereus)
buf$B.cereus<-factor(buf$B.cereus,levels = c('Control','-','±','+','++','+++'))

buf$stain<-factor(ifelse(buf$type=='Normal','Control',buf$`Histopathologic score`),levels=c('Control','0','1','2','3','4'))

buf$group<-ifelse(buf$B.cereus %in% c('-','±'),'negative','positive')
buf$group<-ifelse(buf$type=='Normal','Control',buf$group)
buf$group<-factor(buf$group,levels = c('Control','negative','positive'))

buf$infertility<-ifelse(buf$`Infertility Y/N`=='Y',1,0)
#buf<-buf[buf$B.cereus!='-',]
colnames(buf)


Color4<-c('#FFa5a5','#CB1106')
Color6<-c('#33B211','#FFe2e1',"#FFa5a5", "#FF6D6D", "#FF1D1D",'#CB1106')
Color3<-c('#33B211',"#FFa5a5",'#CB1106')
value<-'size' #size AMH CA-125 CA-19-9
p1<-ggplot(buf[buf$stain %in% c(1,2,3,4),],aes(x=stain, y=get(value),fill=B.cereus,color=B.cereus))+
  geom_boxplot()+theme_pubr(legend='none')+
  labs(x='B.cereus staining intensity',y='Size')+
  scale_fill_manual(values = paste0(Color2,'99'))+scale_color_manual(values = Color2)+
  theme(text=element_text(family='arial',size=12))
p1

p2<-ggplot(buf,aes(x=group, y=get(value), fill=group, color=group))+geom_boxplot()+theme_pubr(legend='none')+
  labs(y=value,x='Histopathologic group')+
  scale_fill_manual(values = paste0(Color3,'99'))+scale_color_manual(values = Color3)+
  theme(text=element_text(family='arial',size=12))+stat_compare_means(label='p',comparisons = list(c('negative','positive')))
p3<-ggplot(buf[buf$type=='Endometriosis',],
           aes(x=`Infertility Y/N`, y=get(value), 
               fill=`Infertility Y/N`, 
               color=`Infertility Y/N`))+
  geom_boxplot()+theme_pubr(legend='none')+
  labs(y=value,x='infertility')+
  scale_fill_manual(values = paste0(Color5,'99'))+scale_color_manual(values = Color5)+
  theme(text=element_text(family='arial',size=12))+stat_compare_means(label='p')

p1+p2+p3



value<-'Gravidity' #Gravidity Parity
p1<-ggplot(buf,aes(x=B.cereus, y=get(value),fill=B.cereus,color=B.cereus))+geom_violin()+theme_pubr(legend='none')+labs(x='Histopathologic score',y=value)+
  scale_fill_manual(values = paste0(Color2,'99'))+scale_color_manual(values = Color2)+
  theme(text=element_text(family='arial',size=12))
p2<-ggplot(buf,aes(x=group, y=get(value), fill=group, color=group))+geom_violin()+theme_pubr(legend='none')+labs(y=value)+
  scale_fill_manual(values = paste0(Color4,'99'))+scale_color_manual(values = Color4)+
  theme(text=element_text(family='arial',size=12))+stat_compare_means(label='p')

p1+p2

value<-'infertility'
p1<-ggplot(buf,aes(fill=`Infertility Y/N`,x=B.cereus))+geom_bar(position='dodge')+labs(y=value)+
  scale_fill_manual(values = paste0(Color,'99'))+theme_pubr()+
  theme(text=element_text(family='arial',size=12))
p2<-ggplot(buf,aes(fill=`Infertility Y/N`,x=group))+geom_bar(position='fill')+labs(y=value)+
  scale_fill_manual(values = paste0(Color,'99'))+theme_pubr()+
  theme(text=element_text(family='arial',size=12))
p1<-ggplot(buf,aes(fill=`Infertility Y/N`,x=group))+geom_bar(position='fill')+labs(y=value)+
  scale_fill_manual(values = paste0(Color,'99'))+theme_pubr()+
  theme(text=element_text(family='arial',size=12))
p1+p2

ggplot(buf,aes(x=`Infertility Y/N`,fill=group))+geom_bar()+labs(y=value)+
  theme_pubr()+
  theme(text=element_text(family='arial',size=12))


#stack bar chart
ggplot(buf,aes(fill=`Infertility Y/N`,x=group))+geom_bar(position='fill')+labs(y=value)+
  scale_fill_manual(values = paste0(Color,'99'))+theme_pubr()+
  theme(text=element_text(family='arial',size=12))
#pie chart
p1<-ggplot(buf,aes(fill=`Infertility Y/N`,x=''))+geom_bar(position='fill')+
  coord_polar('y')+labs(y=value)+facet_grid(group~type)+
  scale_fill_manual(values = paste0(Color,'99'))+theme_pubr()+
  theme(text=element_text(family='arial',size=12),axis.text.x=element_blank())#+
  #geom_text(aes(x=1.65,y=0.01,label=paste0(`Infertility Y/N`, '%')),size=4,position = position_stack(vjust=0.5))
Color5<-c('#002A6C','#C2113A')
p2<-ggplot(buf[buf$type=='Endometriosis',],aes(fill=`Infertility Y/N`,x=''))+geom_bar(position='fill',color='white')+
  coord_polar('y')+labs(y=value)+facet_wrap(~group)+
  scale_fill_manual(values = paste0(Color5,'ee'))+theme_void()+
  theme(text=element_text(family='arial',size=16),
        legend.title = element_text(size=12),
        axis.text.x=element_blank(),
        legend.position = 'bottom')+
  geom_text(data=buf[buf$Sample=='A1',],aes(x=1.6,y=0.09),label='18.3 %',size=4)+
  geom_text(data=buf[buf$Sample=='A1',],aes(x=1.6,y=0.58),label='81.7 %',size=4)

buf1<-aggregate(buf$infertility,by = list(buf$`Histopathologic score`),FUN = sum)
buf1[1,'x']<-0


ggplot(buf1,aes(x=Group.1,y=x,fill='#C2113A'))+geom_bar(stat = 'identity')+theme_pubr()+scale_fill_manual(values='#C2113A')+
  theme(text=element_text(family='arial',size=12),
        legend.position = 'none')+
  labs(y='Count of infertility cases', x='B. cereus staining intensity')


ggplot(buf,aes(x=group,y=AMH,fill=group))+geom_boxplot()+theme_pubr()
ggplot(buf,aes(x=`Infertility Y/N`,y=AMH,fill=`Infertility Y/N`))+geom_boxplot()+theme_pubr()+stat_compare_means()
#table(buf[buf$type=='Endometriosis' & buf$group=='positive','infertility'])
p2+stat_compare_means()
p2/p1+plot_layout(heights = c(1,2))


total.data.s.filter[c('condition','Bacteria_s_Bacillus cereus')]
mean(total.data.s.filter[4:10,'Bacteria_s_Bacillus cereus'])/mean(total.data.s.filter[1:3,'Bacteria_s_Bacillus cereus'])

#corr with xcell----
corr.raw.data<-immune.filter.xcell[4:10,12:ncol(immune.filter.xcell)]
corr.matrix<-rcorr(as.matrix(corr.raw.data),type='spearman')
cor.df.xcell.filter<-data.frame(x='',y='',r=0,p=0)
cor.df.xcell.filter<-cor.df.xcell.filter[-1,]
for (i in colnames(corr.matrix$r)){
  for (j in rownames(corr.matrix$r))
    cor.df.xcell.filter[nrow(cor.df.xcell.filter)+1,]<-c(i,j,corr.matrix$r[i,j],corr.matrix$P[i,j])
}
cor.df.xcell.filter$r<-as.numeric(cor.df.xcell.filter$r)
cor.df.xcell.filter$p<-as.numeric(cor.df.xcell.filter$p)

cor.df.xcell.filter<-cor.df.xcell.filter[order(cor.df.xcell.filter$r),]
buf<-cor.df.xcell.filter
buf<-buf[buf$x=='Bacteria_s_Bacillus cereus',]
buf<-buf[buf$y!='Bacteria_s_Bacillus cereus',]
buf$p.sym<-ifelse(buf$p<0.1,round(buf$p,3),'ns')
buf$p.sym<-ifelse(buf$p<0.1,round(buf$p,3),round(buf$p,3))
buf$type<-ifelse(buf$r>0,'po','ne')
buf$Y<-gsub('_',' ',buf$y)
buf<-buf[order(buf$r,decreasing = F),]
buf$Y<-factor(buf$Y,level=buf$Y)

ggplot(buf, aes(y=Y,x=r,fill=type))+ geom_bar(position="dodge", stat="identity")+scale_fill_manual(values=c('red','blue'))+
  theme_pubr()+theme(text=element_text(family="arial",size=12),legend.position = 'none')+geom_text(aes(label=p.sym,x=-r/abs(r)/12,y=Y))+
  labs(y='')

buf1<-immune.filter.xcell[c(3,5,12:ncol(immune.filter.xcell))]
immune.melt<-melt(buf1,id.vars = c('condition','tissue'))
immune.melt$variable<-gsub('_',' ',immune.melt$variable)
immune.melt<-immune.melt[immune.melt$variable!='Bacteria s Bacillus cereus',]
immune.melt$variable<-factor(immune.melt$variable,levels=buf$Y)
immune.melt$group<-ifelse(immune.melt$tissue=='Endometrium','Patient_Endometrium',immune.melt$condition)
immune.melt$group<-factor(immune.melt$group,levels=c('Control','Patient_Endometrium','Endometriosis'))

ggplot(immune.melt,aes(x=condition,y=value,fill=condition))+geom_boxplot()+theme_pubr()+
  facet_wrap(~variable,scales = 'free',ncol=4)+
  stat_compare_means(label='p',label.y.npc = 0.9)+
  scale_fill_manual(values=Color)
ggplot(immune.melt,aes(x=group,y=value,fill=group))+geom_boxplot()+theme_pubr()+
  facet_wrap(~variable,scales = 'free',ncol=8)+
  stat_compare_means(label='p',label.y.npc = 0.9,comparisons = list(c('Control','Patient_Endometrium')))+
  scale_fill_manual(values=Color3)+
  theme(axis.text.x=element_blank())

#dc markergene----
dc.table<-read.table('/data2/msi/endometriosis/endometriosis/xcell_DC_marker.csv')
dc.marker<-c()
for (i in dc.table$Celltype_Source_ID...of.genes){
  buf<-unique(strsplit(i,'\t')[[1]])
  title<-buf[1]
  gene<-buf[3:(length(buf)-1)]
  dc.marker<-c(dc.marker,setNames(list(gene),title))
}


dc_type<-'cDC' # none, a, c, i, p
dc_list<-dc.marker[names(dc.marker)[startsWith(names(dc.marker),dc_type)]]

library(ggvenn)
ggvenn(dc_list[c('cDC_HPCA_1','cDC_HPCA_2','cDC_HPCA_3')],show_elements = T,show_percentage = F)/

ggvenn(dc_list[c('cDC_NOVERSHTERN_1','cDC_NOVERSHTERN_2','cDC_NOVERSHTERN_3')],show_elements = T,show_percentage = F)

dc_count<-data.frame(row.names = unique(unlist(dc_list)),rep(0,38))
for (i in names(dc_list)){
  for (j in dc_list[[i]]){
    dc_count[j,]<-dc_count[j,]+1
  }
}

#항생제 치료 전 후 비교====
antibiotic.patient<-read_excel('/data2/msi/Rproject/endometriosis/data/EMS-antibiotic-241010.xlsx',sheet='patient')
antibiotic.patient<-antibiotic.patient[antibiotic.patient$...1!='67-B',]
antibiotic.patient$...1[67]<-'67'
antibiotic.control<-read_excel('/data2/msi/Rproject/endometriosis/data/EMS-antibiotic-241010.xlsx',sheet='control')

plotdata.antibiotic<-antibiotic.patient[antibiotic.patient$`Anti 사용여부`==1,]

plotdata.antibiotic<-melt(plotdata.antibiotic,measure.vars = c('pre-anti 통증점수','post-anti 통증점수'))
plotdata.antibiotic$variable<-ifelse(plotdata.antibiotic$variable=='pre-anti 통증점수','pre-anti','post-anti')
colnames(plotdata.antibiotic)
p1<-ggpaired(plotdata.antibiotic, x = "variable", y = "value",
         line.color = "gray", line.size = 0.4, color='variable',fill='variable')+
  stat_compare_means(paired = T,comparisons = list(c('pre-anti','post-anti')))+
  scale_fill_manual(values=c('#ff5555aa','#5555ffaa'))+scale_color_manual(values=c('#ff5555','#5555ff'))+
  labs(x='',y='Pain score')+theme(legend.position = 'none')
p1
buf1<-antibiotic.patient$`pre-anti 통증점수`
buf2<-antibiotic.patient$`post-anti 통증점수`
buf1<-buf1[!is.na(buf2)]
buf2<-buf2[!is.na(buf2)]
wilcox.test(buf1,buf2,paired=T,)
plotdata.antibiotic1<-antibiotic.patient[antibiotic.patient$`Anti 사용여부`==1,]
#plotdata.antibiotic1$diff.size<-plotdata.antibiotic1$`post-anti size`/plotdata.antibiotic1$size

plotdata.antibiotic1<-melt(plotdata.antibiotic1,measure.vars = c('size','post-anti size'))
plotdata.antibiotic1$variable<-ifelse(plotdata.antibiotic1$variable=='size','pre-anti','post-anti')
colnames(plotdata.antibiotic)
p2<-ggpaired(plotdata.antibiotic1, x = "variable", y = "value",
         line.color = "gray", line.size = 0.4, color='variable',fill='variable')+
  stat_compare_means(paired = TRUE,comparisons = list(c('pre-anti','post-anti')))+
  scale_fill_manual(values=c('#ff5555aa','#5555ffaa'))+scale_color_manual(values=c('#ff5555','#5555ff'))+
  labs(x='',y='Size')+theme(legend.position = 'none')

p1+p2
table(plotdata.antibiotic2$diff.size)
##paired 없이----
plotdata.antibiotic1$variable<-factor(plotdata.antibiotic1$variable,levels=c('pre-anti','post-anti'))
table(plotdata.antibiotic1$variable)
ggplot(plotdata.antibiotic1,aes(x=variable,y=value,color=variable,fill=variable))+geom_boxplot()+theme_pubr()+
  stat_compare_means(paired = TRUE,comparisons = list(c('pre-anti','post-anti')))+
  scale_fill_manual(values=c('#ff5555aa','#5555ffaa'))+scale_color_manual(values=c('#ff5555','#5555ff'))+
  labs(x='',y='Size')+theme(legend.position = 'none')+
ggplot(plotdata.antibiotic1,aes(x=variable,y=value,color=variable,fill=variable))+geom_violin()+theme_pubr()+
  stat_compare_means(paired = TRUE,comparisons = list(c('pre-anti','post-anti')))+
  scale_fill_manual(values=c('#ff5555aa','#5555ffaa'))+scale_color_manual(values=c('#ff5555','#5555ff'))+
  labs(x='',y='Size')+theme(legend.position = 'none')
plotdata.antibiotic$`Anti 사용기간(day)`

##size 변화 그래프
plotdata.antibiotic2<-antibiotic.patient[antibiotic.patient$`Anti 사용여부`==1,]
plotdata.antibiotic2$diff.size<-plotdata.antibiotic2$`post-anti size`/plotdata.antibiotic2$size-1.01
plotdata.antibiotic2<-plotdata.antibiotic2[order(plotdata.antibiotic2$diff.size,decreasing = T),]
plotdata.antibiotic2$`...1`<-factor(plotdata.antibiotic2$`...1`,levels=plotdata.antibiotic2$...1)
plotdata.antibiotic2$`Anti 종류`<-factor(plotdata.antibiotic2$`Anti 종류`)

ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size, fill='#ff5555'))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  labs(y='Size change after antibiotic treatment (%)',x='Patients')+
  theme(axis.text.x=element_blank(),legend.position = 'none',axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  scale_y_continuous(breaks = c(-0.01,-0.11,-0.21,-0.31),labels =c(0,-10,-20,-30))+
  geom_hline(yintercept = -0.01,color='#444444',linetype =2)

ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=age))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_gradient(low="red", high="gray")+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=BMI))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_gradient(low="red", high="gray")+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=`EMS stage`))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_gradient(low="red", high="gray")+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=`laterality`))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_gradient(low="red", high="gray")+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=`size`))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_gradient(low="red", high="gray")+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=`AMH`))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_gradient(low="red", high="gray")+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=`CA-125`))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_gradient(low="red", high="gray")+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=`CA-19-9`))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_gradient(low="red", high="gray")+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=`Gravidity`))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_gradient(low="red", high="gray")+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=`Parity`))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_gradient(low="red", high="gray")+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=`Infertility Y/N`))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_manual(values=c("red", "gray"))+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=`Anti 종류`))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=`Anti 사용기간(day)`))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_gradient(low="red", high="gray")+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=`pre-anti 통증점수`))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_gradient(low="red", high="gray")+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=`post-anti 통증점수`))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_gradient(low="red", high="gray")+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=`post-anti size`))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_gradient(low="red", high="gray")+
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=`Histopathologic scores`))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()

plotdata.antibiotic2<-antibiotic.patient[antibiotic.patient$`Anti 사용여부`==1,]
plotdata.antibiotic2$diff.size<-plotdata.antibiotic2$`post-anti size`/plotdata.antibiotic2$size

plotdata.antibiotic2$`Anti 종류`<-factor(plotdata.antibiotic2$`Anti 종류`)
ggplot(plotdata.antibiotic2,aes(x=`...1`,y=diff.size,fill=age))+geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  scale_fill_gradient(low="red", high="gray")+



#size 수정----
buf<-patient.stain[patient.stain$type=='Endometriosis','Histopathologic score']
antibiotic.patient$`Histopathologic scores`<-buf
antibiotic.patient$`Histopathologic scores`<-as.factor(antibiotic.patient$`Histopathologic scores`)
ggplot(antibiotic.patient[antibiotic.patient$`Histopathologic scores`!=0,],
       aes(x=`Histopathologic scores`,y=size,fill=`Histopathologic scores`,col=`Histopathologic scores`))+geom_boxplot()+
  theme_pubr()+scale_color_manual(values=Color2)+scale_fill_manual(values=paste0(Color2,'aa'))+
  theme(legend.position = 'none')
antibiotic.patient$`EMS stage`

