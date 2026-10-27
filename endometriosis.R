#packages----
library(ggplot2)
library(reshape2)
library(ggpubr)
library(patchwork)
library(vegan)
library(gridExtra)
library(ALDEx2)
library(phyloseq)
library(pheatmap)
library(microbiomeMarker)
#library(tidyverse)
library(edgeR)
library(DESeq2)
library(Hmisc)
library(clusterProfiler)
library(org.Hs.eg.db)
#library(enrichplot)
library(ggrepel)
library(readxl)
library(car)

mb.path<-'microbiome_counts_matrix.csv' #mb_matrix
immune.estimated.path<-'immune_estimation_epic.csv' #immune_estimation_epic
supple.path<-'Supplementary_tables.xlsx' #provide in paper supplementary meterial

mb.df<-read.csv(mb.path,row.names = 1)
mb.matrix<-mb.df[,12:(ncol(mb.df))]
#Figure1----
alpha.s<-mb.df[,c('source_name','condition','tissue')]

mb.bacteria<-mb.matrix[,startsWith(colnames(mb.matrix),'Bacteria')]
alpha.s$Chao1<-estimateR(mb.bacteria)[1,]
alpha.s$shannon<-vegan::diversity(mb.bacteria,index="shannon")

##A.alpha-----
M.alpha.s<-melt(alpha.s,id.vars = c("tissue",'source_name','condition'))
M.alpha.s$variable<-factor(M.alpha.s$variable,levels=c('Chao1','shannon'))


p.diversity.chao<-ggplot(M.alpha.s[M.alpha.s$variable=='Chao1',],
                           aes(x=condition,y=value,fill=condition))+geom_boxplot()+theme_pubr()+
  scale_fill_manual(values = c('#44f','#f44'))+
  facet_wrap(~variable,nrow=1,scales = 'free')+
  stat_compare_means(method='wilcox.test',label = 'p',label.x.npc = 0.5,label.y.npc = 0.95, family='arial',)+
  theme(text=element_text(family="arial",size=12), strip.text=element_text(face='bold'),
        axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
        legend.title = element_blank())+
  labs(y='score')
p.diversity.chao
p.diversity.shannon<-ggplot(M.alpha.s[M.alpha.s$variable=='shannon',],
                              aes(x=condition,y=value,fill=condition))+geom_boxplot()+theme_pubr()+
  scale_fill_manual(values = c('#44f','#f44'))+
  facet_wrap(~variable,nrow=1,scales = 'free')+
  stat_compare_means(method='wilcox.test',label = 'p',label.x.npc = 0.5,label.y.npc = 0.9, family='arial',)+
  theme(text=element_text(family="arial",size=12), strip.text=element_text(face='bold'),
        axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
        legend.title = element_blank())+
  labs(y='score')

p.diversity.shannon

##B.beta====
set.seed(6212)
beta.dist<-vegdist(mb.bacteria,method='bray')
mds<-metaMDS(beta.dist)

mds.data<-as.data.frame(mds$points)

mds.data$type<-mb.df$condition
mds.data$type<-factor(ifelse(mds.data$type=='Control','Control','Endometriosis'),levels=c('Control','Endometriosis'))

set.seed(4230)
pv<-adonis2(mb.matrix[,12:(ncol(mb.matrix)-1)]~condition,mb.matrix, method='bray')$`Pr(>F)`[1] #p value calculate (Pr(>F))

p.filter.beta<-ggplot(mds.data, aes(x=MDS1,y=MDS2,color=type))+geom_point()+theme_pubr()+
  theme(text=element_text(family="arial",size=12),
        legend.title=element_blank())+
  scale_color_manual(values = c('#44f','#f44'))+
  annotate("text", x=min(mds.data$MDS1)+0.4, y=max(mds.data$MDS2), label=paste0("p = ",pv),family='arial')
p.filter.beta

##C.lefse====
taxamet<-data.frame(Kingdom="", Phylum='', Class='', Order='', Family='', Genus='',Species='')
taxamet<-taxamet[-1,]


for (i in colnames(mb.bacteria)){
  taxamet[i,1:6]<-'-'
  taxamet[i,'Species']<-strsplit(i,'_')[[1]][3]
}

#total
OTU = otu_table(t(mb.bacteria), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxamet))
SAM = sample_data(data.frame(mb.df[1:11]))
physeq = phyloseq(OTU, TAX, SAM)
sample_names(SAM)
sample_names(OTU)

rownames(physeq.s@otu_table)<-1:nrow(physeq.s@otu_table)
rownames(physeq.s@tax_table)<-1:nrow(physeq.s@tax_table)

mm_test<- run_lefse( physeq,
                     #wilcoxon_cutoff = 0.05,
                     kw_cutoff = 0.06,
                     norm = "CPM", 
                     group = "condition",
                     lda_cutoff = 2.5,
                     taxa_rank = 'all',
                     multigrp_strat = TRUE
)
mm_test@marker_table

p.lefse<-plot_ef_bar(mm_test)+theme_pubr()+scale_fill_manual(values = c('#44f','#f44'))+labs(fill='')
p.lefse

##D.immune estimated score----
#immune.estimated.df<-read.csv('/data1/endometriosis/immune_estimation_epic.csv',row.names = 1)
immune.estimated.df<-read.csv(immune.estimated.path,row.names = 1)
corr.matrix<-rcorr(as.matrix(immune.estimated.df[4:10,12:ncol(immune.estimated.df)]),type='spearman')
corr.df<-data.frame(x='',y='',r=0,p=0)[-1,]
for (i in colnames(corr.matrix$r)){
  for (j in rownames(corr.matrix$r)){
    corr.df[nrow(corr.df)+1,]<-c(i,j,corr.matrix$r[i,j],corr.matrix$P[i,j])
  }
}

corr.df<-corr.df[order(corr.df$r),]
corr.df<-corr.df[corr.df$x=='Bacteria_s_Bacillus.cereus',]
corr.df<-corr.df[corr.df$y!='Bacteria_s_Bacillus.cereus',]
corr.df$r<-as.numeric(corr.df$r)
corr.df$p<-as.numeric(corr.df$p)

corr.df$p.sym<-ifelse(corr.df$p<0.1,round(corr.df$p,3),round(corr.df$p,3))
corr.df$p.sym<-ifelse(corr.df$p.sym==0,'<0.001',corr.df$p.sym)
corr.df<-corr.df[order(corr.df$r,decreasing = F),]
corr.df$y<-factor(corr.df$y,levels=corr.df$y)

p.immune.corr<-ggplot(corr.df, aes(y=y,x=r))+geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values='#eee')+
  theme_pubr()+theme(text=element_text(family="arial"),legend.position = 'none')+
  geom_text(aes(label=p.sym,x=(-(r)/abs(r)/5),y=y))+
  labs(y='',x='R')
p.immune.corr

#Figure2----
##B.pie----
stain.df<-read_xlsx(supple.path, sheet = 'Supplementary Table 3')
colnames(stain.df)<-stain.df[1,]
stain.df<-stain.df[-1,]

pie.df<-data.frame(Control=table(stain.df$`B. cereus staining intensity*`,stain.df$Cohort)[,1],
                Endometriosis=table(stain.df$`B. cereus staining intensity*`,stain.df$Cohort)[,2],
                row.names = c(0,1,2,3,4))
pie.df$lv<-c('0','1','2','3','4')
pie.df.melt<-melt(pie.df,id.vars ='lv')
pie.df.melt$label<-pie.df.melt$value
pie.df.melt[1:5,'label']<-round(pie.df.melt[1:5,'label']/sum(pie.df.melt[1:5,'label']),3)*100
pie.df.melt[6:10,'label']<-round(pie.df.melt[6:10,'label']/sum(pie.df.melt[6:10,'label']),3)*100
pie.df.melt$percent<-ifelse(paste0(pie.df.melt$label,'%')=='0%','',paste0(pie.df.melt$label,'%'))
pie.df.melt$position<-0
pie.df.melt[1:5,'position']<-sum(pie.df.melt[1:5,'value'])
j<-0
for (i in 1:5){
  j <- j + pie.df.melt[i,'value']/2
  pie.df.melt[i,'position']<-pie.df.melt[i,'position']-j
  j <- j + pie.df.melt[i,'value']/2
}
pie.df.melt[6:10,'position']<-sum(pie.df.melt[6:10,'value'])
j<-0
for (i in 6:10){
  j <- j + pie.df.melt[i,'value']/2
  pie.df.melt[i,'position']<-pie.df.melt[i,'position']-j
  j <- j + pie.df.melt[i,'value']/2
}
pie.df.melt
p.ihc.level<-ggplot(pie.df.melt,aes(fill=lv,x='lv',y=value))+geom_bar(stat='identity',color='black')+coord_polar('y')+
  theme_void()+ facet_wrap(~variable,scale='free',ncol=1)+
  geom_text(aes(x=1.7,y=position,label = percent),
            size = 4,) +
  scale_fill_manual(values=c('#EFF3FF','#BDD7E7','#6BAED6','#3182BD','#08519C'))+
  theme(strip.text=element_text(size=12), legend.title=element_blank())
p.ihc.level

##C&D.antibiotic treatment====
treatment.df<-read_xlsx(supple.path,sheet = 'Supplementary Table 3')
colnames(treatment.df)<-treatment.df[1,]
treatment.df<-treatment.df[-1,]
treatment.df<-treatment.df[-nrow(treatment.df),]
treatment.df<-treatment.df[treatment.df$`Antibiotic\r\nuse`=='Yes',]
treatment.df$`Post-antibiotic\r\npain score`<-as.numeric(treatment.df$`Post-antibiotic\r\npain score`)
treatment.df$`Pre-antibiotic\r\npain score`<-as.numeric(treatment.df$`Pre-antibiotic\r\npain score`)
treatment.df.melt<-melt(treatment.df,measure.vars = c('Pre-antibiotic\r\npain score','Post-antibiotic\r\npain score'))


treatment.df.melt$variable<-ifelse(treatment.df.melt$variable=='Pre-antibiotic\r\npain score','pre-anti','post-anti')

p.treatment.pain<-ggpaired(treatment.df.melt, x = "variable", y = "value",
             line.color = "gray", line.size = 0.4,fill='variable')+
              stat_compare_means(paired = TRUE,comparisons = list(c('pre-anti','post-anti')))+
              scale_fill_manual(values=c('#ff5555aa','#5555ffaa'))+
              labs(x='',y='Size')+theme(legend.position = 'none')
p.treatment.pain


treatment.df$`Post-antibiotic\r\nsize (cm)`<-as.numeric(treatment.df$`Post-antibiotic\r\nsize (cm)`)
treatment.df$`Pre-antibiotic\r\nsize (cm)`<-as.numeric(treatment.df$`Pre-antibiotic\r\nsize (cm)`)
treatment.df$diff.size<-treatment.df$`Post-antibiotic\r\nsize (cm)`/treatment.df$`Pre-antibiotic\r\nsize (cm)`-1.01
treatment.df$order<-factor(treatment.df$PID,
                           levels=as.list(treatment.df[order(treatment.df$diff.size,decreasing = T),'PID'])$PID)

p.treatment.size<-ggplot(treatment.df,aes(x=order,y=diff.size,
                                fill=''))+
                  geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
                  labs(y='Size change after antibiotic treatment (%)',x='Patients')+
                  scale_fill_manual(values='#ff5555')+
                  theme(axis.text.x=element_blank(),
                        legend.position = 'none',
                        axis.line.x = element_blank(),
                        axis.ticks.x=element_blank())+
                  scale_y_continuous(breaks = c(-0.01,-0.11,-0.21,-0.31),labels =c(0,-10,-20,-30))+
                  geom_hline(yintercept = -0.01,color='#444444',linetype =2)
p.treatment.size
##E&F.clinical variance====
stain.df.patient<-stain.df[stain.df$Cohort=='Endometriosis',]
stain.df.patient$Infertility<-ifelse(stain.df.patient$Infertility=='Yes',1,0)
stain.df.patient$size<-as.numeric(stain.df.patient$`Pre-antibiotic\r\nsize (cm)`)
stain.df.patient<-stain.df.patient[-nrow(stain.df.patient),]

p.ihc.size<-ggplot(stain.df.patient[stain.df.patient$`B. cereus staining intensity*`!='0',],
                   aes(x=`B. cereus staining intensity*`, y=size))+geom_boxplot()+
            theme_pubr(legend='none')+labs(x='Histopathologic score',y='Size (cm)')+
            theme(text=element_text(family='arial',size=12))
p.ihc.size

p.ihc.infertility<-ggplot(stain.df.patient[stain.df.patient$Infertility==1,],
                          aes(fill=Infertility,x=`B. cereus staining intensity*`))+
                    geom_bar(position='dodge')+labs(y='infertility')+
                    theme_pubr()+
                    scale_x_discrete(limits = c('0', '1', '2', '3', '4'))+
                    theme(text=element_text(family='arial',size=12),
                          legend.position = 'none')
p.ihc.infertility

#Figure3----
##B.pie----
stain.df<-read_xlsx(supple.path, sheet = 'Supplementary Table 5')
colnames(stain.df)<-stain.df[1,]
stain.df<-stain.df[-1,]

pie.df<-data.frame(Endometriosis=table(stain.df$`B. cereus staining intensity*`),
                   row.names = c(0,1,2,3,4))
pie.df$lv<-c('0','1','2','3','4')
pie.df.melt<-pie.df
pie.df.melt$label<-as.numeric(pie.df.melt$Endometriosis.Freq)
pie.df.melt[1:5,'label']<-round(pie.df.melt[1:5,'label']/sum(pie.df.melt[1:5,'label']),3)*100

pie.df.melt$percent<-ifelse(paste0(pie.df.melt$label,'%')=='0%','',paste0(pie.df.melt$label,'%'))
pie.df.melt$position<-0
pie.df.melt[1:5,'position']<-sum(pie.df.melt[1:5,'Endometriosis.Freq'])
j<-0
for (i in 1:5){
  j <- j + pie.df.melt[i,'Endometriosis.Freq']/2
  pie.df.melt[i,'position']<-pie.df.melt[i,'position']-j
  j <- j + pie.df.melt[i,'Endometriosis.Freq']/2
}
pie.df.melt
p.ihc.level1<-ggplot(pie.df.melt,aes(fill=lv,x='lv',y=Endometriosis.Freq))+
  geom_bar(stat='identity',color='black')+coord_polar('y')+
  theme_void()+ 
  geom_text(aes(x=1.7,y=position,label = percent),
            size = 4,) +
  scale_fill_manual(values=c('#EFF3FF','#BDD7E7','#6BAED6','#3182BD','#08519C'))+
  theme(strip.text=element_text(size=12), legend.title=element_blank())
p.ihc.level1

##C&D.antibiotic treatment====
treatment.df<-read_xlsx(supple.path,sheet = 'Supplementary Table 5')
colnames(treatment.df)<-treatment.df[1,]
treatment.df<-treatment.df[-1,]
treatment.df<-treatment.df[-nrow(treatment.df),]
treatment.df<-treatment.df[treatment.df$`Antibiotic\r\nuse`=='Yes',]
treatment.df$`Post-antibiotic\r\npain score`<-as.numeric(treatment.df$`Post-antibiotic\r\npain score`)
treatment.df$`Pre-antibiotic\r\npain score`<-as.numeric(treatment.df$`Pre-antibiotic\r\npain score`)
treatment.df.melt<-melt(treatment.df,measure.vars = c('Pre-antibiotic\r\npain score','Post-antibiotic\r\npain score'))


treatment.df.melt$variable<-ifelse(treatment.df.melt$variable=='Pre-antibiotic\r\npain score','pre-anti','post-anti')

p.treatment.pain1<-ggpaired(treatment.df.melt, x = "variable", y = "value",
                           line.color = "gray", line.size = 0.4,fill='variable')+
  stat_compare_means(paired = TRUE,comparisons = list(c('pre-anti','post-anti')))+
  scale_fill_manual(values=c('#ff5555aa','#5555ffaa'))+
  labs(x='',y='Size')+theme(legend.position = 'none')
p.treatment.pain1

treatment.df$`Post-antibiotic\r\nsize (cm)`<-as.numeric(treatment.df$`Post-antibiotic\r\nsize (cm)`)
treatment.df$`Pre-antibiotic\r\nsize (cm)`<-as.numeric(treatment.df$`Pre-antibiotic\r\nsize (cm)`)
treatment.df$diff.size<-treatment.df$`Post-antibiotic\r\nsize (cm)`/treatment.df$`Pre-antibiotic\r\nsize (cm)`-1.01
treatment.df$order<-factor(treatment.df$PID,
                           levels=as.list(treatment.df[order(treatment.df$diff.size,decreasing = T),'PID'])$PID)

p.treatment.size1<-ggplot(treatment.df,aes(x=order,y=diff.size,
                                          fill=''))+
  geom_bar(stat = "identity", position = "dodge")+theme_pubr()+
  labs(y='Size change after antibiotic treatment (%)',x='Patients')+
  scale_fill_manual(values='#ff5555')+
  theme(axis.text.x=element_blank(),
        legend.position = 'none',
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank())+
  scale_y_continuous(breaks = c(-0.01,-0.11,-0.21,-0.31),labels =c(0,-10,-20,-30))+
  geom_hline(yintercept = -0.01,color='#444444',linetype =2)
p.treatment.size1
##E&F.clinical variance====
stain.df$Infertility<-ifelse(stain.df$Infertility=='Yes',1,0)
stain.df$size<-as.numeric(stain.df$`Pre-antibiotic\r\nsize (cm)`)
stain.df<-stain.df[-nrow(stain.df),]

p.ihc.size1<-ggplot(stain.df[stain.df$`B. cereus staining intensity*`!='0',],
                   aes(x=`B. cereus staining intensity*`, y=size))+geom_boxplot()+
  theme_pubr(legend='none')+labs(x='Histopathologic score',y='Size (cm)')+
  theme(text=element_text(family='arial',size=12))
p.ihc.size1

p.ihc.infertility1<-ggplot(stain.df[stain.df$Infertility==1,],
                          aes(fill=Infertility,x=`B. cereus staining intensity*`))+
  geom_bar(position='dodge')+labs(y='infertility')+
  theme_pubr()+
  scale_x_discrete(limits = c('0', '1', '2', '3', '4'))+
  theme(text=element_text(family='arial',size=12),
        legend.position = 'none')
p.ihc.infertility1
