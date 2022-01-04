# Flow cytometry data analysys

library(ggplot2)
library(tidyverse)
library(lmerTest)
library(rstatix)
library(gridExtra)

# Import the text file dataset from selected directory
data <- read.delim("~/Directory/flow_cytometry_dataset.txt")

# get dataset ready for analisys
# Convert batch effects to a factor variable
data$batch=as.factor(data$batch)

# Convert age groups to a factor variable
data$age_group=as.factor(data$age_group)

# Calculate B cell to T cell ratio
data$bcell_tcell_ratio=data$bcells/data$cd3tcells

# Calculate CD4+ cell to CD8+ T cell ratio
data$cd4tcell_cd8tcell_ratio=data$cd4tcells/data$cd8tcells

# Reoder column to have batch effects in the last column of dataset
data=data[,c(1:15, 17:18, 16)]

# Filter for samples that have social rank information to create the social status dataset
data_status=data %>% filter (social_rank %in% c("HIGH", "MEDIUM", "LOW"))

# Orfer social ranks from high to low
data_status <- within(data_status, social_rank <- factor(social_rank, levels = c("HIGH", "MEDIUM", "LOW")))

#Dataset is now ready to start working

# Create Figure 1A histograms
# Complete datsaset histogram
ggplot(data, aes(x=age)) + geom_histogram(color="black", fill="black", alpha=0.3, binwidth =1.7 ) +
  xlab("Age") + ylab("Count") + theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), legend.position="none", 
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.y = element_text(size = 18), 
        axis.line = element_line(colour = "black"))

# By Sex Histogram
colors=c('dodgerblue4', 'darkorange1')
ggplot(data, aes(x=age, color=sex, fill=sex)) + xlab("Age") + ylab("Count") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none",
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.y = element_text(size = 18), 
        axis.line = element_line(colour = "black")) + facet_wrap(~sex) + 
  geom_histogram(data=subset(data,sex=="MALE"), fill="dodgerblue4", color="dodgerblue4", alpha=.4, binwidth = 2.5) + 
  geom_histogram(data=subset(data,sex=='FEMALE'), fill="darkorange1", color="darkorange1", alpha=.4, binwidth = 2.5)

# Run principal compinent analisys on immune cell typyes 
pca_data<-prcomp(data[,6:17], scale = TRUE)
summary(pca_data)

# Bind principal component one with data and rename PC1 column
data<-cbind(data, pca_data$x[,1])
data=data %>% rename(PC1= `pca_data$x[, 1]`)

# Run linear mixed model on PC1 as a function of sex and age
PC1_lmer=lmer(PC1 ~ sex + age + batch + (1|sample), data=data)
summary(PC1_lmer)

# Create figure 1B plot
# Subtract out partial residuals of batch effects from model to create residuals plot 
data$PC1_wout_batch = data$PC1 - (model.matrix(PC1_lmer)[,c("batch2","batch3")]%*%summary(PC1_lmer)$
                                    coef[c("batch2","batch3"),1])

#PC1_plot
ggplot(data, aes(x = age, y = PC1_wout_batch)) + xlab("Age") + ylab("PC1 (34.6%)") +
  geom_smooth(method="lm", color="black") + geom_point(color="dodgerblue4")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18),  axis.title.y = element_text(size = 18), 
        axis.line = element_line(colour = "black"))

# Evaluate the correlation of every cell type with PC1 (Pairwise comparison)
bcell_cor=cor.test(data$bcells, data$PC1)
cd3tcell_cor=cor.test(data$cd3tcells, data$PC1)
bcell_tcell_cor=cor.test(data$bcell_tcell_ratio, data$PC1)
cd4tcells_cor=cor.test(data$cd4tcells, data$PC1)
cd8tcells_cor=cor.test(data$cd8tcells, data$PC1)
cd4tcell_cd8tcell_cor=cor.test(data$cd4tcell_cd8tcell_ratio, data$PC1)
cd4tregs_cor=cor.test(data$cd4tregulatory, data$PC1)
cd8tregs_cor=cor.test(data$cd8tregulatory, data$PC1)
nkcells_cor=cor.test(data$nkcells, data$PC1)
class_monocytes_cor=cor.test(data$class_monocytes, data$PC1)
nonclass_monocytes_cor=cor.test(data$nonclass_monocytes, data$PC1)
inter_monocytes_cor=cor.test(data$inter_monocytes, data$PC1)

#Create Supplementary Table 
pearson_correlation=c(bcell_cor$estimate, cd3tcell_cor$estimate, bcell_tcell_cor$estimate, 
                      cd4tcells_cor$estimate, cd8tcells_cor$estimate, cd4tcell_cd8tcell_cor$estimate,
                      cd4tregs_cor$estimate, cd8tcells_cor$estimate, nkcells_cor$estimate,
                      class_monocytes_cor$estimate, nonclass_monocytes_cor$estimate, 
                      inter_monocytes_cor$estimate)
pearson_correlation=round(pearson_correlation, digits = 2)

correlation_p_value=c(bcell_cor$p.value, cd3tcell_cor$p.value, bcell_tcell_cor$p.value, 
                      cd4tcells_cor$p.value, cd8tcells_cor$p.value, cd4tcell_cd8tcell_cor$p.value,
                      cd4tregs_cor$p.value, cd8tcells_cor$p.value, nkcells_cor$p.value,
                      class_monocytes_cor$p.value, nonclass_monocytes_cor$p.value, 
                      inter_monocytes_cor$p.value)
correlation_p_value=p_round(correlation_p_value, digits = 2)

cells= c('CD20+ B Cells', 'CD3+ T Cells', 'B Cell/T Cell Ratio', 'CD4+ T Cells', 'CD8+ T Cells', 
         'CD4/CD8 T Cell Ratio', 'CD4+ T Regulatory Cells', 'CD8+ T Regulatory Cells', 
         'NK Cells', 'Classical Monocytes', 'Non-Classical Monocytes', 'Intermediate Monocytes')

correlation_table=data.frame(cells, pearson_correlation, correlation_p_value)  
correlation_table <- correlation_table[order(pearson_correlation),]
View(correlation_table)

# Create Figure 1C 
#PC1_correlations
ggplot(data=correlation_table, aes(x=pearson_correlation, y=reorder(cells,pearson_correlation))) + 
  geom_bar(stat="identity", fill='darkgoldenrod2') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", panel.background = element_blank(), axis.text=element_text(size=14),
        axis.title.x = element_text(size = 18), axis.title.y = element_blank(),
        axis.line = element_line(colour = "black")) + xlab("Correlation with PC1")

# Run linear mixed-effects models on cell ratios and cell types

# CD20+ B cell to CD3+ T cell ratio
bcell_tcell_ratio_lmer=lmer(bcell_tcell_ratio ~ age + sex + batch + (1|sample), data=data)
summary(bcell_tcell_ratio_lmer)

# CD20+ B cells
bcells_lmer=lmer(bcells ~ age + sex + batch + (1|sample), data=data)
summary(bcells_lmer)

# CD3+ T cells
cd3tcells_lmer=lmer(cd3tcells ~ age + sex + batch + (1|sample), data=data)
summary(cd3tcells_lmer)

# CD4+ T cells
cd4tcells_lmer=lmer(cd4tcells ~ age + sex + batch + (1|sample), data=data)
summary(cd4tcells_lmer)

# CD8+ T cells
cd8tcells_lmer=lmer(cd8tcells ~ age + sex + batch + (1|sample), data=data)
summary(cd8tcells_lmer)

# CD4+ to CD8+ T cell ratio
cd4_cd8_ratio_lmer=lmer(cd4tcell_cd8tcell_ratio ~ age + sex + batch + (1|sample), data=data)
summary(cd4_cd8_ratio_lmer)

# CD4+ T regulatory cells
cd4tregulatory_lmer=lmer(cd4tregulatory ~ age + sex + batch + (1|sample), data=data)
summary(cd4tregulatory_lmer)

# CD8+ T regulatory cells
cd8tregulatory_lmer=lmer(cd8tregulatory ~ age + sex + batch + (1|sample), data=data)
summary(cd8tregulatory_lmer)

# Natural Killer cells
nkcells_lmer=lmer(nkcells ~ age + sex + batch + (1|sample), data=data)
summary(nkcells_lmer)

# Classical Monocytes
class_monocytes_lmer=lmer(class_monocytes ~ age + sex + batch + (1|sample), data=data)
summary(class_monocytes_lmer)

# Intermediate Monocytes
inter_monocytes_lmer=lmer(inter_monocytes ~ age + sex + batch + (1|sample), data=data)
summary(inter_monocytes_lmer)

# Run linear mixed-effects models on cell ratios and cell types for age groups supplementary analysys
# B cells
bcells_agegroup_lmer=lmer(bcells ~ age_group + sex + batch + (1|sample), data=data)
summary(bcells_agegroup_lmer)

# CD3+ T cells
cd3tcells_agegroup_lmer=lmer(cd3tcells ~ age_group + sex + batch + (1|sample), data=data)
summary(cd3tcells_agegroup_lmer)

# CD8+ T cells
cd8tcells_agegroup_lmer=lmer(cd8tcells ~ age_group + sex + batch + (1|sample), data=data)
summary(cd8tcells_agegroup_lmer)

# CD4+ T regulatory cells
cd4tregulatory_agegroup_lmer=lmer(cd4tregulatory ~ age_group + sex + batch + (1|sample), data=data)
summary(cd4tregulatory_agegroup_lmer)

# run same models with group 2 as reference group
data$age_group <- relevel( data$age_group, ref="2" )

cd3tcells_agegroup_lmer=lmer(cd3tcells ~ age_group + sex + batch + (1|sample), data=data)
summary(cd3tcells_agegroup_lmer)

# CD8+ T cells
cd8tcells_agegroup_lmer=lmer(cd8tcells ~ age_group + sex + batch + (1|sample), data=data)
summary(cd8tcells_agegroup_lmer)

# CD4+ T regulatory cells
cd4tregulatory_agegroup_lmer=lmer(cd4tregulatory ~ age_group + sex + batch + (1|sample), data=data)
summary(cd4tregulatory_agegroup_lmer)

# run same models with group 3 as reference group
data$age_group <- relevel( data$age_group, ref="3" )

cd3tcells_agegroup_lmer=lmer(cd3tcells ~ age_group + sex + batch + (1|sample), data=data)
summary(cd3tcells_agegroup_lmer)

# CD8+ T cells
cd8tcells_agegroup_lmer=lmer(cd8tcells ~ age_group + sex + batch + (1|sample), data=data)
summary(cd8tcells_agegroup_lmer)

# CD4+ T regulatory cells
cd4tregulatory_agegroup_lmer=lmer(cd4tregulatory ~ age_group + sex + batch + (1|sample), data=data)
summary(cd4tregulatory_agegroup_lmer)

# Susbstract the effect of batch from every model and substitute every cell type
# proportion to work with resudials for better plot visualizations

# CD20+ B cell to CD3+ T cell residuals
data$bcell_tcell_ratio = data$bcell_tcell_ratio - 
  (model.matrix(bcell_tcell_ratio_lmer)[,c("batch2","batch3")]
                            %*%summary(bcell_tcell_ratio_lmer)$coef[c("batch2","batch3"),1]) 

# CD20+ B cells rediduals
data$bcells = data$bcells - (model.matrix(bcells_lmer)[,c("batch2","batch3")]
                                               %*%summary(bcells_lmer)$coef[c("batch2","batch3"),1]) 

# CD8+ T cells rediduals
data$cd8tcells = data$cd8tcells - (model.matrix(cd8tcells_lmer)[,c("batch2","batch3")]
                             %*%summary(cd8tcells_lmer)$coef[c("batch2","batch3"),1]) 

# CD4+ T regulatory cells 
data$cd4tregulatory = data$cd4tregulatory - (model.matrix(cd4tregulatory_lmer)[,c("batch2","batch3")]
                             %*%summary(cd4tregulatory_lmer)$coef[c("batch2","batch3"),1]) 

# CD8+ T regulatory cells residuals
data$cd8tregulatory = data$cd8tregulatory - (model.matrix(cd8tregulatory_lmer)[,c("batch2","batch3")]
                             %*%summary(cd8tregulatory_lmer)$coef[c("batch2","batch3"),1]) 

# Natural Killer cells residuals
data$nkcells = data$nkcells - (model.matrix(nkcells_lmer)[,c("batch2","batch3")]
                             %*%summary(nkcells_lmer)$coef[c("batch2","batch3"),1]) 

# Classical Monocytes resdiauls 
data$class_monocytes = data$class_monocytes - (model.matrix(class_monocytes_lmer)[,c("batch2","batch3")]
                             %*%summary(class_monocytes_lmer)$coef[c("batch2","batch3"),1]) 


# Generate figure 2 plots

# Figure 2A
#bcell_tcell_ratio_plot= ## MMW: same as above
ggplot(data, aes(x = age, y = bcell_tcell_ratio)) + xlab("Age") + 
  ylab("CD20+/CD3+ Ratio") + geom_smooth(method="lm", color="black") + 
  geom_point(color="dodgerblue4", alpha=0.8) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18),  axis.title.y = element_text(size = 18),
        axis.line = element_line(colour = "black"))
  
# Figure 2B
ggplot(data, aes(x = age, y = bcells)) + xlab("Age") + 
  ylab("% CD20+ B Cells") + geom_smooth(method="lm", color="black") + 
  geom_point(color="dodgerblue4", alpha=0.8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18),  axis.title.y = element_text(size = 18),
        axis.line = element_line(colour = "black"))  

# Figure 2C
ggplot(data, aes(x = age, y = cd8tcells)) + xlab("Age") + 
  ylab("% CD8+ T Cells") + geom_smooth(method="lm", color="black") + 
  geom_point(color="dodgerblue4", alpha=0.8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18),  axis.title.y = element_text(size = 18),
        axis.line = element_line(colour = "black"))  

# Figure 2D
ggplot(data, aes(x = age, y = cd4tcell_cd8tcell_ratio)) + xlab("Age") + 
  ylab("CD4+/CD8+ Ratio") + geom_smooth(method="lm", color="black") + 
  geom_point(color="dodgerblue4", alpha=0.8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18),  axis.title.y = element_text(size = 18),
        axis.line = element_line(colour = "black"))  

# Figure 2E
ggplot(data, aes(x = age, y = cd4tregulatory)) + xlab("Age") + 
  ylab("% CD4+ T Regulatory Cells") + geom_smooth(method="lm", color="black") + 
  geom_point(color="dodgerblue4", alpha=0.8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18),  axis.title.y = element_text(size = 18),
        axis.line = element_line(colour = "black"))  

# Figure 2F
ggplot(data, aes(x = age, y = cd8tregulatory)) + xlab("Age") + 
  ylab("% CD8+ T Regulatory Cells") + geom_smooth(method="lm", color="black") + 
  geom_point(color="dodgerblue4", alpha=0.8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18),  axis.title.y = element_text(size = 18),
        axis.line = element_line(colour = "black"))  

# Generate figure 3 plots

# Figure 3A
ggplot(data, aes(x = age, y = class_monocytes)) + xlab("Age") + 
  ylab("% CD14+/CD16-/HLADR+ Classical Monocytes") + geom_smooth(method="lm", color="black") + 
  geom_point(color="dodgerblue4", alpha=0.8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18),  axis.title.y = element_text(size = 15),
        axis.line = element_line(colour = "black")) 

# Figure 3B
ggplot(data, aes(x = age, y = inter_monocytes)) + xlab("Age") + 
  ylab("% CD14+/CD16+/HLADR+ Intermediate Monocytes") + geom_smooth(method="lm", color="black") + 
  geom_point(color="dodgerblue4", alpha=0.8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18),  axis.title.y = element_text(size = 14),
        axis.line = element_line(colour = "black")) 

# Figure 3C
ggplot(data, aes(x = age, y = nkcells)) + xlab("Age") + 
  ylab("% Natural Killer Cells") + geom_smooth(method="lm", color="black") + 
  geom_point(color="dodgerblue4", alpha=0.8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18),  axis.title.y = element_text(size = 18),
        axis.line = element_line(colour = "black")) 

# Eliminate the effect of batch + the effect of age from cell types where sex differences were detected 
# to work with their residuals
# extract age beta and procceed to substract the effect of batch and age from model

# CD4+ T cells
summary(cd4tcells_lmer)$coefficients
data$cd4tcellssex = data$cd4tcells - (model.matrix(cd4tcells_lmer)[,c("batch2","batch3")]
                      %*%summary(cd4tcells_lmer)$coef[c("batch2","batch3"),1]) - (data$age*-0.1013208) 

# CD4+ to CD8+ T cell ratio
summary(cd4_cd8_ratio_lmer)$coefficients
data$cd4tcell_cd8tcell_ratiosex = data$cd4tcell_cd8tcell_ratio - 
  (model.matrix(cd4_cd8_ratio_lmer)[,c("batch2","batch3")]
   %*%summary(cd4_cd8_ratio_lmer)$coef[c("batch2","batch3"),1]) - (data$age*-0.06352201) 

# Generate figure 4 plots

# Figure 4A
ggplot(data, aes(x=sex,y=cd4tcellssex)) +
  geom_boxplot(color=colors, alpha=0.4, fill=colors) + 
  xlab("Sex") + ylab("% CD4+ T Cells") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18),  axis.title.y = element_text(size = 18), 
        axis.line = element_line(colour = "black")) + geom_jitter(color= "dimgray", width=.1, alpha=0.7)

# Figure 4B
ggplot(data, aes(x=sex,y=cd4tcell_cd8tcell_ratiosex)) +
  geom_boxplot(color=colors, alpha=0.4, fill=colors) + 
  xlab("Sex") + ylab("CD4+/CD8+ T Cell Ratio") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18),  axis.title.y = element_text(size = 18), 
        axis.line = element_line(colour = "black")) + geom_jitter(color= "dimgray", width=.1, alpha=0.7)


# Work on the subset of data for which we have social rank data available

# Run linear mixed-effects models on CD4+ T regulatory cells to test the 
# interaction between age and social rank on levels of these cells

# CD4+ T regulatory cells
cd4tregulatory_social_lmer=lmer(cd4tregulatory ~ social_rank * age + sex + batch + (1|sample), 
                                data=data_status)
summary(cd4tregulatory_social_lmer)

# Relevel social ranks to make medium ranking inidividuals the reference group
data_status$social_rank <- relevel(data_status$social_rank, ref="MEDIUM" )

#run liner mixed-effects model
# CD4+ T regulatory cells
cd4tregulatory_social_lmer2=lmer(cd4tregulatory ~ social_rank * age + sex + batch + (1|sample), 
                                data=data_status)
summary(cd4tregulatory_social_lmer2)

# Substract the effect of batch on linear model in order to work with residuals for better plot 
# visualization
summary(cd4tregulatory_social_lmer)$coefficients
data_status$cd4tregulatory = data_status$cd4tregulatory - 
  (model.matrix(cd4tregulatory_social_lmer)[,c("batch2","batch3")]
                   %*%summary(cd4tregulatory_social_lmer)$coef[c("batch2","batch3"),1]) 

# Generate plot for the interaction
colors=c('darkorange1', 'darkseagreen4', 'darkorchid', 'darkorange1', 'darkseagreen4', 'darkorchid')
ggplot(data_status, aes(age, cd4tregulatory)) + 
  geom_point(aes(colour = social_rank)) +
  geom_smooth(method = "lm", color='black') + facet_wrap(~social_rank) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text=element_text(size=11), 
        axis.title.x = element_text(size = 16),  axis.title.y = element_text(size = 16), 
        axis.line = element_line(colour = "black"), legend.position = 'none') + 
  scale_colour_manual(values = colors) + xlab("Age") + ylab("% CD4+ T Regulatory Cells")


# Supplementary figures

# Figure 1
colors=c('cadetblue', 'cornsilk2', 'coral3', 'azure3')
# Make group 1 reference in order to have groups in increasing order
data <- within(data, age_group <- factor(age_group, levels = c("1", "2", "3", "4")))
# B Cells
p1<-ggplot(data, aes(x=age_group,y=bcells)) + geom_boxplot(color=colors, alpha=0.4, fill=colors) + 
  xlab("Age") + ylab("% CD20+ B Cells") +  theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), panel.background = element_blank(), 
                axis.text=element_text(size=14), axis.title.x = element_text(size = 18), 
                axis.title.y = element_text(size = 18), axis.line = element_line(colour = "black")) + 
  geom_jitter(color= "dimgray")
#CD3+ T Cells
p2<-ggplot(data, aes(x=age_group,y=cd3tcells)) + geom_boxplot(color=colors, alpha=0.4, fill=colors) + 
  xlab("Age") + ylab("% CD3+ T Cells") +  theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(), panel.background = element_blank(), 
                   axis.text=element_text(size=14), axis.title.x = element_text(size = 18), 
                     axis.title.y = element_text(size = 18), axis.line = element_line(colour = "black")) + 
  geom_jitter(color= "dimgray")
#CD8+ T Cells
p3<-ggplot(data, aes(x=age_group,y=cd8tcells)) + geom_boxplot(color=colors, alpha=0.4, fill=colors) + 
  xlab("Age") + ylab("% CD8+ T Cells") +  theme(panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), panel.background = element_blank(), 
                            axis.text=element_text(size=14), axis.title.x = element_text(size = 18), 
                  axis.title.y = element_text(size = 18), axis.line = element_line(colour = "black")) + 
  geom_jitter(color= "dimgray")
#CD4+ T Regulatory Cells
p4<-ggplot(data, aes(x=age_group,y=cd4tregulatory)) + geom_boxplot(color=colors, alpha=0.4, fill=colors) + 
  xlab("Age") + ylab("% CD4+ T Regulatory Cells") +  theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), panel.background = element_blank(), 
                          axis.text=element_text(size=14), axis.title.x = element_text(size = 18), 
                        axis.title.y = element_text(size = 18), axis.line = element_line(colour = "black")) + 
  geom_jitter(color= "dimgray")
grid.arrange(p1,p2,p3,p4, nrow=2)

#Supplementary Figure 2
# CD3+ T cells rediduals
data$cd3tcells = data$cd3tcells - (model.matrix(cd3tcells_lmer)[,c("batch2","batch3")]
                                   %*%summary(cd3tcells_lmer)$coef[c("batch2","batch3"),1]) 
ggplot(data, aes(x = age, y = cd3tcells)) + xlab("Age") + 
  ylab("% CD3+ T Cells") + geom_smooth(method="lm", color="black") + 
  geom_point(color="dodgerblue4", alpha=0.8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18),  axis.title.y = element_text(size = 18),
        axis.line = element_line(colour = "black"))  

#Supplementary Figure 3
ggplot(data, aes(x = age, y = cd4tcells)) + xlab("Age") + 
  ylab("% CD4+ T Cells") + geom_smooth(method="lm", color="black") + 
  geom_point(color="dodgerblue4", alpha=0.8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text=element_text(size=14), 
        axis.title.x = element_text(size = 18),  axis.title.y = element_text(size = 18),
        axis.line = element_line(colour = "black"))  




