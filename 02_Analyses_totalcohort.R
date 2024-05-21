# Analyses IBD-prediction project
# Q3-Q4 2023, Femke Prins & Iwan Hidding 
# Script for merged cohort analyses

#### Data input 
Meta_merged # Clinical phenotypes
Tax_merged_baseline # All microbial taxonomic data, not filtered
Merged_species # Only species, filtered
Merged_genera # Only genera, filtered
Merged_CLR_1s # CLR transformed abundance of species
Selected_metabolites_merged # metabolites filtered on 70% presence

#### 6 MONTHS RESPONSE ---- -------------------------------------------
# Species: alpha diversity ----
calculate_alpha <- function(inDF,IDcol="RowNames",
                            metrics=c("shannon","simpson","invsimpson","richness"),
                            DIVlvls=c("taxS")) {
  # select IDs column
  if (IDcol == "RowNames") {
    DIVMatrix <- data.frame(RN=rownames(inDF))
  } else {
    DIVMatrix <- data.frame(IDcol=inDF[[IDcol]])
    colnames(DIVMatrix)[1] <- IDcol
  }
  # iterate over metrics, calculate each
  # NOTE: richness is not implemented in vegan, requires special treatment
  for (l in DIVlvls) {
    if (grepl('^tax.$',l)) {
      toUse <- gsub('^tax','',l) } }
  for (m in metrics) {
    print(paste0('  > calculating ',m,'[',l,']'))
    if (m=="richness") {
      inDFpa <- inDF
      inDFpa[inDFpa > 0] <- 1
      dv <- rowSums(inDFpa)
    } else {
      dv <- diversity(inDF,index = m)
    }
    DIVMatrix <- cbind.data.frame(DIVMatrix,dv)
    colnames(DIVMatrix)[ncol(DIVMatrix)] <- paste0('DIV.',toUse,'.',m)
  }
  return(DIVMatrix)
}
Alpha_metrices <- calculate_alpha(Merged_species, IDcol="RowNames", metrics=c("shannon","simpson","invsimpson","richness"), DIVlvls=c("taxS"))
Data_div_total <- merge(Alpha_metrices, Meta_merged, by.x = "RN", by.y = "Fecal_sample_ID_1")
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL

my_comparisons_all=list(c("yes", "no"))
my_cols <- c("no"="#e6be7f","yes"="#9ec2e2")
custom_labels <- c("yes" = "responder", "no" = "non-responder")

alpha_diversity <- Data_div_total %>%
  ggplot(aes(x = ResponseScoreYesNo, y = DIV.S.shannon, fill = ResponseScoreYesNo, color = ResponseScoreYesNo)) +
  geom_jitter(alpha = 1, width = 0.1) +
  geom_violin(trim = FALSE, position = position_dodge(0.9), alpha = 0.5, linewidth = 0.8) +
  geom_boxplot(alpha = 0.3, width = 0.3, size = 0.8) +
  # Set theme
  theme_light() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black", hjust = 0),
    legend.position = "none",
    strip.text = element_text(size = 12)) +
  ylab("Shannon Diversity index") +
  scale_color_manual(values = c("#ce7e00", "#3d85c6")) +
  scale_fill_manual(values = my_cols) +
  scale_x_discrete(labels = custom_labels) +
  stat_compare_means(comparisons = my_comparisons_all, method = "wilcox.test", paired = FALSE, hide.ns = FALSE)
alpha_diversity

#Calculate means
Data_div_total %>%
  dplyr::group_by(ResponseScoreYesNo) %>%
  dplyr::summarize(mean_DIV_S_shannon = mean(DIV.S.shannon, na.rm = TRUE))

Data_div_total_merged <- Data_div_total


# Species: beta diversity PCoA ----
#Calculating beta diversity
vegdist(Merged_CLR_1s, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- merge(PC, Meta_merged, by.x = 'row.names', by.y = "Fecal_sample_ID_1", all = FALSE)

#Making plots
my_cols <- c("no"="#e6be7f","yes"="#9ec2e2")
ref <- reformulate("ResponseScoreYesNo","cbind(V1,V2,V3,V4,V5)")
centroids <- aggregate(ref,PCoA_meta,mean) #calculate centroids

PCoA1_2 <- PCoA_meta %>% ggplot(aes(x=V1, y=V2, color=ResponseScoreYesNo)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_2 <- PCoA1_2 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V2"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V2",col="ResponseScoreYesNo"),alpha=0.8, colour="black") 
print(PCoA1_2)

PCoA1_3 <- PCoA_meta %>% ggplot(aes(x=V1, y=V3, color=ResponseScoreYesNo)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_3 <- PCoA1_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V3",col="ResponseScoreYesNo"),alpha=0.8, colour="black") 
print(PCoA1_3)

PCoA1_4 <- PCoA_meta %>% ggplot(aes(x=V1, y=V4, color=ResponseScoreYesNo)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo4 (", var_expl[4],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_4 <- PCoA1_4 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V4"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V4",col="ResponseScoreYesNo"),alpha=0.8, colour="black") 
print(PCoA1_4)

PCoA2_3 <- PCoA_meta %>% ggplot(aes(x=V2, y=V3, color=ResponseScoreYesNo)) + 
  xlab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA2_3 <- PCoA2_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V2",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V2",y="V3",col="ResponseScoreYesNo"),alpha=0.8, colour="black") 
print(PCoA2_3)

Combined_PCoAs <- ggarrange(PCoA1_2, PCoA1_3, PCoA1_4, PCoA2_3, ncol = 2, nrow = 2)
print(Combined_PCoAs)

# Test significant difference of centroid
wilcoxVarsTouse <- c("V1","V2","V3","V4")
Wilcox_PCoA_Results <- NULL
for (i in wilcoxVarsTouse[1:4]) {
  #create dataframe with only relevant info
  PCoA_meta %>% dplyr::select(c(i, ResponseScoreYesNo)) -> tmDF
  #do wilcoxon
  pairwise.wilcox.test(tmDF[,1], tmDF$ResponseScoreYesNo, p.adjust.method="none") -> wilcox_df
  reshape2::melt(wilcox_df$p.value) -> wilcoxon_coordinates
  #save info in a dataframe
  Wilcox_all <- data.frame(Coordinate=i,
                           Variable1=wilcoxon_coordinates[,1],
                           Variable2=wilcoxon_coordinates[,2],
                           P_value=wilcoxon_coordinates[,3])
  print(Wilcox_all)
  Wilcox_PCoA_Results <- rbind.data.frame(Wilcox_PCoA_Results,Wilcox_all)
}

#betadisper for homogeneity
df_merged <- merge(Merged_CLR_1s, Meta_merged, by.x = 0, by.y = "Fecal_sample_ID_1")
rownames(df_merged) <- df_merged$Row.names
df_merged$Row.names <- NULL
dis <- vegdist(df_merged[,1:65],method = "euclidean") #create distance matrix

#calculate multivariate dispersions
mod <- betadisper(dis, df_merged$ResponseScoreYesNo)
anova(mod)
boxplot(mod, xlab = "")

# Species: PERMANOVA analysis ----
# adonis multivariate 
df_pheno <- Meta_merged %>% select(Fecal_sample_ID_1, Study_ID, ResponseScoreYesNo, V1_CurrentIBDDiagnosis, V1_Sex, V1_BMI, V1_AgeatFecalSampling, V1_AntibioticsWithin3months, V1_PPI, reads_sample_1, V1_ResectionAny, cohort)
df_pheno[sapply(df_pheno, is.character)] <- lapply(df_pheno[sapply(df_pheno, is.character)], as.factor)
df_pheno$V1_PPI <- as.factor(df_pheno$V1_PPI)
df_pheno <- df_pheno[df_pheno$Fecal_sample_ID_1 != "Ust45_1_199", ] #remove this sample <1 million reads
df_pheno <- df_pheno[df_pheno$Fecal_sample_ID_1 != "Ust31_1_133", ]  #remove this sample <1 million reads
df_pheno <- df_pheno[complete.cases(df_pheno$Fecal_sample_ID_1), ]

All_df <- merge(df_pheno, Merged_CLR_1s, by.x = 'Fecal_sample_ID_1', by.y = 'row.names')
ad_taxa <- All_df[c(13:77)]
ad_taxa_dm <- vegdist(ad_taxa,method = "euclidean")
ad_test <- All_df[c(3:12)]
ad1 <- adonis2(formula = ad_taxa_dm ~ ResponseScoreYesNo + V1_Sex + V1_BMI + V1_AgeatFecalSampling + reads_sample_1 + cohort, data = ad_test, permutations = 1000, by = "margin")

model.frame(ad1) -> adonis_mv
adonis_mv$Significant="No"
adonis_mv$Significant[adonis_mv$`Pr(>F)`<0.05]="Yes"
adonis_mv$name <- rownames(adonis_mv)

# Make a plot
my_colours <- c("Yes" = "darkolivegreen3", "No" = "indianred3")
outline_colors <- c("Yes" = "darkgreen", "No" = "darkred")

adonis_plot <- ggplot(adonis_mv, aes(R2, reorder(name, R2), R2, fill = Significant)) +
  geom_col(width = 0.5, aes(color = Significant)) +
  labs(x = "Explained variance (R2)", y = NULL) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(size = 8, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black", hjust = 0))
adonis_plot

# Model 2 with extra covariates
ad2 <- adonis2(formula = ad_taxa_dm ~ ResponseScoreYesNo + V1_CurrentIBDDiagnosis + V1_Sex + V1_BMI + V1_AgeatFecalSampling + V1_AntibioticsWithin3months + V1_PPI + reads_sample_1 + V1_ResectionAny + cohort, data = ad_test, permutations = 1000, by = "margin")

model.frame(ad2) -> adonis_mv
adonis_mv$Significant="No"
adonis_mv$Significant[adonis_mv$`Pr(>F)`<0.05]="Yes"
adonis_mv$name <- rownames(adonis_mv)

# Make a plot
my_colours <- c("Yes" = "darkolivegreen3", "No" = "indianred3")
outline_colors <- c("Yes" = "darkgreen", "No" = "darkred")

adonis_plot <- ggplot(adonis_mv, aes(R2, reorder(name, R2), R2, fill = Significant)) +
  geom_col(width = 0.5, aes(color = Significant)) +
  labs(x = "Explained variance (R2)", y = NULL) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(size = 8, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black", hjust = 0))
adonis_plot

# Species: differential abundance analysis ----
#metadata <- df_pheno
#ID <- "Fecal_sample_ID_1"
#CLR_transformed_data <- Ustek_CLR_1s
#pheno_list <- "ResponseScoreYesNo"

DDA_taxa <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df$Fecal_sample_ID_1[is.na(df[colnames(df) == pheno]) == F] -> To_keep
      df_pheno_dda = filter(df, Fecal_sample_ID_1 %in% To_keep )
      Model2 = as.formula(paste( c(Bug2,  " ~ V1_Sex + V1_BMI +",pheno2, "+ V1_AgeatFecalSampling + V1_CurrentIBDDiagnosis + V1_AntibioticsWithin3months + V1_PPI + reads_sample_1 + V1_ResectionAny + cohort"), collapse="" ))
      lm(Model2, df_pheno_dda) -> resultmodel2
      as.data.frame(summary(resultmodel2)$coefficients)[4,1:4] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(Bug =Bug, Pheno=pheno) -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$`Pr(>|t|)`, method = "BH")
  
  return(p)
}

df_pheno <- df_pheno %>% filter(!is.na(Fecal_sample_ID_1))
df_pheno <- as.data.frame(df_pheno)
diff_ab_bacteria_Merged <- DDA_taxa(df_pheno, "Fecal_sample_ID_1", Merged_CLR_1s, c("ResponseScoreYesNo"))
p_nominal_sig_Merged_6m <- diff_ab_bacteria_Merged %>% filter(`Pr(>|t|)`<0.05)

# Plot the abundances of nominal significant findings (using CLR data)
species_sig <- p_nominal_sig_Merged_6m$Bug
species_to_plot <- c("ResponseScoreYesNo", species_sig)
All_df_sub <- All_df[, species_to_plot]
All_df_long <- tidyr::pivot_longer(All_df_sub, cols = -ResponseScoreYesNo, names_to = "Species", values_to = "Abundance")
All_df_long$Species <- str_extract(All_df_long$Species, "s__[A-Za-z0-9_]+")
All_df_long$Species <- sub("^s_", "", All_df_long$Species)
All_df_long$Species <- gsub("_", " ", All_df_long$Species)

ggplot(All_df_long, aes(x = ResponseScoreYesNo, y = Abundance)) +
  geom_boxplot() +
  facet_wrap(~ Species, scales = "free_y", ncol = 3) +
  xlab("Responder") +
  ylab("CLR transformed Relative Abundance") +
  ggtitle("Boxplots of nominal significant differentially abundant species")

# Plot the abundances of nominal significant findings (using non transformed data)
#Merged_species_1_sub <- Merged_species[, species_sig]
#Merged_species_1_sub$Fecal_sample_ID_1 <- rownames(Merged_species_1_sub)
#df_to_plot <- Merged_species_1_sub %>%
#  left_join(Meta_merged %>% select(Fecal_sample_ID_1, ResponseScoreYesNo), by = "Fecal_sample_ID_1")
#df_to_plot$Fecal_sample_ID_1 <- NULL
#df_to_plot_long <- tidyr::pivot_longer(df_to_plot, cols = -ResponseScoreYesNo, names_to = "Species", values_to = "Abundance")
#df_to_plot_long$Species <- str_extract(df_to_plot_long$Species, "s__[A-Za-z0-9_]+")

#ggplot(df_to_plot_long, aes(x = ResponseScoreYesNo, y = Abundance)) +
#  geom_boxplot() +
#  facet_wrap(~ Species, scales = "free_y", ncol = 2) +
#  xlab("ResponseScoreYesNo") +
#  ylab("Relative Abundance") +
#  ggtitle("Boxplots of Relative Abundance non transformed")


# Pathways: beta diversity PCoA ----
#Calculating beta diversity
vegdist(PWYs_merged_baseline_filt_clr, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- merge(PC, Meta_merged, by.x = 'row.names', by.y = "Fecal_sample_ID_1", all = FALSE)

#Making plots
my_cols <- c("no"="#e6be7f","yes"="#9ec2e2")
ref <- reformulate("ResponseScoreYesNo","cbind(V1,V2,V3,V4,V5)")
centroids <- aggregate(ref,PCoA_meta,mean) #calculate centroids

PCoA1_2 <- PCoA_meta %>% ggplot(aes(x=V1, y=V2, color=ResponseScoreYesNo)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_2 <- PCoA1_2 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V2"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V2",col="ResponseScoreYesNo"),alpha=0.8, colour="black") 
print(PCoA1_2)

PCoA1_3 <- PCoA_meta %>% ggplot(aes(x=V1, y=V3, color=ResponseScoreYesNo)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_3 <- PCoA1_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V3",col="ResponseScoreYesNo"),alpha=0.8, colour="black") 
print(PCoA1_3)

PCoA1_4 <- PCoA_meta %>% ggplot(aes(x=V1, y=V4, color=ResponseScoreYesNo)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo4 (", var_expl[4],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_4 <- PCoA1_4 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V4"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V4",col="ResponseScoreYesNo"),alpha=0.8, colour="black") 
print(PCoA1_4)

PCoA2_3 <- PCoA_meta %>% ggplot(aes(x=V2, y=V3, color=ResponseScoreYesNo)) + 
  xlab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA2_3 <- PCoA2_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V2",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V2",y="V3",col="ResponseScoreYesNo"),alpha=0.8, colour="black") 
print(PCoA2_3)

Combined_PCoAs <- ggarrange(PCoA1_2, PCoA1_3, PCoA1_4, PCoA2_3, ncol = 2, nrow = 2)
print(Combined_PCoAs)

# Test significant difference of centroid
wilcoxVarsTouse <- c("V1","V2","V3","V4")
Wilcox_PCoA_Results <- NULL
for (i in wilcoxVarsTouse[1:4]) {
  #create dataframe with only relevant info
  PCoA_meta %>% dplyr::select(c(i, ResponseScoreYesNo)) -> tmDF
  #do wilcoxon
  pairwise.wilcox.test(tmDF[,1], tmDF$ResponseScoreYesNo, p.adjust.method="none") -> wilcox_df
  reshape2::melt(wilcox_df$p.value) -> wilcoxon_coordinates
  #save info in a dataframe
  Wilcox_all <- data.frame(Coordinate=i,
                           Variable1=wilcoxon_coordinates[,1],
                           Variable2=wilcoxon_coordinates[,2],
                           P_value=wilcoxon_coordinates[,3])
  print(Wilcox_all)
  Wilcox_PCoA_Results <- rbind.data.frame(Wilcox_PCoA_Results,Wilcox_all)
}

# Pathways: PERMANOVA analysis ----
df_pheno <- Meta_merged %>% select(Fecal_sample_ID_1, Study_ID, ResponseScoreYesNo, V1_CurrentIBDDiagnosis, V1_Sex, V1_BMI, V1_AgeatFecalSampling, V1_AntibioticsWithin3months, V1_PPI, reads_sample_1, V1_ResectionAny, cohort)
df_pheno[sapply(df_pheno, is.character)] <- lapply(df_pheno[sapply(df_pheno, is.character)], as.factor)
df_pheno$V1_PPI <- as.factor(df_pheno$V1_PPI)
df_pheno <- df_pheno[df_pheno$Fecal_sample_ID_1 != "Ust45_1_199", ] #remove this sample <1 million reads
df_pheno <- df_pheno[df_pheno$Fecal_sample_ID_1 != "Ust31_1_133", ]  #remove this sample <1 million reads
df_pheno <- df_pheno[complete.cases(df_pheno$Fecal_sample_ID_1), ]

# adonis multivariate 
All_df <- merge(df_pheno, PWYs_merged_baseline_filt_clr, by.x = 'Fecal_sample_ID_1', by.y = 'row.names')
ad_pwy <- All_df[c(13:205)]
ad_pwy_dm <- vegdist(ad_pwy,method = "euclidean")
ad_test <- All_df[c(3:12)]
ad1 <- adonis2(formula = ad_pwy_dm ~ ResponseScoreYesNo + V1_Sex + V1_BMI + V1_AgeatFecalSampling + reads_sample_1 + cohort, data = ad_test, permutations = 1000, by = "margin")

model.frame(ad1) -> adonis_mv
adonis_mv$Significant="No"
adonis_mv$Significant[adonis_mv$`Pr(>F)`<0.05]="Yes"
adonis_mv$name <- rownames(adonis_mv)

# Make a plot
my_colours <- c("Yes" = "darkolivegreen3", "No" = "indianred3")
outline_colors <- c("Yes" = "darkgreen", "No" = "darkred")

adonis_plot <- ggplot(adonis_mv, aes(R2, reorder(name, R2), R2, fill = Significant)) +
  geom_col(width = 0.5, aes(color = Significant)) +
  labs(x = "Explained variance (R2)", y = NULL) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(size = 8, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black", hjust = 0))
adonis_plot

# Model 2 with extra covariates
ad2 <- adonis2(formula = ad_pwy_dm ~ ResponseScoreYesNo + V1_CurrentIBDDiagnosis + V1_Sex + V1_BMI + V1_AgeatFecalSampling + V1_AntibioticsWithin3months + V1_PPI + reads_sample_1 + V1_ResectionAny + cohort, data = ad_test, permutations = 1000, by = "margin")

model.frame(ad2) -> adonis_mv
adonis_mv$Significant="No"
adonis_mv$Significant[adonis_mv$`Pr(>F)`<0.05]="Yes"
adonis_mv$name <- rownames(adonis_mv)

# Make a plot
my_colours <- c("Yes" = "darkolivegreen3", "No" = "indianred3")
outline_colors <- c("Yes" = "darkgreen", "No" = "darkred")

adonis_plot <- ggplot(adonis_mv, aes(R2, reorder(name, R2), R2, fill = Significant)) +
  geom_col(width = 0.5, aes(color = Significant)) +
  labs(x = "Explained variance (R2)", y = NULL) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(size = 8, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black", hjust = 0))
adonis_plot

# Pathways: differential abundance analysis ----
DDA_pwy <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by="row.names")
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df$Fecal_sample_ID_1[is.na(df[colnames(df) == pheno]) == F] -> To_keep
      df_pheno_dda = filter(df, Fecal_sample_ID_1 %in% To_keep )
      Model2 = as.formula(paste( c(Bug2,  " ~ V1_Sex + V1_BMI +",pheno2, "+ V1_AgeatFecalSampling + V1_CurrentIBDDiagnosis + V1_AntibioticsWithin3months + V1_PPI + reads_sample_1 + V1_ResectionAny + cohort"), collapse="" ))
      lm(Model2, df_pheno_dda) -> resultmodel2
      as.data.frame(summary(resultmodel2)$coefficients)[4,1:4] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(Bug =Bug, Pheno=pheno) -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$`Pr(>|t|)`, method = "BH")
  
  return(p)
}

df_pheno <- as.data.frame(df_pheno)
diff_ab_pathways_Merged <- DDA_pwy(df_pheno, "Fecal_sample_ID_1", PWYs_merged_baseline_filt_clr, c("ResponseScoreYesNo"))
p_pwy_nominal_sig_Merged_6m <- diff_ab_pathways_Merged %>% filter(`Pr(>|t|)`<0.05)

# Plot the abundances of nominal significant findings (using CLR data)
pathways_sig <- p_pwy_nominal_sig_Merged_6m$Bug
pathways_to_plot <- c("ResponseScoreYesNo", pathways_sig)
All_df_sub <- All_df[, pathways_to_plot]
All_df_long <- tidyr::pivot_longer(All_df_sub, cols = -ResponseScoreYesNo, names_to = "Pathways", values_to = "Abundance")

ggplot(All_df_long, aes(x = ResponseScoreYesNo, y = Abundance)) +
  geom_boxplot() +
  facet_wrap(~ Pathways, scales = "free_y", ncol = 4) +
  xlab("Responder") +
  ylab("CLR transformed Relative Abundance") +
  ggtitle("Boxplots of nominal significant differentially abundant pathways")

# Metabolites: differential abundance analysis ----
#metadata <- df_pheno
#CLR_transformed_data <- Selected_metabolites_CLR
#ID <- "Study_ID"
#pheno_list <- "ResponseScoreYesNo"

DDA_metabolites <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Metabolite in Prevalent){
    if (! Metabolite %in% colnames(df)){ next }
    Metabolite2 = paste(c("`",Metabolite, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df$Study_ID[is.na(df[colnames(df) == pheno]) == F] -> To_keep
      df_pheno_dda = filter(df, Study_ID %in% To_keep )
      Model2 = as.formula(paste( c(Metabolite2, "~", pheno2, " + V1_Sex + V1_BMI + V1_AgeatFecalSampling + V1_AntibioticsWithin3months + V1_PPI + V1_CurrentIBDDiagnosis + V1_ResectionAny + cohort"), collapse="" ))
      lm(Model2, df_pheno_dda) -> resultmodel2
      as.data.frame(summary(resultmodel2)$coefficients)[2,1:4] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(Metabolite =Metabolite, Pheno=pheno) -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$`Pr(>|t|)`, method = "BH")
  
  return(p)
}

df_pheno_m <- Meta_merged %>% select(Fecal_sample_ID_1, Study_ID, ResponseScoreYesNo, V1_CurrentIBDDiagnosis, V1_Sex, V1_BMI, V1_AgeatFecalSampling, V1_AntibioticsWithin3months, V1_PPI, V1_ResectionAny, reads_sample_1, cohort)
df_pheno_m[sapply(df_pheno_m, is.character)] <- lapply(df_pheno_m[sapply(df_pheno_m, is.character)], as.factor)
df_pheno_m$V1_PPI <- as.factor(df_pheno_m$V1_PPI)
df_pheno_m <- df_pheno_m[df_pheno_m$Study_ID %in% row.names(Selected_metabolites_merged), ]
df_pheno_m <- as.data.frame(df_pheno_m)
df_pheno_m <- df_pheno_m[complete.cases(df_pheno_m$Fecal_sample_ID_1), ]

metabolites_merged_2 <- Selected_metabolites_merged[rownames(Selected_metabolites_merged) %in% df_pheno_m$Study_ID, ]

diff_metabolites_Merged <- DDA_metabolites(df_pheno_m, "Study_ID", metabolites_merged_2, c("ResponseScoreYesNo"))
p_nominal_sig_mMerged_6m <- diff_metabolites_Merged %>% filter(`Pr(>|t|)`<0.05)
p_FDR_005_sig_mMerged_6m <- diff_metabolites_Merged %>% filter(FDR<0.05)
p_FDR_010_sig_mMerged_6m <- diff_metabolites_Merged %>% filter(FDR<0.10)

# What are these metabolites?
metabolite_key <- read_excel("~/Documents/MDPhD/Hfst_IBD_drugresponse/Metabolomics_data_Arnau.XLSX", sheet = "Chemical Annotation")
metabolites_of_interest <- c(p_FDR_010_sig_mMerged_6m$Metabolite)

df_mtb = filter(metabolite_key, CHEM_ID %in% metabolites_of_interest)

# Merge datasets
df <- df_pheno_m
row.names(df) <- df$Study_ID
df_metabolites_merged <- metabolites_merged_2
colnames(df_metabolites_merged) <- paste0("x_", colnames(df_metabolites_merged))
df<-merge(df, df_metabolites_merged, by='row.names')
df$Row.names <- NULL

target_metabolites <- unique(df_mtb$CHEM_ID)
target_metabolites <- paste0("x_", target_metabolites)

# Check emmeans
create_lm <- function(target_metabolites) {
  lm_formula <- as.formula(paste(target_metabolites, "~ ResponseScoreYesNo + V1_AgeatFecalSampling + 
                                 V1_Sex + V1_BMI + V1_CurrentIBDDiagnosis + V1_AntibioticsWithin3months + V1_PPI +
                                 cohort"))
  return(lm(lm_formula, data = df))
}

# Create a list of linear models
lm_list <- lapply(target_metabolites, create_lm)

# Calculate estimated marginal means
emm_list <- lapply(lm_list, emmeans, ~ResponseScoreYesNo)

# Combine the results into a single data frame
result_df <- bind_rows(
  lapply(1:length(target_metabolites), function(i) {
    emm_df <- as.data.frame(emm_list[[i]])
    emm_df$Metabolite <- target_metabolites[i]
    return(emm_df)
  }))

df_mtb$CHEM_ID <- sub("^", "x_", df_mtb$CHEM_ID)
results_df_mix <- merge(result_df, df_mtb, by.x = "Metabolite", by.y = "CHEM_ID")
results_df_mix <- results_df_mix %>% group_by(Metabolite) %>% mutate(increased_in_responder = ifelse(emmean[ResponseScoreYesNo == "yes"] > emmean[ResponseScoreYesNo == "no"], "yes", "no"))
results_df_mix$increased_in_responder <- as.integer(results_df_mix$increased_in_responder == "yes")

my_cols <- c("no"="#e6be7f","yes"="#9ec2e2")
g <- ggplot(results_df_mix, aes(x = emmean, y = reorder(PLOT_NAME, increased_in_responder), color = ResponseScoreYesNo)) +
  geom_point(aes(x = emmean), size = 3) +
  theme_light() +
  geom_errorbar(aes(xmin = emmean - SE, xmax = emmean + SE), width = 0.2) +
  theme(
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.y = element_text(size = 10, face = "italic", colour = "black"),
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.position = "top",
    legend.text = element_text(size = 8, face = "bold", colour = "black"),
    legend.title = element_text(size = 8, face = "bold", colour = "black")) +
  scale_colour_manual(values = my_cols) +
  labs(x = "EMM (log transformed data)", y = "Metabolites") +
  labs(color = "Responder") #+ 
#facet_grid(SUPER_PATHWAY~., scales = "free", space = "free")
print(g)

# Create a volcano plot
for_volcano <- diff_metabolites_Merged %>% select(Estimate, `Pr(>|t|)`, Metabolite, FDR)
colnames(for_volcano)[colnames(for_volcano) == 'Pr(>|t|)'] <- 'pval'

metabolite_key <- read_excel("~/Documents/MDPhD/Hfst_IBD_drugresponse/Metabolomics_data_Arnau.XLSX", sheet = "Chemical Annotation")
key2 = filter(metabolite_key, CHEM_ID %in% for_volcano$Metabolite)

for_volcano2 <- merge(for_volcano, key2, by.x = "Metabolite", by.y = "CHEM_ID" )
for_volcano2$delabel <- ifelse(for_volcano2$PLOT_NAME %in% head(for_volcano2[order(for_volcano2$pval), "PLOT_NAME"], 4), for_volcano2$PLOT_NAME, NA)

ggplot(for_volcano2, aes(Estimate,-log10(pval), fill=SUPER_PATHWAY, label = delabel)) + 
  theme_bw() + geom_hline(yintercept = -log10(0.05), col="red") + xlim (-2.1,2.1) +
  geom_hline(yintercept = -log10(3.0e-04), col="darkgreen") +
  geom_point(shape=21, size=2.5) + 
  scale_fill_manual(name = NULL, values  = c("firebrick3",  "steelblue2",  "gold3", "limegreen","grey77",  "pink1", "salmon", "turquoise", "black", "azure")) +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  labs(x = "Linear Regression Coefficient", y = "-log10(p-value)") +
  ggtitle('Fecal metabolites in responders versus non-responders') +
  geom_text_repel(size=3)


# Enterotype analysis ----
#Import and create pseudocounts (optie 1)
taxa_g <- Merged_genera 

cc_rd <- Meta_merged
cc_rd <- as.data.frame(cc_rd[complete.cases(cc_rd$Fecal_sample_ID_1), ])
rownames(cc_rd) <- cc_rd$Fecal_sample_ID_1
cc_rd <- cc_rd[rownames(cc_rd) != "Ust45_1_199", ] #remove this sample <1 million reads
cc_rd <- cc_rd[rownames(cc_rd) != "Ust31_1_133", ]  #remove this sample <1 million reads

cc_rd2=merge(cc_rd, taxa_g, by=0)
row.names(cc_rd2)=cc_rd2$Row.names
cc_rd3=cc_rd2[,c(2:124)]
predictors2=cc_rd2[,c(125:172)]
taxa_pred=as.data.frame(predictors2*cc_rd3$reads_sample_1)
count <- as.matrix(taxa_pred)

#Perform DMM
fit <- lapply(1:7, dmn, count = count, verbose=TRUE)
lplc <- base::sapply(fit, DirichletMultinomial::laplace) # AIC / BIC / Laplace
bic <- base::sapply(fit, DirichletMultinomial::BIC) # AIC / BIC / Laplace

plot(lplc, type="b", xlab="Number of Dirichlet Components",
     ylab="Model Fit")
plot(bic, type="b", xlab="Number of Dirichlet Components",
     ylab="Model Fit")

best <- fit[[which.min(lplc)]] #best fit based on lplc

ass <- as.data.frame(apply(mixture(best), 1, which.max)) #assign enterotype to sample
colnames(ass) <- c("enterotype")

#Print enterotypes composition
best <- fit[[which.min(lplc)]]

mixturewt(best)
head(mixture(best), 3)

plots <- list()
for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("Genera", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    mutate(Genera = gsub(".*g__", "g__", Genera)) %>%
    arrange(value) %>%
    mutate(Genera = factor(Genera, levels = unique(Genera))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = Genera, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Community type", k))
  print(p)
  plots[[k]] <- p
}

enterotypes <- wrap_plots(plots, ncol = 3)
print(enterotypes)

#Merge PCoA data and enterotype information
vegdist(Merged_CLR_1s, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- merge(PC, Meta_merged, by.x = 'row.names', by.y = "Fecal_sample_ID_1", all = FALSE)

PCoA_e <- merge(PCoA_meta, ass, by.x = 'Row.names', by.y = 0)
PCoA_e$enterotype <- as.factor(PCoA_e$enterotype)

ref <- reformulate("enterotype","cbind(V1,V2,V3,V4,V5)")
centroids <- aggregate(ref,PCoA_e,mean) #calculate centroids

PCoA1_2 <- PCoA_e %>% ggplot(aes(x=V1, y=V2, color=enterotype)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() +
  theme(legend.title = element_blank())
PCoA1_2 <- PCoA1_2 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V2"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V2",col="enterotype"),alpha=0.8, colour="black") 
print(PCoA1_2)

#Test association response and Bact2
df_entero <- as.data.frame(PCoA_e[, c("enterotype", "ResponseScoreYesNo", "Response_2years", "Response_HBI_SCCAI",
                                      "V1_Sex", "V1_AgeatFecalSampling", "V1_BMI", "V1_SerumCRP", "V1_FecalCalprotectine",
                                      "V1_Currentsmoker", "V1_DiseaseDurationYears", "V1_ResectionAny", "V1_PreviousantiTNF")])

df_entero$Bact2 <- ifelse(df_entero$enterotype == 2, "yes", "no")
df_entero[sapply(df_entero, is.character)] <- lapply(df_entero[sapply(df_entero, is.character)], as.factor)
df_entero$V1_PreviousantiTNF <- as.factor(df_entero$V1_PreviousantiTNF)
str(df_entero)

# Function to perform logistic regression and extract relevant information
logistic_regression_summary <- function(response_variable, predictor_variable, data) {
  formula <- as.formula(paste(response_variable, "~", predictor_variable))
  logit_model <- glm(formula, data = data, family = binomial(link = "logit"))
  
  return(data.frame(
    variable = predictor_variable,
    coefficient = coef(logit_model)[2],
    SE = summary(logit_model)$coefficients[, "Std. Error"][2],
    p_Value = summary(logit_model)$coefficients[, "Pr(>|z|)"][2],
    response_score = response_variable,
    n_obs = nobs(logit_model)
  ))
}

# List of variables
test_var <- c("V1_Sex", "V1_AgeatFecalSampling", "V1_BMI", "V1_SerumCRP",
              "V1_FecalCalprotectine", "V1_Currentsmoker", "V1_DiseaseDurationYears",
              "V1_ResectionAny", "V1_PreviousantiTNF", "Bact2")

# List of response variables
response_vars <- c("ResponseScoreYesNo")

# List to store results in dataframe format
results_df <- do.call(rbind, lapply(response_vars, function(response_var) {
  do.call(rbind, lapply(test_var, function(variable) {
    logistic_regression_summary(response_var, variable, df_entero)
  }))
}))

#Results logit Bact2 enterotype yes/no
rownames(results_df) <- NULL
results_df$FDR<-p.adjust(results_df$p_Value, method = "BH")

#Chi square test (taking into account all three enterotypes, not just bact2 yes/no)
df_table <- table(df_entero$enterotype, df_entero$ResponseScoreYesNo)
test1 <- chisq.test(df_table)

#Calculate richness
calculate_alpha <- function(inDF,IDcol="RowNames",
                            metrics=c("shannon","simpson","invsimpson","richness"),
                            DIVlvls=c("taxS")) {
  # select IDs column
  if (IDcol == "RowNames") {
    DIVMatrix <- data.frame(RN=rownames(inDF))
  } else {
    DIVMatrix <- data.frame(IDcol=inDF[[IDcol]])
    colnames(DIVMatrix)[1] <- IDcol
  }
  # iterate over metrics, calculate each
  # NOTE: richness is not implemented in vegan, requires special treatment
  for (l in DIVlvls) {
    if (grepl('^tax.$',l)) {
      toUse <- gsub('^tax','',l) } }
  for (m in metrics) {
    print(paste0('  > calculating ',m,'[',l,']'))
    if (m=="richness") {
      inDFpa <- inDF
      inDFpa[inDFpa > 0] <- 1
      dv <- rowSums(inDFpa)
    } else {
      dv <- diversity(inDF,index = m)
    }
    DIVMatrix <- cbind.data.frame(DIVMatrix,dv)
    colnames(DIVMatrix)[ncol(DIVMatrix)] <- paste0('DIV.',toUse,'.',m)
  }
  return(DIVMatrix)
}
Alpha_metrices <- calculate_alpha(Merged_species, IDcol="RowNames", metrics=c("shannon","simpson","invsimpson","richness"), DIVlvls=c("taxS"))
Data_div_total <- merge(Alpha_metrices, PCoA_e, by.x = "RN", by.y = "Row.names")
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL

# Group by enterotype and calculate the mean DIV.S.richness for each group
mean_DIV_S_richness <- Data_div_total %>%
  group_by(enterotype) %>%
  summarise(mean_DIV_S_richness = mean(DIV.S.richness, na.rm = TRUE))




#### 2 YEARS RESPONSE ---- -----------------------------------------
# Species: alpha diversity ----
calculate_alpha <- function(inDF,IDcol="RowNames",
                            metrics=c("shannon","simpson","invsimpson","richness"),
                            DIVlvls=c("taxS")) {
  # select IDs column
  if (IDcol == "RowNames") {
    DIVMatrix <- data.frame(RN=rownames(inDF))
  } else {
    DIVMatrix <- data.frame(IDcol=inDF[[IDcol]])
    colnames(DIVMatrix)[1] <- IDcol
  }
  # iterate over metrics, calculate each
  # NOTE: richness is not implemented in vegan, requires special treatment
  for (l in DIVlvls) {
    if (grepl('^tax.$',l)) {
      toUse <- gsub('^tax','',l) } }
  for (m in metrics) {
    print(paste0('  > calculating ',m,'[',l,']'))
    if (m=="richness") {
      inDFpa <- inDF
      inDFpa[inDFpa > 0] <- 1
      dv <- rowSums(inDFpa)
    } else {
      dv <- diversity(inDF,index = m)
    }
    DIVMatrix <- cbind.data.frame(DIVMatrix,dv)
    colnames(DIVMatrix)[ncol(DIVMatrix)] <- paste0('DIV.',toUse,'.',m)
  }
  return(DIVMatrix)
}
Alpha_metrices <- calculate_alpha(Merged_species, IDcol="RowNames", metrics=c("shannon","simpson","invsimpson","richness"), DIVlvls=c("taxS"))
Data_div_total <- merge(Alpha_metrices, Meta_merged, by.x = "RN", by.y = "Fecal_sample_ID_1")
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL
Data_div_total <- Data_div_total[!(Data_div_total$Response_2years %in% c("NA (overleden)", "NA (verhuisd)")), ]


my_comparisons_all=list(c("yes", "no"))
my_cols <- c("no"="#e6be7f","yes"="#9ec2e2")
custom_labels <- c("yes" = "responder", "no" = "non-responder")

alpha_diversity <- Data_div_total %>%
  ggplot(aes(x = Response_2years, y = DIV.S.shannon, fill = Response_2years, color = Response_2years)) +
  geom_jitter(alpha = 1, width = 0.1) +
  geom_violin(trim = FALSE, position = position_dodge(0.9), alpha = 0.5, linewidth = 0.8) +
  geom_boxplot(alpha = 0.3, width = 0.3, size = 0.8) +
  # Set theme
  theme_light() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black", hjust = 0),
    legend.position = "none",
    strip.text = element_text(size = 12)) +
  ylab("Shannon Diversity index") +
  scale_color_manual(values = c("#ce7e00", "#3d85c6")) +
  scale_fill_manual(values = my_cols) +
  scale_x_discrete(labels = custom_labels) +
  stat_compare_means(comparisons = my_comparisons_all, method = "wilcox.test", paired = FALSE, hide.ns = FALSE)
alpha_diversity

#Calculate means
Data_div_total %>%
  dplyr::group_by(Response_2years) %>%
  dplyr::summarize(mean_DIV_S_shannon = mean(DIV.S.shannon, na.rm = TRUE))


# Species: beta diversity PCoA----
#Calculating beta diversity
vegdist(Merged_CLR_1s, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- merge(PC, Meta_merged, by.x = 'row.names', by.y = "Fecal_sample_ID_1", all = FALSE)
PCoA_meta <- PCoA_meta[!(PCoA_meta$Response_2years %in% c("NA (overleden)", "NA (verhuisd)")), ]

#Making plots
my_cols <- c("no"="#e6be7f","yes"="#9ec2e2")
ref <- reformulate("Response_2years","cbind(V1,V2,V3,V4,V5)")
centroids <- aggregate(ref,PCoA_meta,mean) #calculate centroids

PCoA1_2 <- PCoA_meta %>% ggplot(aes(x=V1, y=V2, color=Response_2years)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_2 <- PCoA1_2 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V2"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V2",col="Response_2years"),alpha=0.8, colour="black") 
print(PCoA1_2)

PCoA1_3 <- PCoA_meta %>% ggplot(aes(x=V1, y=V3, color=Response_2years)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_3 <- PCoA1_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V3",col="Response_2years"),alpha=0.8, colour="black") 
print(PCoA1_3)

PCoA1_4 <- PCoA_meta %>% ggplot(aes(x=V1, y=V4, color=Response_2years)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo4 (", var_expl[4],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_4 <- PCoA1_4 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V4"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V4",col="Response_2years"),alpha=0.8, colour="black") 
print(PCoA1_4)

PCoA2_3 <- PCoA_meta %>% ggplot(aes(x=V2, y=V3, color=Response_2years)) + 
  xlab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA2_3 <- PCoA2_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V2",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V2",y="V3",col="Response_2years"),alpha=0.8, colour="black") 
print(PCoA2_3)

Combined_PCoAs <- ggarrange(PCoA1_2, PCoA1_3, PCoA1_4, PCoA2_3, ncol = 2, nrow = 2)
print(Combined_PCoAs)

# Test significant difference of centroid
wilcoxVarsTouse <- c("V1","V2","V3","V4")
Wilcox_PCoA_Results <- NULL
for (i in wilcoxVarsTouse[1:4]) {
  #create dataframe with only relevant info
  PCoA_meta %>% dplyr::select(c(i, Response_2years)) -> tmDF
  #do wilcoxon
  pairwise.wilcox.test(tmDF[,1], tmDF$Response_2years, p.adjust.method="none") -> wilcox_df
  reshape2::melt(wilcox_df$p.value) -> wilcoxon_coordinates
  #save info in a dataframe
  Wilcox_all <- data.frame(Coordinate=i,
                           Variable1=wilcoxon_coordinates[,1],
                           Variable2=wilcoxon_coordinates[,2],
                           P_value=wilcoxon_coordinates[,3])
  print(Wilcox_all)
  Wilcox_PCoA_Results <- rbind.data.frame(Wilcox_PCoA_Results,Wilcox_all)
}

#betadisper for homogeneity
df_merged_m <- merge(Merged_CLR_1s, Meta_merged, by.x = 0, by.y = "Fecal_sample_ID_1")
rownames(df_merged_m) <- df_merged_m$Row.names
df_merged_m$Row.names <- NULL
df_merged_m <- df_merged_m[!(df_merged_m$Response_2years %in% c("NA (overleden)", "NA (verhuisd)")), ]
dis <- vegdist(df_merged_m[,1:65],method = "euclidean") #create distance matrix

#calculate multivariate dispersions
mod <- betadisper(dis, df_merged_m$Response_2years)
anova(mod)
boxplot(mod, xlab = "")

# Species: PERMANOVA analysis ----
# adonis multivariate 
df_pheno <- Meta_merged %>% select(Fecal_sample_ID_1, Study_ID, Response_2years, V1_CurrentIBDDiagnosis, V1_Sex, V1_BMI, V1_AgeatFecalSampling, V1_AntibioticsWithin3months, V1_PPI, reads_sample_1, V1_ResectionAny, cohort)
df_pheno[sapply(df_pheno, is.character)] <- lapply(df_pheno[sapply(df_pheno, is.character)], as.factor)
df_pheno$V1_PPI <- as.factor(df_pheno$V1_PPI)
df_pheno <- df_pheno[df_pheno$Fecal_sample_ID_1 != "Ust45_1_199", ] #remove this sample <1 million reads
df_pheno <- df_pheno[df_pheno$Fecal_sample_ID_1 != "Ust31_1_133", ]  #remove this sample <1 million reads
df_pheno <- df_pheno[complete.cases(df_pheno$Fecal_sample_ID_1), ]
df_pheno <- df_pheno[!(df_pheno$Response_2years %in% c("NA (overleden)", "NA (verhuisd)")), ]

All_df <- merge(df_pheno, Merged_CLR_1s, by.x = 'Fecal_sample_ID_1', by.y = 'row.names')
ad_taxa <- All_df[c(13:77)]
ad_taxa_dm <- vegdist(ad_taxa,method = "euclidean")
ad_test <- All_df[c(3:12)]
ad1 <- adonis2(formula = ad_taxa_dm ~ Response_2years + V1_Sex + V1_BMI + V1_AgeatFecalSampling + reads_sample_1 + cohort, data = ad_test, permutations = 1000, by = "margin")

model.frame(ad1) -> adonis_mv
adonis_mv$Significant="No"
adonis_mv$Significant[adonis_mv$`Pr(>F)`<0.05]="Yes"
adonis_mv$name <- rownames(adonis_mv)

# Make a plot
my_colours <- c("Yes" = "darkolivegreen3", "No" = "indianred3")
outline_colors <- c("Yes" = "darkgreen", "No" = "darkred")

adonis_plot <- ggplot(adonis_mv, aes(R2, reorder(name, R2), R2, fill = Significant)) +
  geom_col(width = 0.5, aes(color = Significant)) +
  labs(x = "Explained variance (R2)", y = NULL) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(size = 8, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black", hjust = 0))
adonis_plot

# Model 2 with extra covariates
ad2 <- adonis2(formula = ad_taxa_dm ~ Response_2years + V1_CurrentIBDDiagnosis + V1_Sex + V1_BMI + V1_AgeatFecalSampling + V1_AntibioticsWithin3months + V1_PPI + reads_sample_1 + V1_ResectionAny + cohort, data = ad_test, permutations = 1000, by = "margin")

model.frame(ad2) -> adonis_mv
adonis_mv$Significant="No"
adonis_mv$Significant[adonis_mv$`Pr(>F)`<0.05]="Yes"
adonis_mv$name <- rownames(adonis_mv)

# Make a plot
my_colours <- c("Yes" = "darkolivegreen3", "No" = "indianred3")
outline_colors <- c("Yes" = "darkgreen", "No" = "darkred")

adonis_plot <- ggplot(adonis_mv, aes(R2, reorder(name, R2), R2, fill = Significant)) +
  geom_col(width = 0.5, aes(color = Significant)) +
  labs(x = "Explained variance (R2)", y = NULL) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(size = 8, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black", hjust = 0))
adonis_plot

# Species: differential abundance analysis ----
#metadata <- df_pheno
#ID <- "Fecal_sample_ID_1"
#CLR_transformed_data <- Ustek_CLR_1s
#pheno_list <- "Response_2years"

DDA_taxa <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df$Fecal_sample_ID_1[is.na(df[colnames(df) == pheno]) == F] -> To_keep
      df_pheno_dda = filter(df, Fecal_sample_ID_1 %in% To_keep )
      Model2 = as.formula(paste( c(Bug2,  " ~ V1_Sex + V1_BMI +",pheno2, "+ V1_AgeatFecalSampling + V1_CurrentIBDDiagnosis + V1_AntibioticsWithin3months + V1_PPI + reads_sample_1 + V1_ResectionAny + cohort"), collapse="" ))
      lm(Model2, df_pheno_dda) -> resultmodel2
      as.data.frame(summary(resultmodel2)$coefficients)[4,1:4] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(Bug =Bug, Pheno=pheno) -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$`Pr(>|t|)`, method = "BH")
  
  return(p)
}

df_pheno <- df_pheno %>% filter(!is.na(Fecal_sample_ID_1))
df_pheno <- as.data.frame(df_pheno)
diff_ab_bacteria_Merged <- DDA_taxa(df_pheno, "Fecal_sample_ID_1", Merged_CLR_1s, c("Response_2years"))
p_nominal_sig_Merged_2y <- diff_ab_bacteria_Merged %>% filter(`Pr(>|t|)`<0.05)

# Plot the abundances of nominal significant findings (using CLR data)
species_sig <- p_nominal_sig_Merged_2y$Bug
species_to_plot <- c("Response_2years", species_sig)
#All_df_sub <- All_df[, species_to_plot]
#All_df_long <- tidyr::pivot_longer(All_df_sub, cols = -Response_2years, names_to = "Species", values_to = "Abundance")
#All_df_long$Species <- str_extract(All_df_long$Species, "s__[A-Za-z0-9_]+")
#All_df_long$Species <- sub("^s_", "", All_df_long$Species)
#All_df_long$Species <- gsub("_", " ", All_df_long$Species)

#ggplot(All_df_long, aes(x = Response_2years, y = Abundance)) +
#  geom_boxplot() +
#  facet_wrap(~ Species, scales = "free_y", ncol = 3) +
#  xlab("Responder") +
#  ylab("CLR transformed Relative Abundance") +
#  ggtitle("Boxplots of nominal significant differentially abundant species")

# Plot the abundances of nominal significant findings (using non transformed data)
#Merged_species_1_sub <- Merged_species[, species_sig]
#Merged_species_1_sub$Fecal_sample_ID_1 <- rownames(Merged_species_1_sub)
#df_to_plot <- Merged_species_1_sub %>%
#  left_join(Meta_merged %>% select(Fecal_sample_ID_1, Response_2years), by = "Fecal_sample_ID_1")
#df_to_plot$Fecal_sample_ID_1 <- NULL
#df_to_plot_long <- tidyr::pivot_longer(df_to_plot, cols = -Response_2years, names_to = "Species", values_to = "Abundance")
#df_to_plot_long$Species <- str_extract(df_to_plot_long$Species, "s__[A-Za-z0-9_]+")

#ggplot(df_to_plot_long, aes(x = Response_2years, y = Abundance)) +
#  geom_boxplot() +
#  facet_wrap(~ Species, scales = "free_y", ncol = 2) +
#  xlab("Response_2years") +
#  ylab("Relative Abundance") +
#  ggtitle("Boxplots of Relative Abundance non transformed")


# Pathways: beta diversity PCoA ----
#Calculating beta diversity
vegdist(PWYs_merged_baseline_filt_clr, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- merge(PC, Meta_merged, by.x = 'row.names', by.y = "Fecal_sample_ID_1", all = FALSE)
PCoA_meta <- PCoA_meta[!(PCoA_meta$Response_2years %in% c("NA (overleden)", "NA (verhuisd)")), ]

#Making plots
my_cols <- c("no"="#e6be7f","yes"="#9ec2e2")
ref <- reformulate("Response_2years","cbind(V1,V2,V3,V4,V5)")
centroids <- aggregate(ref,PCoA_meta,mean) #calculate centroids

PCoA1_2 <- PCoA_meta %>% ggplot(aes(x=V1, y=V2, color=Response_2years)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_2 <- PCoA1_2 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V2"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V2",col="Response_2years"),alpha=0.8, colour="black") 
print(PCoA1_2)

PCoA1_3 <- PCoA_meta %>% ggplot(aes(x=V1, y=V3, color=Response_2years)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_3 <- PCoA1_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V3",col="Response_2years"),alpha=0.8, colour="black") 
print(PCoA1_3)

PCoA1_4 <- PCoA_meta %>% ggplot(aes(x=V1, y=V4, color=Response_2years)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo4 (", var_expl[4],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_4 <- PCoA1_4 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V4"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V4",col="Response_2years"),alpha=0.8, colour="black") 
print(PCoA1_4)

PCoA2_3 <- PCoA_meta %>% ggplot(aes(x=V2, y=V3, color=Response_2years)) + 
  xlab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA2_3 <- PCoA2_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V2",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V2",y="V3",col="Response_2years"),alpha=0.8, colour="black") 
print(PCoA2_3)

Combined_PCoAs <- ggarrange(PCoA1_2, PCoA1_3, PCoA1_4, PCoA2_3, ncol = 2, nrow = 2)
print(Combined_PCoAs)

# Test significant difference of centroid
wilcoxVarsTouse <- c("V1","V2","V3","V4")
Wilcox_PCoA_Results <- NULL
for (i in wilcoxVarsTouse[1:4]) {
  #create dataframe with only relevant info
  PCoA_meta %>% dplyr::select(c(i, Response_2years)) -> tmDF
  #do wilcoxon
  pairwise.wilcox.test(tmDF[,1], tmDF$Response_2years, p.adjust.method="none") -> wilcox_df
  reshape2::melt(wilcox_df$p.value) -> wilcoxon_coordinates
  #save info in a dataframe
  Wilcox_all <- data.frame(Coordinate=i,
                           Variable1=wilcoxon_coordinates[,1],
                           Variable2=wilcoxon_coordinates[,2],
                           P_value=wilcoxon_coordinates[,3])
  print(Wilcox_all)
  Wilcox_PCoA_Results <- rbind.data.frame(Wilcox_PCoA_Results,Wilcox_all)
}

# Pathways: PERMANOVA analysis ----
df_pheno <- Meta_merged %>% select(Fecal_sample_ID_1, Study_ID, Response_2years, V1_CurrentIBDDiagnosis, V1_Sex, V1_BMI, V1_AgeatFecalSampling, V1_AntibioticsWithin3months, V1_PPI, reads_sample_1, V1_ResectionAny, cohort)
df_pheno[sapply(df_pheno, is.character)] <- lapply(df_pheno[sapply(df_pheno, is.character)], as.factor)
df_pheno$V1_PPI <- as.factor(df_pheno$V1_PPI)
df_pheno <- df_pheno[df_pheno$Fecal_sample_ID_1 != "Ust45_1_199", ] #remove this sample <1 million reads
df_pheno <- df_pheno[df_pheno$Fecal_sample_ID_1 != "Ust31_1_133", ]  #remove this sample <1 million reads
df_pheno <- df_pheno[complete.cases(df_pheno$Fecal_sample_ID_1), ]
df_pheno <- df_pheno[!(df_pheno$Response_2years %in% c("NA (overleden)", "NA (verhuisd)")), ]

# adonis multivariate 
All_df <- merge(df_pheno, PWYs_merged_baseline_filt_clr, by.x = 'Fecal_sample_ID_1', by.y = 'row.names')
ad_pwy <- All_df[c(13:205)]
ad_pwy_dm <- vegdist(ad_pwy,method = "euclidean")
ad_test <- All_df[c(3:12)]
ad1 <- adonis2(formula = ad_pwy_dm ~ Response_2years + V1_Sex + V1_BMI + V1_AgeatFecalSampling + reads_sample_1 + cohort, data = ad_test, permutations = 1000, by = "margin")

model.frame(ad1) -> adonis_mv
adonis_mv$Significant="No"
adonis_mv$Significant[adonis_mv$`Pr(>F)`<0.05]="Yes"
adonis_mv$name <- rownames(adonis_mv)

# Make a plot
my_colours <- c("Yes" = "darkolivegreen3", "No" = "indianred3")
outline_colors <- c("Yes" = "darkgreen", "No" = "darkred")

adonis_plot <- ggplot(adonis_mv, aes(R2, reorder(name, R2), R2, fill = Significant)) +
  geom_col(width = 0.5, aes(color = Significant)) +
  labs(x = "Explained variance (R2)", y = NULL) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(size = 8, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black", hjust = 0))
adonis_plot

# Model 2 with extra covariates
ad2 <- adonis2(formula = ad_pwy_dm ~ Response_2years + V1_CurrentIBDDiagnosis + V1_Sex + V1_BMI + V1_AgeatFecalSampling + V1_AntibioticsWithin3months + V1_PPI + reads_sample_1 + V1_ResectionAny + cohort, data = ad_test, permutations = 1000, by = "margin")

model.frame(ad2) -> adonis_mv
adonis_mv$Significant="No"
adonis_mv$Significant[adonis_mv$`Pr(>F)`<0.05]="Yes"
adonis_mv$name <- rownames(adonis_mv)

# Make a plot
my_colours <- c("Yes" = "darkolivegreen3", "No" = "indianred3")
outline_colors <- c("Yes" = "darkgreen", "No" = "darkred")

adonis_plot <- ggplot(adonis_mv, aes(R2, reorder(name, R2), R2, fill = Significant)) +
  geom_col(width = 0.5, aes(color = Significant)) +
  labs(x = "Explained variance (R2)", y = NULL) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(size = 8, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black", hjust = 0))
adonis_plot

# Pathways: differential abundance analysis ----
DDA_pwy <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by="row.names")
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df$Fecal_sample_ID_1[is.na(df[colnames(df) == pheno]) == F] -> To_keep
      df_pheno_dda = filter(df, Fecal_sample_ID_1 %in% To_keep )
      Model2 = as.formula(paste( c(Bug2,  " ~ V1_Sex + V1_BMI +",pheno2, "+ V1_AgeatFecalSampling + V1_CurrentIBDDiagnosis + V1_AntibioticsWithin3months + V1_PPI + reads_sample_1 + V1_ResectionAny + cohort"), collapse="" ))
      lm(Model2, df_pheno_dda) -> resultmodel2
      as.data.frame(summary(resultmodel2)$coefficients)[4,1:4] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(Bug =Bug, Pheno=pheno) -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$`Pr(>|t|)`, method = "BH")
  
  return(p)
}

df_pheno <- as.data.frame(df_pheno)
diff_ab_pathways_Merged <- DDA_pwy(df_pheno, "Fecal_sample_ID_1", PWYs_merged_baseline_filt_clr, c("Response_2years"))
p_pwy_nominal_sig_Merged_2y <- diff_ab_pathways_Merged %>% filter(`Pr(>|t|)`<0.05)

# Plot the abundances of nominal significant findings (using CLR data)
pathways_sig <- p_pwy_nominal_sig_Merged_2y$Bug
pathways_to_plot <- c("Response_2years", pathways_sig)
All_df_sub <- All_df[, pathways_to_plot]
All_df_long <- tidyr::pivot_longer(All_df_sub, cols = -Response_2years, names_to = "Pathways", values_to = "Abundance")

ggplot(All_df_long, aes(x = Response_2years, y = Abundance)) +
  geom_boxplot() +
  facet_wrap(~ Pathways, scales = "free_y", ncol = 4) +
  xlab("Responder") +
  ylab("CLR transformed Relative Abundance") +
  ggtitle("Boxplots of nominal significant differentially abundant pathways")

# Metabolites: differential abundance analysis ----
#metadata <- df_pheno
#CLR_transformed_data <- Selected_metabolites_CLR
#ID <- "Study_ID"
#pheno_list <- "Response_2years"

DDA_metabolites <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Metabolite in Prevalent){
    if (! Metabolite %in% colnames(df)){ next }
    Metabolite2 = paste(c("`",Metabolite, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df$Study_ID[is.na(df[colnames(df) == pheno]) == F] -> To_keep
      df_pheno_dda = filter(df, Study_ID %in% To_keep )
      Model2 = as.formula(paste( c(Metabolite2, "~", pheno2, " + V1_Sex + V1_BMI + V1_AgeatFecalSampling + V1_AntibioticsWithin3months + V1_PPI + V1_CurrentIBDDiagnosis + V1_ResectionAny + cohort"), collapse="" ))
      lm(Model2, df_pheno_dda) -> resultmodel2
      as.data.frame(summary(resultmodel2)$coefficients)[2,1:4] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(Metabolite =Metabolite, Pheno=pheno) -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$`Pr(>|t|)`, method = "BH")
  
  return(p)
}

df_pheno_m <- Meta_merged %>% select(Fecal_sample_ID_1, Study_ID, Response_2years, V1_CurrentIBDDiagnosis, V1_Sex, V1_BMI, V1_AgeatFecalSampling, V1_AntibioticsWithin3months, V1_PPI, V1_ResectionAny, reads_sample_1, cohort)
df_pheno_m[sapply(df_pheno_m, is.character)] <- lapply(df_pheno_m[sapply(df_pheno_m, is.character)], as.factor)
df_pheno_m$V1_PPI <- as.factor(df_pheno_m$V1_PPI)
df_pheno_m <- df_pheno_m[df_pheno_m$Study_ID %in% row.names(Selected_metabolites_merged), ]
df_pheno_m <- as.data.frame(df_pheno_m)
df_pheno_m <- df_pheno_m[!(df_pheno_m$Response_2years %in% c("NA (overleden)", "NA (verhuisd)")), ]

diff_metabolites_Merged <- DDA_metabolites(df_pheno_m, "Study_ID", Selected_metabolites_merged, c("Response_2years"))
p_nominal_sig_mMerged_2y <- diff_metabolites_Merged %>% filter(`Pr(>|t|)`<0.05)
p_FDR005_sig_mMerged_2y <- diff_metabolites_Merged %>% filter(FDR<0.05)
p_FDR010_sig_mMerged_2y <- diff_metabolites_Merged %>% filter(FDR<0.10)

# What are these metabolites?
metabolite_key <- read_excel("~/Documents/MDPhD/Hfst_IBD_drugresponse/Metabolomics_data_Arnau.XLSX", sheet = "Chemical Annotation")
metabolites_of_interest <- c(p_nominal_sig_mMerged_2y$Metabolite)

df_mtb = filter(metabolite_key, CHEM_ID %in% metabolites_of_interest)

# Merge datasets
df <- df_pheno_m
row.names(df) <- df$Study_ID
df_metabolites_merged <- Selected_metabolites_merged
colnames(df_metabolites_merged) <- paste0("x_", colnames(df_metabolites_merged))
df<-merge(df, df_metabolites_merged, by='row.names')
df$Row.names <- NULL

target_metabolites <- unique(df_mtb$CHEM_ID)
target_metabolites <- paste0("x_", target_metabolites)

# Check emmeans
create_lm <- function(target_metabolites) {
  lm_formula <- as.formula(paste(target_metabolites, "~ Response_2years + V1_AgeatFecalSampling + 
                                 V1_Sex + V1_BMI + V1_CurrentIBDDiagnosis + V1_AntibioticsWithin3months + V1_PPI +
                                 cohort"))
  return(lm(lm_formula, data = df))
}

# Create a list of linear models
lm_list <- lapply(target_metabolites, create_lm)

# Calculate estimated marginal means
emm_list <- lapply(lm_list, emmeans, ~Response_2years)

# Combine the results into a single data frame
result_df <- bind_rows(
  lapply(1:length(target_metabolites), function(i) {
    emm_df <- as.data.frame(emm_list[[i]])
    emm_df$Metabolite <- target_metabolites[i]
    return(emm_df)
  }))

df_mtb$CHEM_ID <- sub("^", "x_", df_mtb$CHEM_ID)
results_df_mix <- merge(result_df, df_mtb, by.x = "Metabolite", by.y = "CHEM_ID")
results_df_mix <- results_df_mix %>% group_by(Metabolite) %>% mutate(increased_in_responder = ifelse(emmean[Response_2years == "yes"] > emmean[Response_2years == "no"], "yes", "no"))
results_df_mix$increased_in_responder <- as.integer(results_df_mix$increased_in_responder == "yes")

my_cols <- c("no"="#e6be7f","yes"="#9ec2e2")
g <- ggplot(results_df_mix, aes(x = emmean, y = reorder(PLOT_NAME, increased_in_responder), color = Response_2years)) +
  geom_point(aes(x = emmean), size = 3) +
  theme_classic() +
  geom_errorbar(aes(xmin = emmean - SE, xmax = emmean + SE), width = 0.2) +
  theme(
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.y = element_text(size = 10, face = "italic", colour = "black"),
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.position = "top",
    legend.text = element_text(size = 8, face = "bold", colour = "black"),
    legend.title = element_text(size = 8, face = "bold", colour = "black")) +
  scale_colour_manual(values = my_cols) +
  labs(x = "EMM (log transformed data)", y = "Metabolites") +
  labs(color = "Responder") #+ 
#facet_grid(SUPER_PATHWAY~., scales = "free", space = "free")
print(g)

# Create a volcano plot
for_volcano <- diff_metabolites_Merged %>% select(Estimate, `Pr(>|t|)`, Metabolite, FDR)
colnames(for_volcano)[colnames(for_volcano) == 'Pr(>|t|)'] <- 'pval'

metabolite_key <- read_excel("~/Documents/MDPhD/Hfst_IBD_drugresponse/Metabolomics_data_Arnau.XLSX", sheet = "Chemical Annotation")
key2 = filter(metabolite_key, CHEM_ID %in% for_volcano$Metabolite)

for_volcano2 <- merge(for_volcano, key2, by.x = "Metabolite", by.y = "CHEM_ID" )
for_volcano2$delabel <- ifelse(for_volcano2$PLOT_NAME %in% head(for_volcano2[order(for_volcano2$pval), "PLOT_NAME"], 10), for_volcano2$PLOT_NAME, NA)

ggplot(for_volcano2, aes(Estimate,-log10(pval), fill=SUPER_PATHWAY, label = delabel)) + 
  theme_bw() + geom_hline(yintercept = -log10(0.05), col="red") + xlim (-2.1,2.1) +
  geom_point(shape=21, size=2.5) + 
  scale_fill_manual(name = NULL, values  = c("firebrick3",  "steelblue2",  "gold3", "limegreen","grey77",  "pink1", "salmon", "turquoise", "black", "azure")) +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  labs(x = "Linear Regression Coefficient", y = "-log10(p-value)") +
  ggtitle('Fecal metabolites in responders versus non-responders') +
  geom_text_repel(size=3)


#### CLINICAL SCORE RESPONSE ---- ----------------------------------
# Species: alpha diversity ----
calculate_alpha <- function(inDF,IDcol="RowNames",
                            metrics=c("shannon","simpson","invsimpson","richness"),
                            DIVlvls=c("taxS")) {
  # select IDs column
  if (IDcol == "RowNames") {
    DIVMatrix <- data.frame(RN=rownames(inDF))
  } else {
    DIVMatrix <- data.frame(IDcol=inDF[[IDcol]])
    colnames(DIVMatrix)[1] <- IDcol
  }
  # iterate over metrics, calculate each
  # NOTE: richness is not implemented in vegan, requires special treatment
  for (l in DIVlvls) {
    if (grepl('^tax.$',l)) {
      toUse <- gsub('^tax','',l) } }
  for (m in metrics) {
    print(paste0('  > calculating ',m,'[',l,']'))
    if (m=="richness") {
      inDFpa <- inDF
      inDFpa[inDFpa > 0] <- 1
      dv <- rowSums(inDFpa)
    } else {
      dv <- diversity(inDF,index = m)
    }
    DIVMatrix <- cbind.data.frame(DIVMatrix,dv)
    colnames(DIVMatrix)[ncol(DIVMatrix)] <- paste0('DIV.',toUse,'.',m)
  }
  return(DIVMatrix)
}
Alpha_metrices <- calculate_alpha(Merged_species, IDcol="RowNames", metrics=c("shannon","simpson","invsimpson","richness"), DIVlvls=c("taxS"))
Data_div_total <- merge(Alpha_metrices, Meta_merged, by.x = "RN", by.y = "Fecal_sample_ID_1")
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL
Data_div_total <- Data_div_total[complete.cases(Data_div_total$Response_HBI_SCCAI), ]

my_comparisons_all=list(c("yes", "no"))
my_cols <- c("no"="#e6be7f","yes"="#9ec2e2")
custom_labels <- c("yes" = "responder", "no" = "non-responder")

alpha_diversity <- Data_div_total %>%
  ggplot(aes(x = Response_HBI_SCCAI, y = DIV.S.shannon, fill = Response_HBI_SCCAI, color = Response_HBI_SCCAI)) +
  geom_jitter(alpha = 1, width = 0.1) +
  geom_violin(trim = FALSE, position = position_dodge(0.9), alpha = 0.5, linewidth = 0.8) +
  geom_boxplot(alpha = 0.3, width = 0.3, size = 0.8) +
  # Set theme
  theme_light() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black", hjust = 0),
    legend.position = "none",
    strip.text = element_text(size = 12)) +
  ylab("Shannon Diversity index") +
  scale_color_manual(values = c("#ce7e00", "#3d85c6")) +
  scale_fill_manual(values = my_cols) +
  scale_x_discrete(labels = custom_labels) +
  stat_compare_means(comparisons = my_comparisons_all, method = "wilcox.test", paired = FALSE, hide.ns = FALSE)
alpha_diversity

#Calculate means
Data_div_total %>%
  dplyr::group_by(Response_HBI_SCCAI) %>%
  dplyr::summarize(mean_DIV_S_shannon = mean(DIV.S.shannon, na.rm = TRUE))


# Species: beta diversity PCoA----
#Calculating beta diversity
vegdist(Merged_CLR_1s, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- merge(PC, Meta_merged, by.x = 'row.names', by.y = "Fecal_sample_ID_1", all = FALSE)
PCoA_meta <- PCoA_meta[complete.cases(PCoA_meta$Response_HBI_SCCAI), ]

#Making plots
my_cols <- c("no"="#e6be7f","yes"="#9ec2e2")
ref <- reformulate("Response_HBI_SCCAI","cbind(V1,V2,V3,V4,V5)")
centroids <- aggregate(ref,PCoA_meta,mean) #calculate centroids

PCoA1_2 <- PCoA_meta %>% ggplot(aes(x=V1, y=V2, color=Response_HBI_SCCAI)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_2 <- PCoA1_2 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V2"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V2",col="Response_HBI_SCCAI"),alpha=0.8, colour="black") 
print(PCoA1_2)

PCoA1_3 <- PCoA_meta %>% ggplot(aes(x=V1, y=V3, color=Response_HBI_SCCAI)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_3 <- PCoA1_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V3",col="Response_HBI_SCCAI"),alpha=0.8, colour="black") 
print(PCoA1_3)

PCoA1_4 <- PCoA_meta %>% ggplot(aes(x=V1, y=V4, color=Response_HBI_SCCAI)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo4 (", var_expl[4],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_4 <- PCoA1_4 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V4"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V4",col="Response_HBI_SCCAI"),alpha=0.8, colour="black") 
print(PCoA1_4)

PCoA2_3 <- PCoA_meta %>% ggplot(aes(x=V2, y=V3, color=Response_HBI_SCCAI)) + 
  xlab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA2_3 <- PCoA2_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V2",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V2",y="V3",col="Response_HBI_SCCAI"),alpha=0.8, colour="black") 
print(PCoA2_3)

Combined_PCoAs <- ggarrange(PCoA1_2, PCoA1_3, PCoA1_4, PCoA2_3, ncol = 2, nrow = 2)
print(Combined_PCoAs)

# Test significant difference of centroid
wilcoxVarsTouse <- c("V1","V2","V3","V4")
Wilcox_PCoA_Results <- NULL
for (i in wilcoxVarsTouse[1:4]) {
  #create dataframe with only relevant info
  PCoA_meta %>% dplyr::select(c(i, Response_HBI_SCCAI)) -> tmDF
  #do wilcoxon
  pairwise.wilcox.test(tmDF[,1], tmDF$Response_HBI_SCCAI, p.adjust.method="none") -> wilcox_df
  reshape2::melt(wilcox_df$p.value) -> wilcoxon_coordinates
  #save info in a dataframe
  Wilcox_all <- data.frame(Coordinate=i,
                           Variable1=wilcoxon_coordinates[,1],
                           Variable2=wilcoxon_coordinates[,2],
                           P_value=wilcoxon_coordinates[,3])
  print(Wilcox_all)
  Wilcox_PCoA_Results <- rbind.data.frame(Wilcox_PCoA_Results,Wilcox_all)
}

#betadisper for homogeneity
df_merged <- merge(Merged_CLR_1s, Meta_merged, by.x = 0, by.y = "Fecal_sample_ID_1")
rownames(df_merged) <- df_merged$Row.names
df_merged$Row.names <- NULL
df_merged <- df_merged[complete.cases(df_merged$Response_HBI_SCCAI), ]
dis <- vegdist(df_merged[,1:65],method = "euclidean") #create distance matrix

#calculate multivariate dispersions
mod <- betadisper(dis, df_merged$Response_HBI_SCCAI)
anova(mod)
boxplot(mod, xlab = "")

# Species: PERMANOVA analysis ----
# adonis multivariate 
df_pheno <- Meta_merged %>% select(Fecal_sample_ID_1, Study_ID, Response_HBI_SCCAI, V1_CurrentIBDDiagnosis, V1_Sex, V1_BMI, V1_AgeatFecalSampling, V1_AntibioticsWithin3months, V1_PPI, reads_sample_1, V1_ResectionAny, cohort)
df_pheno[sapply(df_pheno, is.character)] <- lapply(df_pheno[sapply(df_pheno, is.character)], as.factor)
df_pheno$V1_PPI <- as.factor(df_pheno$V1_PPI)
df_pheno <- df_pheno[df_pheno$Fecal_sample_ID_1 != "Ust45_1_199", ] #remove this sample <1 million reads
df_pheno <- df_pheno[df_pheno$Fecal_sample_ID_1 != "Ust31_1_133", ]  #remove this sample <1 million reads
df_pheno <- df_pheno[complete.cases(df_pheno$Fecal_sample_ID_1), ]
df_pheno <- df_pheno[complete.cases(df_pheno$Response_HBI_SCCAI), ]

All_df <- merge(df_pheno, Merged_CLR_1s, by.x = 'Fecal_sample_ID_1', by.y = 'row.names')
ad_taxa <- All_df[c(13:77)]
ad_taxa_dm <- vegdist(ad_taxa,method = "euclidean")
ad_test <- All_df[c(3:12)]
ad1 <- adonis2(formula = ad_taxa_dm ~ Response_HBI_SCCAI + V1_Sex + V1_BMI + V1_AgeatFecalSampling + reads_sample_1 + cohort, data = ad_test, permutations = 1000, by = "margin")

model.frame(ad1) -> adonis_mv
adonis_mv$Significant="No"
adonis_mv$Significant[adonis_mv$`Pr(>F)`<0.05]="Yes"
adonis_mv$name <- rownames(adonis_mv)

# Make a plot
my_colours <- c("Yes" = "darkolivegreen3", "No" = "indianred3")
outline_colors <- c("Yes" = "darkgreen", "No" = "darkred")

adonis_plot <- ggplot(adonis_mv, aes(R2, reorder(name, R2), R2, fill = Significant)) +
  geom_col(width = 0.5, aes(color = Significant)) +
  labs(x = "Explained variance (R2)", y = NULL) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(size = 8, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black", hjust = 0))
adonis_plot

# Model 2 with extra covariates
ad2 <- adonis2(formula = ad_taxa_dm ~ Response_HBI_SCCAI + V1_CurrentIBDDiagnosis + V1_Sex + V1_BMI + V1_AgeatFecalSampling + V1_AntibioticsWithin3months + V1_PPI + reads_sample_1 + V1_ResectionAny + cohort, data = ad_test, permutations = 1000, by = "margin")

model.frame(ad2) -> adonis_mv
adonis_mv$Significant="No"
adonis_mv$Significant[adonis_mv$`Pr(>F)`<0.05]="Yes"
adonis_mv$name <- rownames(adonis_mv)

# Make a plot
my_colours <- c("Yes" = "darkolivegreen3", "No" = "indianred3")
outline_colors <- c("Yes" = "darkgreen", "No" = "darkred")

adonis_plot <- ggplot(adonis_mv, aes(R2, reorder(name, R2), R2, fill = Significant)) +
  geom_col(width = 0.5, aes(color = Significant)) +
  labs(x = "Explained variance (R2)", y = NULL) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(size = 8, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black", hjust = 0))
adonis_plot

# Species: differential abundance analysis ----
#metadata <- df_pheno
#ID <- "Fecal_sample_ID_1"
#CLR_transformed_data <- Ustek_CLR_1s
#pheno_list <- "Response_HBI_SCCAI"

DDA_taxa <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df$Fecal_sample_ID_1[is.na(df[colnames(df) == pheno]) == F] -> To_keep
      df_pheno_dda = filter(df, Fecal_sample_ID_1 %in% To_keep )
      Model2 = as.formula(paste( c(Bug2,  " ~ V1_Sex + V1_BMI +",pheno2, "+ V1_AgeatFecalSampling + V1_CurrentIBDDiagnosis + V1_AntibioticsWithin3months + V1_PPI + reads_sample_1 + V1_ResectionAny + cohort"), collapse="" ))
      lm(Model2, df_pheno_dda) -> resultmodel2
      as.data.frame(summary(resultmodel2)$coefficients)[4,1:4] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(Bug =Bug, Pheno=pheno) -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$`Pr(>|t|)`, method = "BH")
  
  return(p)
}

df_pheno <- df_pheno %>% filter(!is.na(Fecal_sample_ID_1))
df_pheno <- as.data.frame(df_pheno)
diff_ab_bacteria_Merged <- DDA_taxa(df_pheno, "Fecal_sample_ID_1", Merged_CLR_1s, c("Response_HBI_SCCAI"))
p_nominal_sig_Merged_HS <- diff_ab_bacteria_Merged %>% filter(`Pr(>|t|)`<0.05)

# Plot the abundances of nominal significant findings (using CLR data)
species_sig <- p_nominal_sig_Merged_HS$Bug
species_to_plot <- c("Response_HBI_SCCAI", species_sig)
All_df_sub <- All_df[, species_to_plot]
All_df_long <- tidyr::pivot_longer(All_df_sub, cols = -Response_HBI_SCCAI, names_to = "Species", values_to = "Abundance")
All_df_long$Species <- str_extract(All_df_long$Species, "s__[A-Za-z0-9_]+")
All_df_long$Species <- sub("^s_", "", All_df_long$Species)
All_df_long$Species <- gsub("_", " ", All_df_long$Species)

ggplot(All_df_long, aes(x = Response_HBI_SCCAI, y = Abundance)) +
  geom_boxplot() +
  facet_wrap(~ Species, scales = "free_y", ncol = 4) +
  xlab("Responder") +
  ylab("CLR transformed Relative Abundance") +
  ggtitle("Boxplots of nominal significant differentially abundant species")

# Plot the abundances of nominal significant findings (using non transformed data)
#Merged_species_1_sub <- Merged_species[, species_sig]
#Merged_species_1_sub$Fecal_sample_ID_1 <- rownames(Merged_species_1_sub)
#df_to_plot <- Merged_species_1_sub %>%
#  left_join(Meta_merged %>% select(Fecal_sample_ID_1, Response_HBI_SCCAI), by = "Fecal_sample_ID_1")
#df_to_plot$Fecal_sample_ID_1 <- NULL
#df_to_plot_long <- tidyr::pivot_longer(df_to_plot, cols = -Response_HBI_SCCAI, names_to = "Species", values_to = "Abundance")
#df_to_plot_long$Species <- str_extract(df_to_plot_long$Species, "s__[A-Za-z0-9_]+")

#ggplot(df_to_plot_long, aes(x = Response_HBI_SCCAI, y = Abundance)) +
#  geom_boxplot() +
#  facet_wrap(~ Species, scales = "free_y", ncol = 2) +
#  xlab("Response_HBI_SCCAI") +
#  ylab("Relative Abundance") +
#  ggtitle("Boxplots of Relative Abundance non transformed")


# Pathways: beta diversity PCoA ----
#Calculating beta diversity
vegdist(PWYs_merged_baseline_filt_clr, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- merge(PC, Meta_merged, by.x = 'row.names', by.y = "Fecal_sample_ID_1", all = FALSE)
PCoA_meta <- PCoA_meta[complete.cases(PCoA_meta$Response_HBI_SCCAI), ]

#Making plots
my_cols <- c("no"="#e6be7f","yes"="#9ec2e2")
ref <- reformulate("Response_HBI_SCCAI","cbind(V1,V2,V3,V4,V5)")
centroids <- aggregate(ref,PCoA_meta,mean) #calculate centroids

PCoA1_2 <- PCoA_meta %>% ggplot(aes(x=V1, y=V2, color=Response_HBI_SCCAI)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_2 <- PCoA1_2 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V2"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V2",col="Response_HBI_SCCAI"),alpha=0.8, colour="black") 
print(PCoA1_2)

PCoA1_3 <- PCoA_meta %>% ggplot(aes(x=V1, y=V3, color=Response_HBI_SCCAI)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_3 <- PCoA1_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V3",col="Response_HBI_SCCAI"),alpha=0.8, colour="black") 
print(PCoA1_3)

PCoA1_4 <- PCoA_meta %>% ggplot(aes(x=V1, y=V4, color=Response_HBI_SCCAI)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo4 (", var_expl[4],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_4 <- PCoA1_4 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V4"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V4",col="Response_HBI_SCCAI"),alpha=0.8, colour="black") 
print(PCoA1_4)

PCoA2_3 <- PCoA_meta %>% ggplot(aes(x=V2, y=V3, color=Response_HBI_SCCAI)) + 
  xlab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.6) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA2_3 <- PCoA2_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V2",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V2",y="V3",col="Response_HBI_SCCAI"),alpha=0.8, colour="black") 
print(PCoA2_3)

Combined_PCoAs <- ggarrange(PCoA1_2, PCoA1_3, PCoA1_4, PCoA2_3, ncol = 2, nrow = 2)
print(Combined_PCoAs)

# Test significant difference of centroid
wilcoxVarsTouse <- c("V1","V2","V3","V4")
Wilcox_PCoA_Results <- NULL
for (i in wilcoxVarsTouse[1:4]) {
  #create dataframe with only relevant info
  PCoA_meta %>% dplyr::select(c(i, Response_HBI_SCCAI)) -> tmDF
  #do wilcoxon
  pairwise.wilcox.test(tmDF[,1], tmDF$Response_HBI_SCCAI, p.adjust.method="none") -> wilcox_df
  reshape2::melt(wilcox_df$p.value) -> wilcoxon_coordinates
  #save info in a dataframe
  Wilcox_all <- data.frame(Coordinate=i,
                           Variable1=wilcoxon_coordinates[,1],
                           Variable2=wilcoxon_coordinates[,2],
                           P_value=wilcoxon_coordinates[,3])
  print(Wilcox_all)
  Wilcox_PCoA_Results <- rbind.data.frame(Wilcox_PCoA_Results,Wilcox_all)
}

# Pathways: PERMANOVA analysis ----
df_pheno <- Meta_merged %>% select(Fecal_sample_ID_1, Study_ID, Response_HBI_SCCAI, V1_CurrentIBDDiagnosis, V1_Sex, V1_BMI, V1_AgeatFecalSampling, V1_AntibioticsWithin3months, V1_PPI, reads_sample_1, V1_ResectionAny, cohort)
df_pheno[sapply(df_pheno, is.character)] <- lapply(df_pheno[sapply(df_pheno, is.character)], as.factor)
df_pheno$V1_PPI <- as.factor(df_pheno$V1_PPI)
df_pheno <- df_pheno[df_pheno$Fecal_sample_ID_1 != "Ust45_1_199", ] #remove this sample <1 million reads
df_pheno <- df_pheno[df_pheno$Fecal_sample_ID_1 != "Ust31_1_133", ]  #remove this sample <1 million reads
df_pheno <- df_pheno[complete.cases(df_pheno$Fecal_sample_ID_1), ]
df_pheno <- df_pheno[complete.cases(df_pheno$Response_HBI_SCCAI), ]

# adonis multivariate 
All_df <- merge(df_pheno, PWYs_merged_baseline_filt_clr, by.x = 'Fecal_sample_ID_1', by.y = 'row.names')
ad_pwy <- All_df[c(13:205)]
ad_pwy_dm <- vegdist(ad_pwy,method = "euclidean")
ad_test <- All_df[c(3:12)]
ad1 <- adonis2(formula = ad_pwy_dm ~ Response_HBI_SCCAI + V1_Sex + V1_BMI + V1_AgeatFecalSampling + reads_sample_1 + cohort, data = ad_test, permutations = 1000, by = "margin")

model.frame(ad1) -> adonis_mv
adonis_mv$Significant="No"
adonis_mv$Significant[adonis_mv$`Pr(>F)`<0.05]="Yes"
adonis_mv$name <- rownames(adonis_mv)

# Make a plot
my_colours <- c("Yes" = "darkolivegreen3", "No" = "indianred3")
outline_colors <- c("Yes" = "darkgreen", "No" = "darkred")

adonis_plot <- ggplot(adonis_mv, aes(R2, reorder(name, R2), R2, fill = Significant)) +
  geom_col(width = 0.5, aes(color = Significant)) +
  labs(x = "Explained variance (R2)", y = NULL) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(size = 8, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black", hjust = 0))
adonis_plot

# Model 2 with extra covariates
ad2 <- adonis2(formula = ad_pwy_dm ~ Response_HBI_SCCAI + V1_CurrentIBDDiagnosis + V1_Sex + V1_BMI + V1_AgeatFecalSampling + V1_AntibioticsWithin3months + V1_PPI + reads_sample_1 + V1_ResectionAny + cohort, data = ad_test, permutations = 1000, by = "margin")

model.frame(ad2) -> adonis_mv
adonis_mv$Significant="No"
adonis_mv$Significant[adonis_mv$`Pr(>F)`<0.05]="Yes"
adonis_mv$name <- rownames(adonis_mv)

# Make a plot
my_colours <- c("Yes" = "darkolivegreen3", "No" = "indianred3")
outline_colors <- c("Yes" = "darkgreen", "No" = "darkred")

adonis_plot <- ggplot(adonis_mv, aes(R2, reorder(name, R2), R2, fill = Significant)) +
  geom_col(width = 0.5, aes(color = Significant)) +
  labs(x = "Explained variance (R2)", y = NULL) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(size = 8, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black", hjust = 0))
adonis_plot

# Pathways: differential abundance analysis ----
DDA_pwy <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by="row.names")
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df$Fecal_sample_ID_1[is.na(df[colnames(df) == pheno]) == F] -> To_keep
      df_pheno_dda = filter(df, Fecal_sample_ID_1 %in% To_keep )
      Model2 = as.formula(paste( c(Bug2,  " ~ V1_Sex + V1_BMI +",pheno2, "+ V1_AgeatFecalSampling + V1_CurrentIBDDiagnosis + V1_AntibioticsWithin3months + V1_PPI + reads_sample_1 + V1_ResectionAny + cohort"), collapse="" ))
      lm(Model2, df_pheno_dda) -> resultmodel2
      as.data.frame(summary(resultmodel2)$coefficients)[4,1:4] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(Bug =Bug, Pheno=pheno) -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$`Pr(>|t|)`, method = "BH")
  
  return(p)
}

df_pheno <- as.data.frame(df_pheno)
diff_ab_pathways_Merged <- DDA_pwy(df_pheno, "Fecal_sample_ID_1", PWYs_merged_baseline_filt_clr, c("Response_HBI_SCCAI"))
p_pwy_nominal_sig_Merged_HS<- diff_ab_pathways_Merged %>% filter(`Pr(>|t|)`<0.05)

# Plot the abundances of nominal significant findings (using CLR data)
pathways_sig <- p_pwy_nominal_sig_Merged_HS$Bug
pathways_to_plot <- c("Response_HBI_SCCAI", pathways_sig)
All_df_sub <- All_df[, pathways_to_plot]
All_df_long <- tidyr::pivot_longer(All_df_sub, cols = -Response_HBI_SCCAI, names_to = "Pathways", values_to = "Abundance")

ggplot(All_df_long, aes(x = Response_HBI_SCCAI, y = Abundance)) +
  geom_boxplot() +
  facet_wrap(~ Pathways, scales = "free_y", ncol = 4) +
  xlab("Responder") +
  ylab("CLR transformed Relative Abundance") +
  ggtitle("Boxplots of nominal significant differentially abundant pathways")

# Metabolites: differential abundance analysis ----
#metadata <- df_pheno
#CLR_transformed_data <- Selected_metabolites_CLR
#ID <- "Study_ID"
#pheno_list <- "Response_HBI_SCCAI"

DDA_metabolites <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Metabolite in Prevalent){
    if (! Metabolite %in% colnames(df)){ next }
    Metabolite2 = paste(c("`",Metabolite, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df$Study_ID[is.na(df[colnames(df) == pheno]) == F] -> To_keep
      df_pheno_dda = filter(df, Study_ID %in% To_keep )
      Model2 = as.formula(paste( c(Metabolite2, "~", pheno2, " + V1_Sex + V1_BMI + V1_AgeatFecalSampling + V1_AntibioticsWithin3months + V1_PPI + V1_CurrentIBDDiagnosis + V1_ResectionAny + cohort"), collapse="" ))
      lm(Model2, df_pheno_dda) -> resultmodel2
      as.data.frame(summary(resultmodel2)$coefficients)[2,1:4] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(Metabolite =Metabolite, Pheno=pheno) -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$`Pr(>|t|)`, method = "BH")
  
  return(p)
}

df_pheno_m <- Meta_merged %>% select(Fecal_sample_ID_1, Study_ID, Response_HBI_SCCAI, V1_CurrentIBDDiagnosis, V1_Sex, V1_BMI, V1_AgeatFecalSampling, V1_AntibioticsWithin3months, V1_PPI, V1_ResectionAny, reads_sample_1, cohort)
df_pheno_m[sapply(df_pheno_m, is.character)] <- lapply(df_pheno_m[sapply(df_pheno_m, is.character)], as.factor)
df_pheno_m$V1_PPI <- as.factor(df_pheno_m$V1_PPI)
df_pheno_m <- df_pheno_m[df_pheno_m$Study_ID %in% row.names(Selected_metabolites_merged), ]
df_pheno_m <- as.data.frame(df_pheno_m)
df_pheno_m <- df_pheno_m[complete.cases(df_pheno_m$Response_HBI_SCCAI), ]

diff_metabolites_Merged <- DDA_metabolites(df_pheno_m, "Study_ID", Selected_metabolites_merged, c("Response_HBI_SCCAI"))
p_nominal_sig_mMerged_HS <- diff_metabolites_Merged %>% filter(`Pr(>|t|)`<0.05)
p_FDR005_sig_mMerged_HS <- diff_metabolites_Merged %>% filter(FDR<0.05)
p_FDR010_sig_mMerged_HS <- diff_metabolites_Merged %>% filter(FDR<0.10)

# What are these metabolites?
metabolite_key <- read_excel("~/Documents/MDPhD/Hfst_IBD_drugresponse/Metabolomics_data_Arnau.XLSX", sheet = "Chemical Annotation")
metabolites_of_interest <- c(p_nominal_sig_mMerged_HS$Metabolite)

df_mtb = filter(metabolite_key, CHEM_ID %in% metabolites_of_interest)

# Merge datasets
df <- df_pheno_m
row.names(df) <- df$Study_ID
df_metabolites_merged <- Selected_metabolites_merged
colnames(df_metabolites_merged) <- paste0("x_", colnames(df_metabolites_merged))
df<-merge(df, df_metabolites_merged, by='row.names')
df$Row.names <- NULL

target_metabolites <- unique(df_mtb$CHEM_ID)
target_metabolites <- paste0("x_", target_metabolites)

# Check emmeans
create_lm <- function(target_metabolites) {
  lm_formula <- as.formula(paste(target_metabolites, "~ Response_HBI_SCCAI + V1_AgeatFecalSampling + 
                                 V1_Sex + V1_BMI + V1_CurrentIBDDiagnosis + V1_AntibioticsWithin3months + V1_PPI +
                                 cohort"))
  return(lm(lm_formula, data = df))
}

# Create a list of linear models
lm_list <- lapply(target_metabolites, create_lm)

# Calculate estimated marginal means
emm_list <- lapply(lm_list, emmeans, ~Response_HBI_SCCAI)

# Combine the results into a single data frame
result_df <- bind_rows(
  lapply(1:length(target_metabolites), function(i) {
    emm_df <- as.data.frame(emm_list[[i]])
    emm_df$Metabolite <- target_metabolites[i]
    return(emm_df)
  }))

df_mtb$CHEM_ID <- sub("^", "x_", df_mtb$CHEM_ID)
results_df_mix <- merge(result_df, df_mtb, by.x = "Metabolite", by.y = "CHEM_ID")
results_df_mix <- results_df_mix %>% group_by(Metabolite) %>% mutate(increased_in_responder = ifelse(emmean[Response_HBI_SCCAI == "yes"] > emmean[Response_HBI_SCCAI == "no"], "yes", "no"))
results_df_mix$increased_in_responder <- as.integer(results_df_mix$increased_in_responder == "yes")

my_cols <- c("no"="#e6be7f","yes"="#9ec2e2")
g <- ggplot(results_df_mix, aes(x = emmean, y = reorder(PLOT_NAME, increased_in_responder), color = Response_HBI_SCCAI)) +
  geom_point(aes(x = emmean), size = 3) +
  theme_classic() +
  geom_errorbar(aes(xmin = emmean - SE, xmax = emmean + SE), width = 0.2) +
  theme(
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.y = element_text(size = 10, face = "italic", colour = "black"),
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.position = "top",
    legend.text = element_text(size = 8, face = "bold", colour = "black"),
    legend.title = element_text(size = 8, face = "bold", colour = "black")) +
  scale_colour_manual(values = my_cols) +
  labs(x = "EMM (log transformed data)", y = "Metabolites") +
  labs(color = "Responder") #+ 
#facet_grid(SUPER_PATHWAY~., scales = "free", space = "free")
print(g)

# Create a volcano plot
for_volcano <- diff_metabolites_Merged %>% select(Estimate, `Pr(>|t|)`, Metabolite, FDR)
colnames(for_volcano)[colnames(for_volcano) == 'Pr(>|t|)'] <- 'pval'

metabolite_key <- read_excel("~/Documents/MDPhD/Hfst_IBD_drugresponse/Metabolomics_data_Arnau.XLSX", sheet = "Chemical Annotation")
key2 = filter(metabolite_key, CHEM_ID %in% for_volcano$Metabolite)

for_volcano2 <- merge(for_volcano, key2, by.x = "Metabolite", by.y = "CHEM_ID" )
for_volcano2$delabel <- ifelse(for_volcano2$PLOT_NAME %in% head(for_volcano2[order(for_volcano2$pval), "PLOT_NAME"], 10), for_volcano2$PLOT_NAME, NA)

ggplot(for_volcano2, aes(Estimate,-log10(pval), fill=SUPER_PATHWAY, label = delabel)) + 
  theme_bw() + geom_hline(yintercept = -log10(0.05), col="red") + xlim (-2.1,2.1) +
  geom_point(shape=21, size=2.5) + 
  scale_fill_manual(name = NULL, values  = c("firebrick3",  "steelblue2",  "gold3", "limegreen","grey77",  "pink1", "salmon", "turquoise", "black", "azure")) +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  labs(x = "Linear Regression Coefficient", y = "-log10(p-value)") +
  ggtitle('Fecal metabolites in responders versus non-responders') +
  geom_text_repel(size=3)


#### Creating list of differential abundant features ----
# 6 months
s_6m <- p_nominal_sig_Merged_6m[, c("Bug", "Pr(>|t|)", "FDR")]
colnames(s_6m)[colnames(s_6m) == "Bug"] <- "Feature"
colnames(s_6m)[colnames(s_6m) == "Pr(>|t|)"] <- "Pval"

m_6m <- p_nominal_sig_mMerged_6m[, c("Metabolite", "Pr(>|t|)", "FDR")]
colnames(m_6m)[colnames(m_6m) == "Metabolite"] <- "Feature"
colnames(m_6m)[colnames(m_6m) == "Pr(>|t|)"] <- "Pval"

p_6m <- p_pwy_nominal_sig_Merged_6m[, c("Bug", "Pr(>|t|)", "FDR")]
colnames(p_6m)[colnames(p_6m) == "Bug"] <- "Feature"
colnames(p_6m)[colnames(p_6m) == "Pr(>|t|)"] <- "Pval"

dda_6m <- full_join(s_6m,m_6m)
dda_6m <- full_join(dda_6m,p_6m)

# 2 years
s_2y <- p_nominal_sig_Merged_2y[, c("Bug", "Pr(>|t|)", "FDR")]
colnames(s_2y)[colnames(s_2y) == "Bug"] <- "Feature"
colnames(s_2y)[colnames(s_2y) == "Pr(>|t|)"] <- "Pval"

m_2y <- p_nominal_sig_mMerged_2y[, c("Metabolite", "Pr(>|t|)", "FDR")]
colnames(m_2y)[colnames(m_2y) == "Metabolite"] <- "Feature"
colnames(m_2y)[colnames(m_2y) == "Pr(>|t|)"] <- "Pval"

p_2y<- p_pwy_nominal_sig_Merged_2y[, c("Bug", "Pr(>|t|)", "FDR")]
colnames(p_2y)[colnames(p_2y) == "Bug"] <- "Feature"
colnames(p_2y)[colnames(p_2y) == "Pr(>|t|)"] <- "Pval"

dda_2y <- full_join(s_2y,m_2y)
dda_2y <- full_join(dda_2y,p_2y)

# 2 years
s_HS <- p_nominal_sig_Merged_HS[, c("Bug", "Pr(>|t|)", "FDR")]
colnames(s_HS)[colnames(s_HS) == "Bug"] <- "Feature"
colnames(s_HS)[colnames(s_HS) == "Pr(>|t|)"] <- "Pval"

m_HS <- p_nominal_sig_mMerged_HS[, c("Metabolite", "Pr(>|t|)", "FDR")]
colnames(m_HS)[colnames(m_HS) == "Metabolite"] <- "Feature"
colnames(m_HS)[colnames(m_HS) == "Pr(>|t|)"] <- "Pval"

p_HS<- p_pwy_nominal_sig_Merged_HS[, c("Bug", "Pr(>|t|)", "FDR")]
colnames(p_HS)[colnames(p_HS) == "Bug"] <- "Feature"
colnames(p_HS)[colnames(p_HS) == "Pr(>|t|)"] <- "Pval"

dda_HS <- full_join(s_HS,m_HS)
dda_HS <- full_join(dda_HS,p_HS)

