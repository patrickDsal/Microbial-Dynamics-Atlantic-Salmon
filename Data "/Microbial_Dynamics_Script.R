rm(list=ls())
# Input files for package
# ape package used for phylogenetic tree reading

library(viridis)
library(phyloseq)
library(ggplot2)
library(microeco)
library(mecodev)
library(ape)
library(tidyverse)
library(magrittr)
library(picante)
library(dplyr)
library(file2meco)
library(microbiomer)
library(car)
library(rstatix)

setwd("C:/Users/Paddy/Documents/PhD/Chapter River(new)/Assembly/Month/animal_Microbiome")

pseq <- readRDS("phyloseq_file.rds")
dataset <- phyloseq2meco(pseq)

# Print the updated sample data
print(sample_data(pseq))
sample_data <- meta_to_df(pseq)
str(sample_data)

if ("sample_id" %in% colnames(sample_data)) {
  sample_data$sample_id <- NULL
}

# Rename the column
colnames(sample_data)[colnames(sample_data) == "Weigth"] <- "Weight"
rownames(sample_data) <- sample_names(pseq)
# Replace the sample data in the phyloseq object with the updated data
sample_data(pseq) <- sample_data(sample_data)
saveRDS(pseq, file = "phyloseq_file.rds")

#Shannon

# Normality test for each SampleType
normality_results <- sample_data %>%
  group_by(SampleType) %>%
  summarise(normality_p = shapiro.test(Shannon)$p.value)

# Levene’s test for homogeneity of variance
levene_test_result <- leveneTest(Shannon ~ SampleType, data = sample_data)

# Prepare results storage
results_list <- list()
results_list$normality_results <- normality_results
results_list$levene_test_result <- data.frame(Test = "Levene's Test", p_value = levene_test_result$`Pr(>F)`[1])

# Choose test based on assumptions
if (all(normality_results$normality_p > 0.05) && levene_test_result$`Pr(>F)`[1] > 0.05) {
  # Perform ANOVA
  anova_result <- aov(Shannon ~ SampleType, data = sample_data)
  anova_summary <- summary(anova_result)[[1]]
  
  # Store ANOVA results
  results_list$anova_result <- data.frame(
    Test = "ANOVA",
    Df = anova_summary$Df[1],
    Sum_Sq = anova_summary$`Sum Sq`[1],
    Mean_Sq = anova_summary$`Mean Sq`[1],
    F_value = anova_summary$`F value`[1],
    p_value = anova_summary$`Pr(>F)`[1]
  )
  
  # Perform Tukey's HSD for pairwise comparisons
  tukey_result <- TukeyHSD(anova_result, "SampleType")$SampleType
  tukey_results_df <- as.data.frame(tukey_result)
  tukey_results_df$Comparison <- rownames(tukey_results_df)
  colnames(tukey_results_df)[4] <- "p_value"
  
  # Adjust p-values (Benjamini-Hochberg correction)
  tukey_results_df <- tukey_results_df %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
    separate(Comparison, into = c("group1", "group2"), sep = "-")  # Extract group names
  
  results_list$pairwise_results <- tukey_results_df
} else {
  # Perform Kruskal-Wallis test
  kruskal_result <- kruskal_test(Shannon ~ SampleType, data = sample_data)
  results_list$kruskal_result <- data.frame(Test = "Kruskal-Wallis", p_value = kruskal_result$p)
  
  # Perform pairwise Wilcoxon test
  pairwise_results <- sample_data %>%
    pairwise_wilcox_test(Shannon ~ SampleType, p.adjust.method = "BH") 
  
  results_list$pairwise_results <- pairwise_results
}

# Combine all results into one dataframe
combined_results <- bind_rows(
  results_list$normality_results,
  results_list$levene_test_result,
  results_list$anova_result,
  results_list$kruskal_result,
  results_list$pairwise_results
)

# Print all pairwise comparisons
print(combined_results, n=50)
results_list$pairwise_results
# Save results
write.csv(combined_results, "pairwise_significance_SampleType_Shannon_Mix.csv", row.names = FALSE)

niche_colors <- c("Water" = "#18115c", "PC" = "#e9c46a", "Macroinvertebrate" = "#c65b25")

get_significance_stars <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")  # For non-significant results
  }
}

# Filter significant results and assign y positions
significant_results <- results_list$pairwise_results %>%
  filter(p.adj < 0.05) %>%  # Keep only significant comparisons
  arrange(p.adj) %>%
  mutate(
    y.position = seq(8.05, 8.4, length.out = n()),  # Adjust y-values manually
    significance = sapply(p.adj, get_significance_stars)  # Assign stars
  )

# Check the results
print(significant_results)
custom_order2 <- c("PC", "Macroinvertebrate", "Water")

# Create the boxplot
boxplot <- ggplot(sample_data, aes(x = factor(SampleType, levels = custom_order2), y = Shannon)) +  
  geom_boxplot(outlier.shape = NA, aes(fill = SampleType)) +  # Fill boxplots with niche colors
  geom_jitter(width = 0.2, size = 3, aes(color = SampleType), alpha = 0.5) +  
  scale_fill_manual(values = niche_colors) +  # Apply niche colors to boxplot fill
  scale_color_manual(values = niche_colors) +  # Keep niche colors for points
  labs(x = NULL, y = "Shannon Index") +  
  theme_minimal(base_size = 14) +  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12),   
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5),  # Border only around the plot
    panel.grid = element_blank()  # Remove gridlines
  )

# Add significant pairwise comparisons with stars
if (nrow(significant_results) > 0) {
  boxplot <- boxplot + 
    ggsignif::geom_signif(
      data = significant_results, 
      aes(xmin = group1, xmax = group2, annotations = significance, y_position = y.position), 
      vjust = 0.6, 
      textsize = 5, 
      parse = FALSE, 
      tip_length = 0.01, 
      manual = TRUE, 
      inherit.aes = FALSE, 
      size = 0.3
    )
}

# Print the plot
print(boxplot)

# Save the plot
ggsave(filename = "Shannon_Niche_comparison_Mix.svg", plot = boxplot, width = 18, height = 14, dpi = 300, units = "cm")


#Chao1

# Normality test for each SampleType
normality_results <- sample_data %>%
  group_by(SampleType) %>%
  summarise(normality_p = shapiro.test(Chao1)$p.value)

# Levene’s test for homogeneity of variance
levene_test_result <- leveneTest(Chao1 ~ SampleType, data = sample_data)

# Prepare results storage
results_list <- list()
results_list$normality_results <- normality_results
results_list$levene_test_result <- data.frame(Test = "Levene's Test", p_value = levene_test_result$`Pr(>F)`[1])

# Choose test based on assumptions
if (all(normality_results$normality_p > 0.05) && levene_test_result$`Pr(>F)`[1] > 0.05) {
  # Perform ANOVA
  anova_result <- aov(Chao1 ~ SampleType, data = sample_data)
  anova_summary <- summary(anova_result)[[1]]
  
  # Store ANOVA results
  results_list$anova_result <- data.frame(
    Test = "ANOVA",
    Df = anova_summary$Df[1],
    Sum_Sq = anova_summary$`Sum Sq`[1],
    Mean_Sq = anova_summary$`Mean Sq`[1],
    F_value = anova_summary$`F value`[1],
    p_value = anova_summary$`Pr(>F)`[1]
  )
  
  # Perform Tukey's HSD for pairwise comparisons
  tukey_result <- TukeyHSD(anova_result, "SampleType")$SampleType
  tukey_results_df <- as.data.frame(tukey_result)
  tukey_results_df$Comparison <- rownames(tukey_results_df)
  colnames(tukey_results_df)[4] <- "p_value"
  
  # Adjust p-values (Benjamini-Hochberg correction)
  tukey_results_df <- tukey_results_df %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
    separate(Comparison, into = c("group1", "group2"), sep = "-")  # Extract group names
  
  results_list$pairwise_results <- tukey_results_df
} else {
  # Perform Kruskal-Wallis test
  kruskal_result <- kruskal_test(Chao1 ~ SampleType, data = sample_data)
  results_list$kruskal_result <- data.frame(Test = "Kruskal-Wallis", p_value = kruskal_result$p)
  
  # Perform pairwise Wilcoxon test
  pairwise_results <- sample_data %>%
    pairwise_wilcox_test(Chao1 ~ SampleType, p.adjust.method = "BH") 
  
  results_list$pairwise_results <- pairwise_results
}

# Combine all results into one dataframe
combined_results <- bind_rows(
  results_list$normality_results,
  results_list$levene_test_result,
  results_list$anova_result,
  results_list$kruskal_result,
  results_list$pairwise_results
)

# Print all pairwise comparisons
print(combined_results, n=50)
results_list$pairwise_results
# Save results
write.csv(combined_results, "pairwise_significance_SampleType_Chao1_Mix.csv", row.names = FALSE)

niche_colors <- c("Water" = "#18115c", "PC" = "#e9c46a", "Macroinvertebrate" = "#c65b25")

get_significance_stars <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")  # For non-significant results
  }
}

# Filter significant results and assign y positions
significant_results <- results_list$pairwise_results %>%
  filter(p.adj < 0.05) %>%  # Keep only significant comparisons
  arrange(p.adj) %>%
  mutate(
    y.position = seq(7000, 7600, length.out = n()),  # Adjust y-values manually
    significance = sapply(p.adj, get_significance_stars)  # Assign stars
  )

# Check the results
print(significant_results)
custom_order2 <- c("PC", "Macroinvertebrate", "Water")

# Create the boxplot
boxplot <- ggplot(sample_data, aes(x = factor(SampleType, levels = custom_order2), y = Chao1)) +  
  geom_boxplot(outlier.shape = NA, aes(fill = SampleType)) +  # Fill boxplots with niche colors
  geom_jitter(width = 0.2, size = 3, aes(color = SampleType), alpha = 0.5) +  
  scale_fill_manual(values = niche_colors) +  # Apply niche colors to boxplot fill
  scale_color_manual(values = niche_colors) +  # Keep niche colors for points
  labs(x = NULL, y = "Chao1 Richness") +  
  theme_minimal(base_size = 14) +  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12),   
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5),  # Border only around the plot
    panel.grid = element_blank()  # Remove gridlines
  )

# Add significant pairwise comparisons with stars
if (nrow(significant_results) > 0) {
  boxplot <- boxplot + 
    ggsignif::geom_signif(
      data = significant_results, 
      aes(xmin = group1, xmax = group2, annotations = significance, y_position = y.position), 
      vjust = 0.6, 
      textsize = 5, 
      parse = FALSE, 
      tip_length = 0.01, 
      manual = TRUE, 
      inherit.aes = FALSE, 
      size = 0.3
    )
}

# Print the plot
print(boxplot)

# Save the plot
ggsave(filename = "Chao1_Niche_comparison_Mix.svg", plot = boxplot, width = 18, height = 14, dpi = 300, units = "cm")


# -------------------------------------Beta
rm(list=ls())

#load phyloseq object
pseq <- readRDS("phyloseq_file.rds")

#Phyloseq to microeco object
dataset <- phyloseq2meco(pseq)

dataset$cal_betadiv(unifrac = T)
dataset$beta_diversity

t1 <- trans_beta$new(dataset = dataset, group = "SampleType", measure = "wei_unifrac")
# t1$res_ordination is the ordination result list
t1$cal_ordination(ordination = "PCoA")
class(t1$res_ordination)

niche_colors <- c("Water" = "#18115c", "PC" = "#e9c46a", "Macroinvertebrate" = "#c65b25")

b1 <- t1$plot_ordination(plot_color = "SampleType", plot_shape=NULL , plot_type = c("point", "centroid"),color_values = niche_colors,
                         point_size = 4, point_alpha = 1, centroid_segment_alpha = 0.3, centroid_segment_size = 0.3, centroid_segment_linetype = 1,
                         shape_values =c(15,16,17,18,19,20,21),
                         add_sample_label =  )

b1 <- b1 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
  theme(text=element_text(family="Calibri", face="bold", size=11))+
  theme(legend.key.size = unit(0.5, 'cm'))+
  theme(legend.title = element_text(size=9))

b1

ggsave(filename = "Beta_Mix_wei_unifrac.svg", plot = b1, width = 18, height = 18, dpi = 300, units = "cm")

#Ancom
library(BiocManager)
BiocManager::install("microbiome")
BiocManager::install("ANCOMBC")

library(ANCOMBC)
library(microbiome)
library(DT)
library(file2meco)
# from microtable to phyloseq object

pseq <- meco2phyloseq(dataset)
pseq
#alternative load phyloseq object
pseq <- readRDS("phyloseq_file.rds")

#AncomBC2
#create microeco dataset and transform to physeq...

#add otu column
sample_data(pseq)$SampleType <- factor(sample_data(pseq)$SampleType, 
                                       levels = c("PC", "Macroinvertebrate", "Water"))

# Check the levels to ensure the referenc
levels(sample_data(pseq)$SampleType)
unique(sample_data(pseq)$SampleType)
# Check the structure of the SampleType variable
str(sample_data(pseq)$SampleType)

# Print each level individually to see if there are hidden characters
for (level in levels(sample_data(pseq)$SampleType)) {
  cat(paste0("'", level, "'"), "\n")
}

#Run ANCOM
output = ancombc2(
  data = pseq, 
  assay_name = "Environmental", 
  tax_level = "Genus",
  fix_formula = "SampleType", 
  rand_formula = NULL,
  p_adj_method = "holm", 
  pseudo_sens = TRUE,
  prv_cut = 0.05, 
  lib_cut = 1500, 
  s0_perc = 0.05,
  group = "SampleType", 
  struc_zero = FALSE, 
  neg_lb = TRUE,
  alpha = 0.05, 
  n_cl = 6, 
  verbose = TRUE,
  global = TRUE, 
  pairwise = TRUE,  # Enable pairwise comparisons
  dunnet = FALSE,   # Disable Dunnett test (not needed for pairwise)
  trend = FALSE,
  iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
  em_control = list(tol = 1e-5, max_iter = 100),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100)
)

dunn <- output$res_dunn
write.csv(dunn,"AncomBC_dunn_OTU_0.05.csv", row.names = TRUE)

#structural zeros
tab_zero = output$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")

#bias correcte dlog abundances
bias_correct_log_table = output$bias_correct_log_table

#Multiple pairwise comparisons

ancombc_table = output$res
str(ancombc_table)

# Create a dataframe with only significant taxa
ancombc_table <- ancombc_table %>%
  mutate(significant_water = p_SampleTypeWater < 0.05,
         significant_macroinvertebrate = p_SampleTypeMacroinvertebrate < 0.05)

str(ancombc_table)
# Create a dataframe with only significant taxa
significant_taxa <- ancombc_table %>%
  filter(p_SampleTypeWater < 0.05 | p_SampleTypeMacroinvertebrate < 0.05) %>%
  select(taxon, lfc_SampleTypeWater, lfc_SampleTypeMacroinvertebrate, p_SampleTypeWater, p_SampleTypeMacroinvertebrate) %>%
  gather(key = "group", value = "lfc", -taxon, -p_SampleTypeWater, -p_SampleTypeMacroinvertebrate) %>%
  mutate(
    group = case_when(
      group == "lfc_SampleTypeWater" ~ "Water",
      group == "lfc_SampleTypeMacroinvertebrate" ~ "Macroinvertebrate"
    ),
    # Determine if the taxa is significant for each group
    significant = case_when(
      group == "Water" & p_SampleTypeWater < 0.05 ~ "Significant",
      group == "Macroinvertebrate" & p_SampleTypeMacroinvertebrate < 0.05 ~ "Significant",
      TRUE ~ "Not Significant"
    )
  )

# Create Bar Plot
bar_plot <- ggplot(significant_taxa, aes(x = reorder(taxon, lfc), y = lfc, fill = significant)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +  # Flip the axes so taxa names are readable
  facet_wrap(~group) +  # Separate plots for Water and Macroinvertebrate
  labs(x = "Taxon", y = "Log Fold Change", title = "Significant Taxa Enriched in Groups") +
  scale_fill_manual(values = c("grey", "red")) +
  theme_minimal(base_size = 12) +
  theme(axis.text = element_text(size = 8), strip.text = element_text(size = 10))

# Print the Bar Plot
print(bar_plot)



ancombc_table = output$res
str(ancombc_table)

df_fig_pair1 <- ancombc_table %>%
  filter(diff_SampleTypeMacroinvertebrate == TRUE | diff_SampleTypeWater == TRUE) %>%
  mutate(
    # Pairwise log fold changes
    lfc_pc_vs_water = ifelse(diff_SampleTypeWater == TRUE, round(lfc_SampleTypeWater, 2), 0),
    lfc_pc_vs_macro = ifelse(diff_SampleTypeMacroinvertebrate == TRUE, round(lfc_SampleTypeMacroinvertebrate, 2), 0),
    
    # Compute Water vs Macro difference manually (if needed)
    lfc_macro_vs_water = ifelse(diff_SampleTypeWater == TRUE & diff_SampleTypeMacroinvertebrate == TRUE,
                                round(lfc_SampleTypeWater - lfc_SampleTypeMacroinvertebrate, 2), 0)
  ) %>%
  pivot_longer(cols = lfc_pc_vs_water:lfc_macro_vs_water, names_to = "group", values_to = "value") %>%
  arrange(taxon)

# Rename groups for clarity
df_fig_pair1 <- df_fig_pair1 %>%
  mutate(group = case_when(
    group == "lfc_pc_vs_water" ~ "PC - Water",
    group == "lfc_pc_vs_macro" ~ "PC - Macroinvertebrate",
    group == "lfc_macro_vs_water" ~ "Macroinvertebrate - Water",
    TRUE ~ as.character(group)
  ))

# Create a second dataframe for color coding based on significance
df_fig_pair2 <- ancombc_table %>%
  filter(diff_SampleTypeMacroinvertebrate == TRUE | diff_SampleTypeWater == TRUE) %>%
  mutate(
    # Set colors based on passed significance
    lfc1 = ifelse(passed_ss_SampleTypeMacroinvertebrate == TRUE & diff_SampleTypeMacroinvertebrate == TRUE, "aquamarine3", "black"),
    lfc2 = ifelse(passed_ss_SampleTypeWater == TRUE & diff_SampleTypeWater == TRUE, "aquamarine3", "black"),
    lfc3 = ifelse(passed_ss_SampleTypeWater == TRUE & diff_SampleTypeMacroinvertebrate == TRUE, "aquamarine3", "black")
  ) %>%
  pivot_longer(cols = lfc1:lfc3, names_to = "group", values_to = "color") %>%
  arrange(taxon)

# Join both dataframes for the final dataframe
df_fig_pair <- df_fig_pair1 %>%
  left_join(df_fig_pair2, by = c("taxon", "group"))

# Rename groups again for clarity
df_fig_pair <- df_fig_pair %>%
  mutate(group = case_when(
    group == "lfc_pc_vs_water" ~ "PC - Water",
    group == "lfc_pc_vs_macro" ~ "PC - Macroinvertebrate",
    group == "lfc_macro_vs_water" ~ "Macroinvertebrate - Water",
    TRUE ~ as.character(group)
  ))

# Factorize the 'group' variable for ordering in the plot
df_fig_pair$group = factor(df_fig_pair$group, 
                           levels = c("PC - Macroinvertebrate",
                                      "PC - Water", 
                                      "Macroinvertebrate - Water"))

# Remove "g__" from the taxon column for cleaner presentation
df_fig_pair <- df_fig_pair %>%
  mutate(taxon = sub("^g__", "", taxon))

str(df_fig_pair)
# Create a ranking metric based on log fold change and p-values
ranking_metric <- df_fig_pair %>%
  mutate(
    rank_metric = -log10(p_SampleTypeMacroinvertebrate.x + 1e-10) * abs(lfc_SampleTypeMacroinvertebrate.x) +
      -log10(p_SampleTypeWater.x + 1e-10) * abs(lfc_SampleTypeWater.x)
  )

# View the updated dataframe to confirm
head(ranking_metric)

# Select top 60 taxa based on ranking metric (if desired)
top_taxa <- ranking_metric %>%
  arrange(desc(rank_metric)) %>%
  slice(1:60)

# Alternatively, manually select taxa (optional)
df_fig_pair <- df_fig_pair %>%
  filter(taxon %in% c("Flavobacterium", "Methylotenera", "Sphingorhabdus", "Novosphingobium", 
                      "Elsteraceae", "Candidatus_Omnitrophus", "Omnitrophales", "Aquabacterium", 
                      "Schlegelella", "Corynebacterium", "Staphylococcus", "Acidovorax", "Turicella", 
                      "Streptococcus", "Clostridium_sensu_stricto_1", "Vibrio", "Psychrobacter", 
                      "Aliivibrio", "Delftia", "Microbacterium", "Bacillus", "Photobacterium", 
                      "Mycoplasma", "Lactococcus", "Streptococcus", "Pseudoalteromonas"))

group_colors <- c(
  "Water" = "#18115c", 
  "PC" = "#e9c46a", 
  "Macroinvertebrate" = "#c65b25"
)


df_fig_pair <- df_fig_pair %>%
  mutate(genus = gsub(".*g__([A-Za-z]+).*", "\\1", taxon))
str(df_fig_pair)
# Check if we correctly extracted the genus (see a summary of the 'genus' column)
table(df_fig_pair$genus)

df_pc_macro <- df_fig_pair %>%
  filter(group == "PC - Macroinvertebrate") %>%
  mutate(color = ifelse(value < 0, "PC", "Macroinvertebrate"))

# Plot for PC vs. Macroinvertebrate
fig_pc_macro <- df_pc_macro %>%
  ggplot(aes(x = genus, y = value, fill = color)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") + 
  scale_fill_manual(values = c("PC" = "#e9c46a", "Macroinvertebrate" = "#c65b25")) + 
  theme_minimal() +
  coord_flip() +  # Flip the plot so genus names appear on the y-axis
  theme(
    plot.title = element_text(hjust = 0.5, size=12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "italic"),
    panel.border = element_rect(color = "black", fill=NA, size = 1),  # Add border around the plot
    axis.line = element_line(color = "black", size = 1),  # Add axis lines
    legend.position = "none",  # Remove the legend
    panel.grid = element_blank() 
  ) +
  labs(
    title = "PC vs. Macroinvertebrate",
    x = "Genus",
    y = "Log Fold Change"
  )

fig_pc_macro


# Filter data for PC vs. Water comparison
df_pc_water <- df_fig_pair %>%
  filter(group == "PC - Water") %>%
  mutate(color = ifelse(value < 0, "PC", "Water"))  # Corrected color mapping

# Plot for PC vs. Water
fig_pc_water <- df_pc_water %>%
  ggplot(aes(x = genus, y = value, fill = color)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") + 
  scale_fill_manual(values = c("PC" = "#e9c46a", "Water" = "#18115c")) +  # Reverse the color mapping
  theme_minimal() +
  coord_flip() +  # Flip the plot so genus names appear on the y-axis
  theme(
    plot.title = element_text(hjust = 0.5, size=12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "italic"),
    panel.border = element_rect(color = "black", fill=NA, size = 1),  # Add border around the plot
    axis.line = element_line(color = "black", size = 1),  # Add axis lines
    legend.position = "none",
    panel.grid = element_blank()   # Remove the legend
  ) +
  labs(
    title = "PC vs. Water",
    x = "Genus",
    y = "Log Fold Change"
  )

fig_pc_water

# Filter data for Macroinvertebrate vs. Water comparison
df_macro_water <- df_fig_pair %>%
  filter(group == "Macroinvertebrate - Water") %>%
  mutate(color = ifelse(value < 0, "Macroinvertebrate", "Water"))

# Plot for Macroinvertebrate vs. Water
fig_macro_water <- df_macro_water %>%
  ggplot(aes(x = genus, y = value, fill = color)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") + 
  scale_fill_manual(values = c("Macroinvertebrate" = "#c65b25", "Water" = "#18115c")) + 
  theme_minimal() +
  coord_flip() +  # Flip the plot so genus names appear on the y-axis
  theme(
    plot.title = element_text(hjust = 0.5, size=12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "italic"),
    panel.border = element_rect(color = "black", fill=NA, size = 1),  # Add border around the plot
    axis.line = element_line(color = "black", size = 1),  # Add axis lines
    legend.position = "none",
    panel.grid = element_blank()   # Remove the legend
  ) +
  labs(
    title = "Macroinvertebrate vs. Water",
    x = "Genus",
    y = "Log Fold Change"
  )

fig_macro_water



library(patchwork)

fig_pc_water <- fig_pc_water + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
fig_macro_water <- fig_macro_water + theme(axis.title.y = element_blank(), axis.text.y = element_blank())

# Now combine the plots horizontally and keep the shared y-axis (genus)
combined_plot <- (fig_pc_macro + fig_pc_water + fig_macro_water) + 
  plot_layout(ncol = 3) 


# Display the combined plot
combined_plot

ggsave(filename = "Ancombc2_Barplot_Mix_New_nogrid.svg", plot = combined_plot, width = 26, height = 18, dpi = 300, units = "cm")


#-------------------------------------------Time Series Plots----------------------------------------------------------
rm(list=ls())


library(ggplot2)
library(dplyr) 
library(viridis) 
library(Interpol.T) 
library(lubridate) 
library(ggExtra) 
library(tidyr) 
library(tidyverse)
library(timetk)


setwd("YourPath")


#-----------------------------------------------------------------------------------------


data<- read.csv("WaterLevel19.csv",header=T)

data$Day <- as.POSIXct((paste(data$Date, data$Time)), format="%Y-%m-%d %H:%M:%S")

data <- na.omit(data)

# Setup for the plotly charts (# FALSE returns ggplots)
interactive <- FALSE

p1 <-data %>% 
  plot_time_series(Day,WaterLevel,
                   .color_var = year(Day),      # Color applied to Week transformation
                   # Facet formatting
                   .interactive = interactive,
                   .title = NULL,
                   .legend_show = FALSE,
                   .x_lab = "Date (30-min intervals)",
                   .y_lab = "Water level (m) ",)
b1 <- p1 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
  theme(legend.position="none")+
  theme(text=element_text(family="Calibri", face="bold", size=11))

b1

ggsave(filename = "WaterLevel.svg", plot = b1, width = 14, height = 10, dpi = 300, units = "cm")

p1<-data %>% 
  plot_time_series(Day,Flow,
                   .color_var = year(Day),      # Color applied to Week transformation
                   # Facet formatting
                   .facet_scales = "free", 
                   .interactive = interactive,
                   .title = NULL,
                   .legend_show = FALSE,
                   .x_lab = "Date (30-min intervals)",
                   .y_lab = "Flow rate (l/s) ",)
b1 <- p1 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
  theme(legend.position="none")+
  theme(text=element_text(family="Calibri", face="bold", size=11))

b1

ggsave(filename = "Flow.png", plot = b1, width = 14, height = 10, dpi = 300, units = "cm")

#------------------------ conductivity/DO------------------------

data<- read.csv("Conductivity19.csv",header=T)

data$Day <- as.POSIXct((paste(data$Date, data$Time)), format="%Y-%m-%d %H:%M:%S")

data <- na.omit(data)

interactive <- FALSE

b1<-data %>% 
  plot_time_series(Day,Conductivity,
                   .color_var = year(Day),      # Color applied to Week transformation
                   # Facet formatting
                   .facet_ncol = 2, 
                   .facet_scales = "free", 
                   .interactive = interactive,
                   .title = NULL,
                   .legend_show = FALSE,
                   .x_lab = "Date (2-min intervals)",
                   .y_lab = "Conductivity (mS/cm) ",)

b1 <- b1 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
  theme(legend.position="none")+
  theme(text=element_text(family="Calibri", face="bold", size=11))

b1

ggsave(filename = "Conductivity.png", plot = b1, width = 14, height = 10, dpi = 300, units = "cm")

p1 <- data %>% 
  plot_time_series(Day,(DO_mgl),
                   .color_var = year(Day),      # Color applied to Week transformation
                   # Facet formatting
                   .facet_scales = "free", 
                   .interactive = interactive,
                   .title = NULL,
                   .legend_show = FALSE,
                   .x_lab = "Date (2-min intervals)",
                   .y_lab = "Dissolved oxygen (mg/l) ",)
b1 <- p1 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
  theme(legend.position="none")+
  theme(text=element_text(family="Calibri", face="bold", size=11))

b1

ggsave(filename = "DO.png", plot = b1, width = 14, height = 10, dpi = 300, units = "cm")


#--------------------------------------------------------------------------------Alpha Diversity Modelling-------------------------------

rm(list=ls())
setwd("C:/Users/Paddy/Documents/PhD/Chapter River(new)/Assembly/Month/Animal_Microbiome")

library(viridis)
library(phyloseq)
library(ggplot2)
library(microeco)
library(mecodev)
library(ape)
library(tidyverse)
library(magrittr)
library(picante)
library(dplyr)
library(file2meco)
library(microbiomer)
library(car)
library(rstatix)
library(lme4)
library(MASS)
library(ggfortify)
library(boot)
library(car)
library(pscl)
library(ggpmisc)
library(MuMIn)
library(performance)
library(corrplot)
library(ggcorrplot)
library(gridExtra)
library(see)


#load phyloseq object
pseq <- readRDS("phyloseq_file.rds")
#Phyloseq to microeco object
dataset <- phyloseq2meco(pseq)

meta <- meta_to_df(pseq)
#subset the data; here we just choose gut samples
data <- subset(meta, Environment == "Gut")

# Examine summary statistics
summary(data)

#----Chao1
summary(data$Chao1)
data$Chao1 <- round(data$Chao1)

#distribution examples:
hist(rnorm(1000), main = "Histogram of Normal Distribution", xlab = "Value")
hist(rpois(1000, lambda = 2), main = "Histogram of Poisson Distribution", xlab = "Value")
hist(rnbinom(1000, size = 2, mu = 5), main = "Histogram of Negative Binomial Distribution", xlab = "Value")
hist(rexp(1000, rate = 1), main = "Histogram of Exponential Distribution", xlab = "Value")
hist(rgamma(1000, shape = 2, scale = 2), main = "Histogram of Gamma Distribution", xlab = "Value")
hist(rbinom(1000, size = 10, prob = 0.5), main = "Histogram of Binomial Distribution", xlab = "Value")

# Check the distribution of Chao1 alpha diversity
hist(data$Chao1, breaks = 50, main = "Histogram of Chao1 Alpha Diversity", xlab = "Chao1")

#------gamma distribution
fit_gamma <- fitdistr(data$Chao1, "gamma")
shape <- fit_gamma$estimate["shape"]
rate <- fit_gamma$estimate["rate"]

# Overlay the fitted Gamma distribution on the histogram
hist(data$Chao1, breaks = 30, probability = TRUE, main = "Histogram of Chao1 with Fitted Gamma", xlab = "Chao1")
curve(dgamma(x, shape = shape, rate = rate), col = "blue", lwd = 2, add = TRUE)

# Q-Q plot to compare the quantiles of the observed data with the quantiles of the Gamma distribution

qqPlot(data$Chao1, dist = "gamma", shape = shape, rate = rate, main = "Q-Q Plot for Gamma Distribution")

#------fit gamma model
model_gamma_glm <- glm(Chao1 ~ Month + Temp + Sex + Origin + Length + Weight , family = Gamma(link = "log"), data = data)
summary(model_gamma_glm)
check_model(model_gamma_glm)

par(mfrow = c(2, 2))
plot(model_gamma_glm)

#-------Fit normal distribution
fit_normal <- fitdistr(data$Chao1, "normal")

print(fit_normal)

#visualize fit
# Create histogram
hist(data$Chao1, breaks = 20, probability = TRUE, 
     main = "Histogram of Chao1 with Fitted Normal Distribution",
     xlab = "Chao1", col = "lightblue")

# Overlay the fitted normal distribution
curve(dnorm(x, mean = fit_normal$estimate["mean"], sd = fit_normal$estimate["sd"]),
      col = "darkblue", lwd = 2, add = TRUE)

# Generate Q-Q plot
qqnorm(data$Chao1, main = "Q-Q Plot of Chao1")
qqline(data$Chao1, col = "red")

#---- Fit a linear regression model
model_lm <- lm(Chao1 ~ Temp + Sex + Origin + Length + Weight + Month, data = data)

#check Variance Inflation Factor (VIF) for colinearity
vif(model_lm)

# Check normality of residuals
residuals <- resid(model_lm)
hist(residuals, breaks = 20, main = "Histogram of Residuals", xlab = "Residuals")

# Shapiro-Wilk test for normality
shapiro.test(residuals)

#Performance
check_model(model_lm)

#-------------Example for Poisson
hist(data$Chao1, probability = TRUE, main = "Histogram of Chao1 with Poisson Fit", xlab = "Chao1")
curve(dpois(x, lambda = mean(data$Chao1)), col = "red", add = TRUE)

#-------Example for Negative Binomial
hist(data$Chao1, probability = TRUE, main = "Histogram of Chao1 with Negative Binomial Fit", xlab = "Chao1")
fit_nb <- fitdistr(data$Chao1, "Negative Binomial")
curve(dnbinom(x, size = fit_nb$estimate["size"], mu = fit_nb$estimate["mu"]), col = "blue", add = TRUE)


# ----fit other model types

model_poisson <- glm(Chao1 ~ Temp + Sex + Origin + Length + Weight +Month, data = data, family = poisson)
model_glmnb <- glm.nb(Chao1 ~ Temp + Sex + Origin + Length + Weight +Month, data = data, maxit = 1000)


autoplot(model_lm)
autoplot(model_poisson)
autoplot(model_glmnb)
autoplot(model_gamma_glm)

# Compare models using AIC
AIC(model_lm, model_glmnb, model_poisson, model_gamma_glm)
BIC(model_lm, model_glmnb, model_poisson, model_gamma_glm)

#Deviance (Lower deviance values indicate a better fit to the data)
deviance(model_glmnb)
deviance(model_gamma_glm)

#pseud r squared (Pseudo R-squared values, indicating a similar proportion of variance explained)
pR2(model_glmnb)
pR2(model_gamma_glm)

#cross validation (Lower delta values in cross-validation indicate better predictive performance)
cv_glm_nb <- cv.glm(data, model_glmnb, K = 10)
cv_glm_gamma <- cv.glm(data, model_gamma_glm, K = 10)
print(cv_glm_nb$delta)
print(cv_glm_gamma$delta)


setwd("C:/Users/Paddy/Documents/PhD/Chapter River(new)/Assembly/Month/Viva")

# check collinearity of number predictors
# Select relevant columns from your data
library(RColorBrewer)
# Select predictors
# Select your predictors
predictors <- c("Temp", "Waterlevel", "Flow", "Conductivity", "DO")

# Calculate the correlation matrix for the selected predictors
cor_matrix <- cor(data[, predictors])

# Save as SVG with larger sizes
svg("correlation_plot_300dpi.svg", width = 10, height = 8)

# Plot the correlation matrix using corrplot
corrplot(cor_matrix, 
         type = "upper",              # Only show upper triangle
         order = "hclust",            # Hierarchical clustering to order variables
         col = brewer.pal(n = 8, name = "RdYlBu"),  # Color palette
         number.cex = 1.2,            # Size of numbers in the circles
         tl.cex = 1.2,                # Size of variable labels (x and y axis labels)
         tl.col = "black",            # Color of variable labels
         main = NULL,
         main.cex = 2.5,              # Size of the title text
         addCoef.col = "white",       # Add the correlation numbers in the circles
         number.font = 2,             # Font of the numbers
         number.digits = 2,           # Number of digits to display in the correlation numbers
         cex.axis = 1.8)              # Size of the axis labels

# Close the SVG device to save the file
dev.off()


# Create a formula for the linear model
formula <- as.formula(paste("Chao1 ~", paste(predictors, collapse = " + ")))

# Fit a linear model using the predictors
model <- lm(formula, data = data)

check_model(model)
vif_values <- vif(model)
print(vif_values)

# Select the predictor columns from the dataset
predictor_data <- data[predictors]
cor_matrix <- cor(predictor_data, use = "complete.obs")

# Save the correlation plot to a PNG file
png("correlation_plot.png", width = 800, height = 600)
corrplot(cor_matrix, method = "color", addCoef.col = "black", tl.col = "black", tl.srt = 45)
dev.off()

# use MuMIn for best model

options(na.action= "na.fail")

model <- glm(Chao1 ~ Month + Temp + Sex + Origin + Length + Weight , family = Gamma(link = "log"), data = data)
summary(model)
# Create all possible combinations of predictors

MuMIn::dredge(global.model = model)


model.best <- glm(Chao1 ~ Month , family = Gamma(link = "log"), data = data)
summary(model.best)
check_model(model.best)
check_predictions(model.best)

ggsave(filename = "Chao1_checkmodel.svg", plot = check, width = 14, height = 10, dpi = 300, units = "cm")

deviance(model.best)
pR2(model.best)
cv_model.best <- cv.glm(data, model.best, K = 10)
print(cv_model.best$delta)

m1 <- glm(Chao1 ~ Month + Sex + Weight,family = Gamma(link = "log"), data=data) 
m2 <- glm(Chao1 ~ Month + Weight ,family = Gamma(link = "log"), data=data) 
m3 <- glm(Chao1 ~ Month + Sex + Weight,family = Gamma(link = "log"), data=data) 
m4 <- glm(Chao1 ~ Month, family = Gamma(link = "log"), data=data)


#check 4 best models
check_model (m1)

library(see)
compare_performance(m1, m2, m3, m4)
p1 <- plot(compare_performance(m1, m2, m3, m4))
test_performance(m1, m2, m3, m4)
p1


ggsave(filename = "Chao1_indices.svg", plot = p1, width = 14, height = 10, dpi = 300, units = "cm")

#same can be done for Shannon

#--plot
#Define colors
library(ggpubr)

myCol <- c( Diptera ="#460000" , Coleoptera ="#9400d3", Trichoptera ="#8B8000",Ephemeroptera= "#FFC0CB", Plecoptera ="#FFCC00" )
myCol2 <- c( AJan_19="#0066CC",  CMar_19="#79deaa" ,EMay_19="#0e5b05",
             FJun_19="#ede07b", GJul_19="#ffa346", KNov_19="#a71b2b")

bxp <- ggboxplot(data, x = "Month", y = "Chao1", fill = "Month", 
                 palette = myCol2 )
bxp

bxp <- bxp + theme_bw() + theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank()) +
  theme(text=element_text(family="Calibri", face="bold", size=11))+
  theme(legend.position = "none")
#+
#   scale_x_discrete(labels=c("Jan_19" = "Jan19", "Mar_19" = "Mar19",
#                             "May_19" = "May19","Jun_19" = "Jun19", "Jul_19" = "Jul19",
#                            "Nov_19" = "Nov19", "May_19_H" = "May19",
#                           "Jun_19_S" = "Jun19"))

bxp

ggsave(filename = "Chao1_alphadiversity_Month_nostat.svg", plot = bxp, width = 14, height = 10, dpi = 300, units = "cm")


# --------------------------------------------------------------------------------------Beta Diversity

#load phyloseq object
pseq <- readRDS("phyloseq_file.rds")
#Phyloseq to microeco object
dataset <- phyloseq2meco(pseq)

#calculate beta diversity in microeco
dataset$cal_betadiv(unifrac = TRUE)
dataset$beta_diversity

#subset the dataset to your needs
Gut <- clone(dataset)
Gut$sample_table <- subset(Gut$sample_table, SampleType == "PC")
Gut$tidy_dataset()

#or Water
Water<- clone(dataset)
Water$sample_table <- subset(Water$sample_table, SampleType == "Water")
Water$tidy_dataset()

#or Macro
Macroinvertebrate <- clone(dataset)
Macroinvertebrate$sample_table <- subset(Macroinvertebrate$sample_table, SampleType == "Macroinvertebrate")
Macroinvertebrate$tidy_dataset()


myCol2 <- c( AJan_19="#0066CC",  CMar_19="#79deaa" ,EMay_19="#0e5b05",
             FJun_19="#ede07b", GJul_19="#ffa346", KNov_19="#a71b2b")

#you can also compute by Origin or another parameter of choice...
t1 <- trans_beta$new(dataset = Gut, group = "Month2", measure = "bray")

t1$cal_ordination(ordination = "PCoA")

b1 <- t1$plot_ordination(plot_color = "Month", plot_shape= "Origin" , plot_type = c("point", "centroid"),color_values = myCol2,
                         point_size = 4, point_alpha = 1, centroid_segment_alpha = 0.3, centroid_segment_size = 0.3, centroid_segment_linetype = 1,
                         shape_values =c(8,15,17,19),
                         add_sample_label =  )

b1 <- b1 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
  theme(text=element_text(family="Calibri", face="bold", size=11))+
  theme(legend.key.size = unit(0.5, 'cm'))+
  theme(legend.title = element_text(size=9))

b1

ggsave(filename = "Beta_Month_Origin.svg", plot = b1, width = 14, height = 14, dpi = 300, units = "cm")

# manova for all groups when manova_all = TRUE
t1$cal_manova(manova_all = TRUE)
t1$res_manova
# manova for each paired groups
t1$cal_manova(manova_all = FALSE)
t1$res_manova
# manova for specified group set: such as "Group + Type"
t1$cal_manova(manova_set = "Origin + Month")
t1$res_manova
# PERMDISP for the whole comparison and for each paired groups
t1$cal_betadisper()
t1$res_betadisper

#--------------------------------------------------------------------------------------DBRDA---------------------------------------------------------
# Distance Based Redundancy Analysis; here Water as example. Adapt accordingly
library(agricolae)
t1 <- trans_env$new(dataset = Water, env_cols = 39:43)

t1$cal_diff(method = "anova", group = "Month")
# place all the plots into a list
tmp <- list()
for(i in colnames(t1$data_env)){
  tmp[[i]] <- t1$plot_diff(measure = i, add_sig_text_size = 5, xtext_size = 12) + theme(plot.margin = unit(c(0.1, 0, 0, 1), "cm"))
}
plot(gridExtra::arrangeGrob(grobs = tmp, ncol = 3))

t1$cal_ordination(method = "dbRDA", use_measure = "bray")
t1$trans_ordination()
t1$res_ordination
t1$res_ordination_trans

anova(t1$res_ordination, permutations = 999, by = "terms")
vegan::RsquareAdj(t1$res_ordination)

b1 <- t1$plot_ordination(plot_color = "Month", plot_shape = NULL, plot_type = c("point"),color_values = myCol2, 
                         point_alpha = 0.8, centroid_segment_alpha = 0.3, centroid_segment_size = 0.3, centroid_segment_linetype = 1)
b1 <- b1 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
  theme(text=element_text(family="Calibri", face="bold", size=11))+
  theme(legend.key.size = unit(0.5, 'cm'))+
  theme(legend.title = element_text(size=9))

b1

ggsave(filename = "RDA_Month_Water.svg", plot = b1, width = 20, height = 14, dpi = 300, units = "cm")

# ---------------------------------------------------------------------------------------------------Source Tracking-------------------------------------------
# This step also requires some manipulation in excel in order to make the data rdy for source tracking analysis

rm(list=ls())
#load phyloseq object
pseq <- readRDS("phyloseq_file.rds")
dataset <- phyloseq2meco(pseq)
meta <- meta_to_df(pseq)


dataset$cal_abund()
genus_data <- dataset[["taxa_abund"]][["Genus"]]
str(genus_data)
genus_data <- data.frame(Genus = row.names(genus_data), genus_data)

# Split sample_info by Month2
sample_info_split <- split(sample_info, sample_info$Month2)

# Get unique Month2 values in correct order
unique_months <- unique(sample_info$Month2)

# Create an empty list for genus_data subsets
genus_data_split <- list()

# Loop through each unique Month2 value and subset genus_data accordingly
for (month in unique_months) {
  # Get the corresponding SampleIDs for the given Month2
  sample_ids <- sample_info$SampleID[sample_info$Month2 == month]
  
  # Subset genus_data to only include those SampleIDs
  genus_data_subset <- genus_data[, c("Genus", sample_ids), drop = FALSE]
  
  # Store the subsetted data in the list
  genus_data_split[[as.character(month)]] <- genus_data_subset
}

# Ensure file names match Month2 values correctly
for (month in unique_months) {
  correct_filename <- gsub(" ", "_", as.character(month))  # Replace spaces with underscores if needed
  
  write.csv(sample_info_split[[as.character(month)]], 
            paste0("sample_info_", correct_filename, ".csv"), 
            row.names = FALSE)
  
  write.csv(genus_data_split[[as.character(month)]], 
            paste0("genus_data_", correct_filename, ".csv"), 
            row.names = FALSE)
}

##Now FEAST
#after we created all our files we can run FEAST
rm(list=ls())

library(devtools)
devtools::install_github("cozygene/FEAST")
library(FEAST)

setwd("C:/Users/Paddy/Documents/PhD/Assembly/FEAST")

metadata <- Load_metadata(metadata_path = "Meta_Jan19.txt")
otus <- Load_CountMatrix(CountMatrix_path = "genus_data_Jan_19.txt")
otus <- ceiling(otus)

FEAST_output <- FEAST(C = otus, metadata = metadata, different_sources_flag = 0, dir_path = "C:/Users/Paddy/Documents/PhD/Assembly/FEAST",
                      outfile="Jan")

#Mar
metadata <- Load_metadata(metadata_path = "Meta_Mar19.txt")
otus <- Load_CountMatrix(CountMatrix_path = "genus_data_Mar_19.txt")
otus <- ceiling(otus)

FEAST_output <- FEAST(C = otus, metadata = metadata, different_sources_flag = 0, dir_path = "C:/Users/Paddy/Documents/PhD/Assembly/FEAST",
                      outfile="Mar")
#May
metadata <- Load_metadata(metadata_path = "Meta_May19.txt")
otus <- Load_CountMatrix(CountMatrix_path = "genus_data_May_19.txt")
otus <- ceiling(otus)

FEAST_output <- FEAST(C = otus, metadata = metadata, different_sources_flag = 0, dir_path = "C:/Users/Paddy/Documents/PhD/Assembly/FEAST",
                      outfile="May")

#June
metadata <- Load_metadata(metadata_path = "Meta_Jun19.txt")
otus <- Load_CountMatrix(CountMatrix_path = "genus_data_Jun_19.txt")
otus <- ceiling(otus)

FEAST_output <- FEAST(C = otus, metadata = metadata, different_sources_flag = 0, dir_path = "C:/Users/Paddy/Documents/PhD/Assembly/FEAST",
                      outfile="Jun")

#July
metadata <- Load_metadata(metadata_path = "Meta_Jul19.txt")
otus <- Load_CountMatrix(CountMatrix_path = "genus_data_Jul_19.txt")
otus <- ceiling(otus)

FEAST_output <- FEAST(C = otus, metadata = metadata, different_sources_flag = 0, dir_path = "C:/Users/Paddy/Documents/PhD/Assembly/FEAST",
                      outfile="Jul")

#Nov
metadata <- Load_metadata(metadata_path = "Meta_Nov19.txt")
otus <- Load_CountMatrix(CountMatrix_path = "genus_data_Nov_19.txt")
otus <- ceiling(otus)

FEAST_output <- FEAST(C = otus, metadata = metadata, different_sources_flag = 0, dir_path = "C:/Users/Paddy/Documents/PhD/Assembly/FEAST",
                      outfile="Nov")

#Feast result
rm(list=ls())
setwd("C:/Users/Paddy/Documents/PhD/Assembly")

sample_info <- read.csv("metafile_mit_alpha.csv",row.names = 1, header=T, stringsAsFactors = FALSE)
Feast <- read.csv("Feast_Result_Meta.csv",row.names = NULL, header=T, stringsAsFactors = FALSE)

str(Feast)
str(sample_info)

merged_df <- merge(Feast, sample_info[, c("SampleID", "Sex", "Origin", "Shannon", "Observed")], 
                   by.x = "SampleID", by.y = "SampleID", all.x = TRUE)

str(merged_df)


# Check for correlation between taxa contributions and the explanatory variables
cor(merged_df[, c("Diptera", "Coleoptera", "Trichoptera", "Ephemeroptera", "Plecoptera", "Water_top", "Water_mid", "Water_bot", "Unknown")], 
    merged_df[, c("Shannon", "Observed")])

# Check for any unusual or missing values
summary(merged_df)

library(GGally)

# Create pairwise scatter plot matrix to explore relationships between the taxa contributions and other variables
ggpairs(merged_df[, c("Diptera", "Coleoptera", "Trichoptera", "Ephemeroptera", "Plecoptera", "Water_top", "Water_mid", "Water_bot", "Unknown", "Shannon", "Observed")])

# Boxplot for 'Diptera' by Month
boxplot(Diptera ~ Month, data = merged_df, main = "Diptera by Month", ylab = "Contribution")

# Boxplot for 'Coleoptera' by Month
boxplot(Coleoptera ~ Month, data = merged_df, main = "Coleoptera by Month", ylab = "Contribution")

long_df <- merged_df %>%
  select(SampleID, Month, Diptera, Coleoptera, Trichoptera, Ephemeroptera, Plecoptera, Water_top, Water_mid, Water_bot, Unknown) %>%
  pivot_longer(cols = c(Diptera, Coleoptera, Trichoptera, Ephemeroptera, Plecoptera, Water_top, Water_mid, Water_bot, Unknown),
               names_to = "Source", values_to = "Contribution")
str(long_df)

avg_contribution <- long_df %>%
  group_by(Month, Source) %>%
  summarise(Average_Contribution = mean(Contribution, na.rm = TRUE)) %>%
  ungroup()

avg_contribution <- avg_contribution %>%
  group_by(Month) %>%
  mutate(Percentage_Contribution = (Average_Contribution / sum(Average_Contribution)) * 100) %>%
  ungroup()


taxon_colors <- c( Diptera ="#460000" , Coleoptera ="#9400d3", Trichoptera ="#8B8000",Ephemeroptera= "#FFC0CB", Plecoptera ="#FFCC00" )

avg_contribution$Month <- factor(avg_contribution$Month, 
                                 levels = c("Jan", "Mar", "May", "Jun", "Jul", "Nov"))

avg_contribution <- as.data.frame(avg_contribution)
str(avg_contribution)
# Create the stacked bar plot

ggplot(avg_contribution, aes(x = Month, y = Percentage_Contribution, fill = Source)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage_Contribution, 1), "%")), 
            position = position_stack(vjust = 0.5), 
            color = "grey95", size = 2.5) +  # Adjust text size here
  scale_fill_manual(values = taxon_colors) +  # Apply custom colors
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # Rotate x-axis labels
    plot.title = element_blank(),  # Remove the title
    panel.grid.major = element_blank(),  # Optionally remove grid lines
    panel.grid.minor = element_blank(),  # Optionally remove minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, size = 1)  # Add border around the plot
  )

# Check stats

shapiro_test_results <- long_df %>%
  group_by(Source, Month) %>%
  summarise(p_value = shapiro.test(Contribution)$p.value)

# Check p-values for normality
shapiro_test_results

statistical_test <- function(df) {
  # Shapiro test results
  shapiro_results <- df %>%
    group_by(Source, Month) %>%
    summarise(p_value = shapiro.test(Contribution)$p.value) %>%
    ungroup()
  
  # Check if normality assumption is met for each source
  df %>%
    group_by(Source) %>%
    do({
      shapiro_result <- shapiro_results %>%
        filter(Source == unique(.$Source)) %>%
        pull(p_value)
      
      # If normality is met (p > 0.05), run ANOVA
      if (all(shapiro_result > 0.05)) {
        aov_result <- aov(Contribution ~ Month, data = .)
        tidy(aov_result)
      } else {
        # Otherwise, run Kruskal-Wallis
        kruskal_result <- kruskal.test(Contribution ~ Month, data = .)
        tibble(Source = unique(.$Source), p_value = kruskal_result$p.value)
      }
    })
}

# Run the statistical test for each source
statistical_test_results <- statistical_test(long_df)
statistical_test_results

#Diptera
# Subset the data for Diptera
diptera_data <- long_df %>% filter(Source == "Diptera")

# Perform pairwise Wilcoxon test
posthoc_diptera <- pairwise.wilcox.test(diptera_data$Contribution, diptera_data$Month, p.adjust.method = "bonferroni")

# View results
print(posthoc_diptera)


# Origin
long_df <- merged_df %>%
  select(SampleID, Origin, Diptera, Coleoptera, Trichoptera, Ephemeroptera, Plecoptera, Water_top, Water_mid, Water_bot, Unknown) %>%
  pivot_longer(cols = c(Diptera, Coleoptera, Trichoptera, Ephemeroptera, Plecoptera, Water_top, Water_mid, Water_bot, Unknown),
               names_to = "Source", values_to = "Contribution")
str(long_df)

avg_contribution <- long_df %>%
  group_by(Origin, Source) %>%
  summarise(Average_Contribution = mean(Contribution, na.rm = TRUE)) %>%
  ungroup()

avg_contribution <- avg_contribution %>%
  group_by(Origin) %>%
  mutate(Percentage_Contribution = (Average_Contribution / sum(Average_Contribution)) * 100) %>%
  ungroup()


taxon_colors <- c( Diptera ="#460000" , Coleoptera ="#9400d3", Trichoptera ="#8B8000",Ephemeroptera= "#FFC0CB", Plecoptera ="#FFCC00" )

avg_contribution$Origin <- factor(avg_contribution$Origin, 
                                  levels = c("F", "HFF", "HWF", "W"))

avg_contribution <- as.data.frame(avg_contribution)
str(avg_contribution)
# Create the stacked bar plot

ggplot(avg_contribution, aes(x = Origin, y = Percentage_Contribution, fill = Source)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage_Contribution, 1), "%")), 
            position = position_stack(vjust = 0.5), 
            color = "grey95", size = 2.5) +  # Adjust text size here
  scale_fill_manual(values = taxon_colors) +  # Apply custom colors
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # Rotate x-axis labels
    plot.title = element_blank(),  # Remove the title
    panel.grid.major = element_blank(),  # Optionally remove grid lines
    panel.grid.minor = element_blank(),  # Optionally remove minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, size = 1)  # Add border around the plot
  )


# ----------------------------------------------------------------------------------------------------------Neutral model------------------------------------------------------

#load phyloseq object
pseq <- readRDS("phyloseq_file.rds")
map<- meta_to_df(pseq)
otu <- otu_table(pseq)

nReads=10000

otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # relative abundance  
occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
  rownames_to_column('otu')

ggplot(data=occ_abun, aes(x=log10(otu_rel), y=otu_occ)) +
  geom_point(pch=21, fill='white') +
  labs(x="log10(mean relative abundance)", y="Occupancy")


PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
  gather(sample_id, abun, -otu) %>%
  left_join(map, by = 'sample_id') %>%
  group_by(otu, Month) %>%
  summarise(plot_freq=sum(abun>0)/length(abun),        # frequency of detection between time points
            coreSite=ifelse(plot_freq == 1, 1, 0), # 1 only if occupancy 1 with specific genotype, 0 if not
            detect=ifelse(plot_freq > 0, 1, 0)) %>%    # 1 if detected and 0 if not detected with specific genotype
  group_by(otu) %>%
  summarise(sumF=sum(plot_freq),
            sumG=sum(coreSite),
            nS=length(Month)*2,
            Index=(sumF+sumG)/nS) # calculating weighting Index based on number of time points detected and 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by='otu') %>%
  transmute(otu=otu,
            rank=Index) %>%
  arrange(desc(rank))

BCaddition <- NULL


otu_start=otu_ranked$otu[1]
start_matrix <- as.matrix(otu[otu_start,])
#start_matrix <- t(start_matrix)
x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 
BCaddition <- rbind(BCaddition,df_s)


for(i in 2:500){
  otu_add=otu_ranked$otu[i]
  add_matrix <- as.matrix(otu[otu_add,])
  # add_matrix <- t(add_matrix)
  start_matrix <- rbind(start_matrix, add_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_a <- data.frame(x_names,x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))
}


x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names,x)
names(df_full)[2] <- length(rownames(otu))
BCfull <- left_join(BCaddition,df_full, by='x_names')

rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))
Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

#Creating thresholds for core inclusion 

#Method: 
#A) Elbow method (first order difference) (script modified from https://pommevilla.github.io/random/elbows.html)
fo_difference <- function(pos){
  left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
  right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
  return(left - right)
}
BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
elbow <- which.max(BC_ranked$fo_diffs)

#B) Final increase in BC similarity of equal or greater then x%; choose x accordingly 
lastCall <- last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.04)])))
ggplot(BC_ranked[1:400,], aes(x=factor(BC_ranked$rank[1:400], levels=BC_ranked$rank[1:400]))) +
  geom_point(aes(y=proportionBC)) +
  theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
  geom_vline(xintercept=elbow, lty=3, col='red', cex=.5) +
  geom_vline(xintercept=lastCall, lty=3, col='blue', cex=.5) +
  labs(x='ranked OTUs',y='Bray-Curtis similarity') +
  annotate(geom="text", x=elbow+10, y=.15, label=paste("Elbow method"," (",elbow,")", sep=''), color="red")+    
  annotate(geom="text", x=lastCall-4, y=.08, label=paste("Last 2% increase (",lastCall,')',sep=''), color="blue")

#Create a column defining "core" OTUs

occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:lastCall]] <- 'core'

spp=t(otu)
taxon=as.vector(rownames(otu))
pool=t(otu)

#log transform spp test
#spp <- log(spp + 1)
#pool<- log(pool + 1)
#define env as pool
#file.path2="otu_table_Jan_19E.csv"
#pool<-read.table (file.path2, check.names = FALSE, header = TRUE, dec = ".", sep = ",", row.names = 1, comment.char = "")

#pool=t(pool)
#file.path2="tax_table_rar.csv"
#taxon<-read.table(file.path2, check.names = FALSE, header = TRUE, dec = ".", sep = ",", row.names = 1, comment.char = "")


#Models for the whole community

# Burns Model scnm.fit

sncm.fit <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
  require(minpack.lm)
  require(Hmisc)
  require(stats4)
  
  options(warn=-1)
  
  #Calculate the number of individuals per community
  N <- mean(apply(spp, 1, sum))
  
  #Calculate the average relative abundance of each taxa across communities
  if(is.null(pool)){
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  }
  
  #Calculate the occurrence frequency of each taxa across communities
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  #Combine
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  
  #Calculate the limit of detection
  d = 1/N
  
  ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
  m.ci <- confint(m.fit, 'm', level=0.95)
  
  ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
  sncm.LL <- function(m, sigma){
    R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
    R = dnorm(R, 0, sigma)
    -sum(log(R))
  }
  m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
  
  ##Calculate Akaike's Information Criterion (AIC)
  aic.fit <- AIC(m.mle, k=2)
  bic.fit <- BIC(m.mle)
  
  ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
  
  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for binomial model
  bino.LL <- function(mu, sigma){
    R = freq - pbinom(d, N, p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.bino <- AIC(bino.mle, k=2)
  bic.bino <- BIC(bino.mle)
  
  ##Goodness of fit for binomial model
  bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
  Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
  
  bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for Poisson model
  pois.LL <- function(mu, sigma){
    R = freq - ppois(d, N*p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.pois <- AIC(pois.mle, k=2)
  bic.pois <- BIC(pois.mle)
  
  ##Goodness of fit for Poisson model
  pois.pred <- ppois(d, N*p, lower.tail=FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
  
  pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Results
  if(stats==TRUE){
    fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
    fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
    return(fitstats)
  } else {
    A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
    A <- as.data.frame(A)
    colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
    if(is.null(taxon)){
      B <- A[order(A[,1]),]
    } else {
      B <- merge(A, taxon, by=0, all=TRUE)
      row.names(B) <- B[,1]
      B <- B[,-1]
      B <- B[order(B[,1]),]
    }
    return(B)
  }
}

# run model
#file.path3="otu_table_env.csv"
#pool<-read.table (file.path3, check.names = FALSE, header = TRUE, dec = ".", sep = ",", row.names = 1, comment.char = "") 


obs.np=sncm.fit(spp, taxon, stats=FALSE, pool=NULL)
obs.np <- subset(obs.np, select = -c(y))
obs.np <- na.omit(obs.np)
sta.np=sncm.fit(spp, taxon, stats=TRUE, pool=NULL)
sta.np.16S <- sta.np 

above.pred=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness
below.pred=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness

ap = obs.np$freq > (obs.np$pred.upr)
bp = obs.np$freq < (obs.np$pred.lwr)


b1 <- ggplot() +
  geom_point(data=occ_abun[occ_abun$fill=='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='cornsilk', alpha=.3)+
  geom_point(data=occ_abun[occ_abun$fill!='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='darkblue',size=2)+
  geom_line(color='black', data=obs.np, size=1, aes(y=obs.np$freq.pred, x=log10(obs.np$p)), alpha=1) +
  geom_line(color='black', lty='dashed', size=1, data=obs.np, aes(y=obs.np$pred.upr, x=log10(obs.np$p)), alpha=1)+
  geom_line(color='black', lty='dashed', size=1, data=obs.np, aes(y=obs.np$pred.lwr, x=log10(obs.np$p)), alpha=1)+
  labs(x="log10(mean relative abundance)", y="Occupancy")

b1 <- b1 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
  theme(text=element_text(family="Calibri", face="bold", size=11))+
  theme(legend.key.size = unit(0.5, 'cm'))+
  theme(legend.title = element_text(size=9))

b1

ggsave(filename = "Core_Occupancy.svg", plot = b1, width = 20, height = 14, dpi = 300, units = "cm")



# merge dataframes to obtain core taxa

library(tibble)
df <- tibble::rownames_to_column(obs.np, "otu")

df_merge <- merge(df,occ_abun,by="otu",all.x = TRUE, all.y=TRUE)


#find upper and lower OTUs + middle

otu.upper<-df_merge[df_merge$freq>df_merge$pred.upr,]
otu.under<-df_merge[df_merge$freq<df_merge$pred.lwr,]
otu.middle <- df_merge[df_merge$freq >= df_merge$pred.lwr & df_merge$freq <= df_merge$pred.upr, ]

# add upper,neutral and lower tag to the df_merge dataframe

# Create a new column called 'model' in df
df_merge$model <- NA

# Loop through each row in df
for (i in 1:nrow(df_merge)) {
  
  # Check if the otu_ID is in otu.upper
  if (df_merge[i, "otu"] %in% otu.upper$otu) {
    df_merge[i, "model"] <- "upper"
  }
  
  # Check if the otu_ID is in otu.middle
  else if (df_merge[i, "otu"] %in% otu.middle$otu) {
    df_merge[i, "model"] <- "neutral"
  }
  
  # Check if the otu_ID is in otu.under
  else if (df_merge[i, "otu"] %in% otu.under$otu) {
    df_merge[i, "model"] <- "under"
  }
  
}


write.csv(df_merge, "model_stats_complete.csv", row.names = FALSE)

#plot with upper/under
b1 <- ggplot() +
  # geom_point(data=occ_abun[occ_abun$fill=='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='white', alpha=.3)+
  geom_point(data=otu.middle, aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='white',size=1.5)+
  geom_point(data=otu.upper, aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='cornflowerblue',size=1.5)+
  geom_point(data=otu.under, aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='coral',size=1.5)+
  # geom_point(data=occ_abun[occ_abun$fill!='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21,shape = 21,size=2)+
  geom_line(color='black', data=obs.np, size=1, aes(y=obs.np$freq.pred, x=log10(obs.np$p)), alpha=1) +
  geom_line(color='black', lty='dashed', size=1, data=obs.np, aes(y=obs.np$pred.upr, x=log10(obs.np$p)), alpha=1)+
  geom_line(color='black', lty='dashed', size=1, data=obs.np, aes(y=obs.np$pred.lwr, x=log10(obs.np$p)), alpha=1)+
  labs(x="log10(mean relative abundance)", y="Occurance frequency")

b1 <- b1 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
  theme(text=element_text(family="Calibri", face="bold", size=11))+
  theme(legend.key.size = unit(0.5, 'cm'))+
  theme(legend.title = element_text(size=9))

b1



ggsave(filename = "Neutral_model_New.svg", plot = b1, width = 18, height = 14, dpi = 300, units = "cm")


otu.upper %>% count(fill)
otu.under %>% count(fill)
otu.middle %>% count(fill)
occ_abun %>% count(fill)

# viva
otu.under
otu.upper
otu_new <- tibble::rownames_to_column(otu, "otu")

otu_model <- otu_new %>%
  filter(otu %in% otu.upper$otu | otu %in% otu.under$otu)
write.csv(otu_model,"otu_deterministic.csv", row.names = FALSE)

# subset otu.s upper dataframe to have only the core OTUs, in order to get only core otus that are selected by the host
CoreDF <- otu.upper %>% filter(fill=="core")

otu_new <- tibble::rownames_to_column(otu, "otu")

# subset to genereta upper core otu table
otu_core<- otu_new %>%
  filter(otu %in% CoreDF$otu)

write.csv(otu_core,"otu_core_upper_new.csv", row.names = FALSE)

# subset otu.s under dataframe to have only the core OTUs that are under
CoreDF2 <- otu.under %>% filter(fill=="core")

otu_new2 <- tibble::rownames_to_column(otu, "otu")

# subset to generet under core otu table
otu_core2<- otu_new2 %>%
  filter(otu %in% CoreDF2$otu)

write.csv(otu_core2,"otu_core_under_new.csv", row.names = FALSE)

# subset otu.s under dataframe to have only the core OTUs that are nuetral
CoreDF3 <- otu.middle %>% filter(fill=="core")

otu_new3 <- tibble::rownames_to_column(otu, "otu")

# subset to generet neutral core otu table
otu_core3<- otu_new3 %>%
  filter(otu %in% CoreDF3$otu)

write.csv(otu_core3,"otu_core_neutral_new.csv", row.names = FALSE)


#subset only upper OTUs
otu_model<- otu_new %>%
  filter(otu %in% otus.upper$otu)

write.csv(otu_model,"otu_upper.csv", row.names = FALSE)

#save model stats
write.csv(sta.np,"model_Stats_new.csv", row.names = T)
write.csv(df_merge,"otu_stats_new.csv", row.names = T)

#otu all; to include stats in all models by running river_eco and combine it with tax table
write.csv(otu,"otu_all.csv", row.names = T)

### CORE pie

setwd("C:/Users/Paddy/Documents/PhD/Chapter River(new)/Assembly/Month/Core")
pie <- read.csv("core_barchart.csv",row.names = 1, header=T, stringsAsFactors = FALSE)

library(ggplot2)
library(ggrepel)
library(tidyverse)

# Get the positions
df2 <- pie %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

b1 <- ggplot(pie, aes(x = "" , y = value, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "RdBu", direction=-1) +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 6.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Group")) +
  theme_void()

b1 <- b1 + 
  theme(text=element_text(family="Calibri", face="bold", size=11)) 

b1

ggsave(filename = "NoCore_OTU_pie_reverse.svg", plot = b1, width = 14, height = 14, dpi = 300, units = "cm")

