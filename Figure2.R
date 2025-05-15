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

setwd("YourPath")

#load phyloseq object
pseq <- readRDS("phyloseq_file.rds")

#Phyloseq to microeco object
dataset <- phyloseq2meco(pseq)

# Print the updated sample data
print(sample_data(pseq))
sample_data <- meta_to_df(pseq)
str(sample_data)

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



#Distance based redundancy analysis


