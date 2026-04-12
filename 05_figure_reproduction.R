#!/usr/bin/env Rscript
# ============================================================
# Figure Reproduction Script
# Reads all data from a single Source_Data.xlsx
# Sheet names: Fig.1a, Fig.1c, Fig.2a, Fig.2b, Fig.2c, etc.
# ============================================================
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(stringr)
library(scales)

SOURCE_FILE <- "Source_Data.xlsx"
dir.create("./output/", showWarnings = FALSE)

# ============================================================
# FIGURE 1a - NMDS Ordination
# ============================================================
df_1a <- read_excel(SOURCE_FILE, sheet = "Fig.1a")
custom_colors_1a <- c(
  "Plasmid-associated community_A" = "#face97",
  "Plasmid-associated community_B" = "#f9aa47",
  "Plasmid-associated community_C" = "#dd8a23",
  "Total community_A"              = "#97cae1",
  "Total community_B"              = "#4b9bbe",
  "Total community_C"              = "#3e6d8d"
)
df_1a <- df_1a %>%
  mutate(Method_State = paste(Community_Type, State, sep = "_"))

p_1a <- ggplot(df_1a, aes(x = NMDS1, y = NMDS2, color = Method_State, shape = Stage)) +
  geom_point(size = 4.5, alpha = 0.95, stroke = 0.8) +
  scale_color_manual(values = custom_colors_1a) +
  labs(
    x = "NMDS1 (Bray-Curtis Dissimilarity)",
    y = "NMDS2 (Bray-Curtis Dissimilarity)",
    color = "Community Type_State",
    shape = "Sampling Stage"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right")

ggsave("./output/Fig1a_reproduced.pdf", p_1a, width = 10, height = 7)
cat("Saved: Fig1a\n")

# ============================================================
# FIGURE 1b - Shannon Diversity (Boxplot, faceted by Method)
# ============================================================
df_1b <- read_excel(SOURCE_FILE, sheet = "Fig.1b")

if ("Shannon index" %in% names(df_1b)) {
  df_1b <- df_1b %>% rename(Shannon = `Shannon index`)
}

df_1b <- df_1b %>%
  mutate(
    Method = factor(Method, levels = c("Plasmid-associated community", "Total community")),
    Site_Label = factor(Site_Label, levels = c("Inf", "AS", "Eff")),
    City_Label = factor(City_Label, levels = c("A", "B", "C"))
  )

kw_labels_1b <- df_1b %>%
  filter(is.finite(Shannon)) %>%
  group_by(Method) %>%
  summarise(
    p = kruskal.test(Shannon ~ Site_Label)$p.value,
    y_max = max(Shannon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sig = case_when(p < 0.05 ~ "*", TRUE ~ "ns"),
    label = paste0(sig, ", p = ", signif(p, 3)),
    x = 1.05,
    y = y_max + 0.08 * y_max
  )

dodge_w <- 0.75
p_1b <- ggplot(df_1b, aes(x = Site_Label, y = Shannon, fill = City_Label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5,
               position = position_dodge(dodge_w)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_w),
    size = 3, aes(color = City_Label, shape = City_Label)
  ) +
  facet_wrap(~ Method, ncol = 2, scales = "fixed") +
  geom_text(
    data = kw_labels_1b,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 4, fontface = "italic", color = "grey30"
  ) +
  labs(x = NULL, y = "Shannon Index") +
  scale_fill_manual(values = c("A" = "steelblue", "B" = "darkorange", "C" = "darkgreen")) +
  scale_color_manual(values = c("A" = "steelblue", "B" = "darkorange", "C" = "darkgreen")) +
  scale_shape_manual(values = c("A" = 16, "B" = 17, "C" = 15)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  ) +
  guides(fill = guide_legend(title = "State"),
         color = guide_legend(title = "State"),
         shape = guide_legend(title = "State"))

ggsave("./output/Fig1b_reproduced.pdf", p_1b, width = 10, height = 4)
cat("Saved: Fig1b\n")

# ============================================================
# FIGURE 1c - Stacked Bar (Relative Abundance)
# ============================================================
df_1c <- read_excel(SOURCE_FILE, sheet = "Fig.1c")

base_colors_1c <- c(
  "#8ED1E9", "#6EB7D3", "#529DBC", "#396D88", "#786184", "#927B9C", "#AC95B4",
  "#D4685D", "#D97870", "#DE8A84", "#E39D98", "#F4A259", "#F6B16F", "#F8C086",
  "#FAD09C", "#FCD8C0", "#FEE6D5", "#f1e9ae", "#cbdab7", "#89a467"
)
family_levels_1c <- c(setdiff(sort(unique(df_1c$Family)), "Other"), "Other")
color_map_1c <- setNames(
  c(colorRampPalette(base_colors_1c)(length(family_levels_1c) - 1), "#CCCCCC"),
  family_levels_1c
)
label_city <- c("A" = "State A", "B" = "State B", "C" = "State C")

df_1c <- df_1c %>%
  mutate(
    Family = factor(Family, levels = family_levels_1c),
    Stage = factor(Stage, levels = c("Inf", "AS", "Eff")),
    Community_Type = factor(Community_Type,
      levels = c("Plasmid-associated community", "Total community"))
  )

p_1c <- ggplot(df_1c, aes(x = Stage, y = Relative_Abundance, fill = Family)) +
  geom_bar(stat = "identity", width = 0.85) +
  facet_grid(Community_Type ~ Stage, labeller = labeller(State = label_city)) +
  scale_fill_manual(values = color_map_1c) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Treatment stages", y = "Relative Abundance (%)", fill = "Family") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        legend.position = "bottom")

ggsave("./output/Fig1c_reproduced.pdf", p_1c, width = 10, height = 12)
cat("Saved: Fig1c\n")

# ============================================================
# FIGURE 2a-c - Violin Plots
# ============================================================
city_colors <- c("A" = "steelblue", "B" = "darkorange", "C" = "darkgreen")
my_comparisons_state <- list(c("A", "B"), c("A", "C"), c("B", "C"))

source_2a <- read_excel(SOURCE_FILE, sheet = "Fig.2a")
source_2b <- read_excel(SOURCE_FILE, sheet = "Fig.2b")
source_2c <- read_excel(SOURCE_FILE, sheet = "Fig.2c")

make_violin_plot <- function(df, value_col, ylab) {
  ymax <- max(df[[value_col]], na.rm = TRUE)
  y_pos <- c(ymax * 1.10, ymax * 1.22, ymax * 1.34)
  ylim_top <- max(y_pos) * 1.06

  ggplot(df, aes(x = State, y = .data[[value_col]], fill = State)) +
    geom_violin(trim = TRUE, alpha = 0.5) +
    geom_boxplot(width = 0.06, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(aes(color = State), width = 0.15, size = 1, alpha = 0.5) +
    scale_fill_manual(values = city_colors) +
    scale_color_manual(values = city_colors) +
    labs(x = "State", y = ylab) +
    theme_classic(base_size = 18) +
    theme(legend.position = "none", plot.margin = margin(t = 12, r = 8, b = 8, l = 8)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
    coord_cartesian(ylim = c(NA, ylim_top), clip = "off") +
    stat_compare_means(comparisons = my_comparisons_state, method = "wilcox.test",
      p.adjust.method = "BH", label = "p.signif", hide.ns = FALSE,
      step.increase = 0, label.y = y_pos, bracket.size = 0.5)
}

plot_2a <- make_violin_plot(source_2a, "MAGs_per_cluster", "Host range\n(MAGs per plasmid cluster)")
plot_2b <- make_violin_plot(source_2b, "Classes_per_cluster", "Host range\n(classes per plasmid cluster)")
plot_2c <- make_violin_plot(source_2c, "Clusters_per_MAG", "Plasmid burden\n(clusters per MAG)")

fig2 <- ggarrange(plot_2a, plot_2b, plot_2c, ncol = 3, widths = c(1, 1, 1))
ggsave("./output/Fig2_reproduced.pdf", fig2, width = 15, height = 5)
cat("Saved: Fig2\n")

# ============================================================
# FIGURE 3a - Bipartite Network (Host-Plasmid Connections)
# ============================================================
if (!requireNamespace("ggnewscale", quietly = TRUE)) install.packages("ggnewscale", repos = "https://cloud.r-project.org")
library(ggnewscale)

df_3a <- read_excel(SOURCE_FILE, sheet = "Fig.3a")

if (any(df_3a$Stage %in% c("Inf", "AS", "Eff"))) {
  df_3a <- df_3a %>% mutate(Stage = case_when(
    Stage == "Inf" ~ "Influent", Stage == "AS" ~ "Activated sludge",
    Stage == "Eff" ~ "Effluent", TRUE ~ Stage))
}

Y_BOTTOM <- 100; Y_TOP <- 10000
BACT_X <- 1.52; PLAS_X <- 1.83
EDGE_X_START <- 1.55; EDGE_X_END <- 1.80
BACT_TILE_W <- 0.06; PLAS_TILE_W <- 0.04

all_classes_3a <- sort(unique(df_3a$Bacterial_Class))
custom_palette_3a <- c(
  "#d53e4f", "#f46d43", "#fdae61", "#f2d4b7",
  "#fee08b", "#e6f598", "#abdda4", "#66c2a5",
  "#5c8c9c", "#3288bd", "#c0c4c9", "#674794"
)
class_colors_3a <- setNames(rep_len(custom_palette_3a, length(all_classes_3a)), all_classes_3a)
location_colors_3a <- c("Influent" = "#66C2A5", "Activated sludge" = "#FC8D62", "Effluent" = "#8DA0CB")

plot_bipartite_site <- function(site_data, site_label) {
  bacterial_nodes <- site_data %>%
    distinct(Bacterial_Family, Bacterial_Class) %>%
    arrange(Bacterial_Class, Bacterial_Family) %>%
    mutate(n = n(),
           Position = Y_BOTTOM + (row_number() - 1) / (n - 1) * (Y_TOP - Y_BOTTOM),
           Color = class_colors_3a[Bacterial_Class]) %>% select(-n)
  bact_tile_h <- (Y_TOP - Y_BOTTOM) / nrow(bacterial_nodes) * 0.85

  plasmid_nodes <- site_data %>%
    distinct(Plasmid_Cluster, Sample, Stage) %>%
    arrange(Stage, Plasmid_Cluster) %>%
    mutate(n = n(),
           Position = Y_BOTTOM + (row_number() - 1) / (n - 1) * (Y_TOP - Y_BOTTOM),
           Color = location_colors_3a[Stage]) %>% select(-n)
  plas_tile_h <- (Y_TOP - Y_BOTTOM) / nrow(plasmid_nodes) * 0.85

  edges <- site_data %>%
    left_join(bacterial_nodes %>% select(Bacterial_Family, Position), by = "Bacterial_Family") %>%
    rename(From_Position = Position) %>%
    left_join(plasmid_nodes %>% select(Plasmid_Cluster, Position, Sample, Color),
              by = c("Plasmid_Cluster", "Sample")) %>%
    rename(To_Position = Position)

  ggplot() +
    geom_segment(data = edges, aes(x = EDGE_X_START, xend = EDGE_X_END,
      y = From_Position, yend = To_Position, color = Color), size = 0.2, alpha = 0.5) +
    geom_tile(data = bacterial_nodes, aes(x = BACT_X, y = Position, fill = Bacterial_Class),
      width = BACT_TILE_W, height = bact_tile_h) +
    scale_fill_manual(values = class_colors_3a, name = "Bacterial Class", drop = TRUE) +
    new_scale_fill() +
    geom_tile(data = plasmid_nodes, aes(x = PLAS_X, y = Position, fill = Stage),
      width = PLAS_TILE_W, height = plas_tile_h) +
    scale_fill_manual(values = location_colors_3a, name = "Stage",
      limits = c("Influent", "Activated sludge", "Effluent")) +
    scale_color_identity(name = "Plasmid Links") +
    scale_x_continuous(limits = c(1.4, 1.95), breaks = c(1.52, 1.83),
                       labels = c("Bacteria", "Plasmids")) +
    scale_y_continuous(limits = c(Y_BOTTOM - 500, Y_TOP + 500), expand = c(0, 0)) +
    labs(title = paste("Host-Plasmid Connections -", site_label), x = NULL, y = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.line = element_blank(),
          plot.margin = margin(10, 10, 10, 10))
}

for (state_name in unique(df_3a$State)) {
  site_data <- df_3a %>% filter(State == state_name)
  p_3a <- plot_bipartite_site(site_data, state_name)
  state_short <- gsub("State ", "", state_name)
  ggsave(paste0("./output/Fig3a_", state_short, "_reproduced.pdf"), p_3a, width = 15, height = 14)
}
cat("Saved: Fig3a (one per state)\n")

# ============================================================
# FIGURE 3b - Dot Matrix
# ============================================================
df_3b <- read_excel(SOURCE_FILE, sheet = "Fig.3b")
df_3b$Stage <- factor(df_3b$Stage, levels = c("Inf", "AS", "Eff"))

class_order_3b <- df_3b %>%
  distinct(State, Bacterial_Class, Stages_Detected) %>%
  arrange(desc(Stages_Detected)) %>% pull(Bacterial_Class) %>% unique()
df_3b$Bacterial_Class <- factor(df_3b$Bacterial_Class, levels = class_order_3b)

p_3b <- ggplot(df_3b, aes(x = Stage, y = reorder(Bacterial_Class, Stages_Detected), group = Bacterial_Class)) +
  geom_line(color = "gray70", linewidth = 0.5) +
  geom_point(aes(fill = as.factor(Stages_Detected)), shape = 21, size = 2, color = "black") +
  scale_fill_manual(values = c("1" = "#fdae61", "2" = "#4575b4", "3" = "#313695"), name = "# of Stages") +
  facet_grid(. ~ State) +
  labs(x = "Stages", y = "Bacterial Class") +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 9), axis.title = element_text(face = "bold"),
        panel.grid = element_blank(), panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom", legend.title = element_text(face = "bold"))

ggsave("./output/Fig3b_reproduced.pdf", p_3b, width = 6, height = 6)
cat("Saved: Fig3b\n")

# ============================================================
# FIGURE 3c - Heatmap
# ============================================================
df_3c <- read_excel(SOURCE_FILE, sheet = "Fig.3c")
df_3c <- df_3c %>%
  mutate(State = factor(State, levels = c("State A", "State B", "State C")),
         Stage = factor(Stage, levels = c("Influent", "AS", "Eff")),
         Bin_Category = factor(Bin_Category, levels = c("1 Bin", "2 Bins", "3-5 Bins", ">5 Bins")),
         Class_Category = factor(Class_Category, levels = c("1 Class", "2-3 Classes", ">3 Classes")))

p_3c <- ggplot(df_3c, aes(x = Bin_Category, y = Class_Category, fill = Plasmid_Count)) +
  geom_tile(color = "white") + geom_text(aes(label = Plasmid_Count), size = 4) +
  scale_fill_gradient(low = "#deebf7", high = "#3182bd", name = "# Plasmid clusters") +
  facet_grid(Stage ~ State) +
  labs(x = "Number of bins (MAG) per plasmid cluster", y = "Number of classes") +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom", strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_hline(yintercept = c(3.6, 7.2), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(4.6, 9.2), linetype = "dashed", color = "gray40")

ggsave("./output/Fig3c_reproduced.pdf", p_3c, width = 8, height = 8)
cat("Saved: Fig3c\n")

# ============================================================
# FIGURE 4a-b - Boxplots by Stage
# ============================================================
stage_levels <- c("Influent", "Activated Sludge", "Effluent")
stage_colors <- c("Influent" = "#1b9e77", "Activated Sludge" = "#d95f02", "Effluent" = "#7570b3")
my_comparisons_stage <- list(c("Influent", "Activated Sludge"), c("Influent", "Effluent"), c("Activated Sludge", "Effluent"))

fig4a_data <- read_excel(SOURCE_FILE, sheet = "Fig.4a")
fig4b_data <- read_excel(SOURCE_FILE, sheet = "Fig.4b")
fig4a_data$Stage <- factor(fig4a_data$Stage, levels = stage_levels)
fig4b_data$Stage <- factor(fig4b_data$Stage, levels = stage_levels)

kw_a <- kruskal.test(Unique_Bacterial_Classes ~ Stage, data = fig4a_data)
p4a <- ggplot(fig4a_data, aes(x = Stage, y = Unique_Bacterial_Classes, fill = Stage)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.18, alpha = 0.35, size = 1) +
  scale_fill_manual(values = stage_colors) +
  stat_compare_means(comparisons = my_comparisons_stage, method = "wilcox.test",
    label = "p.signif", hide.ns = FALSE, step.increase = 0.08, bracket.size = 0.5) +
  annotate("text", x = 1, y = max(fig4a_data$Unique_Bacterial_Classes, na.rm = TRUE) * 1.05,
    label = paste0("p = ", signif(kw_a$p.value, 2)), hjust = 0, size = 4, fontface = "italic") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
  labs(tag = "a", x = NULL, y = "Number of Classes") +
  theme_classic(base_size = 14) + theme(legend.position = "none", plot.tag = element_text(face = "bold", size = 18))

kw_b <- kruskal.test(Unique_MAGs ~ Stage, data = fig4b_data)
p4b <- ggplot(fig4b_data, aes(x = Stage, y = Unique_MAGs, fill = Stage)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.18, alpha = 0.35, size = 1) +
  scale_fill_manual(values = stage_colors) +
  stat_compare_means(comparisons = my_comparisons_stage, method = "wilcox.test",
    label = "p.signif", hide.ns = FALSE, step.increase = 0.08, bracket.size = 0.5) +
  annotate("text", x = 1, y = max(fig4b_data$Unique_MAGs, na.rm = TRUE) * 1.05,
    label = paste0("p = ", signif(kw_b$p.value, 2)), hjust = 0, size = 4, fontface = "italic") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
  labs(tag = "b", x = NULL, y = "Number of MAGs") +
  theme_classic(base_size = 14) + theme(legend.position = "none", plot.tag = element_text(face = "bold", size = 18))

fig4ab <- plot_grid(p4a, p4b, ncol = 2, align = "hv")
ggsave("./output/Fig4ab_reproduced.pdf", fig4ab, width = 12, height = 6)
cat("Saved: Fig4ab\n")

# ============================================================
# FIGURE 4c - Scatter Plot (Host Range vs Burden with Correlation)
# ============================================================
df_4c <- read_excel(SOURCE_FILE, sheet = "Fig.4c")
stage_colors_full <- c("Influent" = "#1b9e77", "Activated Sludge" = "#d95f02", "Effluent" = "#7570b3")

cor_stats_4c <- df_4c %>%
  group_by(State) %>%
  summarise(cor_result = list(cor.test(Unique_Hosts_MAGs, Unique_Bacterial_Classes, method = "spearman")),
            .groups = "drop") %>%
  mutate(
    rho = sapply(cor_result, function(x) round(x$estimate, 2)),
    p_value = sapply(cor_result, function(x) {
      if (x$p.value < 0.001) { "<0.001" } else { format(round(x$p.value, 3), nsmall = 3) }
    }),
    label = paste0("r = ", rho, ", p ", ifelse(p_value == "<0.001", "<0.001", paste0("= ", p_value)))
  )

label_positions_4c <- df_4c %>%
  group_by(State) %>%
  summarise(
    x = max(Unique_Hosts_MAGs, na.rm = TRUE) * 0.6,
    y = max(Unique_Bacterial_Classes, na.rm = TRUE) * 0.99,
    .groups = "drop"
  ) %>%
  left_join(cor_stats_4c, by = "State")

p_4c <- ggplot(df_4c, aes(x = Unique_Hosts_MAGs, y = Unique_Bacterial_Classes)) +
  geom_point(aes(color = Stage, shape = Plasmid_Type), alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "grey80", linewidth = 0.7, alpha = 0.15) +
  geom_text(data = label_positions_4c, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 6, fontface = "bold") +
  facet_grid(. ~ State, labeller = labeller(State = c("A" = "State A", "B" = "State B", "C" = "State C"))) +
  scale_color_manual(values = stage_colors_full) +
  scale_shape_manual(values = c("AMR Cluster" = 16, "Non-AMR Cluster" = 17)) +
  labs(
    x = "Number of Unique Hosts (MAGs)",
    y = "Number of Unique Bacterial Classes",
    color = "Stage",
    shape = "Plasmid Type"
  ) +
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
    strip.text = element_text(face = "bold", size = 16, color = "black"),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 16),
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  )

ggsave("./output/Fig4c_reproduced.pdf", p_4c, width = 14, height = 6)
cat("Saved: Fig4c\n")

# ============================================================
# FIGURE 5d - ARG-Host Dot Plot (with subclass-based gene ordering)
# ============================================================
df_5d <- read_excel(SOURCE_FILE, sheet = "Fig.5d")

df_5d <- df_5d %>%
  mutate(
    Concordance = factor(Concordance,
      levels = c("Concordant", "Non-concordant")),
    alpha_group = ifelse(Concordance == "Non-concordant", 0.22, 0.75)
  )

fam_levels_5d <- df_5d %>%
  group_by(Host_Family) %>%
  summarise(total = sum(Plasmid_Clusters), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(Host_Family)

gene_order_5d <- df_5d %>%
  group_by(ARG_Gene, ARG_Subclass) %>%
  summarise(total = sum(Plasmid_Clusters), .groups = "drop") %>%
  mutate(subclass_key = ifelse(ARG_Subclass %in% c("Unknown", "Mixed"),
                               paste0("ZZZ_", ARG_Subclass), ARG_Subclass)) %>%
  arrange(subclass_key, desc(total), ARG_Gene) %>%
  pull(ARG_Gene)

df_5d <- df_5d %>%
  mutate(
    State = factor(State, levels = c("A", "B", "C"), labels = c("State A", "State B", "State C")),
    Host_Family = factor(Host_Family, levels = fam_levels_5d),
    ARG_Gene = factor(ARG_Gene, levels = gene_order_5d),
    ARG_Subclass = factor(ARG_Subclass)
  )

p_5d <- ggplot(df_5d, aes(x = ARG_Gene, y = Host_Family)) +
  geom_point(aes(size = Plasmid_Clusters, color = ARG_Subclass,
                 shape = Concordance, alpha = alpha_group)) +
  facet_grid(. ~ State, scales = "fixed") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_alpha_identity(guide = "none") +
  scale_size_continuous(
    name = "# plasmid clusters",
    breaks = c(1, 2, 4, 8, 10),
    limits = c(1, max(10, max(df_5d$Plasmid_Clusters, na.rm = TRUE))),
    range = c(1.2, 6)
  ) +
  labs(x = "ARG", y = "Bacterial Host (Family)", color = "ARG subclass", shape = "Concordance") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.6, "cm"))

ggsave("./output/Fig5d_reproduced.pdf", p_5d, width = 18, height = 11)
cat("Saved: Fig5d\n")

# ============================================================
# SI FIGURE S1 - Plasmid Burden per Family (Boxplot)
# ============================================================
if (!requireNamespace("ggtext", quietly = TRUE)) install.packages("ggtext", repos = "https://cloud.r-project.org")
library(ggtext)

df_S1 <- read_excel(SOURCE_FILE, sheet = "Fig.S1")

class_totals_S1 <- df_S1 %>% count(Bacterial_Class, name = "Total") %>% arrange(desc(Total))
top10_S1 <- class_totals_S1 %>% slice_head(n = 10) %>% pull(Bacterial_Class)
df_S1 <- df_S1 %>% mutate(Class_Grouped = ifelse(Bacterial_Class %in% top10_S1, Bacterial_Class, "Other Class")) %>%
  filter(Class_Grouped != "Other Class")

uq_S1 <- sort(unique(df_S1$Class_Grouped))
pal_S1 <- colorRampPalette(c("#1f77b4","#ff7f0e","#2ca02c","#9467bd","#e377c2","#8c564b","#d62728","#17becf"))(length(uq_S1))
cmap_S1 <- setNames(pal_S1, uq_S1)

df_S1 <- df_S1 %>% mutate(Family_Styled = paste0("<span style='color:", cmap_S1[Class_Grouped], "'>", Bacterial_Family, "</span>"))
forder_S1 <- df_S1 %>% distinct(Family_Styled, Class_Grouped) %>% arrange(Class_Grouped, Family_Styled) %>% pull(Family_Styled)
df_S1$Family_Styled <- factor(df_S1$Family_Styled, levels = forder_S1)
df_S1$State <- factor(df_S1$State, levels = c("State A", "State B", "State C"))

pS1 <- ggplot(df_S1, aes(x = Plasmid_Count, y = Family_Styled)) +
  geom_boxplot(aes(fill = Class_Grouped), outlier.shape = NA, alpha = 0.4, width = 0.6, color = "black") +
  geom_jitter(aes(color = Class_Grouped), position = position_jitter(height = 0.2, width = 0), size = 1.2, alpha = 0.6) +
  facet_wrap(~ State, ncol = 3) +
  scale_fill_manual(values = cmap_S1, name = "Class") +
  scale_color_manual(values = cmap_S1, guide = "none") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "Number of plasmid clusters", y = "Plasmid host (Family level)", fill = "Class") +
  theme_classic(base_size = 13) +
  theme(axis.text.y = ggtext::element_markdown(size = 9), axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom", strip.text = element_text(size = 13, face = "bold")) +
  geom_hline(data = df_S1 %>% group_by(Class_Grouped) %>% summarise(mf = max(as.numeric(Family_Styled)), .groups = "drop"),
             aes(yintercept = mf + 0.5), linetype = "dashed", color = "gray50", inherit.aes = FALSE)

ggsave("./output/FigS1_reproduced.pdf", pS1, width = 12, height = 10)
cat("Saved: FigS1\n")

# ============================================================
# SI FIGURE S2a-b - Boxplots by Stage per State
# ============================================================
figS2a_data <- read_excel(SOURCE_FILE, sheet = "Fig.S2a")
figS2b_data <- read_excel(SOURCE_FILE, sheet = "Fig.S2b")
stage_colors_short <- c("Inf" = "#1b9e77", "AS" = "#d95f02", "Eff" = "#7570b3")
comps_short <- list(c("Inf", "AS"), c("Inf", "Eff"), c("AS", "Eff"))

figS2a_data$Stage <- factor(figS2a_data$Stage, levels = c("Inf", "AS", "Eff"))
figS2b_data$Stage <- factor(figS2b_data$Stage, levels = c("Inf", "AS", "Eff"))

pS2a <- ggplot(figS2a_data, aes(x = Stage, y = Unique_Bacterial_Classes, fill = Stage)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) + geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  scale_fill_manual(values = stage_colors_short) +
  stat_compare_means(method = "kruskal.test", label.y = max(figS2a_data$Unique_Bacterial_Classes) * 1.4, label = "p.format") +
  stat_compare_means(comparisons = comps_short, method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  facet_wrap(~ State, ncol = 3) + labs(x = "Stages", y = "Number of Classes") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black"))

pS2b <- ggplot(figS2b_data, aes(x = Stage, y = Unique_MAGs, fill = Stage)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) + geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  scale_fill_manual(values = stage_colors_short) +
  stat_compare_means(method = "kruskal.test", label.y = max(figS2b_data$Unique_MAGs) * 1.4, label = "p.format") +
  stat_compare_means(comparisons = comps_short, method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  facet_wrap(~ State, ncol = 3) + labs(x = "Stages", y = "Number of MAGs") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black"))

figS2 <- ggarrange(pS2a, pS2b, ncol = 2, labels = NULL, align = "hv")
ggsave("./output/FigS2_reproduced.pdf", figS2, width = 14, height = 6)
cat("Saved: FigS2\n")

# ============================================================
# SI FIGURE S3a - Stacked Bar (Mobility by State x Stage)
# ============================================================
df_S3a <- read_excel(SOURCE_FILE, sheet = "Fig.S3a")
df_S3a <- df_S3a %>% mutate(Stage = factor(Stage, levels = c("Inf", "AS", "Eff")))

pS3a <- ggplot(df_S3a, aes(x = Stage, y = Percentage, fill = Mobility_Category)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 2.5) +
  facet_wrap(~ State, labeller = labeller(State = function(x) paste("State", x))) +
  scale_fill_manual(
    values = c("conjugative" = "#3182bd", "mobilizable/conjugative" = "#9ecae1", "no_marker" = "#bdbdbd"),
    name = "Mobility category") +
  labs(x = "Treatment Stage", y = "Percentage (%)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

ggsave("./output/FigS3a_reproduced.pdf", pS3a, width = 8, height = 5)
cat("Saved: FigS3a\n")

# ============================================================
# SI FIGURE S3b - MOB/MPF Marker Heatmap
# ============================================================
df_S3b <- read_excel(SOURCE_FILE, sheet = "Fig.S3b")
top_S3b <- df_S3b %>% group_by(MOB_MPF_Marker) %>%
  summarise(tot = sum(Cluster_Count), .groups = "drop") %>%
  arrange(desc(tot)) %>% head(10) %>% pull(MOB_MPF_Marker)

pS3b <- df_S3b %>% filter(MOB_MPF_Marker %in% top_S3b) %>%
  mutate(MOB_MPF_Marker = factor(MOB_MPF_Marker, levels = rev(top_S3b))) %>%
  ggplot(aes(x = State, y = MOB_MPF_Marker, fill = Cluster_Count)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Cluster_Count), size = 3.5) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", name = "# clusters") +
  labs(x = "State", y = "MOB / MPF marker") +
  theme_minimal(base_size = 12) + theme(panel.grid = element_blank())

ggsave("./output/FigS3b_reproduced.pdf", pS3b, width = 5, height = 5)
cat("Saved: FigS3b\n")

# ============================================================
# SI FIGURE S3c - Conjugation Gene Heatmap
# ============================================================
df_S3c <- read_excel(SOURCE_FILE, sheet = "Fig.S3c")
top_S3c <- df_S3c %>% group_by(Conjugation_Gene) %>%
  summarise(tot = sum(Cluster_Count), .groups = "drop") %>%
  arrange(desc(tot)) %>% head(20) %>% pull(Conjugation_Gene)

pS3c <- df_S3c %>% filter(Conjugation_Gene %in% top_S3c) %>%
  mutate(Conjugation_Gene = factor(Conjugation_Gene, levels = rev(top_S3c))) %>%
  ggplot(aes(x = State, y = Conjugation_Gene, fill = Cluster_Count)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Cluster_Count), size = 2.8) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", name = "# clusters") +
  labs(x = "State", y = "Conjugation gene") +
  theme_minimal(base_size = 11) + theme(panel.grid = element_blank())

ggsave("./output/FigS3c_reproduced.pdf", pS3c, width = 5.5, height = 7)
cat("Saved: FigS3c\n")

# ============================================================
# SI FIGURE S4 - Key Host Genus Abundance (Line + Error Bars)
# ============================================================
df_S4 <- read_excel(SOURCE_FILE, sheet = "Fig.S4")
df_S4 <- df_S4 %>%
  mutate(Stage = factor(Stage, levels = c("Inf", "AS", "Eff")))

pS4 <- ggplot(df_S4, aes(x = Stage, y = Mean_Abundance_Pct, color = Host_Genus, group = Host_Genus)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.8) +
  geom_errorbar(aes(ymin = pmax(0, Mean_Abundance_Pct - SE),
                    ymax = Mean_Abundance_Pct + SE),
                width = 0.15, linewidth = 0.4) +
  facet_wrap(~ State, labeller = labeller(State = c("A" = "State A", "B" = "State B", "C" = "State C"))) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Treatment Stage", y = "Relative Abundance (%)", color = "Host Genus") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8, face = "italic"),
        strip.text = element_text(size = 11)) +
  guides(color = guide_legend(nrow = 2))

ggsave("./output/FigS4_reproduced.pdf", pS4, width = 9, height = 5)
cat("Saved: FigS4\n")

# ============================================================
# SI FIGURE S6a - Reads After QC and Mapped to MAGs
# ============================================================
df_S6a <- read_excel(SOURCE_FILE, sheet = "Fig.S6a")
state_cols_S6 <- c("A" = "#1b9e77", "B" = "#d95f02", "C" = "#7570b3")

s6a_long <- df_S6a %>%
  pivot_longer(cols = c(Reads_After_QC, Reads_Mapped_to_MAGs),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
    Reads_After_QC = "Reads after QC (R1+R2)",
    Reads_Mapped_to_MAGs = "Reads mapped to MAGs"))

pS6a <- ggplot(s6a_long, aes(x = State, y = Value, fill = State)) +
  geom_violin(trim = TRUE, alpha = 0.45) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95) +
  geom_jitter(width = 0.12, alpha = 0.75, size = 1.6) +
  facet_wrap(~ Metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = state_cols_S6) +
  stat_compare_means(method = "kruskal.test", label.x = 1, label.y.npc = 1, hjust = 0.4, vjust = -6, size = 3.5) +
  stat_compare_means(comparisons = list(c("A","B"), c("A","C"), c("B","C")),
    method = "wilcox.test", label = "p.signif", hide.ns = F, size = 3.5) +
  labs(x = "State", y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "bottom")

ggsave("./output/FigS6a_reproduced.pdf", pS6a, width = 10, height = 5)
cat("Saved: FigS6a\n")

# ============================================================
# SI FIGURE S6c - MAG Counts per Sample
# ============================================================
df_S6c <- read_excel(SOURCE_FILE, sheet = "Fig.S6c")

s6c_long <- df_S6c %>%
  pivot_longer(cols = c(MAGs_All, MAGs_MQHQ), names_to = "Set", values_to = "nMAG") %>%
  mutate(Set = recode(Set, MAGs_All = "All MAGs", MAGs_MQHQ = "MQ/HQ only"))

pS6c <- ggplot(s6c_long, aes(x = State, y = nMAG, fill = State)) +
  geom_violin(trim = TRUE, alpha = 0.45) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95) +
  geom_jitter(width = 0.12, alpha = 0.75, size = 1.6) +
  facet_wrap(~ Set, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = state_cols_S6) +
  stat_compare_means(method = "kruskal.test", label.x = 1, label.y.npc = 1, hjust = 0.4, vjust = -7, size = 3.5) +
  stat_compare_means(comparisons = list(c("A","B"), c("A","C"), c("B","C")),
    method = "wilcox.test", label = "p.signif", hide.ns = F, size = 3.5) +
  labs(x = "State", y = "# MAGs per sample") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "none")

ggsave("./output/FigS6c_reproduced.pdf", pS6c, width = 8, height = 5)
cat("Saved: FigS6c\n")

# ============================================================
# SI FIGURE S6e - Completeness vs Contamination Scatter
# ============================================================
df_S6e <- read_excel(SOURCE_FILE, sheet = "Fig.S6e")
stage_cols_S6 <- c("Inf" = "#1f78b4", "AS" = "#e31a1c", "Eff" = "#33a02c")

ann_S6e <- df_S6e %>%
  group_by(State) %>%
  summarise(n_total = n(), n_mq = sum(Is_MQHQ, na.rm = TRUE),
            pct_mq = 100 * n_mq / n_total, .groups = "drop") %>%
  mutate(label = sprintf("MQ/HQ: %d/%d (%.1f%%)", n_mq, n_total, pct_mq))

pS6e <- ggplot(df_S6e, aes(x = Contamination, y = Completeness, color = Stage)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_vline(xintercept = 10, linetype = "dashed") +
  geom_hline(yintercept = 50, linetype = "dashed") +
  facet_wrap(~ State, nrow = 1) +
  geom_text(data = ann_S6e, aes(x = Inf, y = -Inf, label = label),
            inherit.aes = FALSE, hjust = 1.05, vjust = -0.35, size = 5) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = stage_cols_S6) +
  labs(x = "Contamination (%)", y = "Completeness (%)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave("./output/FigS6e_reproduced.pdf", pS6e, width = 10, height = 5)
cat("Saved: FigS6e\n")

# ============================================================
# SI FIGURE S6f - Genome Size Distribution
# ============================================================
df_S6f <- read_excel(SOURCE_FILE, sheet = "Fig.S6f")

s6f_dat <- bind_rows(
  df_S6f %>% mutate(Set = "All MAGs"),
  df_S6f %>% filter(Is_MQHQ == TRUE) %>% mutate(Set = "MQ/HQ only")
) %>% mutate(Set = factor(Set, levels = c("All MAGs", "MQ/HQ only")))

pS6f <- ggplot(s6f_dat, aes(x = State, y = Genome_Size_bp, color = Stage)) +
  geom_violin(aes(group = State), fill = "grey90", color = NA, trim = TRUE) +
  geom_jitter(width = 0.12, alpha = 0.35, size = 1.2) +
  facet_wrap(~ Set, nrow = 1, scales = "free_y") +
  scale_color_manual(values = stage_cols_S6) +
  scale_y_log10() +
  stat_compare_means(aes(group = State), method = "kruskal.test", label.x = 1, label.y.npc = 1, hjust = 0.37, vjust = -8, size = 3) +
  stat_compare_means(aes(group = State), comparisons = list(c("A","B"), c("A","C"), c("B","C")),
    method = "wilcox.test", label = "p.signif", hide.ns = F, size = 3.5) +
  labs(x = "State", y = "Genome size (bp, log10)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave("./output/FigS6f_reproduced.pdf", pS6f, width = 8, height = 5)
cat("Saved: FigS6f\n")

# ============================================================
# SI FIGURE S8 - Plasmid Cluster Length Distribution (Violin)
# ============================================================
df_S8 <- read_excel(SOURCE_FILE, sheet = "Fig.S8")
site_colors_S8 <- c("Inf" = "steelblue", "AS" = "darkorange", "Eff" = "darkgreen")

pS8 <- ggplot(df_S8, aes(x = State, y = Total_Contig_Length)) +
  geom_violin(aes(fill = State), alpha = 0.5, trim = TRUE, color = NA) +
  geom_jitter(aes(color = Site), width = 0.2, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = c("A" = "steelblue", "B" = "darkorange", "C" = "darkgreen")) +
  scale_color_manual(values = site_colors_S8, limits = c("Inf", "AS", "Eff")) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "State", y = "Plasmid Length (bp)", fill = "State", color = "Sampling Site") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 13)
  )

ggsave("./output/FigS8_reproduced.pdf", pS8, width = 8, height = 6)
cat("Saved: FigS8\n")

# ============================================================
cat("\n== All figures reproduced successfully ==\n")
