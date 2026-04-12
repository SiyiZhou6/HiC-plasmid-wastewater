#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(readr)
  library(optparse)
  suppressWarnings({
    quietly_loaded <- suppressMessages(require(ggplot2, quietly = TRUE) &&
                                       require(ggnetwork, quietly = TRUE) &&
                                       require(RColorBrewer, quietly = TRUE) &&
                                       require(grid, quietly = TRUE))
  })
})
opt <- OptionParser(option_list=list(
  make_option("--final_remove", help="Filtered Hi-C edge list (TSV with columns: col1, col2, type_col1, type_col2, V3)"),
  make_option("--min_weight", type="integer", default=15, help="Minimum Hi-C contact count for edge retention (default: 15)"),
  make_option("--inclusive", action="store_true", default=FALSE, help="Use >= threshold instead of >."),
  make_option("--steps", type="integer", default=10, help="Walktrap random walk steps (default: 10)"),
  make_option("--amr", help="(optional) AMRFinderPlus output TSV to highlight AMR-linked clusters"),
  make_option("--out_csv", help="Output CSV path for cluster membership (Node, Cluster)"),
  make_option("--plot_prefix", help="(optional) Prefix for PNG network plots")
))
opt <- parse_args(opt)
stopifnot(!is.null(opt$final_remove), !is.null(opt$out_csv))
df <- read.delim(opt$final_remove, header = TRUE, check.names = FALSE)
if (!all(c("type_col1","type_col2","col1","col2","V3") %in% names(df))) {
  stop("final_remove file missing required columns: type_col1, type_col2, col1, col2, V3")
}
if (opt$inclusive) {
  edges <- df %>% dplyr::filter(type_col1 == "plasmid", type_col2 == "plasmid", V3 >= opt$min_weight)
} else {
  edges <- df %>% dplyr::filter(type_col1 == "plasmid", type_col2 == "plasmid", V3 >  opt$min_weight - 1)
}
edges <- edges %>% dplyr::select(col1, col2, V3) %>% dplyr::distinct()
if (nrow(edges) == 0) {
  writeLines("No edges passed the filter; writing an empty membership file.")
  readr::write_csv(tibble::tibble(Node=character(), Cluster=integer()), opt$out_csv)
  quit(save="no", status=0)
}
g <- igraph::graph_from_data_frame(edges[, c("col1","col2")], directed = FALSE)
igraph::E(g)$weight <- edges$V3
wt <- igraph::cluster_walktrap(g, weights = igraph::E(g)$weight, steps = opt$steps)
membership_df <- tibble::tibble(Node = igraph::V(g)$name, Cluster = igraph::membership(wt))
readr::write_csv(membership_df, opt$out_csv)
message("Wrote membership: ", opt$out_csv)
if (!is.null(opt$plot_prefix) && exists("ggnetwork")) {
  amr_clusters <- integer(0)
  if (!is.null(opt$amr) && file.exists(opt$amr)) {
    amr <- tryCatch({
      read.delim(opt$amr, header = TRUE, check.names = FALSE)
    }, error = function(e) {
      suppressWarnings(read.delim(opt$amr, header = FALSE, skip = 1, fill = TRUE, check.names = FALSE))
    })
    cn <- colnames(amr)
    es <- cn[grepl("^Element(\.| )?subtype$", cn, ignore.case = TRUE)][1]
    cid <- cn[grepl("^Contig(\.| )?id$", cn, ignore.case = TRUE)][1]
    if (!is.na(es) && !is.na(cid)) {
      amr_contigs <- amr %>% dplyr::filter(.data[[es]] == "AMR") %>% dplyr::pull(dplyr::all_of(cid)) %>% unique()
      ctoc <- membership_df
      amr_clusters <- ctoc %>% dplyr::filter(Node %in% amr_contigs) %>% dplyr::pull(Cluster) %>% unique()
    }
  }
  n <- ggnetwork::ggnetwork(g, arrow.gap = 0)
  membership_mapping <- membership_df %>% dplyr::mutate(Cluster = as.factor(Cluster))
  n <- n %>% dplyr::left_join(membership_mapping, by = setNames("Node", "name")) %>%
    dplyr::mutate(Cluster = ifelse(is.na(Cluster), "unknown", Cluster),
           highlight_color = ifelse(as.integer(as.character(Cluster)) %in% amr_clusters, "tomato2", "grey"))
  num_clusters <- length(unique(n$Cluster)); if (num_clusters < 3) num_clusters <- 3
  mypal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(num_clusters)
  clustered_plot <- ggplot2::ggplot(n, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggnetwork::geom_edges(color = "white", alpha = 0.4) +
    ggnetwork::geom_nodes(ggplot2::aes(color = Cluster), size = 3) +
    ggplot2::scale_color_manual(values = mypal) +
    ggplot2::theme_void() +
    ggplot2::ggtitle("Clustered Graph View") +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5))
  amr_plot <- ggplot2::ggplot(n, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggnetwork::geom_edges(color = "white", alpha = 0.4) +
    ggnetwork::geom_nodes(ggplot2::aes(fill = highlight_color), color = "black", shape = 21, size = 3) +
    ggplot2::scale_fill_identity() +
    ggplot2::theme_void() +
    ggplot2::ggtitle("Graph with AMR-Linked Clusters Highlighted") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5))
  ggplot2::ggsave(paste0(opt$plot_prefix, "_clustered.png"), clustered_plot, width = 8, height = 6, dpi = 300)
  ggplot2::ggsave(paste0(opt$plot_prefix, "_amr_highlight.png"), amr_plot, width = 8, height = 6, dpi = 300)
  message("Wrote plots with prefix: ", opt$plot_prefix)
}
