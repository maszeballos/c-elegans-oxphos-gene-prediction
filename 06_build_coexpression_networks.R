library(WGCNA)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()


DATA_DIR <- "data"
OUTPUT_DIR <- file.path("results", "coexpression_networks")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


NETWORK_CONFIGS <- list(
  list(
    name = "adult",
    expression_file = file.path(DATA_DIR, "adult_expression_time_adjusted_tpm.csv"),
    soft_power = 2,
    network_output = file.path(OUTPUT_DIR, "06.1_adult_coexpression_network.RData"),
    gene_module_output = file.path(OUTPUT_DIR, "06.3_adult_gene_modules.csv"),
    soft_power_plot_output = file.path(OUTPUT_DIR, "06.5_adult_soft_power_diagnostics.png")
  ),
  list(
    name = "embryo",
    expression_file = file.path(DATA_DIR, "embryo_expression_time_tpm.csv"),
    soft_power = 3,
    network_output = file.path(OUTPUT_DIR, "06.2_embryo_coexpression_network.RData"),
    gene_module_output = file.path(OUTPUT_DIR, "06.4_embryo_gene_modules.csv"),
    soft_power_plot_output = file.path(OUTPUT_DIR, "06.6_embryo_soft_power_diagnostics.png")
  )
)


load_expression_matrix <- function(expression_file) {
  if (!file.exists(expression_file)) {
    stop("Expression input file does not exist: ", expression_file)
  }

  expression_df <- read.csv(
    file = expression_file,
    header = TRUE,
    check.names = FALSE
  )

  if (!"cell_type" %in% colnames(expression_df)) {
    stop("Expression input file must contain a 'cell_type' column: ", expression_file)
  }

  rownames(expression_df) <- expression_df$cell_type
  expression_df$cell_type <- NULL

  expression_df[] <- lapply(expression_df, as.numeric)
  expression_df
}


filter_expression_matrix <- function(expression_data) {
  sample_gene_check <- goodSamplesGenes(expression_data, verbose = 3)

  if (!sample_gene_check$allOK) {
    if (sum(!sample_gene_check$goodGenes) > 0) {
      printFlush(paste(
        "Removing genes:",
        paste(names(expression_data)[!sample_gene_check$goodGenes], collapse = ", ")
      ))
    }

    if (sum(!sample_gene_check$goodSamples) > 0) {
      printFlush(paste(
        "Removing samples:",
        paste(rownames(expression_data)[!sample_gene_check$goodSamples], collapse = ", ")
      ))
    }

    expression_data <- expression_data[
      sample_gene_check$goodSamples,
      sample_gene_check$goodGenes
    ]
  }

  expression_data
}


plot_soft_power_diagnostics <- function(expression_data, output_file) {
  powers <- c(1:10, seq(12, 20, by = 2))
  soft_power_diagnostics <- pickSoftThreshold(
    expression_data,
    powerVector = powers,
    verbose = 5
  )

  png(output_file, width = 2400, height = 1200, res = 300)
  par(mfrow = c(1, 2))

  plot(
    soft_power_diagnostics$fitIndices[, 1],
    soft_power_diagnostics$fitIndices[, 2],
    xlab = "Soft threshold (power)",
    ylab = "Scale-free topology model fit, signed R^2",
    type = "n",
    main = "Scale independence"
  )
  text(
    soft_power_diagnostics$fitIndices[, 1],
    soft_power_diagnostics$fitIndices[, 2],
    labels = soft_power_diagnostics$fitIndices[, 1],
    col = "red"
  )
  abline(h = 0.80, col = "red")

  plot(
    soft_power_diagnostics$fitIndices[, 1],
    soft_power_diagnostics$fitIndices[, 5],
    xlab = "Soft threshold (power)",
    ylab = "Mean connectivity",
    type = "n",
    main = "Mean connectivity"
  )
  text(
    soft_power_diagnostics$fitIndices[, 1],
    soft_power_diagnostics$fitIndices[, 5],
    labels = soft_power_diagnostics$fitIndices[, 1],
    col = "red"
  )

  dev.off()

  soft_power_diagnostics
}


build_coexpression_network <- function(config) {
  message("Building co-expression network: ", config$name)

  expression_data <- load_expression_matrix(config$expression_file)
  expression_data <- filter_expression_matrix(expression_data)

  soft_power_diagnostics <- plot_soft_power_diagnostics(
    expression_data,
    config$soft_power_plot_output
  )
  soft_power <- config$soft_power

  adjacency_matrix <- adjacency(expression_data, power = soft_power)
  tom_similarity <- TOMsimilarity(adjacency_matrix)
  tom_dissimilarity <- 1 - tom_similarity

  gene_tree <- hclust(as.dist(tom_dissimilarity), method = "average")

  dynamic_modules <- cutreeDynamic(
    dendro = gene_tree,
    distM = tom_dissimilarity,
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    minClusterSize = 50
  )

  module_colors <- labels2colors(dynamic_modules)
  names(module_colors) <- colnames(expression_data)

  module_eigengenes <- moduleEigengenes(
    expression_data,
    colors = module_colors
  )$eigengenes

  module_eigengenes <- module_eigengenes[
    ,
    colSums(is.na(module_eigengenes)) == 0,
    drop = FALSE
  ]

  merged_modules <- mergeCloseModules(
    expression_data,
    module_colors,
    cutHeight = 0.25,
    verbose = 3
  )

  merged_colors <- merged_modules$colors
  names(merged_colors) <- colnames(expression_data)

  merged_module_eigengenes <- merged_modules$newMEs

  gene_modules <- data.frame(
    gene_id = names(merged_colors),
    module = unname(merged_colors),
    stringsAsFactors = FALSE
  )

  write.csv(gene_modules, config$gene_module_output, row.names = FALSE)

  save(
    expression_data,
    adjacency_matrix,
    tom_similarity,
    tom_dissimilarity,
    gene_tree,
    dynamic_modules,
    module_colors,
    module_eigengenes,
    merged_colors,
    merged_module_eigengenes,
    gene_modules,
    soft_power,
    soft_power_diagnostics,
    file = config$network_output
  )

  message("Saved network: ", config$network_output)
  message("Saved gene-module table: ", config$gene_module_output)

  invisible(gene_modules)
}


for (network_config in NETWORK_CONFIGS) {
  build_coexpression_network(network_config)
}
