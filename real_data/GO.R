library(data.table)

file <- "D:/github/GL-neuronized-svg/real_data/bio/GO_overlap.tsv"

# Read the file as text
txt <- readLines(file, warn = FALSE)

# Locate the two table headers
go_header <- grep("^Gene Set Name\\t", txt)
matrix_header <- grep("^Entrez Gene Id\\t", txt)
matrix_marker <- grep("^Gene/Gene Set Overlap Matrix$", txt)

# ============================================================
# 1. Read GO enrichment results
# ============================================================

go_lines <- txt[go_header:(matrix_marker - 1L)]
go_lines <- go_lines[nzchar(go_lines)]

go_result <- fread(
  text = paste(go_lines, collapse = "\n"),
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

# ============================================================
# 2. Read the gene--GO overlap matrix
# ============================================================

gene_lines <- txt[matrix_header:length(txt)]

gene_go_matrix <- fread(
  text = paste(gene_lines, collapse = "\n"),
  sep = "\t",
  header = TRUE,
  fill = TRUE,
  na.strings = c("", "NA"),
  check.names = FALSE
)

# Check dimensions
dim(go_result)       # expected: 100 x 7
dim(gene_go_matrix)  # expected: 341 x 103

library(data.table)
library(ggplot2)

# ============================================================
# 1. Output directory
# ============================================================

out_dir <- paste0("D:/github/GL-neuronized-svg/real_data/bio")

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 2. Basic checks
# ===========================================================

# First three columns contain gene information
gene_info_cols <- names(gene_go_matrix)[1:3]

# The remaining 100 columns correspond to GO terms
go_term_names <- names(gene_go_matrix)[4:103]

# ============================================================
# 3. Match the 100 matrix columns to the GO enrichment table
# ============================================================

term_match <- match(
  go_term_names,
  go_result[["Gene Set Name"]]
)


# Reorder terms by FDR q-value and then p-value
term_order <- order(
  as.numeric(go_result[["FDR q-value"]][term_match]),
  as.numeric(go_result[["p-value"]][term_match])
)

go_term_names <- go_term_names[term_order]
term_match <- term_match[term_order]

# Assign shortened identifiers
gene_set_ids <- seq_along(go_term_names)

# ============================================================
# 4. Convert the gene--GO matrix into a 0/1 matrix
# ============================================================

binary_matrix <- vapply(
  go_term_names,
  function(term_name) {
    x <- gene_go_matrix[[term_name]]
    
    as.integer(
      !is.na(x) &
        nzchar(trimws(as.character(x)))
    )
  },
  integer(nrow(gene_go_matrix))
)

binary_matrix <- as.matrix(binary_matrix)

colnames(binary_matrix) <- gene_set_ids

# Use the second column, Gene Symbol, as row names
gene_symbols <- make.unique(
  as.character(gene_go_matrix[[2]])
)

rownames(binary_matrix) <- gene_symbols

dim(binary_matrix)
# 341 × 100

n_enriched_go_terms <- rowSums(binary_matrix)

# Supplementary Table S1
S1_svg_list <- data.frame(
  `Entrez Gene ID` = gene_go_matrix[[1]],
  `Gene symbol` = gene_symbols,
  `Gene description` = gene_go_matrix[[3]],
  `Number of enriched GO terms` = n_enriched_go_terms,
  check.names = FALSE
)

write.csv(
  S1_svg_list,
  file = file.path(
    out_dir,
    "Supplementary_Table_S1_341_SVG_list.csv"
  ),
  row.names = FALSE,
  na = ""
)

S2_GO_terms <- as.data.frame(
  go_result[term_match, ],
  check.names = FALSE
)

S2_GO_terms <- data.frame(
  `Gene Set ID` = gene_set_ids,
  S2_GO_terms,
  check.names = FALSE
)


write.csv(
  S2_GO_terms,
  file = file.path(
    out_dir,
    "Supplementary_Table_S2_GO_term_key.csv"
  ),
  row.names = FALSE,
  na = ""
)

# ============================================================
# 5. Order genes by the number of associated GO terms
# ============================================================

gene_membership_count <- rowSums(binary_matrix)

binary_plot <- binary_matrix

# ============================================================
# 6. Convert to long format
# ============================================================

plot_data <- as.data.table(
  binary_plot,
  keep.rownames = "Gene symbol"
)

plot_data <- melt(
  plot_data,
  id.vars = "Gene symbol",
  variable.name = "Gene set id",
  value.name = "Membership"
)

plot_data[, `Gene set id` := factor(
  `Gene set id`,
  levels = gene_set_ids
)]

plot_data[, `Gene symbol` := factor(
  `Gene symbol`,
  levels = rownames(binary_plot),#rev(rownames(binary_plot))
)]

# ============================================================
# 7. Plot all 341 genes
# ============================================================

p_all <- ggplot(
  plot_data,
  aes(
    y = `Gene set id`,
    x = `Gene symbol`,
    fill = factor(Membership)
  )
) +
  geom_tile(
    linewidth = 0.08,
    color = "grey90"
  ) +
  scale_fill_manual(
    values = c(
      "0" = "white",
      "1" = "#173A6E"
    ),
    guide = "none"
  ) +
  labs(
    y = "Gene Set ID", #"Ranked GO term",
    x = "Gene"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 5
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(
      fill = NA,
      linewidth = 0.4
    ),
    plot.margin = margin(
      5,
      5,
      5,
      5
    )
  )

ggsave(
  filename = file.path(
    out_dir,
    "Gene_GO_binary_matrix_all_341_genes.pdf"
  ),
  plot = p_all,
  width = 9,
  height = 12
)

ggsave(
  filename = file.path(
    out_dir,
    "Gene_GO_binary_matrix_all_341_genes.png"
  ),
  plot = p_all,
  width = 9,
  height = 14,
  dpi = 300
)
