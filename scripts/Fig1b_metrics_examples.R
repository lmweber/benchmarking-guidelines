##############################################################################
# Using data from collection of benchmarking datasets on IMLS penticton server
# - http://imlspenticton.uzh.ch/robinson_lab/benchmark_collection/
# - section 'Differential gene expression, human'
#
# Lukas Weber, November 2018
##############################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)


# ---------
# load data
# ---------

# load example data
dir_data <- "iCOBRA_example_data/NB_625_625"

cobradata <- COBRAData_from_text(
  truth_file = file.path(dir_data, "NB_625_625_5spc_repl1.truth.txt"),
  result_files = file.path(dir_data, "NB_625_625_5spc_repl1.merged_results.txt"),
  feature_id = "gene"
)

# subset and repeat genes to create more 'typical' TPR-FDR curves
keep_true <- cobradata@truth$gene[cobradata@truth$differential.expression == 1 & cobradata@truth$truelog2foldchanges > 1.6]
set.seed(101)
keep_false <- cobradata@truth$gene[sample(1001:12499, 6000)]
genes_keep <- c(rep(keep_true, 5), keep_false)

cobradata@pval <- cobradata@pval[genes_keep, ]
cobradata@padj <- cobradata@padj[genes_keep, ]
cobradata@truth <- cobradata@truth[genes_keep, ]

# subset to keep only some methods
cobradata@pval <- cobradata@pval[, c("edgeRGLM", "voomlimma")]
cobradata@padj <- cobradata@padj[, c("edgeRGLM", "voomlimma")]

# change method names
colnames(cobradata@pval) <- c("method_1", "method_2")
colnames(cobradata@padj) <- c("method_1", "method_2")


# ------------------
# calculate measures
# ------------------

# calculate performance measures
cobradata <- calculate_adjp(cobradata)

cobraperf <- calculate_performance(
  cobradata,
  binary_truth = "differential.expression",
  cont_truth = "truelog2foldchanges"
)

slotNames(cobraperf)


# --------------
# generate plots
# --------------

# generate plots
cobraplot <- prepare_data_for_plot(
  cobraperf, 
  colorscheme = c("forestgreen", "darkorange", "black"), 
  facetted = TRUE
)

# ROC curves
p_ROC <- 
  plot_roc(cobraplot) + 
  coord_fixed() + 
  theme(
    strip.background = element_blank(), 
    strip.text.x = element_blank(), 
    axis.text.x = element_text(size = 10), 
    axis.text.y = element_text(size = 10), 
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
    legend.position = "none"
  )

# TPR-FDR curves
p_TPRFDR <- 
  plot_fdrtprcurve(cobraplot) + 
  coord_fixed() + 
  theme(
    strip.background = element_blank(), 
    strip.text.x = element_blank(), 
    axis.text.x = element_text(size = 10), 
    axis.text.y = element_text(size = 10), 
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
    legend.position = "none"
  )

# PR curves
# (note: using 'plot_roc' function)
# (note: can't include some methods using this setup due to NAs, e.g. DESeq2)
precision <- cobraplot@fdrtprcurve$TP / (cobraplot@fdrtprcurve$TP + cobraplot@fdrtprcurve$FP)
recall <- cobraplot@fdrtprcurve$TPR
stopifnot(nrow(cobraplot@roc) == nrow(cobraplot@fdrtprcurve))
cobraplot_PR <- cobraplot
cobraplot_PR@roc$FPR <- recall
cobraplot_PR@roc$TPR <- precision

p_PR <- 
  plot_roc(cobraplot_PR) + 
  xlab("recall") + 
  ylab("precision") + 
  coord_fixed() + 
  theme(
    strip.background = element_blank(), 
    strip.text.x = element_blank(), 
    axis.text.x = element_text(size = 10), 
    axis.text.y = element_text(size = 10), 
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
    legend.text = element_text(size = 12), 
    legend.key.size = unit(8, "mm")
  )


# ----------
# save plots
# ----------

# save plots (multiple panels using cowplot)
plot_grid(p_ROC, p_TPRFDR, p_PR, ncol = 3, rel_widths = c(1, 1, 1.375), labels = "b", label_size = 22)

ggsave("Fig1b_metrics_examples.pdf", height = 3.25, width = 11.25)
ggsave("Fig1b_metrics_examples.png", height = 3.25, width = 11.25)


