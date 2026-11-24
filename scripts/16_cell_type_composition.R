library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

obj <- readRDS(
  "data/processed/GSE171213_HC_PD_annotated.rds"
)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

meta <- obj@meta.data
meta$cell <- rownames(meta)

donors <- meta %>%
  distinct(donor, condition)

cell_types <- sort(unique(meta$cell_type))

observed_counts <- meta %>%
  count(donor, condition, cell_type, name = "n_cells")

# Add zero-count donor/cell-type combinations so every donor is represented.
composition <- crossing(
  donors,
  cell_type = cell_types
) %>%
  left_join(
    observed_counts,
    by = c("donor", "condition", "cell_type")
  ) %>%
  mutate(n_cells = replace_na(n_cells, 0L)) %>%
  group_by(donor, condition) %>%
  mutate(
    total_cells = sum(n_cells),
    proportion = n_cells / total_cells
  ) %>%
  ungroup()

write.csv(
  composition,
  "results/tables/cell_type_proportions_by_donor.csv",
  row.names = FALSE
)

p <- ggplot(
  composition,
  aes(x = condition, y = proportion * 100, color = condition)
) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 2.5
  ) +
  facet_wrap(~ cell_type, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Cell-type proportions by donor",
    x = NULL,
    y = "Percent of retained cells"
  ) +
  theme(legend.position = "none")

ggsave(
  "results/figures/cell_type_proportions_by_donor.png",
  p,
  width = 11,
  height = 8,
  dpi = 300
)

message("Saved donor proportions table and figure.")

composition_tests <- composition %>%
  group_by(cell_type) %>%
  group_modify(~ {
    healthy <- .x$proportion[.x$condition == "Healthy"]
    disease <- .x$proportion[.x$condition == "Periodontitis"]

    tibble(
      healthy_median_percent = median(healthy) * 100,
      periodontitis_median_percent = median(disease) * 100,
      difference_percentage_points =
        (median(disease) - median(healthy)) * 100,
      p_value = wilcox.test(
        healthy, disease,
        exact = FALSE
      )$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(FDR_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(FDR_BH)

write.csv(
  composition_tests,
  "results/tables/cell_type_proportion_tests.csv",
  row.names = FALSE
)

print(composition_tests)