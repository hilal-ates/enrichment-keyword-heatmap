library(readxl); library(tidyr); library(dplyr); library(ggplot2); library(tools)

keywords <- c("apoptosis","immune response", "stress","protein","DNA","autophagy")

count_one_keyword <- function(enr_df, kw) {
  if (is.null(enr_df) || nrow(enr_df) == 0) return(0)
  
  if (!"term_name" %in% names(enr_df)) return(0)
  
  term_vec <- as.character(enr_df$term_name)
  term_vec <- term_vec[!is.na(term_vec)]
  
  matches <- grepl(kw, term_vec, ignore.case = TRUE, fixed = FALSE)
  
  sum(matches)
}

##Enrichment Result
out_dir="C:/project1/enrichment/data"
enr_files <- list.files(out_dir, pattern = "^enrichment_.*\\.xlsx$", full.names = TRUE)

rows <- list()
for (f in enr_files) {
  pid <- gsub("enrichment_|\\.xlsx", "", basename(f))
  
  all_df   <- tryCatch(read_excel(f, sheet = "enrichment_allGenes"),   error = function(e) NULL)
  ## If Up/Down Enrichment Results are not available
  #up_df   <- tryCatch(read_excel(f, sheet = "enrichment_up"),   error = function(e) NULL)
  #down_df <- tryCatch(read_excel(f, sheet = "enrichment_down"), error = function(e) NULL)
  
  for (direction in c("all", "up", "down")) {
    df <- switch(direction,
                 #up   = up_df,
                 #down = down_df,
                 all  = all_df)
    
    for (kw in keywords) {
      rows[[length(rows) + 1]] <- data.frame(
        Patient = pid, 
        Direction = direction, 
        Keyword = kw,
        Count = count_one_keyword(df, kw),
        stringsAsFactors = FALSE
      )
    }
  }
}

long_df <- do.call(rbind, rows)

wide_df <- long_df %>%
  pivot_wider(
    names_from = Keyword, 
    values_from = Count,
    values_fill = 0)

"C:/project1/enrichment/output"
write.csv(wide_df, file.path(out_dir, "enrichment_keyword_term_counts.csv"), row.names = FALSE)

patients <- c("Pat1","Pat2","Pat3")

df <- read.csv(file.path(out_dir, "enrichment_keyword_term_counts.csv"), header = TRUE)

df_long <- df %>%
  pivot_longer(
    cols = -c(Patient, Direction), 
    names_to = "Category",
    values_to = "Count")

df_cleaned <- df_long %>%
  mutate(Category = toTitleCase(gsub("\\.", " ", trimws(as.character(Category)))),
         Patient_ID = Patient,
         Level = case_when(
           grepl("", Patient) ~ "Gene",)
  ) %>%
  filter(!is.na(Category) & Category != "" & Category %in% c("Apoptosis", "Immune Response", "Stress",
                                                             "Protein","DNA","Autophagy")) %>%
  complete(Patient_ID, Level, Category, fill = list(Count = 0))

## If Up/Down Enrichment Results are not available
df_plot <- df_cleaned |>
  dplyr::group_by(Patient_ID, Category, Level) |>
  dplyr::summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop")

ggplot(df_plot, aes(x = Level, y = Category, fill = Count)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = Count, color = ifelse(is.na(Count), FALSE, Count > 10)),
    fontface = "bold",
    size = 5.5
  ) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
  facet_wrap(~ Patient_ID, ncol = 3) +
  scale_fill_gradientn(
    colors = c("#f7f7f7", "#92c9df", "#4393c3", "#2166ac", "#053061"),
    name = ""
  ) +
  theme_minimal(base_size = 12) +
  labs(title = "", subtitle = "", x = "", y = "") +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 13, color = "gray40", margin = margin(b = 10)),
    strip.background = element_rect(fill = "#e5e7eb", color = NA),
    strip.text = element_text(face = "bold", size = 13),
    axis.text = element_text(face = "bold", color = "gray20"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "enrTerms_PLevel.png"), width = 12, height = 7, dpi = 300)