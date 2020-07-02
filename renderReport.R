#' ---
#' title: "Limma - multi-group (Global F-test / ANOVA): script report"
#' ---

#'## Dataset/experiment code: `r exper_code`
#'
#'## Grouping information summary: 
#+ echo = FALSE, message = FALSE
kable(group_by(annot_dat, Group = Group) %>% summarise(N = n()))%>%
          kable_styling(bootstrap_options = "striped", full_width = F, position = "left")
#'
#'
#'## Number of proteins affected by the grouping variable: `r n_significant`
#'
#'## Boxplots of the top `r n_top_hits` hits after testing 
#'
#+ fig.width=11.7, fig.height= 7.1, echo = FALSE, message = FALSE, warning = FALSE
print(boxplots)
print(sessionInfo())