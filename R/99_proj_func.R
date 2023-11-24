volcano_plot <- function(data, condition_test){
  plt <- data |>
    ggplot(aes(x = estimate,
                 y = -log10(p.value),
                 color = dif_exp)) +
    geom_point(alpha = 0.4) +
      # geom_text_repel(aes(label = sign),
      #                 max.overlaps = Inf,
      #                 size = 2,
      #                 box.padding = 0.2) +
      # geom_hline(yintercept = 0) +
    labs(title = paste0("Diferentially expressed proteins in the test: ", condition_test, " vs. Non-", condition_test),
         subtitle = "Proteins highlighted either red or turquoise were significant after multiple test correction",
         x = "Estimates", 
         y = "-log10(p)") +
    scale_color_manual(values = c("blue", "grey", "red")) +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) 
  return(plt)
  
}

DEA_proteins <- function(data_in, condition_test){
  
  # Needed for model fitting
  col_name <- deparse(substitute(condition_test))
  
  # First select only the columns corresponding to protein expression data 
  # and the column that corresponds to the condition that is to be tested
  data_long <- data_in |>
    dplyr::select(matches("^NP"),
                  matches("^XP"),
                  matches("^YP"),
                  {{ condition_test }}) |>
    
    pivot_longer(cols = -{{ condition_test }},
                 names_to = "Protein",
                 values_to = "log2_iTRAQ")
  
  # Then create a nested data structure, where data for each protein is stored 
  # as a nested data frame within a new column
  data_long_nested <- data_long |>
    group_by(Protein) |>
    nest() |>
    ungroup()
  
  
  # Add a new variable named model_object computed for every protein
  data_w_model <- data_long_nested |>
    group_by(Protein) |>
    mutate(model_object = map(.x = data,
                              .f = ~lm(formula = paste0("log2_iTRAQ ~", col_name) ,
                                       data = .x)))

  # Tidy the models using functions from the package broom by creating a
  # new variable model_object_tidy
  data_w_model <- data_w_model |>
    mutate(model_object_tidy = map(.x = model_object,
                                   .f = ~tidy(.x,
                                              conf.int = TRUE,
                                              conf.level = 0.95)))
  
  # Unnest model_object_tidy and filter to get the slope term corresponding to
  # to the condition being tested. Select the variables and ungroup.
  estimates <- data_w_model |>
    unnest(model_object_tidy) |>
    filter(term == col_name) |>
    ungroup() |>
    dplyr::select(Protein, p.value, estimate, conf.low, conf.high) |>


    # Finally calculate the adjusted p-value and save as new variable q.value, and
    # based on that value create a new variable called is_significant
    mutate(q.value = p.adjust(p.value)) |>
    mutate(dif_exp = case_when(q.value <= 0.05 & estimate > 0 ~ 'Up',
                               q.value <= 0.05 & estimate < 0 ~ 'Down',
                               q.value > 0.05 ~ 'NS'))

  # Call volcano_plot function
  plt_volcano <- volcano_plot(estimates, col_name)
  return(list(estimates=estimates, plt_volcano=plt_volcano))
}
