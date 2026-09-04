# Penetrance analysis====
#' Penetrance Analysis with Stratification and Plot Generation
#'
#' This function performs a penetrance analysis using survival curves, optionally stratified by a secondary variable.
#' It generates survival plots and calculates pairwise comparisons with optional p-value correction.
#'
#' @param dataframe Dataframe containing the data.
#' @param time_var Name of the time variable.
#' @param event_var Name of the event variable (1 if event occurred, 0 otherwise).
#' @param group_var Name of the grouping variable for survival curves.
#' @param palette Vector of colors for the groups.
#' @param strata_var Optional name of the stratification variable. If provided, the analysis will be performed per stratum.
#' @param plot_title Optional custom title for the plot. If not provided, it will be generated dynamically.
#' @param output_file Optional path to save the plot as a PNG file.
#' @param pairwise_correction Method for p-value correction in pairwise comparisons (passed to p.adjust.method), or NULL for no correction.
#' @param xlim_vals Optional x-axis limits for the plot.
#' @param ylim_vals Optional y-axis limits for the plot.
#'
#' @return A list containing survival plots and pairwise comparisons for each stratum.
#' @export
#' @import survival
#' @importFrom ggplot2 ggsave
#' @importFrom dplyr filter group_split
#' @importFrom purrr map2
#' @importFrom ggpubr ggarrange
#' @importFrom rlang .data sym

#' Penetrance Analysis with Stratification and Plot Generation
#'
#' This function performs a penetrance analysis using survival curves, optionally stratified by a secondary variable.
#' It generates survival plots and calculates pairwise comparisons with optional p-value correction.
#'
#' @param dataframe Dataframe containing the data.
#' @param time_var Name of the time variable.
#' @param event_var Name of the event variable (1 if event occurred, 0 otherwise).
#' @param group_var Name of the grouping variable for survival curves.
#' @param palette Vector of colors for the groups.
#' @param strata_var Optional name of the stratification variable. If provided, the analysis will be performed per stratum.
#' @param plot_title Optional custom title for the plot. If not provided, it will be generated dynamically.
#' @param output_file Optional path to save the plot as a PNG file.
#' @param pairwise_correction Method for p-value correction in pairwise comparisons (passed to p.adjust.method), or NULL for no correction.
#' @param xlim_vals Optional x-axis limits for the plot.
#' @param ylim_vals Optional y-axis limits for the plot.
#' @param break_by Optional spacing between x-axis breaks (in the same units as `time_var`).
#'   If `NULL`, a default value of **10** is used.
#'
#' @return A list containing survival plots and pairwise comparisons for each stratum.
#' @export
penetrance_analysis <- function(dataframe, time_var, event_var, group_var = NULL, palette = NULL,
                                strata_var = NULL, plot_title = NULL, output_file = NULL,
                                pairwise_correction = NULL, xlim_vals = NULL, ylim_vals = NULL,
                                break_by = NULL) {

  # Validar que las variables existan
  if (!(time_var %in% colnames(dataframe))) stop(paste("Variable", time_var, "not found in the dataset."))
  if (!(event_var %in% colnames(dataframe))) stop(paste("Variable", event_var, "not found in the dataset."))
  if (!is.null(group_var) && !(group_var %in% colnames(dataframe))) stop(paste("Variable", group_var, "not found in the dataset."))
  if (!is.null(strata_var) && !(strata_var %in% colnames(dataframe))) stop(paste("Stratification variable", strata_var, "not found in the dataset."))

  # Limpiar datos de NA
  dataframe <- dataframe %>%
    dplyr::filter(!is.na(.data[[time_var]]) & !is.na(.data[[event_var]]))

  if (!is.null(group_var)) {
    dataframe <- dataframe %>% dplyr::filter(!is.na(.data[[group_var]]))
  }

  # Estratificación si se pide
  if (!is.null(strata_var)) {
    dataframe <- dataframe %>% dplyr::filter(!is.na(.data[[strata_var]]))
    stratified_results <- dataframe %>%
      dplyr::group_split(!!rlang::sym(strata_var)) %>%
      purrr::map2(
        .x = .,
        .y = purrr::map_chr(., ~ as.character(unique(.x[[strata_var]]))),
        ~ {
          cat("\n### Análisis Estratificado para", .y, "###\n")
          penetrance_analysis(
            dataframe = .x,
            time_var = time_var,
            event_var = event_var,
            group_var = group_var,
            palette = palette,
            plot_title = paste("Penetrance -", .y),
            output_file = NULL,
            pairwise_correction = pairwise_correction,
            xlim_vals = xlim_vals,
            ylim_vals = ylim_vals,
            break_by = break_by   # <- propagar el argumento
          )
        }
      )
    return(invisible(stratified_results))
  }

  # Preparar datos para survfit
  temp_data <- dataframe %>%
    dplyr::select(Time = dplyr::all_of(time_var), Event = dplyr::all_of(event_var)) %>%
    dplyr::mutate(Group = if (!is.null(group_var)) dataframe[[group_var]] else factor("All Data"))

  if (nrow(temp_data) == 0) stop("No valid data available after removing missing values.")

  cat("\n### Análisis de:", group_var, "###\n")

  # Estadísticas descriptivas
  stats_table <- temp_data %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      Mediana = stats::median(Time, na.rm = TRUE),
      IQR = stats::IQR(Time, na.rm = TRUE),
      `% Eventos` = round(100 * mean(Event == 1, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>% as.data.frame()
  cat("\n### Estadísticas Descriptivas por Grupo ###\n")
  print(stats_table, row.names = FALSE)

  # Ajuste de supervivencia
  fit_penetrance <- survival::survfit(survival::Surv(Time, Event == 1) ~ Group, data = temp_data)

  # P-valor global
  logrank_test <- survival::survdiff(survival::Surv(Time, Event == 1) ~ Group, data = temp_data)
  p_global <- 1 - stats::pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
  cat("\n### P-valor Global (Log-Rank Test):", format.pval(p_global, digits = 4, eps = 0.0001), "###\n")

  # Comparaciones pareadas (opcional)
  if (!is.null(group_var)) {
    adj_method <- if (is.null(pairwise_correction)) "none" else pairwise_correction
    cat(sprintf("\n### Comparaciones Pareadas (método de ajuste: %s) ###\n", toupper(adj_method)))
    penetrance_pairwise <- survminer::pairwise_survdiff(
      survival::Surv(Time, Event == 1) ~ Group,
      data = temp_data,
      p.adjust.method = adj_method
    )
    print(penetrance_pairwise)
  } else {
    penetrance_pairwise <- NULL
  }

  # Límites de ejes
  if (is.null(xlim_vals)) xlim_vals <- range(dataframe[[time_var]], na.rm = TRUE)
  if (is.null(ylim_vals)) ylim_vals <- c(0, 1)

  # Valor por defecto del espaciado de cortes (10) si no se especifica
  if (is.null(break_by)) break_by <- 10

  # Coordenadas p-valor
  pval_x_pos <- min(xlim_vals) + 1
  pval_y_pos <- max(ylim_vals) - 0.05

  # Gráfico de penetrancia
  penetrance_plot <- survminer::ggsurvplot(
    fit_penetrance,
    data = temp_data,
    risk.table = TRUE,
    palette = if (!is.null(palette)) palette else "black",
    pval = TRUE,
    pval.coord = c(pval_x_pos, pval_y_pos),
    conf.int = FALSE,
    xlim = xlim_vals,
    ylim = ylim_vals,
    xlab = "Age at Diagnosis (years)",
    ylab = "Cumulative Probability",
    break.time.by = break_by,     # <- aquí usamos el nuevo parámetro
    fun = "event",
    title = plot_title,
    ggtheme = ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title.x = ggplot2::element_text(size = 12),
        axis.title.y = ggplot2::element_text(size = 12),
        legend.text = ggplot2::element_text(size = 12)
      ),
    risk.table.y.text.col = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    ncensor.plot = FALSE,
    ncensor.plot.height = 0.25,
    conf.int.style = "step",
    surv.median.line = "hv",
    legend.title = if (is.null(group_var)) NULL else "",
    legend.labs = if (is.null(group_var)) NULL else levels(temp_data$Group)
  )

  # Combinar gráfico y tabla
  penetrance_combined <- ggpubr::ggarrange(
    penetrance_plot$plot,
    penetrance_plot$table,
    ncol = 1,
    heights = c(3, 1)
  )
  print(penetrance_combined)

  # Guardar si necesario
  if (!is.null(output_file)) {
    ggplot2::ggsave(filename = output_file,
                    plot = penetrance_combined,
                    dpi = 300, height = 7, width = 7, units = "in")
  }

  invisible(list(
    fit = fit_penetrance,
    p_global = p_global,
    pairwise_comparisons = penetrance_pairwise,
    statistics = stats_table,
    plot = penetrance_combined
  ))
}

#' Survival Analysis with Stratification and Customized Plots
#'
#' @param dataframe Dataframe containing the data.
#' @param time_var Name of the time variable.
#' @param event_var Name of the event variable (1 if event occurred, 0 otherwise).
#' @param group_var Name of the grouping variable for survival curves.
#' @param palette Vector of colors for the groups.
#' @param strata_var Optional stratification variable (analyze per stratum).
#' @param plot_title Optional custom title.
#' @param output_file Optional path to save the combined plot (PNG).
#' @param pairwise_correction p.adjust.method for pairwise tests (or NULL = none).
#' @param xlim_vals Optional x-axis limits.
#' @param ylim_vals Optional y-axis limits.
#' @param break_by Spacing between x-axis breaks (same units as time_var).
#'   Default = 10 if NULL.
#'
#' @return list(fit, pairwise_comparisons, plot)
#' @export
survival_analysis <- function(dataframe, time_var, event_var, group_var = NULL, palette,
                              strata_var = NULL, plot_title = NULL, output_file = NULL,
                              pairwise_correction = NULL, xlim_vals = NULL, ylim_vals = NULL,
                              break_by = NULL) {

  # cheques básicos
  if (!(time_var %in% colnames(dataframe))) stop(paste("Variable", time_var, "not found in the dataset."))
  if (!(event_var %in% colnames(dataframe))) stop(paste("Variable", event_var, "not found in the dataset."))
  if (!is.null(group_var) && !(group_var %in% colnames(dataframe))) stop(paste("Variable", group_var, "not found in the dataset."))
  if (!is.null(strata_var) && !(strata_var %in% colnames(dataframe))) stop(paste("Stratification variable", strata_var, "not found in the dataset."))

  # estratificación (si aplica): propagar break_by
  if (!is.null(strata_var)) {
    stratified_results <- dataframe %>%
      dplyr::group_split(!!rlang::sym(strata_var)) %>%
      purrr::map2(
        .x = .,
        .y = purrr::map_chr(., ~ as.character(unique(.x[[strata_var]]))),
        ~ survival_analysis(
          dataframe = .x,
          time_var = time_var,
          event_var = event_var,
          group_var = group_var,
          palette = palette,
          plot_title = paste("Survival Curve -", .y),
          output_file = NULL,
          pairwise_correction = pairwise_correction,
          xlim_vals = xlim_vals,
          ylim_vals = ylim_vals,
          break_by = break_by
        )
      )
    return(stratified_results)
  }

  # datos limpios
  temp_data <- dataframe %>%
    dplyr::select(Time = dplyr::all_of(time_var), Event = dplyr::all_of(event_var)) %>%
    dplyr::mutate(Group = if (!is.null(group_var)) dataframe[[group_var]] else "All Data") %>%
    dplyr::filter(!is.na(Time) & !is.na(Event))
  if (nrow(temp_data) == 0) stop("No valid data available after removing missing values.")

  if (is.null(plot_title)) plot_title <- "Survival Analysis"
  temp_data$Group <- factor(temp_data$Group)

  # modelo KM
  fit_survival <- survival::survfit(survival::Surv(Time, Event == 1) ~ Group, data = temp_data)

  # comparaciones por pares (opcional)
  if (!is.null(group_var)) {
    adj_method <- if (is.null(pairwise_correction)) "none" else pairwise_correction
    cat(sprintf("\n### Pairwise Comparisons (adjustment: %s) ###\n", toupper(adj_method)))
    survival_pairwise <- survminer::pairwise_survdiff(
      survival::Surv(Time, Event == 1) ~ Group,
      data = temp_data,
      p.adjust.method = adj_method
    )
    print(survival_pairwise)
  } else {
    survival_pairwise <- NULL
  }

  # límites y cortes del eje
  if (is.null(xlim_vals)) xlim_vals <- range(dataframe[[time_var]], na.rm = TRUE)
  if (is.null(ylim_vals)) ylim_vals <- c(0, 1)
  if (is.null(break_by)) break_by <- 10  # <<< por defecto 10

  # posición del p-valor del plot
  pval_x_pos <- min(xlim_vals) + 1
  pval_y_pos <- max(ylim_vals) - 0.05

  # gráfico KM
  survival_plot <- survminer::ggsurvplot(
    fit_survival,
    data = temp_data,
    risk.table = TRUE,
    palette = palette,
    pval = TRUE,
    print = FALSE,
    pval.coord = c(pval_x_pos, pval_y_pos),
    conf.int = FALSE,
    xlim = xlim_vals,
    ylim = ylim_vals,
    xlab = "Time",
    ylab = "Survival Probability",
    break.time.by = break_by,  # <<< nuevo parámetro
    title = plot_title,
    ggtheme = ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title.x = ggplot2::element_text(size = 12),
        axis.title.y = ggplot2::element_text(size = 12),
        legend.text = ggplot2::element_text(size = 12)
      ),
    risk.table.y.text.col = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    ncensor.plot = FALSE,
    ncensor.plot.height = 0.25,
    conf.int.style = "step",
    surv.median.line = "hv",
    legend.title = "",
    legend.labs = if (!is.null(group_var)) levels(temp_data$Group) else "All Data"
  )

  # combinar plot + risk table
  survival_combined <- ggpubr::ggarrange(
    survival_plot$plot,
    survival_plot$table,
    ncol = 1,
    heights = c(3, 1)
  )

  # imprimir y/o guardar
  print(fit_survival)
  if (!is.null(survival_pairwise)) print(survival_pairwise)
  print(survival_combined)

  if (!is.null(output_file)) {
    ggplot2::ggsave(
      filename = output_file,
      plot = survival_combined,
      dpi = 300, height = 7, width = 7, units = "in"
    )
  }

  list(
    fit = fit_survival,
    pairwise_comparisons = survival_pairwise,
    plot = survival_combined
  )
}


# Violin Plots ===========================================================
#' Violin Plot for Distribution Visualization with Optional Stratification
#'
#' This function generates a violin plot to visualize the distribution of a numeric variable
#' grouped by a categorical variable, with optional stratified analysis.
#'
#' @param dataframe Dataframe containing the data.
#' @param variable Name of the numeric variable to analyze.
#' @param group_var Name of the categorical variable used for grouping.
#' @param palette Vector of colors for each group.
#' @param strata_var Optional name of the categorical variable for stratification. If provided,
#'                   violin plots will be generated separately for each level of the stratification variable.
#' @param plot_title Optional title for the plot.
#' @param output_file Optional path to save the plot as a PNG file.
#'
#' @return If `strata_var` is provided, returns a list of violin plots for each stratum.
#'         Otherwise, returns a single violin plot.
#' @export
#' @importFrom dplyr filter group_split
#' @importFrom ggplot2 ggplot aes geom_violin geom_jitter stat_summary labs theme_bw theme element_text scale_fill_manual scale_color_manual
#' @importFrom ggplot2 ggsave
#' @importFrom purrr map2 walk
#' @importFrom rlang .data sym

violin_plot <- function(dataframe, variable, group_var, palette, strata_var = NULL, plot_title = NULL, output_file = NULL) {
  # Validar que las variables existan en el dataframe
  if (!(variable %in% colnames(dataframe))) stop(paste("Variable", variable, "not found in the dataset."))
  if (!(group_var %in% colnames(dataframe))) stop(paste("Grouping variable", group_var, "not found in the dataset."))
  if (!is.null(strata_var) && !(strata_var %in% colnames(dataframe))) {
    stop(paste("Stratification variable", strata_var, "not found in the dataset."))
  }

  # Función auxiliar para crear el gráfico
  create_violin_plot <- function(data, title) {
    ggplot(data) +
      geom_violin(
        aes(x = .data[[group_var]], y = .data[[variable]], fill = .data[[group_var]]),
        alpha = 0.5, color = "black", show.legend = FALSE, scale = "width"
      ) +
      geom_jitter(
        aes(x = .data[[group_var]], y = .data[[variable]], color = .data[[group_var]]),
        width = 0.2, alpha = 0.2, show.legend = FALSE
      ) +
      stat_summary(
        aes(x = .data[[group_var]], y = .data[[variable]]),
        fun = mean, geom = "crossbar", width = 0.3, color = "black", fatten = 0.5
      ) +
      stat_summary(
        aes(x = .data[[group_var]], y = .data[[variable]]),
        fun = mean, geom = "point", shape = 21, size = 3, fill = "black", color = "black"
      ) +
      labs(
        title = title,
        x = "Group",
        y = variable
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette) +
      scale_y_continuous(limits = c(5, 45))  # Fijar la escala del eje Y
  }

  # Si se define la estratificación, analizar por estratos
  if (!is.null(strata_var)) {
    dataframe <- dataframe %>%
      filter(!is.na(.data[[variable]]) & !is.na(.data[[group_var]]) & !is.na(.data[[strata_var]]))

    stratified_results <- dataframe %>%
      group_split(!!sym(strata_var)) %>%
      map2(
        .x = .,
        .y = map_chr(., ~ as.character(unique(.x[[strata_var]]))),
        ~ {
          cat("\n### Stratified Analysis for", .y, "###\n")
          plot_title_stratum <- paste("Violin Plot of", variable, "by", group_var, "-", .y)
          # Crear y retornar el gráfico
          create_violin_plot(.x, plot_title_stratum)
        }
      )

    cat("\nPrinting violin plots for each stratum...\n")
    # Imprimir gráficos sin retornar objetos gráficos adicionales
    walk(stratified_results, ~ print(.x))

    # Retornar solo los objetos gráficos sin impresión automática
    return(invisible(stratified_results))
  }

  # Si no hay estratificación, crear el gráfico directamente
  clean_data <- dataframe %>%
    filter(!is.na(.data[[variable]]) & !is.na(.data[[group_var]]))

  if (is.null(plot_title)) {
    plot_title <- paste("Violin Plot of", variable, "by", group_var)
  }

  final_plot <- create_violin_plot(clean_data, plot_title)

  # Guardar el gráfico si se proporciona un archivo de salida
  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot = final_plot,
      dpi = 300,
      height = 7,
      width = 7,
      units = "in"
    )
  }

  # Retornar el gráfico sin impresión automática
  return(invisible(final_plot))
}







# Continuous Variables Analysis ===========================================================
#' Means Comparison with Stratification and Visualization
#'
#' This function performs mean comparisons (ANOVA, Tukey, Kruskal-Wallis, Dunn tests) and creates boxplots and Tukey/Dunn plots.
#' It supports stratification for more granular analysis.
#'
#' @param dataframe Dataframe containing the data.
#' @param variable Name of the numeric variable to analyze.
#' @param group_var Name of the categorical variable used for grouping.
#' @param palette Vector of colors for the groups.
#' @param strata_var Optional stratification variable. If provided, the analysis will be done per stratum.
#' @param plot_title Optional title for the plots.
#'
#' @return A list containing analysis summaries and plots.
#' @export
#' @importFrom ggplot2 ggplot aes geom_boxplot stat_summary labs theme_bw theme element_text scale_fill_manual ggsave
#' @importFrom dplyr filter group_split summarise group_by mutate
#' @importFrom purrr map2
#' @importFrom stats aov TukeyHSD kruskal.test p.adjust
#' @importFrom FSA dunnTest
#' @importFrom utils combn
#' @importFrom rlang .data sym

means_analysis <- function(dataframe, variable, group_var, palette, strata_var = NULL, plot_title = NULL, output_file = NULL) {
  # Validate that variables exist in the dataset
  if (!(variable %in% colnames(dataframe))) stop(paste("Variable", variable, "not found in the dataset."))
  if (!(group_var %in% colnames(dataframe))) stop(paste("Grouping variable", group_var, "not found in the dataset."))
  if (!is.null(strata_var) && !(strata_var %in% colnames(dataframe))) {
    stop(paste("Stratification variable", strata_var, "not found in the dataset."))
  }

  # If strata_var is provided, perform stratified analysis
  if (!is.null(strata_var)) {
    stratified_results <- dataframe %>%
      group_split(!!sym(strata_var)) %>%
      map2(
        .x = .,
        .y = map_chr(., ~ as.character(unique(.x[[strata_var]]))),
        ~ {
          cat("\n### Stratified Analysis for", .y, "###\n")
          means_analysis(
            dataframe = .x,
            variable = variable,
            group_var = group_var,
            palette = palette,
            plot_title = paste("Boxplot of", variable, "by", group_var, "-", .y),
            output_file = NULL  # Avoid saving automatically inside iterations
          )
        }
      )
    return(stratified_results)  # Return results for each stratum
  }

  # Clean data by filtering NA values
  clean_data <- dataframe %>%
    filter(!is.na(.data[[variable]]) & !is.na(.data[[group_var]]))

  # Create the plot title if not provided
  if (is.null(plot_title)) {
    plot_title <- paste("Boxplot of", variable, "by", group_var)
  }

  # Overall statistical summary
  summary_total <- clean_data %>%
    summarise(
      n = n(),Mean = round(mean(.data[[variable]], na.rm = TRUE), 2),
      Median = median(.data[[variable]], na.rm = TRUE),
      SD = round(sd(.data[[variable]], na.rm = TRUE), 2),
      Range = paste0(min(.data[[variable]], na.rm = TRUE), " - ", max(.data[[variable]], na.rm = TRUE)),
      Q1 = quantile(.data[[variable]], 0.25, na.rm = TRUE),
      Q3 = quantile(.data[[variable]], 0.75, na.rm = TRUE)
    )

  # Statistical summary by group
  summary_by_group <- clean_data %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      n = n(),
      Mean = round(mean(.data[[variable]], na.rm = TRUE), 2),
      Median = median(.data[[variable]], na.rm = TRUE),
      SD = round(sd(.data[[variable]], na.rm = TRUE), 2),
      Range = paste0(min(.data[[variable]], na.rm = TRUE), " - ", max(.data[[variable]], na.rm = TRUE)),
      Q1 = quantile(.data[[variable]], 0.25, na.rm = TRUE),
      Q3 = quantile(.data[[variable]], 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  # Print summaries
  print("Overall Statistical Summary:")
  print(summary_total)
  print("Statistical Summary by Group:")
  print(summary_by_group)

  # Perform ANOVA test
  anova_result <- aov(dataframe[[variable]] ~ dataframe[[group_var]], data = dataframe)
  anova_summary <- summary(anova_result)
  print("ANOVA Result:")
  print(anova_summary)

  tukey_result <- NULL  # Inicializamos como NULL por si ANOVA no es significativo
  if (anova_summary[[1]]$`Pr(>F)`[1] < 0.05) {
    tukey_result <- TukeyHSD(anova_result)
  }

  # Perform Kruskal-Wallis test
  kruskal_result <- kruskal.test(dataframe[[variable]] ~ dataframe[[group_var]], data = dataframe)
  print("Kruskal-Wallis Test Result:")
  print(kruskal_result)

  # Perform Dunn's test safely
  dunn_result <- NULL
  if (requireNamespace("FSA", quietly = TRUE)) {
    dunn_result <- FSA::dunnTest(dataframe[[variable]] ~ dataframe[[group_var]], data = dataframe, method = "bonferroni")
    print("Dunn Test Result:")
    print(dunn_result)
  } else {
    print("Warning: FSA package not installed. Skipping Dunn test.")
  }

  # Create boxplot
  boxplot <- ggplot(clean_data, aes(x = .data[[group_var]], y = .data[[variable]], fill = .data[[group_var]])) +
    geom_boxplot(alpha = 0.5, color = "black", outlier.shape = 16, outlier.size = 2) +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 3, fill = "black", color = "black") +
    labs(
      title = plot_title,
      x = group_var,
      y = variable
    ) +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Display the plot
  print(boxplot)

  # Save the plot if output_file is specified
  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot = boxplot,
      dpi = 300,
      height = 7,
      width = 7,
      units = "in"
    )
  }

  # Return results as a list
  return(list(
    summary_total = summary_total,
    summary_by_group = summary_by_group,
    anova_result = anova_result,
    tukey_result = tukey_result,
    kruskal_result = kruskal_result,
    dunn_result = dunn_result,
    boxplot = boxplot
  ))
}





# Contingency Tables ===========================================================
#' Contingency Table Analysis and Proportion Plot with Stratification
#'
#' This function performs a contingency table analysis, calculates statistical tests
#' (Chi-square and Fisher tests), and generates proportion plots for a specified category
#' within a variable across different groups. Optionally, it can stratify the analysis
#' by another variable.
#'
#' @param dataframe A data frame containing the data.
#' @param variable A string specifying the name of the categorical variable to analyze.
#' @param group_var A string specifying the name of the categorical variable used for grouping.
#' @param category_filter A value representing the specific category of interest within the variable.
#' For example, for a binary variable like `sex`, this could be `"Male"`, or for numeric binary variables, it could be `1`.
#' @param palette A vector of colors for the groups, used in the plot.
#' @param plot_title Optional. A string specifying the title for the plot. If NULL, a default title is generated.
#' @param strata_var Optional. A string specifying the name of the variable to stratify the analysis by. If provided,
#' the analysis will be repeated for each level of this variable.
#'
#' @return A list containing the following:
#' \item{chi2_test}{Result of the Chi-square test.}
#' \item{fisher_test}{Result of the Fisher's exact test (if applicable).}
#' \item{contingency_table}{The contingency table with frequencies.}
#' \item{proportion_plot}{A ggplot object representing the bar plot of proportions.}
#' If stratification is provided, returns a list of results for each stratum.
#' @export
#' @importFrom dplyr filter group_by summarise mutate
#' @importFrom stats chisq.test fisher.test p.adjust
#' @importFrom utils combn
#' @importFrom scales percent_format
#' @importFrom ggplot2 ggplot aes geom_bar geom_text labs theme_bw theme element_text scale_fill_manual position_dodge scale_y_continuous
#' @importFrom rlang .data

contingency_analysis <- function(dataframe, variable, group_var, category_filter, palette, plot_title = NULL, pairwise_comparisons = FALSE) {

  # Check if the variables exist in the dataframe
  if (!(variable %in% colnames(dataframe))) stop(paste("Variable", variable, "not found in the dataset."))
  if (!(group_var %in% colnames(dataframe))) stop(paste("Grouping variable", group_var, "not found in the dataset."))

  # Remove NA values from the main variables
  clean_data <- dataframe %>%
    dplyr::filter(!is.na(.data[[variable]]) & !is.na(.data[[group_var]]))

  # Ensure variables are factors for comparison purposes
  clean_data[[variable]] <- as.factor(clean_data[[variable]])
  clean_data[[group_var]] <- as.factor(clean_data[[group_var]])

  # Create contingency table
  contingency_table <- table(clean_data[[group_var]], clean_data[[variable]])
  print("Contingency table (frequencies):")
  print(contingency_table)

  # Chi-square test
  chi2_test <- chisq.test(contingency_table)
  print("Chi² test result:")
  print(chi2_test)

  # Fix for pairwise comparisons
  pairwise_result <- NULL
  if (pairwise_comparisons) {
    # Verificar que la categoría existe y extraer los conteos correctamente
    category_index <- which(colnames(contingency_table) == as.character(category_filter))

    if (length(category_index) == 0) {
      stop("Error: category_filter does not match any column in the contingency table.")
    }

    success_counts <- contingency_table[, category_index]
    group_sizes <- rowSums(contingency_table)

    # Verificar que haya suficientes observaciones para realizar comparaciones
    if (all(success_counts > 0) && length(success_counts) > 1) {
      pairwise_result <- pairwise.prop.test(success_counts, group_sizes, p.adjust.method = "bonferroni")
      print("Pairwise comparisons of proportions (adjusted p-values):")
      print(pairwise_result)
    } else {
      print("Warning: Insufficient counts or only one group available for pairwise comparisons.")
    }
  }

  # Proportion data for the plot
  proportions_data <- clean_data %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(
      proportion = mean(.data[[variable]] == category_filter, na.rm = TRUE),
      n = n()
    ) %>%
    dplyr::mutate(
      label = paste0(round(proportion * 100, 1), "% (n=", n, ")")
    )

  # Set the default plot title if not provided
  if (is.null(plot_title)) {
    plot_title <- paste("Proportion of", category_filter, "in", variable, "by", group_var)
  }

  # Create the plot
  barplot <- ggplot(proportions_data, aes(x = .data[[group_var]], y = proportion, fill = .data[[group_var]])) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8, color = "black") +
    geom_text(aes(label = label), position = position_dodge(width = 0.9), vjust = -0.5, size = 3.5) +
    labs(
      title = plot_title,
      x = group_var,
      y = "Proportion (% of the group)",
      fill = "Group"
    ) +
    scale_y_continuous(limits = c(0, 0.35), labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Display the plot
  print(barplot)

  # Return results
  return(list(
    chi2_test = chi2_test,
    pairwise_result = pairwise_result,
    contingency_table = contingency_table,
    proportion_plot = barplot
  ))
}


scatter_plot <- function(dataframe, x, y, group_var, palette, strata_var = NULL, plot_title = NULL, output_file = NULL) {


  # Validar que las variables existan en el dataframe
  if (!(x %in% colnames(dataframe))) stop(paste("Variable", x, "not found in the dataset."))
  if (!(y %in% colnames(dataframe))) stop(paste("Variable", y, "not found in the dataset."))
  if (!(group_var %in% colnames(dataframe))) stop(paste("Grouping variable", group_var, "not found in the dataset."))
  if (!is.null(strata_var) && !(strata_var %in% colnames(dataframe))) {
    stop(paste("Stratification variable", strata_var, "not found in the dataset."))
  }

  # Función auxiliar para crear el scatter plot y calcular R2
  create_scatter_plot <- function(data, title) {
    # Ajustar el modelo de regresión lineal y calcular R2
    model <- lm(as.formula(paste(y, "~", x)), data = data)
    r2 <- summary(model)$r.squared
    cat("\nR-squared for", title, "=", round(r2, 3), "\n")

    # Crear el scatter plot con la línea de regresión
    plot <- ggplot(data) +
      geom_point(aes(x = .data[[x]], y = .data[[y]], color = .data[[group_var]]), size = 2, alpha = 0.7) +
      geom_smooth(aes(x = .data[[x]], y = .data[[y]]), method = "lm", se = FALSE, color = "black", linetype = "dashed") +
      labs(
        title = title,
        x = x,
        y = y
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      scale_color_manual(values = palette) +
      scale_y_continuous(limits = c(5, 45))  # Fijar la escala del eje Y

    # Imprimir el gráfico en la pantalla
    print(plot)

    # Retornar el gráfico
    return(plot)
  }

  # Si se define la estratificación, analizar por estratos
  if (!is.null(strata_var)) {
    dataframe <- dataframe %>%
      filter(!is.na(.data[[x]]) & !is.na(.data[[y]]) & !is.na(.data[[group_var]]) & !is.na(.data[[strata_var]]))

    stratified_results <- dataframe %>%
      group_split(!!sym(strata_var)) %>%
      map2(
        .x = .,
        .y = map_chr(., ~ as.character(unique(.x[[strata_var]]))),
        ~ {
          plot_title_stratum <- paste("Scatter Plot of", y, "vs", x, "by", group_var, "-", .y)
          create_scatter_plot(.x, plot_title_stratum)
        }
      )

    # Retornar solo los objetos gráficos sin impresión automática
    return(invisible(stratified_results))
  }

  # Si no hay estratificación, crear el gráfico directamente
  clean_data <- dataframe %>%
    filter(!is.na(.data[[x]]) & !is.na(.data[[y]]) & !is.na(.data[[group_var]]))

  if (is.null(plot_title)) {
    plot_title <- paste("Scatter Plot of", y, "vs", x, "by", group_var)
  }

  final_plot <- create_scatter_plot(clean_data, plot_title)

  # Guardar el gráfico si se proporciona un archivo de salida
  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot = final_plot,
      dpi = 300,
      height = 7,
      width = 7,
      units = "in"
    )
  }

  # Retornar el gráfico sin impresión automática
  return(invisible(final_plot))
}

# Normality Analysis ========================================================
#' Normality Analysis with Kernel Density, QQ-Plot, and Statistical Tests
#'
#' This function performs a normality analysis using visual methods (density plot, QQ-plot)
#' and statistical tests (Anderson-Darling, Shapiro-Wilk, and Kolmogorov-Smirnov).
#' Optionally, it can perform separate analyses for each group.
#'
#' @param dataframe Dataframe containing the data.
#' @param variable Name of the numeric variable to analyze.
#' @param group_var Optional grouping variable for stratified analysis.
#'
#' @return A list containing plots and test results.
#' @export
#' @importFrom ggplot2 ggplot aes geom_density stat_function labs theme_bw stat_qq stat_qq_line
#' @importFrom dplyr filter group_by summarise mutate
#' @importFrom nortest ad.test
#' @importFrom stats ks.test shapiro.test
normality_analysis <- function(dataframe, variable, group_var = NULL) {
  # Validar que las variables existan en el dataframe
  if (!(variable %in% colnames(dataframe))) stop(paste("Variable", variable, "no encontrada en el dataset."))
  if (!is.null(group_var) && !(group_var %in% colnames(dataframe))) stop(paste("Variable de grupo", group_var, "no encontrada en el dataset."))

  # Filtrar filas con valores finitos en la variable de inter?s
  clean_data <- dataframe %>% dplyr::filter(is.finite(.data[[variable]]))

  # Kernel Density vs Normal Curve
  density_plot <- ggplot(clean_data, aes(x = .data[[variable]])) +
    geom_density(fill = "lightblue", alpha = 0.5) +
    stat_function(
      fun = dnorm,
      args = list(
        mean = mean(clean_data[[variable]], na.rm = TRUE),
        sd = sd(clean_data[[variable]], na.rm = TRUE)
      ),
      color = "red", size = 1.2
    ) +
    labs(
      title = paste("Kernel Density vs Normal Curve:", variable),
      x = variable,
      y = "Density"
    ) +
    theme_bw()

  # Mostrar el gr?fico de densidad
  print(density_plot)

  # QQ-Plot
  qq_plot <- ggplot(clean_data, aes(sample = .data[[variable]])) +
    stat_qq(color = "blue") +
    stat_qq_line(color = "red", linetype = "dashed") +
    labs(
      title = paste("QQ Plot:", variable),
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_bw()

  # Mostrar el QQ-Plot
  print(qq_plot)

  # Prueba de Anderson-Darling
  ad_result <- ad.test(clean_data[[variable]])
  print("Resultado de la prueba Anderson-Darling para la variable completa:")
  print(ad_result)

  # Resultados de Anderson-Darling por grupo (si group_var no es NULL)
  if (!is.null(group_var)) {
    ad_results <- clean_data %>%
      dplyr::group_by(.data[[group_var]]) %>%
      dplyr::summarise(
        p_value = ad.test(.data[[variable]])$p.value,
        A_statistic = ad.test(.data[[variable]])$statistic,
        .groups = "drop"
      )
    print("Resultados de la prueba Anderson-Darling por grupo:")
    print(ad_results)
  } else {
    ad_results <- NULL
  }
}

#Odds Ratio================

#' Calculate odds ratio, confidence interval, p-value, etiological fraction, and optional penetrance.
#'
#' Unified OR function using Woolf approximation with optional Haldane correction (+0.5 to all
#' cells) when zeros are present. Fisher's exact test is selected when any expected cell count
#' falls below \code{fisher_threshold}; otherwise chi-square without continuity correction.
#' The etiological fraction CI is bounded so lower <= upper regardless of OR direction.
#' Penetrance is estimated as \code{baseline_risk * (case_rate / control_rate)} using
#' Wilson CIs for the bounds.
#'
#' @param ac_case    Number of cases carrying the allele/exposure.
#' @param n_case     Total number of cases.
#' @param ac_control Number of controls carrying the allele/exposure.
#' @param n_control  Total number of controls.
#' @param baseline_risk Optional population baseline risk for penetrance estimation. Default NULL.
#' @param conf.level Confidence level for all interval estimates. Default 0.95.
#' @param fisher_threshold Minimum *expected* cell count below which Fisher's exact test is used
#'   instead of chi-square. Default 5 (standard epidemiological convention).
#' @param haldane One of \code{"auto"} (apply when any observed cell is 0), \code{TRUE}
#'   (always apply), or \code{FALSE} (never apply, will produce Inf/NaN with zero cells).
#'   Default \code{"auto"}.
#'
#' @return A one-row tibble with columns:
#'   \code{ac_case}, \code{n_case}, \code{ac_control}, \code{n_control},
#'   \code{case_rate}, \code{control_rate},
#'   \code{OR}, \code{OR_lower}, \code{OR_upper},
#'   \code{p_value}, \code{p_method},
#'   \code{EF}, \code{EF_lower}, \code{EF_upper},
#'   \code{penetrance}, \code{penetrance_lower}, \code{penetrance_upper},
#'   \code{haldane_applied}, \code{invalid_input}, \code{note}, \code{baseline_risk}.
#' @export
#' @importFrom stats chisq.test fisher.test qnorm
#' @importFrom tibble tibble
#' @aliases or_jpo ORjpo
compute_or <- function(ac_case,
                       n_case,
                       ac_control,
                       n_control,
                       baseline_risk    = NULL,
                       conf.level       = 0.95,
                       fisher_threshold = 5,
                       haldane          = "auto") {

  # --- Internal: Wilson CI for a proportion ---
  wilson_ci <- function(x, n) {
    p_hat  <- x / n
    z      <- stats::qnorm(1 - (1 - conf.level) / 2)
    denom  <- 1 + z^2 / n
    center <- p_hat + z^2 / (2 * n)
    margin <- z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
    list(lower = (center - margin) / denom,
         upper = (center + margin) / denom,
         mean  = p_hat)
  }

  # Save baseline_risk before any modification (fixes ifelse bug in prior version)
  br_input <- if (is.null(baseline_risk)) NA_real_ else as.numeric(baseline_risk)

  # --- NA tibble helper ---
  na_tibble <- function(note) {
    tibble::tibble(
      ac_case          = as.numeric(ac_case),
      n_case           = as.numeric(n_case),
      ac_control       = as.numeric(ac_control),
      n_control        = as.numeric(n_control),
      case_rate        = NA_real_,
      control_rate     = NA_real_,
      OR               = NA_real_,
      OR_lower         = NA_real_,
      OR_upper         = NA_real_,
      p_value          = NA_real_,
      p_method         = NA_character_,
      EF               = NA_real_,
      EF_lower         = NA_real_,
      EF_upper         = NA_real_,
      penetrance       = NA_real_,
      penetrance_lower = NA_real_,
      penetrance_upper = NA_real_,
      haldane_applied  = NA,
      invalid_input    = TRUE,
      note             = note,
      baseline_risk    = br_input
    )
  }

  # --- Input validation ---
  missing_input  <- any(is.na(c(ac_case, n_case, ac_control, n_control)))
  invalid_counts <- any(c(n_case, n_control) <= 0,
                        ac_case    < 0, ac_control    < 0,
                        ac_case > n_case, ac_control > n_control)

  if (missing_input)  return(na_tibble("Missing input value"))
  if (invalid_counts) return(na_tibble("Invalid case/control counts"))

  # --- 2x2 cells (observed) ---
  a <- ac_case
  b <- n_case    - ac_case
  c <- ac_control
  d <- n_control - ac_control

  # --- Haldane decision ---
  has_zeros <- any(c(a, b, c, d) == 0)
  apply_haldane <- switch(
    as.character(haldane),
    "auto"  = has_zeros,
    "TRUE"  = TRUE,
    "FALSE" = FALSE,
    stop("haldane must be 'auto', TRUE, or FALSE")
  )

  delta <- if (apply_haldane) 0.5 else 0
  a_adj <- a + delta
  b_adj <- b + delta
  c_adj <- c + delta
  d_adj <- d + delta

  # --- Crude rates (from observed, not adjusted) ---
  case_rate    <- a / n_case
  control_rate <- c / n_control

  # --- OR via Woolf (on adjusted cells) ---
  or_estimate <- (a_adj / b_adj) / (c_adj / d_adj)
  log_or_se   <- sqrt(1/a_adj + 1/b_adj + 1/c_adj + 1/d_adj)
  z           <- stats::qnorm(1 - (1 - conf.level) / 2)
  or_ci_lower <- exp(log(or_estimate) - z * log_or_se)
  or_ci_upper <- exp(log(or_estimate) + z * log_or_se)

  # --- p-value: Fisher when min(expected) < fisher_threshold ---
  ct       <- matrix(c(a, b, c, d), nrow = 2)
  n_total  <- a + b + c + d
  expected <- outer(c(a + b, c + d), c(a + c, b + d)) / n_total

  if (min(expected) < fisher_threshold) {
    p_test   <- stats::fisher.test(ct)
    p_method <- "Fisher exact"
  } else {
    p_test   <- stats::chisq.test(ct, correct = FALSE)
    p_method <- "Chi-square"
  }
  p_value <- p_test$p.value

  # --- Etiological Fraction: (OR - 1) / OR, bounds enforced lower <= upper ---
  ef_from_or <- function(or) (or - 1) / or
  EF       <- ef_from_or(or_estimate)
  EF_lower <- pmin(ef_from_or(or_ci_lower), ef_from_or(or_ci_upper))
  EF_upper <- pmax(ef_from_or(or_ci_lower), ef_from_or(or_ci_upper))

  # --- Penetrance ---
  if (is.na(br_input)) {
    penetrance       <- NA_real_
    penetrance_lower <- NA_real_
    penetrance_upper <- NA_real_
  } else {
    if (br_input <= 0) stop("baseline_risk must be a positive numeric value.")

    # Wilson CI on Haldane-adjusted cells for consistency with OR
    n_case_adj    <- n_case    + 2 * delta
    n_control_adj <- n_control + 2 * delta
    case_ci    <- wilson_ci(a_adj, n_case_adj)
    control_ci <- wilson_ci(c_adj, n_control_adj)

    penetrance       <- br_input * (case_ci$mean  / control_ci$mean)
    penetrance_lower <- if (control_ci$upper > 0) br_input * (case_ci$lower / control_ci$upper) else NA_real_
    penetrance_upper <- if (control_ci$lower > 0) br_input * (case_ci$upper / control_ci$lower) else NA_real_
  }

  # --- Note field ---
  notes <- character(0)
  if (apply_haldane)           notes <- c(notes, "Haldane correction applied (+0.5 to all cells)")
  if (has_zeros & !apply_haldane) notes <- c(notes, "Zero cells present; Haldane suppressed by user")
  note_out <- if (length(notes) == 0L) NA_character_ else paste(notes, collapse = "; ")

  tibble::tibble(
    ac_case          = ac_case,
    n_case           = n_case,
    ac_control       = ac_control,
    n_control        = n_control,
    case_rate        = case_rate,
    control_rate     = control_rate,
    OR               = or_estimate,
    OR_lower         = or_ci_lower,
    OR_upper         = or_ci_upper,
    p_value          = p_value,
    p_method         = p_method,
    EF               = EF,
    EF_lower         = EF_lower,
    EF_upper         = EF_upper,
    penetrance       = penetrance,
    penetrance_lower = penetrance_lower,
    penetrance_upper = penetrance_upper,
    haldane_applied  = apply_haldane,
    invalid_input    = FALSE,
    note             = note_out,
    baseline_risk    = br_input
  )
}

# Backward-compatibility aliases
or_jpo <- compute_or
ORjpo  <- compute_or

#Genomic helpers================

#' Impute missing AN values from positional neighbors
#'
#' For variants absent from a dataset (gnomAD or HiC), AN is NA but the locus
#' was still sequenced. This function estimates AN as the median of the nearest
#' non-NA values before and after the variant by genomic position.
#'
#' Locus format expected: "chrN:pos:ref:alt" (colon-separated, position second).
#'
#' @param loci     Character vector of locus strings.
#' @param an_values Numeric vector of AN values (same length as \code{loci}); NAs are imputed.
#'
#' @return Numeric vector with NAs replaced by positional-neighbor medians.
#'   Variants with no available neighbors on either side receive \code{NA}.
#' @export
impute_an_neighbors <- function(loci, an_values) {
  positions <- as.integer(regmatches(loci, regexpr("(?<=:)\\d+", loci, perl = TRUE)))
  result    <- an_values
  na_idx    <- which(is.na(an_values))
  avail_idx <- which(!is.na(an_values))

  if (length(avail_idx) == 0 || length(na_idx) == 0) return(result)

  for (i in na_idx) {
    pos_i  <- positions[i]
    before <- avail_idx[positions[avail_idx] < pos_i]
    after  <- avail_idx[positions[avail_idx] > pos_i]

    neighbors <- c(
      if (length(before) > 0) an_values[utils::tail(before, 1)] else NA_real_,
      if (length(after)  > 0) an_values[utils::head(after,  1)] else NA_real_
    )
    result[i] <- stats::median(neighbors, na.rm = TRUE)
  }
  result
}

#' Add unified gnomAD PASS flag and conflict flag to a variant table
#'
#' Applies a source-aware PASS logic that accounts for platform biases between
#' gnomAD exomes and genomes. Adds two columns to \code{df}:
#' \itemize{
#'   \item \code{gnomad_pass_unified} — logical; TRUE if the variant passes QC
#'     in at least one platform under the rules below.
#'   \item \code{gnomad_filter_conflict} — logical; TRUE if one platform passes
#'     and the other has an explicit quality failure (not AC0).
#' }
#'
#' Rules (applied in order):
#' \enumerate{
#'   \item Both platforms pass → PASS
#'   \item Intronic/non-coding consequence, genomes pass → PASS (trust genomes)
#'   \item Exonic, exomes pass, genomes AC0 (not captured) → PASS
#'   \item Exonic, genomes pass, exomes AC0 → PASS
#'   \item Otherwise → FAIL (conservative)
#' }
#'
#' gnomAD convention: \code{NA} in filter column means PASS; \code{"AC0"} means
#' not observed in that platform (treated as absent, not as failure).
#'
#' @param df A data frame with columns: \code{gnomad_ac}, \code{gnomad_filter_genomes},
#'   \code{gnomad_filter_exomes}, \code{consequence}.
#'
#' @return \code{df} with two additional logical columns appended.
#' @export
#' @importFrom dplyr mutate case_when
gnomad_pass_unified <- function(df) {
  intronic_consequences <- c(
    "intron_variant", "intergenic_variant",
    "3_prime_UTR_variant", "5_prime_UTR_variant",
    "non_coding_transcript_exon_variant", "non_coding_transcript_variant",
    "downstream_gene_variant", "upstream_gene_variant"
  )

  df |>
    dplyr::mutate(
      .is_intronic = grepl(paste(intronic_consequences, collapse = "|"), consequence),
      .in_gnomad   = !is.na(gnomad_ac),
      .gn_pass     = .in_gnomad & is.na(gnomad_filter_genomes),
      .ex_pass     = .in_gnomad & is.na(gnomad_filter_exomes),
      .gn_ac0      = !is.na(gnomad_filter_genomes) & gnomad_filter_genomes == "AC0",
      .ex_ac0      = !is.na(gnomad_filter_exomes)  & gnomad_filter_exomes  == "AC0",
      .gn_fail     = .in_gnomad & !is.na(gnomad_filter_genomes) & gnomad_filter_genomes != "AC0",
      .ex_fail     = .in_gnomad & !is.na(gnomad_filter_exomes)  & gnomad_filter_exomes  != "AC0",

      gnomad_pass_unified = dplyr::case_when(
        !.in_gnomad                          ~ FALSE,
        .gn_pass & .ex_pass                  ~ TRUE,
        .is_intronic & .gn_pass              ~ TRUE,
        !.is_intronic & .ex_pass & .gn_ac0  ~ TRUE,
        !.is_intronic & .gn_pass & .ex_ac0  ~ TRUE,
        TRUE                                 ~ FALSE
      ),
      gnomad_filter_conflict = .in_gnomad & ((.gn_pass & .ex_fail) | (.ex_pass & .gn_fail))
    ) |>
    dplyr::select(-dplyr::starts_with("."))
}

#' Compute per-variant OR enrichment across phenotypes
#'
#' For each variant in \code{cohort} and each phenotype in \code{phenotypes},
#' computes OR (Woolf + Haldane) against three comparators: internal HiC
#' controls, gnomAD NFE, and gnomAD total. AN values are imputed from
#' positional neighbors when missing.
#'
#' Haldane override: \code{ac_case == 0} → OR = NA (no case signal;
#' Haldane correction would produce an uninformative upward pull).
#' \code{ac_ctrl == 0, ac_case > 0} → Haldane applied normally.
#'
#' AN units: HiC uses individuals; gnomAD allele counts are divided by 2
#' before the OR call to convert to individuals.
#'
#' @param cohort          Data frame with columns \code{locus} and per-phenotype
#'   count columns named \code{{pheno}_ac_case}, \code{{pheno}_nac_case},
#'   \code{{pheno}_ac_ctrl}, \code{{pheno}_nac_ctrl}.
#' @param phenotypes      Character vector of phenotype prefixes (lowercase).
#' @param gnomad_ac_nfe   Numeric vector (length = nrow(cohort)): gnomAD NFE AC,
#'   already NA-replaced with 0 for absent variants.
#' @param gnomad_an_nfe   Numeric vector: gnomAD NFE AN in alleles (will be /2).
#'   NAs imputed internally from positional neighbors.
#' @param gnomad_ac_total Numeric vector: gnomAD total AC (NA → 0).
#' @param gnomad_an_total Numeric vector: gnomAD total AN in alleles (NAs imputed).
#'
#' @return A tibble with one row per variant and columns prefixed by phenotype:
#'   \code{{p}_ac_case}, \code{{p}_an_case}, \code{{p}_an_case_imputed},
#'   \code{{p}_ac_ctrl_int}, \code{{p}_an_ctrl_int}, \code{{p}_an_ctrl_int_imputed},
#'   \code{{p}_or_int}, \code{{p}_ci_low_int}, \code{{p}_ci_high_int}, \code{{p}_p_int},
#'   and equivalent \code{_nfe} and \code{_total} blocks.
#' @export
#' @importFrom purrr map pmap_dfr
#' @importFrom dplyr bind_cols select
#' @importFrom tibble as_tibble tibble
compute_enrichment <- function(cohort,
                               phenotypes,
                               gnomad_ac_nfe,
                               gnomad_an_nfe,
                               gnomad_ac_total,
                               gnomad_an_total) {

  na_or <- tibble::tibble(OR = NA_real_, OR_lower = NA_real_,
                          OR_upper = NA_real_, p_value = NA_real_)

  safe_or <- function(ac_case, an_case, ac_ctrl, an_ctrl) {
    if (any(is.na(c(ac_case, an_case, ac_ctrl, an_ctrl)))) return(na_or)
    if (an_case == 0 || an_ctrl == 0)                       return(na_or)
    if (ac_case == 0)                                        return(na_or)
    tryCatch(
      compute_or(ac_case, an_case, ac_ctrl, an_ctrl) |>
        dplyr::select(OR, OR_lower, OR_upper, p_value),
      error = function(e) na_or
    )
  }

  # Impute and convert gnomAD ANs (alleles → individuals)
  an_nfe_imp   <- impute_an_neighbors(cohort$locus, gnomad_an_nfe)
  an_total_imp <- impute_an_neighbors(cohort$locus, gnomad_an_total)
  an_nfe_ind   <- an_nfe_imp   / 2
  an_total_ind <- an_total_imp / 2

  blocks <- purrr::map(phenotypes, function(p) {
    ac_case_raw  <- cohort[[paste0(p, "_ac_case")]]
    nac_case_raw <- cohort[[paste0(p, "_nac_case")]]
    ac_ctrl_raw  <- cohort[[paste0(p, "_ac_ctrl")]]
    nac_ctrl_raw <- cohort[[paste0(p, "_nac_ctrl")]]

    ac_case <- replace(ac_case_raw, is.na(ac_case_raw), 0)
    ac_ctrl <- replace(ac_ctrl_raw, is.na(ac_ctrl_raw), 0)

    an_case_raw <- ac_case_raw + nac_case_raw
    an_ctrl_raw <- ac_ctrl_raw + nac_ctrl_raw
    an_case     <- impute_an_neighbors(cohort$locus, an_case_raw)
    an_ctrl     <- impute_an_neighbors(cohort$locus, an_ctrl_raw)

    or_int   <- purrr::pmap_dfr(list(ac_case, an_case, ac_ctrl,       an_ctrl),    safe_or)
    or_nfe   <- purrr::pmap_dfr(list(ac_case, an_case, gnomad_ac_nfe, an_nfe_ind), safe_or)
    or_total <- purrr::pmap_dfr(list(ac_case, an_case, gnomad_ac_total, an_total_ind), safe_or)

    block <- list(
      ac_case              = ac_case,
      an_case              = an_case,
      an_case_imputed      = is.na(an_case_raw),

      ac_ctrl_int          = ac_ctrl,
      an_ctrl_int          = an_ctrl,
      an_ctrl_int_imputed  = is.na(an_ctrl_raw),
      or_int               = or_int$OR,
      ci_low_int           = or_int$OR_lower,
      ci_high_int          = or_int$OR_upper,
      p_int                = or_int$p_value,

      ac_ctrl_nfe          = gnomad_ac_nfe,
      an_ctrl_nfe          = an_nfe_ind,
      an_ctrl_nfe_imputed  = is.na(gnomad_an_nfe),
      or_nfe               = or_nfe$OR,
      ci_low_nfe           = or_nfe$OR_lower,
      ci_high_nfe          = or_nfe$OR_upper,
      p_nfe                = or_nfe$p_value,

      ac_ctrl_total        = gnomad_ac_total,
      an_ctrl_total        = an_total_ind,
      an_ctrl_total_imputed = is.na(gnomad_an_total),
      or_total             = or_total$OR,
      ci_low_total         = or_total$OR_lower,
      ci_high_total        = or_total$OR_upper,
      p_total              = or_total$p_value
    )
    names(block) <- paste0(p, "_", names(block))
    tibble::as_tibble(block)
  }) |>
    dplyr::bind_cols()

  blocks
}
