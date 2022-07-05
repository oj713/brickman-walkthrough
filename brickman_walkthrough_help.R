suppressPackageStartupMessages({
  library(tidymodels) # also loads ggplot, dplyr, and purrr
  library(stars) # also loads sf
  library(brickman) # for brickman data access
  library(viridis) # for plotting 
})

#' Retrieves predictions from a model for a given brickman scenario by month. 
#' 
#' @param wkf fitted workflow, the workflow to use when generating predictions
#' @param brickman_vars chr, the list of brickman covariates used by the workflow
#' @param year int, the desired year for predictions -- NA if for present day
#' @param scenario chr, the desired conservation scenario for predictions
#' @param augment_preds boolean, should covariates be bound to the output predictions?
#' @param downsample int, the desired resolution of predictions. 0 is original resolution, 
#'    1 is high resolution, 2 is medium resolution, and 3+ is low resolution. 
#'    Downsampling will significantly reduce file size of plots and prediction tables. 
#' @verbose boolean, display some measure of progress?
#' @return list of 12 dfs named by month, where each df contains predictions. 
#'   df columns include lat, lon, month, .pred_0, .pred_1, 
#'   .pred_class, and optional covariates
get_predictions <- function(wkf, 
                            brickman_vars,
                            year = c(NA, 2055, 2075)[1], 
                            scenario = c("PRESENT", "RCP45", "RCP85")[1], 
                            augment_preds = FALSE,
                            verbose = FALSE,
                            downsample = c(0, 1, 2, 3)[2]) {
  
  if (verbose) {message("Beginning brickman variable retrieval")}
  
  # since bathymetry doesn't have a monthly component and doesn't change 
  # signficantly with time, it must be retrieved separately from other 
  # covariates if it exists
  bathy <- NULL
  not_bathymetry <- !brickman_vars %in% "Bathy_depth"
  
  if (!all(not_bathymetry)) {
    bathy <- brickman::read_brickman(scenario = "PRESENT", 
                                     vars = "Bathy_depth", 
                                     interval = "ann", 
                                     form = "stars")
  }
  
  # creating new list of variables omitting bathymetry
  mon_vars <- brickman_vars[not_bathymetry] 
  
  # any future env covariates are expressed as difference from present day, 
  #   so a copy of present day brickman data must be kept on hand 
  present_data <- brickman::read_brickman(scenario="PRESENT", 
                                          vars = mon_vars, 
                                          interval = "mon", 
                                          form = "stars")
  
  # retrieving climate data for desired scenario
  if (scenario == "PRESENT") {
    brickman_data <- present_data
  } else {
    brickman_data <- brickman::read_brickman(scenario = scenario, 
                                             year = year, 
                                             vars = mon_vars, 
                                             add = present_data,
                                             interval = "mon",
                                             form = "stars")
  }
  
  if (verbose) {message("Finished brickman variable retrieval")}
  
  process_month <- function(mon) {
    if (verbose) {message(paste("Processing month:", mon))}
    
    # mon_pred is a new input set of data 
    # taking brickman data and slicing it to only represent desired month
    # merging this monthly data with bathymetry information
    mon_data <- c(brickman_data |> dplyr::slice(month, mon),
                  bathy) |>
      # downsampling, if desired
      stars::st_downsample(n = downsample) |>
      # converting to tibble and making minor name/type changes 
      as_tibble() |>
      rename(lon = x, lat = y) |>
      mutate(MONTH = as.factor(mon)) |>
      na.omit()
    
    # determining whether to return augmented data or simply probability
    non_augment <- function(wkf, data) {
      predict(wkf, data, type = "prob") |>
        bind_cols(select(data, lon, lat, MONTH)) |>  
        mutate(.pred_class = (.pred_1 > .5) |> as.numeric() |> as.factor())
    }
    
    #depending on whether we are augmenting or not
    predict_method <- ifelse(augment_preds, augment, non_augment)
    
    # using input data to generate predictions
    predict_method(wkf, mon_data)
  }
  
  # returning prediction list
  list("Jan" = 1, "Feb" = 2, "Mar" = 3, "Apr" = 4, 
       "May" = 5, "Jun" = 6, "Jul" = 7, "Aug" = 8,
       "Sep" = 9, "Oct" = 10, "Nov" = 11, "Dec" = 12) |>
    lapply(process_month)
}

#' Converts a list of predictions into probability /change in probability plots
#' 
#' @param preds_list list of dfs named by month, used to generate plots
#' @param title chr, desired plot title. Month will be appended. 
#' @param pt_size dbl, the point size of the plot. A suitable value for 
#'   this will depend on the downsampling value for predictions and any bounds. 
#' @param xlim dbl, a vector of length 2 with x limits for plot. NULL if none. 
#' @param ylim dbl, a vector of length 2 with y limits for plot. NULL if none. 
#' @param comparison_list, list of dfs named by month used as base for comparisons. 
#'   Note that months included and downsample value must match preds_list. If
#'   NULL, return raw probability predictions
#' @return a list of ggplot objects representing their corresponding predictions
get_value_plots <- function(preds_list,
                            title = "Predicted Presence Probability",
                            pt_size = .3,
                            xlim = NULL, 
                            ylim = NULL,
                            comparison_list = NULL) {
  
  # determining whether this is a comparison plot 
  is_comparison <- !is.null(comparison_list)
  
  # helper function that generates a plot for a dataframe
  plot_month <- function (preds, mon_name) {
    # generating plot
    plot_base <- ggplot(preds, aes(x = lon, y = lat, col = .pred_1)) +
      geom_point(cex = pt_size, pch = 15) +
      coord_quickmap(xlim = xlim,
                     ylim = ylim) + 
      # aesthetic changes - all optional
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            legend.position = "bottom") +
      labs(x = "Latitude", y = "Longtitude", color = NULL) +
      ggtitle(paste(title, "-", mon_name))
    
    # desired color scheme depends on whether this is raw or comparison plot
    if (is_comparison) {
      plot <- plot_base +
        scale_color_gradientn(limits = c(-.7, .7),
                              colors = c("midnightblue",
                                         "dodgerblue3",
                                         "deepskyblue1",
                                         "gray94", 
                                         "goldenrod1",
                                         "darkorange1", 
                                         "orangered3"),
                              na.value = "white")
    } else {
      plot <- plot_base +
        scale_color_viridis(option = "inferno", limits = c(0, 1))
    }
    
    plot
  }
  
  # for comparison plots, subtract comparison data from original
  if (is_comparison) {
    plot_data_list <- 
      purrr::map2(preds_list, 
                  comparison_list, 
                  ~mutate(.x, .pred_1 = .pred_1 - .y$.pred_1))
  } else {
    plot_data_list <- preds_list
  }
  
  # generating plot list from prediction list
  plot_data_list |>
    purrr::imap(~plot_month(.x, .y))
}

#' Converts a list of predictions into change in presence plots that show 
#'   how predictions shifted relative to a desired threshold. 
#' 
#' @param preds_list list of dfs named by month, used to generate plots
#' @param comparison_list, list of dfs named by month used as base for comparisons. 
#'   Note that months included and downsample value must match preds_list.
#' @param threshold dbl, the desired threshold to distinguish between low and 
#'   high predictive values. 
#' @param title chr, desired plot title. Month will be appended. 
#' @param pt_size dbl, the point size of the plot. A suitable value for 
#'   this will depend on the downsampling value for predictions and any bounds. 
#' @param xlim dbl, a vector of length 2 with x limits for plot. NULL if none. 
#' @param ylim dbl, a vector of length 2 with y limits for plot. NULL if none. 
#' @return a list of ggplot objects representing their corresponding predictions
get_threshold_plots <- function(preds_list,
                                comparison_list,
                                threshold = .5,
                                title = "Change in Presence (Threshold: .5)",
                                pt_size = .3,
                                xlim = NULL,
                                ylim = NULL) {
  
  # naming factor levels 
  # COMPARISONPRESENCE_PREDSPRESENCE
  feedstatus <- list(FALSE_FALSE = "No Presence",
                     FALSE_TRUE = "New Presence", 
                     TRUE_FALSE = "Lost Presence", 
                     TRUE_TRUE = "Presence")
  # palette for colors 
  pal <- c(FALSE_FALSE = "gray90",
           FALSE_TRUE = "red3", 
           TRUE_FALSE = "dodgerblue2", 
           TRUE_TRUE = "goldenrod1")
  
  # helper to generate plot
  plot_month <- function (preds, mon_name) {
    # generating plot
    plot <- ggplot(preds, aes(x = lon, y = lat, col = REL_PRESENCE)) +
      geom_point(cex = pt_size, pch = 15) +
      scale_color_manual(labels = feedstatus,
                         values = pal) + 
      coord_quickmap(xlim = xlim,
                     ylim = ylim) + 
      # aesthetic changes - all optional
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            legend.position = "bottom") +
      guides(colour = guide_legend(override.aes = list(size=2))) +
      labs(x = "Latitude", y = "Longtitude", color = "Presence Status") +
      ggtitle(paste(title, "-", mon_name))
  }
  
  # adding column to predictions data to represent threshold change
  plot_data_list <-  
    purrr::map2(preds_list, 
                comparison_list, 
                ~mutate(.x, 
                        REL_PRESENCE = paste0(.y$.pred_1 > threshold,
                                              "_",
                                              .x$.pred_1 > threshold)))

  # generating plot list from prediction list
  plot_data_list |>
    purrr::imap(~plot_month(.x, .y))
}
