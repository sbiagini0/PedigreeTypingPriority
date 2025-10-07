#' Run typing-priority simulations and export figures (pedigree + power plot)
#' Builds selection sets from `members_selected` (optionally capped by `k_max`),
#' runs `forrel::MPPsims()` with the missing person labeled `"MP"`, generates a
#' type-3 `forrel::powerPlot()`, saves both figures into temporary subfolders
#' under `Output/`, composes a side-by-side JPEG in `output/`, and removes the
#' temporary subfolders.
#' @param pedigree A `ped` object to be analyzed.
#' @param members_selected Character vector of candidate IDs to type.
#' @param baseline Logical; forwarded to `MPPsims(addBaseline)`.
#' @param nsim Integer; number of LR simulations (`MPPsims(lrSims)`).
#' @param nprofiles Integer; simulated profiles per LR simulation (`MPPsims(nProfiles)`).
#' @param threshold Numeric; inclusion threshold for IP (`MPPsims(thresholdIP)`).
#' @param mutations Logical; forwarded to `MPPsims(disableMutations)`.
#' @param ncores Integer; CPU cores for parallel computation (`MPPsims(numCores)`).
#' @param seeds Integer; random seed for reproducibility.
#' @param k_max Optional integer; maximum size of candidate combinations.
#' @return Invisibly, a named list with:
#'   \item{image}{`magick-image` composite (pedigree + power plot).}
#'   \item{path}{Absolute path to the composite JPEG saved in `output/`.}
#'   \item{simData}{The object returned by `forrel::MPPsims()`.}
typingPriority <- function(pedigree, members_selected, baseline, nsim, nprofiles, threshold,
                           mutations, ncores, seeds, k_max = NULL) {
  
  stopifnot(pedtools::is.ped(pedigree))
  if (!is.character(members_selected) || length(members_selected) == 0)
    stop("members_selected must be a non-empty character vector")
  
  set.seed(seeds)
  
  # ---- Build selections ----
  typed_now <- pedtools::typedMembers(pedigree)
  selections <- generate_combinations(
    members_selected = members_selected,
    typedMembers     = typed_now,
    k_max            = k_max
  )
  if (length(selections) == 0) stop("No selections generated")
  
  # ---- Print selections ----
  n_sets <- length(selections)
  hdr <- c(
    sprintf("Total selection sets : %d", n_sets),
    sprintf("Max combo size (kmax): %s", ifelse(is.null(k_max), "all", as.character(k_max)))
  )
  show_idx <- seq_len(n_sets)
  lines_sets <- vapply(show_idx, function(i) {
    sprintf("%s: %s", names(selections)[i], paste(selections[[i]], collapse = ", "))
  }, character(1))
  lines <- c(hdr, "", lines_sets)
  width <- max(nchar(lines), 1L)
  cat("\n",
      paste0("┌", strrep("─", width + 2), "┐"), "\n",
      paste0("│ ", sprintf(paste0("%-", width, "s"), lines), " │", collapse = "\n"), "\n",
      paste0("└", strrep("─", width + 2), "┘"), "\n\n", sep = "")
  
  # ---- Run simulations ----
  simData <- forrel::MPPsims(
    reference        = pedigree,
    missing          = "MP",
    selections       = selections,
    ep               = TRUE,
    ip               = TRUE,
    addBaseline      = baseline,
    lrSims           = nsim,
    nProfiles        = nprofiles,
    thresholdIP      = threshold,
    disableMutations = mutations,
    numCores         = ncores,
    seed             = seeds,
    verbose          = TRUE
  )
  
  # ---- Paths ----
  fam <- pedtools::famid(pedigree)
  output_dir      <- file.path(getwd(), "output")
  pedigree_dir    <- file.path(output_dir, "Pedigree plot")
  simulation_dir  <- file.path(output_dir, "Simulation plot")
  
  if (!dir.exists(output_dir))     dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(pedigree_dir))   dir.create(pedigree_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(simulation_dir)) dir.create(simulation_dir, recursive = TRUE, showWarnings = FALSE)
  
  powerplot_path <- file.path(simulation_dir, paste0(fam, "_powerplot.jpeg"))
  ped_path       <- file.path(pedigree_dir,   paste0(fam, "_pedigree.jpeg"))
  final_path     <- file.path(output_dir,     paste0(fam, "_PedigreeTypingPriority.jpeg"))
  
  # ---- Power plot ----
  powerplot <- forrel::powerPlot(simData, type = 3) +
    ggplot2::labs(
      subtitle = paste0(
        "N° Simulations: ", nsim,
        "  |  Profiles: ", nprofiles,
        "  |  Seed: ", seeds,
        "  |  Date: ", format(Sys.Date(), "%Y-%m-%d")
      )
    ) +
    ggplot2::theme(plot.title.position = "plot")
  
  ggplot2::ggsave(powerplot_path, plot = powerplot, width = 2160, height = 2160, dpi = 300, units = "px")
  
  # ---- Pedigree plot ----
  grDevices::jpeg(ped_path, width = 3840, height = 2160, res = 300)
  dev_id <- grDevices::dev.cur()
  on.exit({
    if (!is.null(grDevices::dev.list()) && dev_id %in% grDevices::dev.list()) {
      try(grDevices::dev.off(dev_id), silent = TRUE)
    }
  }, add = TRUE)
  
  nm <- pedtools::nMarkers(pedigree)
  
  tryCatch({
    plot(
      pedigree,
      hatched  = pedtools::typedMembers(pedigree),
      title    = fam,
      carrier  = members_selected,
      col      = list(red = "MP", black = members_selected),
      margins  = c(5.1, 8, 5.1, 8),
      keep.par = TRUE
    )
    mtext(sprintf("nMarkers: %d", nm), side = 3, line = 0.5, cex = 0.9, font = 3)
    grDevices::dev.off(dev_id)
  }, error = function(e) {
    if (!is.null(grDevices::dev.list()) && dev_id %in% grDevices::dev.list()) {
      try(grDevices::dev.off(dev_id), silent = TRUE)
    }
    stop(sprintf("Failed to create pedigree plot: %s", e$message), call. = FALSE)
  })
  
  # ---- Compose final and clean temp subfolders ----
  ped_img   <- magick::image_read(ped_path)
  power_img <- magick::image_read(powerplot_path)
  img <- magick::image_append(c(ped_img, power_img), stack = FALSE)
  magick::image_write(img, final_path)
  
  # --- Remove temporary subfolders ---
  unlink(pedigree_dir,   recursive = TRUE, force = TRUE)
  unlink(simulation_dir, recursive = TRUE, force = TRUE)
  
  # --- Return plot ---
  return(knitr::include_graphics(final_path))
  
  invisible(list(image = img, path = final_path, simData = simData))
  
  gc()
}
