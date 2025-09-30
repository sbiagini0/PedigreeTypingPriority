#' Load and (if missing) install R packages
#' Installs any packages not yet present and then loads all requested packages.
#' @param pkgs Character vector with package names to ensure installed and loaded.
#' @return Invisibly, a list with the results of `library()` for each package.
load_packages <- function(pkgs) {
  installed_packages <- pkgs %in% rownames(installed.packages())
  if (any(!installed_packages)) {
    install.packages(pkgs[!installed_packages])
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

################################################################################
#' Read a Familias .fam file (from a path or by scanning a folder)
#' If `file` is provided, reads that file. Otherwise scans `path` for a single
#' `.fam` file. Errors if none or more than one are found.
#' @param file Character path to a `.fam` file (optional).
#' @param path Folder to search for `.fam` if `file` is NULL. Default: "data".
#' @return A ped/pedlist on success; otherwise NULL.
read_famfile <- function(file = NULL, path = "data") {
  if (is.null(file)) {
    fams <- list.files(path = path, pattern = "\\.fam$", full.names = TRUE)
    if (length(fams) == 0) {
      stop(sprintf("No .fam files found in '%s'. Provide `file` or place one .fam there.", path))
    }
    if (length(fams) > 1) {
      stop(sprintf(
        "Multiple .fam files found in '%s':\n- %s\nSpecify `file` explicitly.",
        path, paste(basename(fams), collapse = "\n- ")
      ))
    }
    fam_file <- fams[1]
  } else {
    fam_file <- file
  }
  
  tryCatch({
    pedFamilias::readFam(
      fam_file,
      useDVI      = TRUE,
      prefixAdded = "EXTRA",
      simplify1   = TRUE,
      deduplicate = TRUE,
      verbose     = FALSE
    )
  }, error = function(e) {
    message("Error reading Familias file: ", e$message)
    NULL
  })
}

################################################################################
#' Format a pedigree list for analysis
#' Harmonizes markers across pedigrees, removes an internal `_comp1` wrapper if present,
#' reorders individuals, sorts genotypes, renames "Missing person" to "MP",
#' and sets each pedigree's `famid` from the list names.
#' @param mpi A list of `ped` objects
#' @return A cleaned list of `ped` objects with consistent markers, labeling, and `famid`.
pedFormat <- function(mpi) {
  if (is.null(mpi) || length(mpi) == 0)
    return(mpi)
  
  mpi <- pedtools::harmoniseMarkers(mpi)
  
  mpi <- lapply(mpi, function(p) {
    if (is.list(p) && !is.null(p[["_comp1"]])) {
      p <- p[["_comp1"]]
    }
    p
  })
  
  mpi <- lapply(mpi, function(p) {
    if (!pedtools::is.ped(p))
      return(p)
    p <- pedtools::reorderPed(p)
    p <- pedtools::sortGenotypes(p)
    p <- pedtools::relabel(p, new = "MP", old = "Missing person")
  })
  
  for (i in seq_along(mpi)) { 
    if (pedtools::is.ped(mpi[[i]])) { 
      pedtools::famid(mpi[[i]]) <- names(mpi)[i] 
      } }
  
  mpi
}

################################################################################
#' Collect simulation parameters via pop-ups
#' Opens simple dialog boxes to gather and validate simulation settings, then prints
#' a compact summary in the console.
#' @return Named list with:
#'   \item{nsim}{integer, number of simulations (>= 1)}
#'   \item{nprofiles}{integer, profiles per simulation (>= 1)}
#'   \item{baseline}{logical, whether to simulate baseline pedigree}
#'   \item{mutations}{logical, whether to deactivate mutations}
#'   \item{threshold}{numeric, LR threshold (>= 0)}
#'   \item{ncores}{integer, CPU cores used (1..detectCores)}
#'   \item{seed}{integer, random seed}
simulation_parameters <- function() {
  if (!requireNamespace("svDialogs", quietly = TRUE))
    stop("Package 'svDialogs' is required. Install it with install.packages('svDialogs').")
  
  max_cores <- max(1L, parallel::detectCores(logical = TRUE))
  defaults <- list(
    nsim      = 1000L,
    nprofiles = 10L,
    baseline  = "TRUE",
    mutations = "FALSE",
    threshold = 10000,
    ncores    = min(10L, max_cores),
    seed      = as.integer(round(as.numeric(Sys.time()))) %% .Machine$integer.max
  )
  
  ask_num <- function(msg, def, integer = FALSE, min = -Inf, max = Inf) {
    repeat {
      res <- svDialogs::dlgInput(message = paste0(msg, " [Default: ", def, "]"),
                                 default = as.character(def))$res
      if (is.null(res) || res == "") res <- as.character(def)
      x <- suppressWarnings(if (integer) as.integer(res) else as.numeric(res))
      if (!is.na(x) && is.finite(x) && x >= min && x <= max) return(x)
      svDialogs::dlgMessage("Invalid value. Please try again.", type = "ok")
    }
  }
  ask_logical <- function(msg, def = "TRUE") {
    pick <- svDialogs::dlgList(choices = c("TRUE","FALSE"),
                               preselect = def, multiple = FALSE,
                               title = msg)$res
    if (is.null(pick) || pick == "") pick <- def
    identical(pick, "TRUE")
  }
  
  nsim      <- ask_num("Enter the number of simulations", defaults$nsim, integer = TRUE, min = 1)
  nprofiles <- ask_num("Enter the number of profiles to simulate", defaults$nprofiles, integer = TRUE, min = 1)
  baseline  <- ask_logical("Baseline pedigree simulation?", def = defaults$baseline)
  mutations <- ask_logical("Deactivate mutations?", def = defaults$mutations)
  threshold <- ask_num("Enter the LR threshold", defaults$threshold, integer = FALSE, min = 0)
  ncores    <- ask_num(paste0("Enter the number of cores (max: ", max_cores, ")"),
                       defaults$ncores, integer = TRUE, min = 1, max = max_cores)
  seed      <- ask_num("Enter the seed", defaults$seed, integer = TRUE)
  
  params <- list(
    nsim = as.integer(nsim),
    nprofiles = as.integer(nprofiles),
    baseline = baseline,
    mutations = mutations,
    threshold = as.numeric(threshold),
    ncores = as.integer(ncores),
    seed = as.integer(seed)
  )
  
  lines <- c(
    sprintf("Number of simulations       : %d", params$nsim),
    sprintf("Profiles per simulation     : %d", params$nprofiles),
    sprintf("Baseline                    : %s", params$baseline),
    sprintf("Mutations deactivated       : %s", params$mutations),
    sprintf("LR threshold                : %g", params$threshold),
    sprintf("Cores                       : %d", params$ncores),
    sprintf("Seed                        : %d", params$seed)
  )
  width <- max(nchar(lines))
  top    <- paste0("┌", strrep("─", width + 2), "┐")
  bottom <- paste0("└", strrep("─", width + 2), "┘")
  body   <- paste0("│ ", sprintf(paste0("%-", width, "s"), lines), " │")
  cat("\n", top, "\n", paste(body, collapse = "\n"), "\n", bottom, "\n\n", sep = "")
  
  params
}

################################################################################
#' Select a single pedigree from an ped list.
#' Opens a dialog to choose one pedigree from `mpi`, prints a short summary box,
#' and returns the selected `ped` object.
#' @param mpi Named list of `ped` objects
#' @return A `ped` object corresponding to the selected pedigree.
select_pedigree <- function(mpi) {
  pedigrees <- select.list(choices = names(mpi), multiple = FALSE, graphics = TRUE, title = "Select pedigree:")
  
  if (is.null(pedigrees) || pedigrees == "") stop("Error")
  
  pedigree <- mpi[[pedigrees]]
  
  line  <- sprintf("Pedigree selected: %s", pedigree$FAMID)
  width <- nchar(line)
  top    <- paste0("┌", strrep("─", width + 2), "┐")
  bottom <- paste0("└", strrep("─", width + 2), "┘")
  body   <- paste0("│ ", sprintf(paste0("%-", width, "s"), line), " │")
  cat("\n", top, "\n", body, "\n", bottom, "\n\n", sep = "")
  
  return(pedigree)
}

################################################################################
#' Interactively edit a pedigree (add parents or a child)
#' Shows the pedigree plot, then iteratively lets you either add parents to a
#' selected child or add a new child. A special option
#' "Not listed (type a new one)" allows entering a new father/mother if the
#' desired individual is not in the current pedigree. The plot is refreshed
#' after each change.
#' @param pedigree A `ped` object to edit interactively.
#' @return The updated `ped` object.
edit_pedigree <- function(pedigree) {
  plot(pedigree, hatched = pedtools::typedMembers(pedigree), title = paste0("Before: ", pedigree$FAMID), arrows = FALSE)
  repeat {
    plot(pedigree, hatched = pedtools::typedMembers(pedigree), title = paste0("Before: ", pedigree$FAMID), arrows = FALSE)
    action <- select.list(
      choices  = c("Add parents", "Add child", "Finish"),
      multiple = FALSE, graphics = TRUE, title = "What do you want to do?"
    )
    
    if (is.null(action) || action == "" || action == "Finish") break
    #---Parents---
    if (action == "Add parents") {
      id_son <- select.list(
        choices  = pedtools::founders(pedigree),
        multiple = TRUE, graphics = TRUE,
        title    = "Select a child to add parents:"
      )
      if (is.null(id_son) || length(id_son) == 0) id_son <- NULL
      
      father <- svDialogs::dlgInput(message = "Father name:")[["res"]]
      mother <- svDialogs::dlgInput(message = "Mother name:")[["res"]]
      
      pedigree <- pedtools::addParents(
        pedigree, father = father, mother = mother, id = id_son, verbose = TRUE
      )
      
      plot(pedigree, hatched = pedtools::typedMembers(pedigree), title = paste0("After: ", pedigree$FAMID), arrows = FALSE)
    }
    
    #---Child---
    if (action == "Add child") {
      sentinel <- "Not listed (type a new one)"
      
      # Father
      choices_f <- c(pedtools::males(pedigree), sentinel)
      id_father <- select.list(
        choices  = choices_f,
        multiple = FALSE, graphics = TRUE,
        title    = "Select father to add a child:"
      )
      if (identical(id_father, sentinel)) {
        tmp <- svDialogs::dlgInput(message = "Enter NEW father name:")[["res"]]
        id_father <- if (is.null(tmp) || tmp == "") NULL else tmp
      } else if (is.null(id_father) || id_father == "") {
        tmp <- svDialogs::dlgInput(message = "Add father name:")[["res"]]
        id_father <- if (is.null(tmp) || tmp == "") NULL else tmp
      }
      
      # Mother
      choices_m <- c(pedtools::females(pedigree), sentinel)
      id_mother <- select.list(
        choices  = choices_m,
        multiple = FALSE, graphics = TRUE,
        title    = "Select mother to add a child:"
      )
      if (identical(id_mother, sentinel)) {
        tmp <- svDialogs::dlgInput(message = "Enter NEW mother name:")[["res"]]
        id_mother <- if (is.null(tmp) || tmp == "") NULL else tmp
      } else if (is.null(id_mother) || id_mother == "") {
        tmp <- svDialogs::dlgInput(message = "Add mother name:")[["res"]]
        id_mother <- if (is.null(tmp) || tmp == "") NULL else tmp
      }
      
      children <- svDialogs::dlgInput(message = "Child name:")[["res"]]
      sex_children <- svDialogs::dlgInput(message = "Child sex (unknown=0, male=1, female=2):")[["res"]]
      sex_children <- as.numeric(sex_children)
      
      pedigree <- pedtools::addChildren(
        pedigree, father = id_father, mother = id_mother,
        nch = 1, sex = sex_children, ids = children, verbose = TRUE
      )
      
      plot(pedigree, hatched = pedtools::typedMembers(pedigree), title = paste0("After: ", pedigree$FAMID), arrows = FALSE)
      
    }
  }
  
  return(pedigree)
}

################################################################################
#' Select untyped pedigree members.
#' Opens a GUI list of untyped members to pick one or more IDs for simulation; 
#' stops on cancel and prints a short summary box.
#' @param pedigree A `ped` or pedlist from which untyped members are derived.
#' @return Character vector with the selected member IDs.
select_members_for_simulation <- function(pedigree) {
  typed   <- pedtools::typedMembers(pedigree)
  untyped <- pedtools::untypedMembers(pedigree)
  
  members_selected <- select.list(
    choices  = untyped,
    multiple = TRUE,
    graphics = TRUE,
    title    = "Select pedigree members to perform simulation:"
  )
  
  if (is.null(members_selected) || length(members_selected) == 0 || identical(members_selected, "")) stop("Error")
  
  line  <- sprintf("Pedigree members selected: %s", paste(members_selected, collapse = ", "))
  width <- nchar(line)
  top    <- paste0("┌", strrep("─", width + 2), "┐")
  bottom <- paste0("└", strrep("─", width + 2), "┘")
  body   <- paste0("│ ", sprintf(paste0("%-", width, "s"), line), " │")
  cat("\n", top, "\n", body, "\n", bottom, "\n\n", sep = "")
  
  return(members_selected)
}

################################################################################
#' Build selection sets from candidate and typed members
#' Creates all combinations (up to `k_max`) of candidate members to be newly typed
#' and, for each combination, appends all currently typed members not in that combo.
#' @param members_selected Character vector of candidate IDs proposed to type.
#' @param typedMembers Character vector of IDs already typed.
#' @param k_max Optional integer limiting the maximum size of candidate combinations.
#' @return Named list of selections; each element is a character vector:
#'         c(<candidate combo>, <typed not in combo>).
generate_combinations <- function(members_selected, typedMembers, k_max = NULL) {
  n <- length(members_selected)
  if (n == 0) return(list())
  
  ks <- 1:n
  if (!is.null(k_max)) ks <- ks[ks <= k_max]
  
  combos <- unlist(lapply(ks, function(k) combn(members_selected, k, simplify = FALSE)), recursive = FALSE)
  
  selections <- vector("list", length(combos))
  names(selections) <- vapply(combos, function(z) paste(sort(z), collapse = "_"), character(1))
  
  for (i in seq_along(combos)) {
    combo <- combos[[i]]
    typed_not_in_combo <- setdiff(typedMembers, combo)
    selections[[i]] <- c(combo, typed_not_in_combo)
  }
  
  selections
}
