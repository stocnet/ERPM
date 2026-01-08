# ======================================================================================
# Fichier : cub_test.R
# Fonction : scénario de tests ERPM/ERGM
# Utilité  : valider un effet via summary/ergm/erpm (partition et biparti)
# ======================================================================================

# macOS/Linux : UTF-8 propre
Sys.setenv(LANG="fr_FR.UTF-8")
invisible(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"))

# ======================================== INIT ========================================
suppressMessages(library(ergm))

if (!exists(".__debug_loaded", envir = .GlobalEnv)) {
  source("scripts/local_source/debug.R", local = FALSE)
}

# init global (settings + launcher + logging)
source("scripts/local_source/init.R", local = FALSE)
# wrapper ERPM (chargé à part)
source("R/erpm_wrapper.R", local = FALSE)
# debug_source_file("R/erpm_wrapper.R")

init_erpm(selftest = FALSE, verbose = FALSE, install_verbose = FALSE)

if (exists("log_msg")) log_msg("INFO", "Démarrage du script pour ERPM") else message("Démarrage du script pour ERPM")

# Patch ERGM optionnel
if (exists("ergm_patch_enable")) ergm_patch_enable(verbose = VERBOSE)

# ======================================== RUN =========================================

# ---------- Scénario 1 ----------
part1 <- c(1,1, 2,2,2, 3,3,3,3, 4)  # tailles: 2,3,4,1  => N=10 # nolint # nolint: commas_linter.
noeuds1 <- data.frame(
  nom      = LETTERS[1:length(part1)],
  sexe     = c("F","F",  "H","H","H",  "F","F","F","F",  "H"),
  service  = c("RH","RH","Vent","Vent","Vent", "RH","RH","Info","Info","Info"),
  promo    = c(2020,2020, 2019,2019,2019, 2021,2021,2021,2021, 2018),
  stringsAsFactors = FALSE
)

launch_model(engine="summary",
             effects = "cov_fullmatch('sexe')",
             partition = part1, nodes = noeuds1,
             estimate = "CD", eval_loglik = FALSE, verbose = TRUE)

launch_model(engine="summary",
             effects = "cov_fullmatch('sexe', size = 2:3)",
             partition = part1, nodes = noeuds1,
             estimate = "CD", eval_loglik = FALSE, verbose = TRUE)

launch_model(engine="summary",
             effects = "cov_fullmatch('sexe', category = 'F')",
             partition = part1, nodes = noeuds1,
             estimate = "CD", eval_loglik = FALSE, verbose = TRUE)

launch_model(engine="summary",
             effects = "cov_fullmatch('service') + cov_fullmatch('service', category = 'Info', size = c(3,4))",
             partition = part1, nodes = noeuds1,
             estimate = "CD", eval_loglik = FALSE, verbose = TRUE)

launch_model(engine="summary",
             effects = quote(cov_fullmatch(cov = noeuds1$promo) + cov_fullmatch(cov = noeuds1$promo, size = 3:4)),
             partition = part1, nodes = noeuds1,
             estimate = "CD", eval_loglik = FALSE, verbose = TRUE)


# ---------- Scénario 2 ----------
part2 <- c(1,1,1, 2,2, 3,3, 4,4,4,4, 5)  # tailles: 3,2,2,4,1 => N=12
noeuds2 <- data.frame(
  nom      = paste0("V", seq_along(part2)),
  sexe     = c("H","H","H",  "F","F",  "H","F",  "F","F","F","F",  "H"),
  service  = c("Rech","Rech","Rech","RH","RH","Vent","Vent","Vent","Vent","Vent","Vent","RH"),
  promo    = c(2022,2022,2022, 2020,2020, 2019,2021, 2018,2018,2018,2018, 2023),
  stringsAsFactors = FALSE
)

launch_model(engine="summary",
             effects = "cov_fullmatch('sexe') + cov_fullmatch('sexe', category = 'H', size = c(2,3,4))",
             partition = part2, nodes = noeuds2,
             estimate = "CD", eval_loglik = FALSE, verbose = TRUE)

launch_model(engine="summary",
             effects = "cov_fullmatch('service', size = 4) + cov_fullmatch('service', category = 'Vent')",
             partition = part2, nodes = noeuds2,
             estimate = "CD", eval_loglik = FALSE, verbose = TRUE)

launch_model(engine="summary",
             effects = quote(cov_fullmatch(cov = noeuds2$promo, category = '2018')),
             partition = part2, nodes = noeuds2,
             estimate = "CD", eval_loglik = FALSE, verbose = TRUE)


# # ---------- Scénario 3 ----------
part3 <- c(1, 2,2, 3,3,3, 4,4, 5,5,5,5, 6)  # tailles: 1,2,3,2,4,1 => N=13
noeuds3 <- data.frame(
  nom      = paste0("N", seq_along(part3)),
  sexe     = c("F",    "H","H",       "F","F","H",    "F","H",        "F","F","F","F",              "H"),
  service  = c("Info", "Info","Info", "RH","RH","RH", "Rech","Rech",  "Rech","Rech","Rech","Rech",  "Vent"),
  promo    = c(2017, 2020,2020, 2021,2021,2021, 2019,2019, 2018,2018,2018,2018, 2022),
  stringsAsFactors = FALSE
)

launch_model(engine="summary",
             effects = "cov_fullmatch('sexe') + cov_fullmatch('sexe', size = c(1,4))",
             partition = part3, nodes = noeuds3,
             estimate = "CD", eval_loglik = FALSE, verbose = TRUE)

launch_model(engine="summary",
             effects = "cov_fullmatch('service', category = 'Rech') + cov_fullmatch('service', category = 'Info', size = 2)",
             partition = part3, nodes = noeuds3,
             estimate = "CD", eval_loglik = FALSE, verbose = TRUE)

launch_model(engine="summary",
             effects = quote(cov_fullmatch(cov = noeuds3$promo) + cov_fullmatch(cov = noeuds3$promo, category = '2018')),
             partition = part3, nodes = noeuds3,
             estimate = "CD", eval_loglik = FALSE, verbose = TRUE)

# ---------- Paramètres communs ERPM rapides ----------
ctrl_fast <- control.ergm(
  init.method   = "CD",
  CD.samplesize = 50000,
  CD.nsteps     = 5,
  MCMC.burnin   = 4000,
  MCMC.interval = 1000,
  MCMLE.maxit   = 1
)

ctrl_cd_only <- control.ergm(
  init.method      = "CD",
  CD.samplesize    = 100000,
  CD.nsteps        = 15,
  MCMC.burnin      = 20000,
  MCMC.interval    = 2000,
  force.main       = TRUE,
  parallel         = 0,
  seed             = 1
)
tmo <- 10

# ========== ERPM: cas "size=1" et "size=c(1,3)" ==========
p_sz <- c(1,    2,2,      3,3,3)
v_sz <- c("X",  "Y","Y",  "Z","Z","Z")
nodes_sz <- data.frame(val = v_sz, stringsAsFactors = FALSE)

# size = 1
err <- NULL; res <- NULL
tryCatch({
  res <- launch_model(engine="erpm",
                      effects = "cov_fullmatch('val', size = c(1))",
                      partition = p_sz, nodes = nodes_sz,
                      estimate = "MLE", eval_loglik = FALSE,
                      control = ctrl_fast, timeout = tmo, verbose = FALSE)
}, error = function(e) err <<- e$message)
if (is.null(err)) { log_msg("INFO","[ERPM] size=1 uniquement -> OK") } else { log_msg("ERROR",err); stop(err) }

# size = c(1,3)
err <- NULL; res <- NULL
tryCatch({
  res <- launch_model(engine="erpm",
                      effects = "cov_fullmatch('val', size = c(1,3))",
                      partition = p_sz, nodes = nodes_sz,
                      estimate = "CD", eval_loglik = FALSE,
                      control = ctrl_cd_only, timeout = tmo, verbose = FALSE)
}, error = function(e) err <<- e$message)
if (is.null(err)) { log_msg("INFO","[ERPM] size=1,3 -> OK") } else { log_msg("ERROR",err); stop(err) }

# ========== ERPM: cas tie 2v2 et filtre size=4 ==========
p_tie <- c(1,1,1,1, 2,2,2)
v_tie <- c("A","A","B","B", "Z","Z","Z")
nodes_tie <- data.frame(val = v_tie, stringsAsFactors = FALSE)

# sans filtre
err <- NULL; res <- NULL
tryCatch({
  res <- launch_model(engine="erpm",
                      effects = "cov_fullmatch('val')",
                      partition = p_tie, nodes = nodes_tie,
                      estimate = "CD", eval_loglik = FALSE,
                      control = ctrl_fast, timeout = tmo, verbose = FALSE)
}, error = function(e) err <<- e$message)
if (is.null(err)) { log_msg("INFO","[ERPM] tie 2v2 neutralisé -> OK") } else { log_msg("ERROR",err); stop(err) }

# size = 4
err <- NULL; res <- NULL
tryCatch({
  res <- launch_model(engine="erpm",
                      effects = "cov_fullmatch('val', size = c(4))",
                      partition = p_tie, nodes = nodes_tie,
                      estimate = "CD", eval_loglik = FALSE,
                      control = ctrl_fast, timeout = tmo, verbose = FALSE)
}, error = function(e) err <<- e$message)
if (is.null(err)) { log_msg("INFO","[ERPM] tie 2v2 size=4 -> OK") } else { log_msg("ERROR",err); stop(err) }

# ========== ERPM: numérique dense ==========
set.seed(42)
p_dense <- c(1,1,1, 2,2,2, 3,3,3, 4, 5)
val_dense <- c(1,2,3, 4,5,6, 7,8,9, 10, 11)
nodes_dense <- data.frame(val = val_dense, stringsAsFactors = FALSE)

# run sans filtre
err <- NULL; res <- NULL
tryCatch({
  res <- launch_model(engine="erpm",
                      effects = "cov_fullmatch('val')",
                      partition = p_dense, nodes = nodes_dense,
                      estimate = "CD", eval_loglik = FALSE,
                      control = ctrl_fast, timeout = tmo, verbose = FALSE)
}, error = function(e) err <<- e$message)
if (is.null(err)) { log_msg("INFO","[ERPM] dense sans filtre -> OK") } else { log_msg("ERROR",err); stop(err) }

# size = 3
err <- NULL; res <- NULL
tryCatch({
  res <- launch_model(engine="erpm",
                      effects = "cov_fullmatch('val', size = c(3))",
                      partition = p_dense, nodes = nodes_dense,
                      estimate = "CD", eval_loglik = FALSE,
                      control = ctrl_fast, timeout = tmo, verbose = FALSE)
}, error = function(e) err <<- e$message)
if (is.null(err)) { log_msg("INFO","[ERPM] dense size=3 -> OK") } else { log_msg("ERROR",err); stop(err) }

# ======================================== CLEAN =======================================
if (exists("ergm_patch_disable")) ergm_patch_disable(verbose = VERBOSE)
flush.console()
if (exists("log_msg")) log_msg("INFO", "Fin du programme -- Nettoyage de l'environnement global")
clean_global_env(verbose = VERBOSE)