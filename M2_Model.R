rm(list = ls())

library(dplyr)
library(foreach)
library(doParallel)

source("utils.R")

# GLOBAL PARAMETERS ------------------------------------------------------ START
# dataset_file <- "Disease"
# dataset_file <- "EYT_1"
# dataset_file <- "EYT_2"
# dataset_file <- "EYT_3"
# dataset_file <- "Groundnut"
# dataset_file <- "Indica"
# dataset_file <- "Japonica"
dataset_file <- "Maize"
predictors <- c("Line")
predictors_text <- paste0(predictors, collapse = " + ")
iterations_number <- 10000
burn_in <- 5000
base_results_dir <- "results"
type <- "M2"
cores_num <- 3
# GLOBAL PARAMETERS -------------------------------------------------------- END

load(sprintf("data/%s.RData", dataset_file), verbose = TRUE)

OriginalPheno <- Pheno
all_envs <- unique(Pheno$Env)
sorted_lines <- sort(rownames(Geno))
Geno <- Geno[sorted_lines, sorted_lines]

# Prepare clusters for parallelization -----------------------------------------
cluster <- makeCluster(
  cores_num,
  outfile = sprintf("%s.out", data_info$name)
)
registerDoParallel(cluster)

for (i in seq_along(data_info$traits)) {
  SKM::echo("* Trait %s/%s", i, length(data_info$traits))
  trait <- data_info$traits[i]

  Pheno <- prepare_pheno(OriginalPheno, trait)
  folds <- cv0(Pheno, all_envs)
  y <- Pheno$Diff

  Line <- dummy_matrix(Pheno$Line)
  G <- Line %*% Geno %*% t(Line)
  ETA <- list(list(K = G, model = "RKHS"))

  trait_results_dir <- file.path(
    base_results_dir,
    type,
    predictors_text,
    data_info$name,
    trait
  )
  SKM::mkdir(trait_results_dir)

  AllPredictions <- foreach(
    fold = folds,
    .combine = rbind,
    .packages = "dplyr"
  ) %dopar% {
    SKM::echo("\t-Fold")

    y_na <- replace(y, fold$testing, NA)

    model <- BGLR::BGLR(
      y = y_na,
      ETA = ETA,
      response_type = "gaussian",
      verbose = FALSE,
      saveAt = tempdir(),
      nIter = iterations_number,
      burnIn = burn_in
    )

    predictions <- model$yHat[fold$testing]

    # Adjust the predictions
    AdjustedPredictions <- Pheno %>%
      slice(fold$testing) %>%
      mutate(Predicted = predictions) %>%
      adjust_predictions(fold$testing_env)

    # Merge the adjusted predictions with the original pheno
    FoldPredictions <- OriginalPheno %>%
      filter(Env == fold$testing_env) %>%
      rename(Observed = trait) %>%
      merge(AdjustedPredictions, by = "Line") %>%
      select(Line, Env, Observed, Predicted)

    FoldPredictions
  }

  EnvSummary <- AllPredictions %>%
    mutate(Fold = Env) %>%
    SKM::gs_summaries()
  EnvSummary <- EnvSummary$env %>%
    mutate(
      Type = type,
      Predictor = predictors_text,
      Dataset = data_info$name,
      Trait = trait
    ) %>%
    relocate(Type, Predictor, Dataset, Trait, 1)

  Global <- EnvSummary %>%
    slice(n()) %>%
    select(-Env)

  SKM::write_csv(EnvSummary, file.path(trait_results_dir, "env_summary.csv"))
  SKM::write_csv(Global, file.path(trait_results_dir, "global_summary.csv"))
  SKM::write_csv(
    AllPredictions,
    file.path(trait_results_dir, "predictions.csv")
  )
}

stopCluster(cluster)
