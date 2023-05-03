# Utils for all models --------------------------------------------------- START
dummy_matrix <- function(values, remove_level = NULL) {
  data <- data.frame(Value = values)
  Matrix <- model.matrix(~0 + Value, data)
  colnames(Matrix) <- gsub("Value", "", colnames(Matrix))

  if (!is.null(remove_level)) {
    if (is.numeric(remove_level)) {
      remove_col <- remove_level
    } else {
      remove_col <- which(colnames(Matrix) == remove_level)
    }

    if (SKM::is_empty(remove_col)) {
      warning(remove_level, " does not exist in the provided values")
    } else {
      Matrix <- Matrix[, -remove_col, drop = FALSE]
    }
  }

  return(Matrix)
}

prepare_eta <- function(Env, GenoLine, LinexEnv, predictors) {
  predictors <- tolower(predictors)
  ETA <- list()

  if ("env" %in% predictors) {
    ETA$Env <- list(X = Env, model = "FIXED")
  }

  if ("line" %in% predictors) {
    ETA$Line <- list(K = GenoLine, model = "RKHS")
  }

  if ("linexenv" %in% predictors) {
    ETA$LinexEnv <- list(K = LinexEnv, model = "RKHS")
  }

  return(ETA)
}
# Utils for all models ------------------------------------------------------END

# Non stationary model --------------------------------------------------- START
# Given a vector and a number of elements, this function returns a list with
# all the possible combinations of the given number of elements.
combinations <- function(values, elements_num) {
  combs <- combn(values, elements_num)

  # convert combs to list
  combs <- lapply(1:ncol(combs), function(i) combs[, i])

  return(combs)
}

# Prepare pheno for the non-stationary model. This function computes all the
# possible combinations of environments taken by two elements and then
# calculates the difference between the responses of the two elements in each
# combination. This value will be used for training the model and adjust the
# predictions.
prepare_pheno <- function(Pheno, trait_name) {
  Pheno <- Pheno %>%
    select(all_of(c("Line", "Env", trait_name)))

  # Generate all possible combinations of environments taken by two elements
  envs_combinations <- combinations(unique(Pheno$Env), 2)

  FinalPheno <- data.frame()

  for (combination in envs_combinations) {
    env1 <- combination[1]
    env2 <- combination[2]

    PhenoEnv1 <- Pheno %>%
      filter(Env == env1) %>%
      rename(y1 = trait_name)

    PhenoEnv2 <- Pheno %>%
      filter(Env == env2) %>%
      rename(y2 = trait_name)

    # Merge the two data frames by Line and calculate the difference between the
    # two response variables
    PhenoMerged <- merge(
        PhenoEnv1,
        PhenoEnv2,
        by = "Line",
        suffixes = seq(2)
      ) %>%
      mutate(Diff = y1 - y2)

    FinalPheno <- rbind(FinalPheno, PhenoMerged)
  }

  # Order FinalPheno by Env1 and Line
  FinalPheno <- FinalPheno %>%
    arrange(Env1, Line)

  # Pheno final structure:
  #         Line Env1   y1 Env2   y2  Diff
  # 1  CKDHL0008  EBU 6.88  KAK 4.83  2.05
  # 2  CKDHL0008  EBU 6.88  KTI 4.84  2.04
  # 3  CKDHL0039  EBU 6.85  KAK 5.23  1.62

  return(FinalPheno)
}

# Get the folds for CV0. At each fold, one environment is used for testing, that
# is, all the records where the Env1 or Env2 has the testing environment, the
# rest of the records are used for training.
cv0 <- function(Pheno, all_envs) {
  folds <- list()
  all_indices <- seq(nrow(Pheno))

  for (testing_env in all_envs) {
    fold <- list(
      testing_env = testing_env,
      training = which(Pheno$Env1 != testing_env & Pheno$Env2 != testing_env)
    )

    # Indices that are not in the training set are in the testing set
    fold$testing <- setdiff(all_indices, fold$training)
    folds <- append(folds, list(fold))
  }

  return(folds)
}

# Pheno is the testing pheno with the Predicted column
adjust_predictions <- function(Pheno, testing_env, by_mean = FALSE) {
  if (by_mean) {
    env1_means <- Pheno %>%
      group_by(Env1) %>%
      summarise(y1 = mean(y1, na.rm = TRUE), y2 = mean(y2, na.rm = TRUE)) %>%
      as.data.frame()
    rownames(env1_means) <- env1_means$Env

    env2_means <- Pheno %>%
      group_by(Env2) %>%
      summarise(y1 = mean(y1, na.rm = TRUE), y2 = mean(y2, na.rm = TRUE)) %>%
      as.data.frame()
    rownames(env2_means) <- env2_means$Env
  }

  adjust_one <- function(row) {
    if (row$Env1 == testing_env) {
      # Env1 and Env2 refer to the position in the right hand side of the
      # following equations:
      # D12 = yE1 - yE2
      # D13 = yE1 - yE3
      # D23 = yE2 - yE3
      # Use D23 as training and D12 and D13 as tst sets.
      # Then  predict:
      # yE1 = average([yE2 + D12^hat] + [yE3 + D13^hat])

      # row$y1 = row$y2 + row$Predicted
      if (by_mean) {
        return(env2_means[row$Env2, "y2"] + row$Predicted)
      } else {
        return(row$y2 + row$Predicted)
      }
    } else {
      # Env1 and Env2 refer to the position in the right hand side of the
      # following equations:
      # Step 1: Calculate the difference
      # D12 = yE1 - yE2
      # D13 = yE1 - yE3
      # D23 = yE2 - yE3
      # Use D12 as training and D13 and D23 as tst sets.
      # Then  predict:
      # yE3 = average([yE1 - D13^hat] + [yE2 - D23^hat])

      # row$y2 = row$y1 - row$Predicted
      if (by_mean) {
        return(env1_means[row$Env1, "y1"] - row$Predicted)
      } else {
        return(row$y1 - row$Predicted)
      }
    }
  }

  # Adjust predicted value row by row
  for (i in seq(nrow(Pheno))) {
    row <- Pheno[i, ]
    Pheno[i, "Predicted"] <- adjust_one(row)
  }

  # In testing Pheno there are all the combinations of environments for the
  # testing environment. We need to average the predictions for each line.
  Pheno <- Pheno %>%
    group_by(Line) %>%
    summarise(Predicted = mean(Predicted, na.rm = TRUE))

  return(Pheno)
}
# Non stationary model ----------------------------------------------------- END

# LEOEO and ONE LEOEO ---------------------------------------------------- START

cv_loeo <- function(envs) {
  unique_envs <- unique(envs)
  all_indices <- seq_along(envs)
  folds <- list()

  for (env in unique_envs) {
    fold <- list(testing = which(envs == env))
    fold$training <- all_indices[-fold$testing]

    fold$testing_env <- env
    fold$training_env <- "all"

    folds <- append(folds, list(fold))
  }

  class(folds) <- "cv_loeo"

  return(folds)
}

print.cv_loeo <- function(folds, envs) {
  i <- 1
  for (fold in folds) {
    SKM::echo("*** Fold %s***", i)
    SKM::echo("\t+ Testing: %s", fold$testing_env)
    print(table(envs[fold$testing]))

    SKM::echo("\t+ Trainig: %s", fold$training_env)
    print(table(envs[fold$training]))

    i <- i + 1
  }
}
# LEOEO and ONE LEOEO ------------------------------------------------------ END
