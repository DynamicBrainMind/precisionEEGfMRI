# ---- Package Load Test ----

packages <- c(
  "arrow",
  "ggridges",
  "tidyverse",
  "eegUtils",
  "glue",
  "data.table",
  "ggpubr",
  "here",
  "fs",
  "psych",
  "scales",
  "RColorBrewer",
  "reticulate"
)

failed <- list()

for (pkg in packages) {
  message("\n--- Checking: ", pkg, " ---")
  tryCatch(
    {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      message("✅ Successfully loaded ", pkg)
    },
    error = function(e) {
      message("❌ Failed to load ", pkg)
      failed[[pkg]] <<- e  # record the error
    }
  )
}

if (length(failed) == 0) {
  message("\n🎉 All packages loaded successfully!")
} else {
  message("\n⚠️ The following packages failed to load:\n")
  print(names(failed))
  message("\nTracebacks for each failed package:\n")
  for (pkg in names(failed)) {
    cat("\n--- ", pkg, " ---\n", sep = "")
    print(failed[[pkg]])
  }
  quit(status = 1)
}

