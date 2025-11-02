#!/usr/bin/env Rscript
# Run all manual tests

cat("==========================================\n")
cat("Running All Manual Tests\n")
cat("==========================================\n\n")

# Get all test files
test_files <- list.files(
  path = "tests/manual",
  pattern = "^test_.*\\.R$",
  full.names = TRUE
)

# Track results
results <- data.frame(
  test = character(),
  status = character(),
  time = numeric(),
  stringsAsFactors = FALSE
)

# Run each test
for (test_file in test_files) {
  test_name <- basename(test_file)
  cat(sprintf("\n=== Running %s ===\n", test_name))

  start_time <- Sys.time()
  status <- tryCatch({
    source(test_file, local = TRUE)
    "PASS"
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    "FAIL"
  })
  end_time <- Sys.time()

  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

  results <- rbind(results, data.frame(
    test = test_name,
    status = status,
    time = elapsed,
    stringsAsFactors = FALSE
  ))
}

# Print summary
cat("\n\n==========================================\n")
cat("Test Summary\n")
cat("==========================================\n\n")

print(results, row.names = FALSE)

cat(sprintf("\nTotal tests: %d\n", nrow(results)))
cat(sprintf("Passed: %d\n", sum(results$status == "PASS")))
cat(sprintf("Failed: %d\n", sum(results$status == "FAIL")))
cat(sprintf("Total time: %.1f seconds\n", sum(results$time)))

# Exit with appropriate code
if (any(results$status == "FAIL")) {
  quit(status = 1)
} else {
  cat("\n✓ All tests passed!\n")
  quit(status = 0)
}
