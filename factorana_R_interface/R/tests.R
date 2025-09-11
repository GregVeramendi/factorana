# Load your package / source files
# devtools::load_all(".")  # if in a package
# source("R/define_model_component.R")
# source("R/estimate_model.R")
# source("R/define_model_system.R")
# source("R/define_estimation_control.R")

set.seed(123)

# Fake data
df <- data.frame(
  y = rnorm(50),
  x1 = rnorm(50),
  x2 = rnorm(50)
)

# Fake factor model (simplified placeholder)
mock_factor <- list(nfactors = 1, n_types = 2)
class(mock_factor) <- "factor_model"

# Test 1: define_model_component
cat("---- Testing define_model_component ----\n")
mc1 <- define_model_component(
  name = "linear_comp",
  data = df,
  outcome = "y",
  factor = mock_factor,
  covariates = c("x1", "x2"),
  model_type = "linear"
)
print(mc1)

# Test 2: define_model_system
cat("\n---- Testing define_model_system ----\n")
ms <- define_model_system(
  factor = mock_factor,
  components = list(mc1)
)
print(ms)

# Test 3: estimation_control
cat("\n---- Testing define_estimation_control ----\n")
ctrl <- define_estimation_control(num_cores = 1)
print(ctrl)

# Test 4: estimate_model (initialization only)
cat("\n---- Testing estimate_model ----\n")
estimate_model(ms, ctrl)
