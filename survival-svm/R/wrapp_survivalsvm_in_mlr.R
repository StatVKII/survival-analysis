library(survivalsvm)
library(mlr)
library(survival)
#--- creation of an mlr leaner for survivalsvm
makeRLearner.surv.survivalsvm = function() {
  makeRLearnerSurv(
    cl = "surv.survivalsvm",
    package = "survivalsvm",
    par.set = makeParamSet(
      makeDiscreteLearnerParam(id = "type",  default = "regression",
                               values = c("regression", "vanbelle1", "vanbelle2",
                                          "hybrid")),
      makeDiscreteLearnerParam(id = "diff.meth",  default = "makediff3",
                               values = c("makediff1", "makediff2", "makediff3")),
      makeNumericVectorLearnerParam(id = "gamma.mu", 
                                    tunable = TRUE, lower = 2^-5, upper = 2^5),
      makeDiscreteLearnerParam(id = "opt.meth", default = "quadprog",
                               values = c("quadprog", "ipop")),
      makeDiscreteLearnerParam(id = "kernel", default = "lin_kernel",
                               values = c("lin_kernel", "add_kernel", "rbf_kernel",
                                          "rbf4_kernel", "poly_kernel")),
      makeNumericLearnerParam(id = "kernel.pars", tunable = TRUE),
      makeNumericLearnerParam(id = "sgf.sv", default = 5, tunable = FALSE),
      makeNumericLearnerParam(id = "sigf", default = 7, tunable = FALSE),
      makeNumericLearnerParam(id = "maxiter", default = 20, tunable = FALSE),
      makeNumericLearnerParam(id = "margin", default = 0.05, tunable = FALSE),
      makeNumericLearnerParam(id = "bound", default = 10, tunable = FALSE),
      makeNumericLearnerParam(id = "eig.tol", default = 1e-06,
                              tunable = FALSE),
      makeNumericLearnerParam(id = "conv.tol", default = 1e-07,
                              tunable = FALSE),
      makeNumericLearnerParam(id = "posd.tol", default = 1e-08,
                              tunable = FALSE)
    ),
    properties = c("missings", "numerics", "factors", "weights", "prob",
                   "rcens"),
    name = "survival support vector machines",
    short.name = "survivalsvm",
    note = "survivalsvm in mlr"
  )
}

#--- creation of trainer for survivalsvm
trainLearner.surv.survivalsvm = function(.learner, .task, .subset, ...) {
  f <-  getTaskFormula(.task)
  data  <-  getTaskData(.task, subset = .subset)
  mod <-  survivalsvm::survivalsvm(formula = f, data = data, ...)
  return(mod)
}

#--- creation of predictor for survivalsvm
predictLearner.surv.survivalsvm = function(.learner, .model, .newdata, ...) {
  if (.learner$predict.type == "response") {
    predict(object = .model$learner.model, newdata = .newdata, ...)$predicted[1,]
  }
}


## Application on veteran data set
data(veteran, package = "survival")

## Some preprocessing steps
veteran[, "prior"] <- as.factor(veteran$prior == 10)
veteran[, "trt"] <- as.factor(veteran$trt == 2)

# Split data in train and test data sets
set.seed(123)
n <- nrow(veteran)
train_index <- sample(1:n, 120, replace = FALSE)
test_index <- setdiff(1:n, train_index)

## Create a survival task
veteran_task <- makeSurvTask(data = veteran[train_index, ],
                             target = c("time", "status"))

## Classicaly define am mlr based learner
survivalsvm_lrn <- makeLearner(cl = "surv.survivalsvm",
                               par.vals = list(type = "regression",
                                               gamma.mu = 1,
                                               opt.meth = "quadprog",
                                               kernel = "add_kernel"))
model_survivalsvm <- train(survivalsvm_lrn, veteran_task)

## Let us do predictions for new test data set
predict_survivalsvm <- predictLearner(survivalsvm_lrn, model_survivalsvm,
                               veteran[test_index, ])


