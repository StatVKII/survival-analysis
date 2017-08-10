require(Hmisc)
## First create our on performance measure to further analysis
## Let us define a c-index function as performance measure
#' my.ci.fun computes the C-index for "survivalsvmprediction" objects.
#'
#' @param task survival task used to fit the model.
#'        This will be later used by "makeMeasure" from the pakage "mlr".
#' @param model model of class "survivalsvm".
#'        This will be later used by "makeMeasure" from the pakage "mlr".
#' @param pred The predictions.
#' @param feats will be later used by "makeMeasure" from the pakage "mlr".
#' @param extra.args will be later used by "makeMeasure" from the pakage "mlr".
#'
#' @return The C-Index
#'
my.ci.fun <- function(task, model, pred, feats, extra.args) {
  myci = rcorr.cens(x = getPredictionResponse(pred),
                    S = getPredictionTruth(pred))
  return(myci["C Index"])
}
# constructs the C-index measure for survivalsvmprediction objects.
c.i <- makeMeasure(
  id = "ci", name = "C-Index",
  properties = c("surv"),
  minimize = FALSE, best = 1, worst = 0,
  fun = my.ci.fun
)

#--- Wrapper for scale data when required
makePreprocWrapperScale = function(learner, center = TRUE, scale = TRUE) {
  trainfun = function(data, target, args = list(center, scale)) {
    cns = colnames(data)
    nums = setdiff(cns[sapply(data, is.numeric)], target)
    x = as.matrix(data[, nums, drop = FALSE])
    x = scale(x, center = args$center, scale = args$scale)
    control = args
    if (is.logical(control$center) && control$center)
      control$center = attr(x, "scaled:center")
    if (is.logical(control$scale) && control$scale)
      control$scale = attr(x, "scaled:scale")
    data = data[, setdiff(cns, nums), drop = FALSE]
    data = cbind(data, as.data.frame(x))
    return(list(data = data, control = control))
  }
  predictfun = function(data, target, args, control) {
    cns = colnames(data)
    nums = cns[sapply(data, is.numeric)]
    x = as.matrix(data[, nums, drop = FALSE])
    x = scale(x, center = control$center, scale = control$scale)
    data = data[, setdiff(cns, nums), drop = FALSE]
    data = cbind(data, as.data.frame(x))
    return(data)
  }
  if(!("surv.glmboost" %in% class(learner)))
    makePreprocWrapper(
      learner,
      train = trainfun,
      predict = predictfun,
      par.set = makeParamSet(
        makeLogicalLearnerParam("center"),
        makeLogicalLearnerParam("scale")
      ),
      par.vals = list(center = center, scale = scale)
    )
  else
    makePreprocWrapper(
      learner,
      train = trainfun,
      predict = predictfun,
      par.set = makeParamSet(),
      par.vals = list(center = center, scale = scale)
    )
}
#--------- TuneWrapper for survivalsvm ---------------------
tuneWrapperGammaMu <- function(type, kernel, ...,
                               preproc = TRUE, tune.scale = TRUE, center = TRUE,
                               scale = TRUE, resolution = 10L, method = "CV",
                               lower = -10L, upper = 10L, iters.rep = 9L) {
  #task <- makeSurvTask(data = data, target = target)
  lrn <- makePreprocWrapperScale(makeLearner(cl = "surv.survivalsvm",
                                             type = type, kernel = kernel,
                                             ...),
                                 center = center, scale = scale)
  configureMlr(on.learner.error = "warn")
  if(kernel %in% c("lin_kernel", "add_kernel")) {
    # Parameters set for linear and additiv kernel
    if (type != "hybrid") {
      discrete_ps <- if(tune.scale) {
        makeParamSet(
          makeDiscreteParam("gamma.mu", values = 2^(lower:upper)),
          makeLogicalLearnerParam("center"),
          makeLogicalLearnerParam("scale")
        )
      } else {
        makeParamSet(
          makeDiscreteParam("gamma.mu",
                            values = 2^(lower:upper))
        )
      }
    } else {
      discrete_ps <- if(tune.scale) {
        makeParamSet(
          makeNumericVectorParam("gamma.mu", len = 2L, lower = lower,
                                 upper = upper,
                                 trafo = function(x) 2^x),
          makeLogicalLearnerParam("center"),
          makeLogicalLearnerParam("scale")
        )
      } else {
        makeParamSet(
          makeNumericVectorParam("gamma.mu", len = 2L, lower = lower,
                                 upper = upper,
                                 trafo = function(x) 2^x)
        )
      }
    }
  } else {
    # Parameters set for RBF kernel
    if (type != "hybrid") {
      discrete_ps <- if(tune.scale) {
        makeParamSet(
          makeDiscreteParam("gamma.mu", values = 2^(lower:upper)),
          makeDiscreteParam("kernel.pars", values = 2^(-5:5)),
          makeLogicalLearnerParam("center"),
          makeLogicalLearnerParam("scale")
        )
      } else {
        makeParamSet(
          makeDiscreteParam("gamma.mu",
                            values = 2^(lower:upper)),
          makeDiscreteParam("kernel.pars", values = 2^(-5:5))
        )
      }
    } else {
      discrete_ps <- if(tune.scale) {
        makeParamSet(
          makeNumericVectorParam("gamma.mu", len = 2L, lower = lower,
                                 upper = upper,
                                 trafo = function(x) 2^x),
          makeLogicalLearnerParam("center"),
          makeLogicalLearnerParam("scale"),
          makeDiscreteParam("kernel.pars", values = 2^(-5:5))
        )
      } else {
        makeParamSet(
          makeNumericVectorParam("gamma.mu", len = 2L, lower = lower,
                                 upper = upper,
                                 trafo = function(x) 2^x),
          makeDiscreteParam("kernel.pars", values = 2^(-5:5))
        )
      }
    }
    
  }
  
  
  ctrl  <-  makeTuneControlGrid(resolution = resolution)
  inner  <-  makeResampleDesc(method = method, iters = iters.rep)
  survivalsvm.tuned  <-  makeTuneWrapper(lrn,  resampling = inner,
                                         par.set = discrete_ps, control = ctrl,
                                         measures = c.i)
  return(survivalsvm.tuned)
}


## Set up parameter for a nested cross validation mlr
outer <- makeResampleDesc("CV", iters = 5L)

#---- survivalsvm regression.tuned--------------------------------------------
set.seed(123)
tunwrp.gm.reg <- tuneWrapperGammaMu(type = "regression",
                                    kernel = "lin_kernel",
                                    opt.meth = "quadprog", method = "CV",
                                    iters.rep = 5)
bench.vet.ssvm.reg <- benchmark(learners = tunwrp.gm.reg, tasks = veteran_task,
                                resamplings = outer, measures = list(c.i))
save(bench.vet.ssvm.reg, file = "benchVetSsvmReg.rda")
lapply(bench.vet.ssvm.reg$results$veteran.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })