## Changes from Version 0.1.0 to 0.1.1

* LassoShooting.fit: handle of NAs in coefficient vector (Bug 1)
* rlasso: removed bug in construction of model.matrix and change in predict-method (analog for rlassologit)
* rlassoIVselectX / rlassoIV: change in man pages (exogenous variables are automatically used as instruments)
* rlassologit.fit: remove "double scaling"
* predict.rlassologit: force as data.frame
* 22/03/2016: added significance test, R^2, adjusted R^2 to summary.rlasso (description in vignette)
* 22/03/2016: vignette: changes in treatment effect example
* 22/03/2016: rlasso with option model
* 03/05/2016: revise model.matrix and predict for rlasso and rlassologit
* 03/05/2016: change intercept handling for rlasso and rlassologit
* 03/05/2016: add Citation
* 03/05/2016: change plot function (test!)
* 03/05/2016: change of penalty choice (mm=1 half penalty and change old rule)