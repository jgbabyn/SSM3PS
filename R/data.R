#'Data readin and used
#'
#' Data required to run SPAM. Includes surveys, a list containing all surveys and catch for use with datSetup.
#' fstock_wt, which are FPCA smoothed backdated weights, fcomm_wt which are fpca smoothed commerical weights,
#' midy_wt which was the provided midy_wt file,mat which was the provided maturity file,stock_wt original provided stock wt,
#' DFO_stock_wt,catch data on it's own, geac and french survey indices,dfo survey indices
#'
#'
"rDat"

#'Result of running the model 100 times with every possible user configurable option
#'
#' This is a table containing the AIC,BIC, time taken, likelihood, convergence message of
#' running SPAM 100 times on a specfic set parameters, this was done with every possible option
#' of correlation type for F and CRL along with the different recruitment curve options
"modelOpts"

#'The top chosen model, residuals, data, etc.
#'
"topMod"

#'The final chosen model
#'
"finMod"

#'Retrospective plots and retro stuff
#'
"retL"
