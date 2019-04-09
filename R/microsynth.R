#' @title
#' Synthetic control methods for disaggregated, micro-level data.
#'
#' @description
#' Implements the synthetic control method for micro-level data as outlined in
#' Robbins, Saunders, and Kilmer (2017).  \code{microsynth} is designed for use
#' in assessment of the effect of an intervention using longitudinal data.
#' However, it may also be used to calculate propensity score-type weights in
#' cross-sectional data. \code{microsynth} is a generalization
#' of \code{Synth} (see Abadie and Gardeazabal (2003) and Abadie, Diamond,
#' Hainmueller (2010, 2011, 2014)) that is designed for data at a more granular
#' level (e.g., micro-level). For more details see the help vignette:
#' \code{vignette('microsynth', package = 'microsynth')}.
#'
#' \code{microsynth} develops a synthetic control group by searching for weights
#' that exactly match a treatment group to a synthetic control group across
#' a number of variables while also minimizing the discrepancy between the
#' synthetic control group and the treatment group across a set second set of
#' variables.  \code{microsynth} works in two primary steps: 1) calculation of
#' weights and 2) calculation of results.  Time series plots of treatment
#' vs. synthetic control for pertinent outcomes may be performed using the
#' function \code{plot.microsynth()}.
#'
#' The time range over which data are observed is segmented into pre- and
#' post-intervention periods.  Treatment is matched to synthetic control
#' across the pre-intervention period, and the effect of the intervention
#' is assessed across the post-intervention (or evaluation) period.  The input
#' \code{end.pre} (which gives the last pre-intervention time period) is used to
#' delineate between pre- and post-intervention.  Note that if the intervention
#' is not believed to have an instantaneous effect, \code{end.pre} should indicate
#' the time of the intervention.
#'
#' Variables are categorized as outcomes (which are time-variant) and covariates
#' (which are time-invariant).  Using the respective inputs \code{match.covar}
#' and \code{match.out}, the user specifies across which covariates and outcomes
#' (and which pre-intervention time points of the outcomes) treatment is to be
#' exactly matched to synthetic control.  The inputs \code{match.covar.min} and
#' \code{match.out.min} are similar but instead specify variables across which
#' treatment is to be matched to synthetic control as closely as possible.  If
#' there are no variables specified in \code{match.covar.min} and
#' \code{match.out.min}, the function \code{calibrate()} from the \code{survey}
#' package is used to calculate weights.  Otherwise, the function
#' \code{LowRankQP()} from the package of the same name is used.  In the event
#' that the model specified by \code{match.covar} and \code{match.out} is not
#' feasible (i.e., weights do not exist that exactly match treatment and
#' synthetic control subject to the given constraints), a less restrictive
#' backup model is used.
#'
#' \code{microsynth} has the capability to perform
#' statistical inference using Taylor series linearization, a jackknife and
#' permutation methods.  Several sets of weights are calculated.  A set of main
#' weights is calculated that is used to determine a point estimate of the
#' intervention effect.  The main weights can also be used to perform inferences
#' on the point estimator via Taylor series linearization.  If a jackknife is to
#' be used, one set of weights is calculated for each jackknife replication
#' group, and if permutation methods are to be used, one set of weights is
#' calculated for each permutation group.  If treatment and synthetic control
#' are not easily matched based upon the model outlined in \code{match.covar}
#' and \code{match.out} (i.e., an exact solution is infeasible or nearly
#' infeasible), it is recommended that the jackknife not be used for inference.
#'
#' The software provides the user the option to output overall findings in an Excel
#' file.  For each outcome variable, the results list the estimated treatment
#' effect, as well as confidence intervals of the effect and p-values of a
#' hypothesis test that assesses whether the effect is zero.   Such results are
#' produced as needed for each of the three methods of statistical inference
#' noted above.  \code{microsynth} can also apply an omnibus test that examines
#' the presence of a treatment effect jointly across several outcomes.
#'
#' @details
#' \code{microsynth} requires specification of the following inputs:
#' \code{data}, \code{idvar}, \code{intvar}.  \code{data} is a longitudinal data
#' frame; \code{idvar} and \code{intvar} are character strings that specific
#' pertinent columns of \code{data}.  In longitudinal data, \code{timevar}
#' should be specified.  Furthermore, specification of \code{match.out} and
#' \code{match.covar} is recommended.
#'
#' \code{microsynth} can also be used to calculate propensity score-type weights
#' in cross sectional data (in which case \code{timevar} does not need to be
#' specified) as proposed by Hainmueller (2012).
#'
#' @param data A data frame.  If longitudinal, the data must be entered in tall
#'   format (e.g., at the case/time-level with one row for each time period for
#'   each case).  Missingness is not allowed.  All individuals must have non-NA
#'   values of all variables at all time points.
#'
#' @param idvar A character string that gives the variable in \code{data} that
#'   identifies multiple records from the same case.
#'
#' @param intvar A character string that gives the variable in \code{data} that
#'   corresponds to the intervention variable.  The intervention variable
#'   indicates which cases (and times) have received the intervention.  The
#'   variable should be binary, with a 1 indicating treated and 0 indicating
#'   untreated.  If \code{end.pre}
#'   is specified, a case is considered treated if there is 1 or more non-zero
#'   entries in the column indicated by \code{intvar} for that case (at any time
#'   point).  If \code{end.pre} is not specified, an attempt will be made to
#'   use \code{intvar} to determine which time periods will be considered
#'   post-intervention (i.e., the times contained in the evaluation period).
#'   In this case, the evaluation period is considered to begin at the time of
#'   the first non-zero entry in \code{intvar}).
#'
#' @param timevar A character string that gives the variable in
#'   \code{data} that differentiates multiple records from the same case.  Can be
#'   set to \code{NULL} only when used with cross-sectional data (i.e., with one
#'   observation per entry in \code{idvar}).
#'
#' @param w A \code{microsynth} object or a list of the form as returned
#'   by a prior application of \code{microsynth}.
#'   If \code{w = NULL}, weights are calculated from scratch.
#'   Entering a \code{non-NULL} value affords the user the ability to
#'   use previously calculated  weights.
#'
#' @param end.pre An integer that gives the final time point of the
#'   pre-intervention period.  That is, \code{end.pre} is the last time at
#'   which treatment and synthetic control will be matched to one another.
#'   All time points
#'   following \code{end.pre} are considered to be post-intervention and the
#'   behavior of outcomes will be compared between the treatment and synthetic
#'   control groups across those time periods.
#'   Setting \code{end.pre = NULL} will begin the post-intervention period
#'   at the time
#'   that corresponds to the first non-zero entry in the column indicated by
#'   \code{intvar}.
#'
#' @param end.post An integer that gives the maximum post-intervention time that
#'   is taken into when compiling results.  That is, the treatment and synthetic
#'   control groups are compared across the outcomes listed in \code{result.var}
#'   from the first time following the intervention up to \code{end.post}.  Can
#'   be a vector (ordered, increasing) giving multiple values of
#'   \code{end.post}.  In this case, the results will be compiled for each entry
#'   in \code{end.post}.  When \code{end.post = NULL} (the default), it is reset
#'   to the maximum time that appears in the column given by \code{timevar}.
#'
#' @param match.out Either A) logical, B) a vector of variable names that
#'   indicates across which time-varying variables treatment is to be exactly matched
#'   to synthetic control pre-intervention, or C) a
#'   list consisting of variable names and timespans over which variables should
#'   be aggregated before matching.  Note that outcome variables and time-varying
#'   covariates should be included in \code{match.out}.
#'
#'   If \code{match.out = TRUE} (the default), it is set equal to
#'   \code{result.var}; if \code{match.out = NULL} or \code{match.out = FALSE},
#'   no outcome variables are factored into the calculation of weights. If
#'   \code{match.out} is passed a vector of variable
#'   names, then weights are calculated to match treatment and synthetic control
#'   for the value of each variable that appears in \code{match.out} at each
#'   time point from \code{start.pre} to \code{end.pre}. Otherwise, to allow
#'   more flexibility, \code{match.out} may also be a list that gives
#'   an outcome-based model outlining more specific constraints that are to be
#'   exactly satisfied within calibration weighting.  In this case, each entry
#'   of \code{match.out} is a vector of integers, and the names of entries of
#'   \code{match.out} are the outcome variables to which the vectors correspond.
#'   Each element of the vectors gives a number of time points that are to be
#'   aggregated for the respective outcome, with the first element indicating
#'   time points immediately prior the beginning of the post-intervention
#'   period.  The sum of
#'   the elements in each vector should not exceed the number of
#'   pre-intervention time periods in the data.
#'
#'   The following examples show the proper formatting of \code{match.out} as a
#'   list.  Assume that there are two outcomes, Y1 and Y2 (across which
#'   treatment is to be matched to synthetic control), and \code{end.pre = 10}
#'   (i.e., the post-intervention period begins at time 11).
#'   Let \code{match.out = list('Y1' = c(1, 3, 3), 'Y2'=
#'   c(2,5,1))}.  According to this specification, treatment is to be matched to
#'   synthetic control across: a) The value of Y1 at time 10; b) the sum of Y1
#'   across times 7, 8 and 9; c) the sum of Y1 across times 4, 5 and 6; e) The
#'   sum of Y2 across times time 9 and 10; e) the sum of Y2 across times 4, 5,
#'   6, 7, and 8; f) the value of Y2 at time 3.  Likewise, if \code{match.out =
#'   list('Y1' = 10, 'Y2'= rep(1,10))}, Y1 is matched to synthetic control
#'   the entire aggregated pre-intervention time range, and Y2 is matched
#'   at each pre-intervention time point individually.
#'
#' @param match.out.min A vector or list of the same format as \code{match.out}
#'   that is used to specify additional time-varying variables to match
#'   on, but which need not be matched exactly. Weights are calculated so the
#'   distance is minimized between treatment and synthetic control across these
#'   variables. If \code{match.out.min = NULL}, no outcome-based constraints
#'   beyond those indicated by \code{match.out} are imposed (i.e., all outcome
#'   variables will be matched on exactly).
#'
#' @param match.covar Either a logical or a vector of variable names that
#'   indicates which time invariant covariates
#'   are to be used for weighting.  Weights are
#'   calculated so that treatment and synthetic control exactly match across
#'   these variables.  If \code{match.covar = TRUE}, it is set equal to a vector
#'   of variable names corresponding to the time invariant variables that
#'   appear in \code{data}.  If \code{match.covar = FALSE}, it is set to
#'   \code{NULL} (in which case no time-invariant variables are used for matching
#'   when calculating weights).
#'
#' @param match.covar.min A vector of variable names that indicates supplemental
#'   time invariant variables that are to be used for weighting, for which exact
#'   matches are not required. Weights are calculated so the distance is
#'   minimized between treatment and synthetic control across these variables.
#'
#' @param perm An integer giving the number of permutation groups that are used.
#'   If \code{perm = 0}, no permutation groups are generated, permutation
#'   weights are not calculated, and permutations do not factor into the
#'   reported results.  \code{perm} is set to the number of possible permutation
#'   groups if the former exceeds the latter.
#'
#' @param jack An integer giving the number of replication groups that are used
#'   for the jackknife.  \code{jack} can also be a logical indicator.  If
#'   \code{jack = 0} or \code{jack = FALSE}, no jackknife replication groups are
#'   generated, jackknife weights are not calculated, and the jackknife is not
#'   considered when reporting results.  If \code{jack = TRUE}, it is reset to
#'   being equal to the minimum between the number of total cases in the
#'   treatment group and the total number of cases in the control group.
#'   \code{jack} is also reset to that minimum if it, as entered, exceeds that
#'   minimum.
#'
#' @param check.feas A logical indicator of whether or not the feasibility of
#'   the model specified by \code{match.out} is evaluated prior to calculation
#'   of weights.  If \code{check.feas = TRUE}, feasibility is assessed. If
#'   \code{match.out} is found to not specify a feasible model, a less
#'   restrictive feasible backup model will be applied to calculate the main
#'   weights and for jackknife and permutation methods.
#'
#' @param use.backup A logical variable that, when true, indicates whether a
#'   backup model should be used whenever the model specified by
#'   \code{match.out} yields unsatisfactory weights.  Weights are deemed to be
#'   unsatisfactory if they do not sufficiently satisfy the constraints imposed
#'   by \code{match.out} and \code{match.covar}.  Different backup models may be
#'   used for each of the main, jackknife or permutation weights as needed.
#'
#' @param max.mse The maximum error (given as mean-squared error) permissible
#'   for constraints that are to be exactly satisfied.  If \code{max.mse} is not
#'   satisfied by these constraints, and either \code{check.feas = TRUE} or
#'   \code{use.backup = TRUE}, then back-up models are used.
#'
#' @param result.var A vector of variable names giving the outcome
#'   variables for which results will be reported.  Time-varying covariates
#'   should be excluded from \code{result.var}.  If \code{result.var = TRUE}
#'   (the default), \code{result.var} is set as being equal to all time-varying
#'   variables that appear in \code{data}.  If \code{result.var = NULL} or
#'   \code{result.var = FALSE}, results are not tabulated.
#'
#' @param omnibus.var A vector of variable names that indicates the outcome
#'   variables that are to be used within the calculation of the omnibus
#'   statistic.  Can also be a logical indicator.  When \code{omnibus.var =
#'   TRUE}, it is reset as being equal to \code{result.var}.  When
#'   \code{omnibus.var = NULL} or \code{omnibus = FALSE}, no omnibus statistic
#'   is calculated.
#'
#' @param period An integer that gives the granularity of the data that will be
#'   used for plotting and compiling results.  If \code{match.out} and
#'   \code{match.out.min} are provided a vector of variable names, it will also
#'   affect the calculation of weights used for matching. In this case, matching
#'   of treatment and synthetic control is performed at a temporal granularity
#'   defined by \code{period}. For instance, if monthly data are provided and
#'   \code{period = 3}, data are aggregated to quarters for plots and results
#'   (and weighting unless otherwise specified). If \code{match.out} and
#'   \code{match.out.min} are provided a list, \code{period} only affects plots
#'   and how results are displayed.
#'
#'   Note that plotting is performed with
#'   \code{plot.microsynth()}; however, a \code{microsynth} object is required as
#'   input for that function and \code{period} should be specified in the creation
#'   of that object.
#'
#' @param cut.mse The maximum error (given as mean-squared error) permissible
#'   for permutation groups.  Permutation groups with a larger than permissible
#'   error are dropped when calculating results.  The mean-squared error is only
#'   calculated over constraints that are to be exactly satisfied.
#'
#' @param test The type of hypothesis test (one-sided lower, one-sided upper, or
#'   two-sided) that is used when calculating p-values.  Entries of
#'   \code{'lower'}, \code{'upper'}, and \code{'twosided'} are recognized.
#'
#' @param result.file A character string giving the name of a file that will be
#'   created in the home directory containing results.  If \code{result.file =
#'   NULL} (the default), no file is created.  If \code{end.post} has length 1,
#'   a \code{.csv} file is created.  If \code{end.post} has length greater than
#'   one, a formatted \code{.xlsx} file is created with one tab for each element
#'   of \code{end.post}.  If \code{result.file} has a \code{.xlsx} (or
#'   \code{.xls}) extension (e.g., the last five characters of result.file are
#'   '.xlsx'), an \code{.xlsx} file is created regardless of the length of
#'   \code{end.post}.
#'
#' @param use.survey If \code{use.survey = TRUE}, Taylor series linearization is
#'   applied to the estimated treatment effect within each permutation group.
#'   Setting \code{use.survey = TRUE} makes for better inference but increases
#'   computation time substantially.  Confidence intervals for permutation
#'   groups are calculated only when \code{use.survey = TRUE}.
#'
#' @param confidence The level of confidence for confidence intervals.
#'
#' @param start.pre An integer indicating the time point that corresponds to the
#'   beginning of the pre-intervention period used for
#'   matching.  When \code{start.pre = NULL} (default), it is reset to the
#'   minimum time appearing in the column given by \code{timevar}.  If
#'   \code{match.out} (and \code{match.out.min}) are given in list format,
#'   \code{start.pre} is ignored except for plotting.
#'
#' @param scale.var  A variable name.  When comparing the treatment group to all
#'   cases, the latter is scaled to the size of the former with respect to the
#'   variable indicated by \code{scale.var}.  Defaults to the number of units
#'   receiving treatment (i.e., the intercept).
#'
#' @param maxit The maximum number of iterations used within the calibration
#'   routine (\code{calibrate()} from the \code{survey} package) for
#'   calculating weights.
#'
#' @param cal.epsilon The tolerance used within the calibration routine
#'   (\code{calibrate()} from the \code{survey} package) for calculating
#'   weights.
#'
#' @param calfun The calibration function used within the calibration routine
#'   (\code{calibrate()} from the \code{survey} package) for calculating
#'   weights.
#'
#' @param bounds Bounds for calibration weighting (fed into the
#'   \code{calibrate()} from the \code{survey} package).
#'
#' @param printFlag If TRUE, \code{microsynth} will print history on console. Use
#'   \code{printFlag = FALSE} for silent computation.
#'
#' @param n.cores The number of CPU cores to use for parallelization. If
#'   \code{n.cores} is not specified by the user, it is guessed using the
#'   \code{detectCores} function in the parallel package.  If \code{TRUE}
#'   (the default), it is set as \code{detectCores()}.  If \code{NULL}, it is set as
#'   \code{detectCores() - 1}.  If \code{FALSE}, it is set as \code{1}, in which case
#'   parallelization is not invoked.  Note that the
#'   documentation for \code{detectCores} makes clear that it is not failsafe and
#'   could return a spurious number of available cores.
#'
#' @details \code{microsynth} calculates weights using
#'   \code{survey::calibrate()} from the \code{survey} package in circumstances
#'   where a feasible solution exists for all constraints, whereas
#'   \code{LowRankQP::LowRankQP()} is used to assess feasibility and to
#'   calculate weights in the event that a feasible solution to all constraints
#'   does not exist.  The \code{LowRankQP} routine is memory-intensive and can
#'   run quite slowly in data that have a large number of cases.  To prevent
#'   \code{LowRankQP} from being used, set \code{match.out.min = NULL},
#'   \code{match.covar.min= NULL}, \code{check.feas = FALSE}, and
#'   \code{use.backup = FALSE}.
#'
#' @return \code{microsynth} returns a list with up to five elements: a)
#'   \code{w}, b) \code{Results}, c) \code{svyglm.stats},
#'   and d) \code{Plot.Stats}, and e) \code{info}.
#'
#'   \code{w} is a list with six elements: a) \code{Weights}, b)
#'   \code{Intervention},
#'   c) \code{MSE}, d) \code{Model}, e) \code{Summary}, and f) \code{keep.groups}.
#'   Assume there are
#'   C total sets of weights calculated, where C = \code{1 + jack + perm}, and
#'   there are N total cases across the treatment and control groups.
#'   \code{w$Weights} is an N x C matrix, where each column provides a set of
#'   weights.  \code{w$Intervention} is an N x C matrix made of logical
#'   indicators that indicate whether or not the case in the respective row is
#'   considered treated (at any point in time) for the respective column.
#'   Entries of \code{NA} are to be dropped for the respective jackknife
#'   replication group (\code{NA}s only appear in jackknife weights).
#'   \code{w$MSE} is a 6 x C matrix that give the MSEs for each set of weights.
#'   MSEs are listed for the primary and secondary constraints for the first,
#'   second, and third models.  \code{w$Model} is a length-C vector that
#'   indicates whether backup models were used in the calculation of each set of
#'   weights.  \code{w$keep.groups} is a logical vector indicating which groups
#'   are to be used in analysis (groups that are not used have pre-intervention
#'   MSE greater than \code{cut.mse}.  \code{w$Summary} is a three-column matrix
#'   that (for treatment,
#'   synthetic control, and the full dataset), shows aggregate values
#'   of the variables across which treatment and synthetic control are matched.
#'   The summary, which is tabulated only for the primary weights, is also
#'   printed by \code{microsynth} while weights are being calculated.
#'
#'   Further, \code{Results} is a list where each element gives the final
#'   results for each value of \code{end.post}.  Each element of \code{Results}
#'   is itself a matrix with each row corresponding to an outcome variable (and
#'   a row for the omnibus test, if used) and each column denotes estimates of
#'   the intervention effects and p-values, upper, and lower bounds of
#'   confidence intervals as found using Taylor series linearization (Linear),
#'   jackknife (jack), and permutation (perm) methods where needed.
#'
#'   In addition, \code{svyglm.stats} is a list where each element is a
#'   matrix that includes the output from the regression models run using the
#'   \code{svyglm()} function to estimate the treatment effect.  The list has one
#'   element for each value of \code{end.post}, and the matrices each have
#'   one row per variable in \code{result.var}.
#'
#'   Next, \code{Plot.Stats} contains the data that are displayed in the
#'   plots which may be generated using \code{plot.microsynth()}.
#'   \code{Plot.Stats} is a list with four elements (Treatment, Control,
#'   All, Difference).  The first three elements are matrices with one row per
#'   outcome variable and one column per time point.  The last element (which
#'   gives the treatment minus control values) is an array that contains data
#'   for each permutation group in addition to the true treatment area.
#'   Specifically, \code{Plot.Stats$Difference[,,1]} contains the time series of
#'   treatment minus control for the true intervention group;
#'   \code{Plot.Stats$Difference[,,i+1]} contains the time series of treatment
#'   minus control for the i^th permutation group.
#'
#' 	 Lastly, \code{info} documents some input parameters for display by
#'   \code{print()}. A summary of weighted matching variables and of results
#'   can be viewed using \code{\link{summary}}
#'
#' @references Abadie A, Diamond A, Hainmueller J (2010). Synthetic control
#'   methods for comparative case studies: Estimating the effect of California's
#'   tobacco control program.? \emph{Journal of the American Statistical
#'   Association}, 105(490), 493-505.
#'
#'   Abadie A, Diamond A, Hainmueller J (2011). Synth: An R Package for
#'   Synthetic Control Methods in Comparative Case Studies.? \emph{Journal
#'   of Statistical Software}, 42(13), 1-17.
#'
#'   Abadie A, Diamond A, Hainmueller J (2015). Comparative politics and the
#'   synthetic control method. \emph{American Journal of Political Science},
#'   59(2), 495-510.
#'
#'   Abadie A, Gardeazabal J (2003). The economic costs of conflict: A case
#'   study of the Basque Country.? \emph{American Economic Review}, pp. 113-132.
#'
#'   Hainmueller, J. (2012), Entropy Balancing for Causal Effects: A
#'   Multivariate Reweighting Method to Produce Balanced Samples in
#'   Observational Studies,? \emph{Political Analysis}, 20, 25-46.
#'
#'   Robbins MW, Saunders J, Kilmer B (2017). A framework for synthetic control
#'   methods with high-dimensional, micro-level data: Evaluating a neighborhood-
#'   specific crime intervention,? \emph{Journal of the American Statistical
#'   Association}, 112(517), 109-126.
#'
#' @examples
#'
#' # Use seattledmi, block-level panel data, to evaluate a crime intervention.
#'
#' # Declare time-variant (outcome) and time-invariant variables for matching
#' cov.var <- c('TotalPop', 'BLACK', 'HISPANIC', 'Males_1521',
#'        'HOUSEHOLDS', 'FAMILYHOUS', 'FEMALE_HOU', 'RENTER_HOU', 'VACANT_HOU')
#'
#' match.out <- c('i_felony', 'i_misdemea', 'i_drugs', 'any_crime')
#' set.seed(99199) # for reproducibility
#'
#'
#'
#' # Perform matching and estimation, without permutations or jackknife
#' # runtime: < 1 min
#'
#' sea1 <- microsynth(seattledmi,
#'                   idvar="ID", timevar="time", intvar="Intervention",
#'                   start.pre=1, end.pre=12, end.post=16,
#'                   match.out=match.out, match.covar=cov.var,
#'                   result.var=match.out, omnibus.var=match.out,
#'                   test="lower",
#'                   n.cores = min(parallel::detectCores(), 2))
#'
#' # View results
#' summary(sea1)
#' plot_microsynth(sea1)
#'
#' \donttest{
#'
#' # Repeat matching and estimation, with permutations and jackknife
#' # Set permutations and jack-knife to very few groups (2) for
#' # quick demonstration only.
#' # runtime: ~30 min
#' sea2 <- microsynth(seattledmi,
#'                      idvar="ID", timevar="time", intvar="Intervention",
#'                      start.pre=1, end.pre=12, end.post=c(14, 16),
#'                      match.out=match.out, match.covar=cov.var,
#'                      result.var=match.out, omnibus.var=match.out,
#'                      test="lower",
#'                      perm=250, jack=TRUE,
#'                      result.file=file.path(tempdir(), 'ExResults2.xlsx'),
#'                      n.cores = min(parallel::detectCores(), 2))
#'
#' # View results
#' summary(sea2)
#' plot_microsynth(sea2)
#'
#' # Specify additional outcome variables for matching, which makes
#' # matching harder.
#' match.out <- c('i_robbery','i_aggassau','i_burglary','i_larceny',
#'        'i_felony','i_misdemea','i_drugsale','i_drugposs','any_crime')
#'
#' # Perform matching, setting check.feas = T and use.backup = T
#' # to ensure model feasibility
#' # runtime: ~40 minutes
#' sea3 <- microsynth(seattledmi,
#'                    idvar="ID", timevar="time", intvar="Intervention",
#'                    end.pre=12,
#'                    match.out=match.out, match.covar=cov.var,
#'                    result.var=match.out, perm=250, jack=0,
#'                    test="lower", check.feas=TRUE, use.backup = TRUE,
#'                    result.file=file.path(tempdir(), 'ExResults3.xlsx'),
#'                    n.cores = min(parallel::detectCores(), 2))
#'
#'
#' # Aggregate outcome variables before matching, to boost model feasibility
#' match.out <- list( 'i_robbery'=rep(2, 6), 'i_aggassau'=rep(2, 6),
#'          'i_burglary'=rep(1, 12), 'i_larceny'=rep(1, 12),
#'          'i_felony'=rep(2, 6), 'i_misdemea'=rep(2, 6),
#'          'i_drugsale'=rep(4, 3), 'i_drugposs'=rep(4, 3),
#'          'any_crime'=rep(1, 12))
#'
#' # After aggregation, use.backup and cheack.feas no longer needed
#' # runtime: ~40 minutes
#' sea4 <- microsynth(seattledmi, idvar='ID', timevar='time',
#'          intvar='Intervention', match.out=match.out, match.covar=cov.var,
#'          start.pre=1, end.pre=12, end.post=16,
#'          result.var=names(match.out), omnibus.var=names(match.out),
#'          perm=250, jack = TRUE, test='lower',
#'          result.file=file.path(tempdir(), 'ExResults4.xlsx'),
#'          n.cores = min(parallel::detectCores(), 2))
#'
#' # View results
#' summary(sea4)
#' plot_microsynth(sea4)
#'
#'
#' # Generate weights only (for four variables)
#' match.out <- c('i_felony', 'i_misdemea', 'i_drugs', 'any_crime')
#'
#' # runtime: ~ 20 minutes
#' sea5 <- microsynth(seattledmi,  idvar='ID', timevar='time',
#'          intvar='Intervention', match.out=match.out, match.covar=cov.var,
#'          start.pre=1, end.pre=12, end.post=16,
#'          result.var=FALSE, perm=250, jack=TRUE,
#'          n.cores = min(parallel::detectCores(), 2))
#'
#' # View weights
#' summary(sea5)
#'
#' # Generate results only
#' sea6 <- microsynth(seattledmi, idvar='ID', timevar='time',
#'           intvar='Intervention',
#'           start.pre=1, end.pre=12, end.post=c(14, 16),
#'           result.var=match.out, test='lower',
#'           w=sea5, result.file=file.path(tempdir(), 'ExResults6.xlsx'),
#'           n.cores = min(parallel::detectCores(), 2))
#'
#' # View results (including previously-found weights)
#' summary(sea6)
#'
#' # Generate plots only
#' plot_microsynth(sea6, plot.var=match.out[1:2])
#'
#' # Apply microsynth in the traditional setting of Synth
#' # Create macro-level (small n) data, with 1 treatment unit
#' set.seed(86872)
#' ids.t <- names(table(seattledmi$ID[seattledmi$Intervention==1]))
#' ids.c <- names(table(seattledmi$ID[seattledmi$Intervention==0]))
#' ids.synth <- c(base::sample(ids.t, 1), base::sample(ids.c, 100))
#' seattledmi.one <- seattledmi[is.element(seattledmi$ID,
#'            as.numeric(ids.synth)), ]
#'
#' # Apply microsynth to the new macro-level data
#' # runtime: < 5 minutes
#' sea8 <- microsynth(seattledmi.one, idvar='ID', timevar='time',
#'            intvar='Intervention',
#'            start.pre=1, end.pre=12, end.post=16,
#'            match.out=match.out[4],
#'            match.covar=cov.var, result.var=match.out[4],
#'            test='lower', perm=250, jack=FALSE,
#'            check.feas=TRUE, use.backup=TRUE,
#'            n.cores = min(parallel::detectCores(), 2))
#'
#' # View results
#' summary(sea8)
#' plot_microsynth(sea8)
#'
#' # Use microsynth to calculate propensity score-type weights
#' # Prepare cross-sectional data at time of intervention
#' seattledmi.cross <- seattledmi[seattledmi$time==16, colnames(seattledmi)!="time"]#'
#'
#' # Apply microsynth to find propensity score-type weights
#' # runtime: ~5 minutes
#' sea9 <- microsynth(seattledmi.cross, idvar='ID', intvar='Intervention',
#'              match.out=FALSE, match.covar=cov.var, result.var=match.out,
#'              test='lower', perm=250, jack=TRUE,
#'              n.cores = min(parallel::detectCores(), 2))
#'
#' # View results
#' summary(sea9)
#' }
#'
#' @export

microsynth <- function (data, idvar, intvar, timevar = NULL, start.pre = NULL,
                        end.pre = NULL, end.post = NULL, match.out = TRUE, match.covar = TRUE,
                        match.out.min = NULL, match.covar.min = NULL, result.var = TRUE,
                        omnibus.var = result.var, period = 1, scale.var = "Intercept",
                        confidence = 0.9, test = "twosided", perm = 0, jack = 0,
                        use.survey = TRUE, cut.mse = Inf, check.feas = FALSE, use.backup = FALSE,
                        w = NULL, max.mse = 0.01, maxit = 250, cal.epsilon = 1e-04,
                        calfun = "linear", bounds = c(0, Inf), result.file = NULL,
                        printFlag = TRUE, n.cores = TRUE)
{

  # Determine the number of cores to be used. CRAN tops at 2.
  n.cores <- msCluster(n.cores)

  # Declare metrics for print() call (1 of 3)
  info <- list()
  info$match <- unique(match.out)
  info$match.min <- unique(match.out.min)
  info$covar <- unique(match.covar)
  info$covar.min <- unique(match.covar.min)
  info$start.pre <- start.pre
  info$end.pre <- end.pre
  info$end.post <- end.post

  all.tmp <- proc.time()
  if (length(timevar) == 0) {
    if (length(table(data[, idvar])) < NROW(data)) {
      stop("Data are not cross-sectional.  Please specify timevar.")
    }
    else {
      data$Time <- 1
      timevar <- "Time"
      end.pre <- 1
    }
  }
  time.tmp <- data[,timevar]
  time.names <- names(table(time.tmp))
  data[,timevar] <- match(as.character(time.tmp), time.names)

  # Declare more metrics for print() call (2 of 3)
  info$nUnits <- length(unique(data[[idvar]])) # num units
  info$nTreatment <- length(unique(data[idvar][data[intvar]==1]))
  info$nControl <- info$nUnits - info$nTreatment

  # Set up dummy variables for use in determining values of indirectly stated inputs
  if (length(start.pre) > 0 & !is.logical(start.pre)) {
    start.pre <- match(as.character(start.pre), time.names)
  }
  if (length(end.pre) > 0) {
    end.pre <- match(as.character(end.pre), time.names)
  }
  if (length(end.post) > 0 & !is.logical(end.post)) {
    end.post <- match(as.character(end.post), time.names)
  }
  twosided <- TRUE
  if (test == "lower" | test == "upper") {
    twosided <- FALSE
  }
  if (is.logical(match.out)) {
    if (!match.out) {
      reset.match.out <- TRUE
    }
    else {
      reset.match.out <- FALSE
    }
    match.out <- NULL
  }
  else if (length(match.out) == 0) {
    reset.match.out <- TRUE
  }
  else {
    reset.match.out <- FALSE
  }
  if (is.logical(match.covar)) {
    if (!match.covar) {
      reset.match.covar <- TRUE
    }
    else {
      reset.match.covar <- FALSE
    }
    match.covar <- NULL
  }
  else if (length(match.covar) == 0) {
    reset.match.covar <- TRUE
  }
  else {
    reset.match.covar <- FALSE
  }
  if (is.logical(result.var)) {
    if (!result.var) {
      reset.result.var <- TRUE
    }
    else {
      reset.result.var <- FALSE
    }
    result.var <- NULL
  }
  else if (length(result.var) == 0) {
    reset.result.var <- TRUE
  }
  else {
    reset.result.var <- FALSE
  }
  match.out <- remove.vars(match.out, dimnames(data)[[2]],
                           "match.out", printFlag = printFlag)
  match.out.min <- remove.vars(match.out.min, dimnames(data)[[2]],
                               "match.out.min", printFlag = printFlag)
  match.covar <- remove.vars(match.covar, dimnames(data)[[2]],
                             "match.covar", printFlag = printFlag)
  match.covar.min <- remove.vars(match.covar.min, dimnames(data)[[2]],
                                 "match.covar.min", printFlag = printFlag)
  result.var <- remove.vars(result.var, dimnames(data)[[2]],
                            "out.covar", printFlag = printFlag)
  if (!is.logical(omnibus.var)) {
    omnibus.var <- remove.vars(omnibus.var, dimnames(data)[[2]],
                               "omnibus.var", printFlag = printFlag)
  }
  nv.names <- union(match.covar, match.covar.min)
  v.names <- result.var
  if (length(match.out) > 0) {
    v.names <- union(v.names, names(match.out))
  }
  if (length(match.out.min) > 0) {
    if(is.list(match.out.min)) {
      v.names <- union(v.names, names(match.out.min))
    } else {
      v.names <- union(v.names, match.out.min)
    }
  }
  ok.cl <- c("numeric", "integer", "logical")
  classes <- NA
  for (i in 1:NCOL(data)) {
    classes[i] <- class(data[, i])
  }
  ok.col <- is.element(classes, ok.cl)
  rm.col <- colnames(data)[!ok.col]
  rm.col <- setdiff(rm.col, c(idvar, timevar, intvar))
  data <- data[, !is.element(colnames(data), rm.col)]
  v.names <- setdiff(v.names, rm.col)
  nv.names <- setdiff(nv.names, rm.col)

  # Shape panel data into 3D array, generate intervention matrix, etc.
  data <- newreshape(data, nv.names = nv.names, v.names = v.names,
                     timevar = timevar, idvar = idvar, intvar = intvar)
  if (length(result.var) == 0) {
    result.var <- data[[4]]
    if (!reset.result.var) {
      if(printFlag){message("result.var = TRUE.  Resetting: \n", appendLF = FALSE)}
      if(printFlag){message("result.var = c(\"", paste(result.var, collapse = "\",\"",
                                                       sep = ""), "\")\n\n", sep = "", appendLF = FALSE)}
    }
  }
  if (length(match.covar) == 0) {
    match.covar <- data[[3]]
    if (!reset.match.covar) {
      if(printFlag){message("match.covar = TRUE.  Resetting: \n", appendLF = FALSE)}
      if(printFlag){message("match.covar = c(\"", paste(match.covar,
                                                        collapse = "\",\"", sep = ""), "\")\n\n", sep = "",
                            appendLF = FALSE)}
    }
  }
  Intervention <- data[[2]]
  data <- data[[1]]
  times <- as.numeric(colnames(Intervention))
  if (length(end.pre) == 0) {
    eval.times <- which(colSums(Intervention) != 0)
    if (length(eval.times) == 0) {
      stop("There are no intervention cases.\n")
    }
    end.pre <- min(eval.times) - 1
    if (end.pre == 0) {
      stop("There are no pre-intervention time points.")
    }
    end.pre <- times[end.pre]
    if(printFlag){message("Setting end.pre = ", time.names[end.pre], ".\n\n", sep = "",
                          appendLF = FALSE)}
  }
  Intervention[] <- as.integer(rowSums(Intervention) != 0)
  if (length(start.pre) == 0) {
    start.pre <- min(times)
  }
  if (length(end.post) == 0) {
    end.post <- max(times)
  }
  if (length(match.out) == 0) {
    match.out <- result.var
    if (!reset.match.out) {
      if(printFlag){message("match.out = TRUE.  Resetting: \n", appendLF = FALSE)}
      if(printFlag){message("match.out = c(\"", paste(match.out, collapse = "\",\"",
                                                      sep = ""), "\")\n\n", sep = "", appendLF = FALSE)}
    }
  }
  if (!is.list(match.out) & length(match.out) > 0) {
    match.out.tmp <- match.out
    match.out <- list()
    for (i in 1:length(match.out.tmp)) {
      match.out[[i]] <- rep(period, (end.pre - start.pre +
                                       1)%/%period)
    }
    names(match.out) <- match.out.tmp
  } else if (length(match.out) > 0) {
    match.out <- check.matchout(match.out, end.pre - min(times) + 1)
    if (match.out[[2]]) {
      if(printFlag){message("WARNING: match.out calls on time periods that are beyond the data range.\n",
                            sep ="", appendLF = FALSE)}
      if(printFlag){message("match.out is being reset accordingly.\n\n", sep="", appendLF = FALSE)}
    }
    match.out <- match.out[[1]]
  }
  if (!is.list(match.out.min) & length(match.out.min) > 0) {
    match.out.tmp1 <- match.out.min
    match.out.min <- list()
    for (i in 1:length(match.out.tmp1)) {
      match.out.min[[i]] <- rep(period, (end.pre - start.pre +
                                           1)%/%period)
    }
    names(match.out.min) <- match.out.tmp1
  } else if (length(match.out.min) > 0) {
    match.out.min <- check.matchout(match.out.min, end.pre - min(times) + 1)
    if (match.out.min[[2]]) {
      if(printFlag){message("WARNING: match.out.min calls on time periods that are beyond the data range. \n",
                            sep ="", appendLF = FALSE)}
      if(printFlag){message("match.out.min is being reset accordingly.\n\n", sep="", appendLF = FALSE)}
    }
    match.out.min <- match.out.min[[1]]
  }
  if (is.logical(omnibus.var)) {
    if (omnibus.var) {
      if (length(match.out) == 1) {
        omnibus.var <- NULL
      }
      else {
        omnibus.var <- result.var
        if(printFlag){message("omnibus.var = TRUE.  Resetting: \n",
                              appendLF = FALSE)}
        if(printFlag){message("omnibus.var = c(\"", paste(omnibus.var,
                                                          collapse = "\",\"", sep = ""), "\")\n\n", sep = "",
                              appendLF = FALSE)}
      }
    }
    else {
      omnibus.var <- NULL
    }
  }
  int.num <- 1
  dum <- max(colSums(Intervention == 1))
  dum <- min(NROW(Intervention) - dum, dum)
  dum1 <- ((max(end.post) - (max(end.post) - end.pre)%%period -
              end.pre)/period)
  if (dum <= dum1 + 1) {
    if(printFlag){message("WARNING: There is a low number (", dum, ") of cases in the treatment or intervention group.\n",
                          sep = "", appendLF = FALSE)}
    if (jack > 0) {
      jack <- 0
      if(printFlag){message("setting jack = 0.\n", appendLF = FALSE)}
    }
    if (use.survey) {
      use.survey <- FALSE
      if(printFlag){message("Setting use.survey = FALSE.\n", appendLF = FALSE)}
    }
    if(printFlag){message("Be cautious of results involving linearization or confidence intervals.\n\n",
                          appendLF = FALSE)}
  }
  if (reset.match.covar) {
    match.covar <- NULL
  }
  if (reset.match.out) {
    match.out <- NULL
  }
  if (length(rm.col) > 0) {
    for (i in 1:length(rm.col)) {
      if (is.element(rm.col[i], result.var)) {
        if(printFlag){message("WARNING: ", rm.col[i], " is a non-numeric variable.  It will be removed from result.var. \n",
                              sep = "", appendLF = FALSE)}
        result.var <- setdiff(result.var, rm.col[i])
      }
      if (is.element(rm.col[i], omnibus.var)) {
        if(printFlag){message("WARNING: ", rm.col[i], " is a non-numeric variable.  It will be removed from omnibus.var. \n",
                              sep = "", appendLF = FALSE)}
        omnibus.var <- setdiff(omnibus.var, rm.col[i])
      }
      if (is.element(rm.col[i], match.covar)) {
        if(printFlag){message("WARNING: ", rm.col[i], " is a non-numeric variable.  It will be removed from match.covar. \n",
                              sep = "", appendLF = FALSE)}
        match.covar <- setdiff(match.covar, rm.col[i])
      }
      if (is.element(rm.col[i], match.covar.min)) {
        if(printFlag){message("WARNING: ", rm.col[i], " is a non-numeric variable.  It will be removed from match.covar.min. \n",
                              sep = "", appendLF = FALSE)}
        match.covar.min <- setdiff(match.covar.min, rm.col[i])
      }
      if (is.element(rm.col[i], names(match.out))) {
        if(printFlag){message("WARNING: ", rm.col[i], " is a non-numeric variable.  It will be removed from match.out. \n",
                              sep = "", appendLF = FALSE)}
        rm.li <- which(is.element(names(match.out), rm.col[i]))
        match.out <- match.out[-rm.li]
      }
      if (is.element(rm.col[i], names(match.out.min))) {
        if(printFlag){message("WARNING: ", rm.col[i], " is a non-numeric variable.  It will be removed from match.out.min. \n",
                              sep = "", appendLF = FALSE)}
        rm.li1 <- which(is.element(names(match.out.min),
                                   rm.col[i]))
        match.out.min <- match.out.min[-rm.li1]
      }
    }
  }

  # Establish synthetic control weights
  if (length(w) == 0) {
    tmp <- proc.time()
    if(printFlag){message("Calculating weights...", "\n", appendLF = FALSE)}
    w <- get.w(data, match.covar, match.covar.min, match.out,
               match.out.min, boot = perm, jack = jack, Int = Intervention[,
                                                                           as.character(end.pre)], int.val = int.num, trim = NULL,
               end.pre = end.pre, cal.epsilon = cal.epsilon, maxit = maxit,
               bounds = bounds, calfun = calfun, check.feas = check.feas,
               scale.var = scale.var, cut.mse = max.mse, use.backup = use.backup, time.names = time.names,
               printFlag = printFlag, n.cores = n.cores)
    tmp <- proc.time() - tmp
    if(printFlag){message("Calculation of weights complete: Total time = ",
                          round(tmp[3], 2), "\n\n", sep = "", appendLF = FALSE)}
  }
  else {
    if(printFlag){message("Weights have been provided.  Will not calculate weights.\n",
                          appendLF = FALSE)}
    if (class(w) == "microsynth") {
      w <- w$w
    }
    is.correct.w <- is.list(w)
    if (is.correct.w) {
      is.correct.w <- is.correct.w & sum(names(w) != c("Weights",
                                                       "Intervention", "MSE", "Model", "Summary", "keep.groups")) ==
        0
    }
    if (is.correct.w) {
      is.correct.w <- is.correct.w & dim(w$Weights)[1] ==
        dim(data)[1]
    }
    if (is.correct.w) {
      is.correct.w <- is.correct.w & dim(w$Intervention)[1] ==
        dim(data)[1]
    }
    if (!is.correct.w) {
      stop("w is not formatted correctly.")
    }
    perm <- sum(grepl("Perm", colnames(w$Weights)))
    jack <- sum(grepl("Jack", colnames(w$Weights)))
    if(printFlag){message("Setting jack = ", jack, "\n", sep = "", appendLF = FALSE)}
    if(printFlag){message("Setting perm = ", perm, "\n\n", sep = "", appendLF = FALSE)}
  }
  not.jack <- !grepl("Jack", colnames(w$Weights))
  end.post <- end.post - (end.post - end.pre)%%period
  if (!reset.result.var) {
    stats <- list()
    stats1 <- list()
    stats2 <- list()
    delta.out <- list()
    dof <- list()
    out.coefs <- list()
    results <- list()
    for (i in 1:length(end.post)) {
      tmp <- proc.time()
      is.graph <- "Calculating basic statistics"
      is.graph1 <- "Completed calculation of basic statistics"
      if(printFlag){message(is.graph, " for end.post = ", time.names[end.post[i]],
                            "...", "\n", sep = "", appendLF = FALSE)}

      # Calculate basic statistics
      stats[[i]] <- get.stats(data, w$Weights, w$Intervention,
                              w$keep.groups, result.var, end.pre = end.pre,
                              period = period, end.post = end.post[i],
                              omnibus.var = omnibus.var,
                              start.pre = start.pre, cut.mse = cut.mse,
                              twosided = twosided, time.names = time.names)
      if (i == which.max(end.post)) {
        plot.stats <- stats[[i]][[5]]
      }
      stats[[i]] <- stats[[i]][-5]
      tmp <- proc.time() - tmp
      if(printFlag){message(is.graph1, " for end.post = ", time.names[end.post[i]],
                            ".  Time = ", round(tmp[3], 2), "\n\n", sep = "",
                            appendLF = FALSE)}
      if (!reset.result.var) {
        w.tmp <- w$Weights
        Inter.tmp <- w$Intervention
        mse.tmp <- w$keep.groups
        if (!use.survey) {
          keep.surv <- !grepl("Perm", colnames(w.tmp))
          w.tmp <- w.tmp[, keep.surv, drop = FALSE]
          Inter.tmp <- Inter.tmp[, keep.surv, drop = FALSE]
          mse.tmp <- mse.tmp[keep.surv]
        }
        tmp <- proc.time()
        if(printFlag){message("Calculating survey statistics for end.post = ",
                              time.names[end.post[i]], "...", "\n", sep = "", appendLF = FALSE)}

        # Calculate complex (i.e., survey) statistics
        stats.tmp <- get.stats1(data, w.tmp, Inter.tmp,
                                mse.tmp, result.var, end.pre = end.pre, period = period,
                                end.post = end.post[i], omnibus.var = omnibus.var,
                                twosided = twosided, printFlag = printFlag, n.cores = n.cores)
        stats1[[i]] <- stats.tmp[[1]]
        stats2[[i]] <- stats.tmp[[2]]
        delta.out[[i]] <- stats.tmp[[3]]
        dof[[i]] <- stats.tmp[[4]]
        out.coefs[[i]] <- stats.tmp[[5]]
        tmp <- proc.time() - tmp
        if(printFlag){message("Completed calculation of survey statistics for end.post = ",
                              time.names[end.post[i]], ".  Time = ", round(tmp[3], 2),
                              "\n\n", sep = "", appendLF = FALSE)}
        Pct.Chng <- cbind(Pct.Chng = stats[[i]][[2]][1,
                                                     ])
        if (is.element("Omnibus", rownames(Pct.Chng))) {
          Pct.Chng["Omnibus", ] <- NA
        }
        Trt <- cbind(Trt = stats[[i]][[3]][1, ])
        Con <- cbind(Con = stats[[i]][[4]][1, ])
        synth.stats1 <- stats1[[i]]
        synth.stats2 <- stats2[[i]]
        synth.delta.out <- delta.out[[i]]
        synth.dof <- dof[[i]]
        if (test == "lower") {
          Linear.pVal <- cbind(Linear.pVal = stats::pnorm(synth.stats1[1,
                                                                       ], lower.tail = TRUE))
        }
        else if (test == "upper") {
          Linear.pVal <- cbind(Linear.pVal = stats::pnorm(synth.stats1[1,
                                                                       ], lower.tail = FALSE))
        }
        else {
          Linear.pVal <- cbind(Linear.pVal = 2 * stats::pnorm(abs(synth.stats1[1,
                                                                               ]), lower.tail = FALSE))
          if (is.element("Omnibus", rownames(Linear.pVal))) {
            Linear.pVal["Omnibus", ] <- stats::pchisq(synth.stats1[1,
                                                                   "Omnibus"], df = synth.dof, lower.tail = FALSE)
          }
        }
        Linear.CI <- make.ci(synth.stats2[1, ], sqrt(synth.delta.out[1,
                                                                     ]), alpha = 1 - confidence)
        colnames(Linear.CI) <- paste("Linear.", colnames(Linear.CI),
                                     sep = "")
        Jack.pVal <- Jack.CI <- Perm.pVal <- Perm.CI <- NULL
        if (jack > 0) {
          jack.stats1 <- synth.stats1[2, ]
          jack.stats2 <- synth.stats2[2, ]
          jack.delta.out <- synth.delta.out[2, ]
          synth.stats1 <- synth.stats1[-2, , drop = FALSE]
          synth.stats2 <- synth.stats2[-2, , drop = FALSE]
          synth.delta.out <- synth.delta.out[-2, , drop = FALSE]
          if (test == "lower") {
            Jack.pVal <- cbind(Jack.pVal = stats::pnorm(jack.stats1,
                                                        lower.tail = TRUE))
          }
          else if (test == "upper") {
            Jack.pVal <- cbind(Jack.pVal = stats::pnorm(jack.stats1,
                                                        lower.tail = FALSE))
          }
          else {
            Jack.pVal <- cbind(Jack.pVal = 2 * stats::pnorm(abs(jack.stats1),
                                                            lower.tail = FALSE))
            if (is.element("Omnibus", rownames(Jack.pVal))) {
              Jack.pVal["Omnibus", ] <- stats::pchisq(c(jack.stats1)["Omnibus"],
                                                      df = synth.dof, lower.tail = FALSE)
            }
          }
          Jack.CI <- make.ci(jack.stats2, sqrt(jack.delta.out),
                             alpha = 1 - confidence)
          colnames(Jack.CI) <- paste("Jack.", colnames(Jack.CI),
                                     sep = "")
        }
        if (perm > 0 & use.survey) {
          perm.stats1 <- get.pval(list(synth.stats1),
                                  ret.na = TRUE)
          dum <- perm.stats1[[2]]
          perm.stats1 <- perm.stats1[[1]]
          if (test == "lower") {
            Perm.pVal <- cbind(Perm.pVal = c(perm.stats1))
          }
          else if (test == "upper") {
            Perm.pVal <- cbind(Perm.pVal = 1 - c(perm.stats1))
          }
          else {
            Perm.pVal <- cbind(Perm.pVal = 2 * apply(cbind(c(perm.stats1),
                                                           1 - c(perm.stats1)), 1, min))
            if (is.element("Omnibus", rownames(Perm.pVal))) {
              Perm.pVal["Omnibus", ] <- 1 - c(perm.stats1)["Omnibus"]
            }
          }
          Perm.CI <- make.ci2(synth.stats2, synth.delta.out,
                              alpha = 1 - confidence)
          colnames(Perm.CI) <- paste("Perm.", colnames(Perm.CI),
                                     sep = "")
        }
        else if (perm > 0) {
          perm.stats = get.pval(stats[[i]], ret.na = TRUE)
          dum <- perm.stats[[2]]
          perm.stats <- (perm.stats[[1]])[c(1), , drop = FALSE]
          if (test == "lower") {
            Perm.pVal <- cbind(Perm.pVal = c(perm.stats[1,
                                                        ]))
          }
          else if (test == "upper") {
            Perm.pVal <- cbind(Perm.pVal = 1 - c(perm.stats[1,
                                                            ]))
          }
          else {
            Perm.pVal <- cbind(Perm.pVal = 2 * pmin(c(perm.stats[1,
                                                                 ]), 1 - c(perm.stats[1, ])))
            if (is.element("Omnibus", rownames(Perm.pVal))) {
              Perm.pVal["Omnibus", ] <- 1 - c(perm.stats[1,
                                                         ])["Omnibus"]
            }
          }
        }
        results[[i]] <- cbind(Trt, Con, Pct.Chng, Linear.pVal,
                              Linear.CI, Jack.pVal, Jack.CI, Perm.pVal, Perm.CI)
        if (NROW(results[[i]]) == 1) {
          rownames(results[[i]]) <- result.var[1]
        }
        names(results)[i] <- names(out.coefs)[i] <- time.names[end.post[i]]
      }
    }
  }
  if (reset.result.var) {
    if(printFlag){message("No outcome variables specified (e.g., result.var = NULL).\n",
                          appendLF = FALSE)}
    if(printFlag){message("Results will not be tabulated.\n", appendLF = FALSE)}
  }
  if (reset.result.var) {
    if(printFlag){message("Returning weights only.\n", appendLF = FALSE)}
  }

  # Tabulate results for output.
  out <- list()
  i <- 1
  out[[i]] <- w
  names(out)[i] <- "w"
  i <- i + 1
  if (!reset.result.var) {
    out[[i]] <- results
    names(out)[i] <- "Results"
    out[[i + 1]] <- out.coefs
    names(out)[i + 1] <- "svyglm.stats"
    i <- i + 2
  }
  if (!reset.result.var) {
    out[[i]] <- plot.stats
    names(out)[i] <- "Plot.Stats"
    i <- i + 1
  }
  ret.stats <- FALSE
  if (ret.stats) {
    out[[i]] <- stats
    names(out)[i] <- "stats"
    out[[i + 1]] <- stats1
    out[[i + 2]] <- stats2
    out[[i + 3]] <- delta.out
    names(out)[i:(i + 3)] <- c("stats", "stats1", "stats2",
                               "delta.out")
    i <- i + 4
  }
  if (length(result.file) > 0 & !reset.result.var) {
    out.results(results, end.pre, end.post = names(results), result.file, printFlag = printFlag)
  }
  all.tmp <- proc.time() - all.tmp
  if(printFlag){message("microsynth complete: Overall time = ", round(all.tmp[3],
                                                                      2), "\n\n", sep = "", appendLF = FALSE)}

  # Declare final output for print (3 of 3)
  info$nConstraints <- nrow(out$w$Summary) - 1
  info$num.constr <- out$w$num.constr
  out$w$num.constr <- NULL

  # Add descriptive stats to output for print() call
  out$info <- info

  out <- makemicrosynth(out)
  return(out)
}




# Sub-function of microsynth(); calculate p-values when permutation is used
get.pval <- function(stats, p = NCOL(stats[[1]]), k = length(stats), ret.na = FALSE) {
  nams <- colnames(stats[[1]])
  n.boot <- NROW(stats[[1]]) - 1
  out <- matrix(NA, k, p)
  colnames(out) <- nams
  rownames(out) <- paste("Stat", 1:k, sep = "")

  dum <- list()
  for (i in 1:k) {
    synth <- stats[[i]][1, ]
    synth.tmp <- t(matrix(synth, p, n.boot))
    boot <- stats[[i]][-1, , drop = FALSE]
    boot[boot == Inf | boot == -Inf] <- NA
    dum[[i]] <- colSums(is.na(boot))
    out[i, ] <- colMeans(synth.tmp > boot, na.rm = TRUE)
  }
  if (ret.na) {
    return(list(out, dum))
  } else {
    return(out)
  }
}



# Sub-function of microsynth(); reshape data to make 3D array with other necessary outputs
newreshape <- function(data, timevar, idvar, intvar, v.names = NULL, nv.names = NULL) {
  times = names(table(data[, timevar]))
  if (length(v.names) == 0) {
    v.names1 <- v.names <- setdiff(colnames(data), c(nv.names, timevar, idvar))
  } else {
    v.names1 <- c(intvar, v.names)
  }
  keep.nams <- c(idvar, timevar, intvar, nv.names, v.names)
  keep.nams <- names(table(keep.nams))
  data <- data[, keep.nams]

  if (length(times) > 1) {
    newdat <- stats::reshape(data, timevar = timevar, idvar = idvar, v.names = v.names1, direction = "wide")
  } else if (length(times) == 1) {
    newdat <- data[, colnames(data) != timevar]
    append.cols <- setdiff(colnames(newdat), c(idvar, nv.names))
    append.cols <- which(is.element(colnames(newdat), append.cols))
    colnames(newdat)[append.cols] <- paste(colnames(newdat)[append.cols], ".", times, sep = "")
  }
  intcols <- substr(colnames(newdat), 1, nchar(intvar)) == intvar
  intcols <- intcols & (substr(colnames(newdat), nchar(intvar) + 1, nchar(intvar) + 1) == "" | substr(colnames(newdat), nchar(intvar) +
                                                                                                        1, nchar(intvar) + 1) == ".")
  intcols <- colnames(newdat)[intcols]
  int <- newdat[, intcols, drop = FALSE]
  colnames(int) <- gsub(paste(intvar, ".", sep = ""), "", colnames(int))
  intcols1 <- as.numeric(colnames(int))
  int <- int[, order(intcols1, decreasing = FALSE), drop = FALSE]
  rownames(int) <- newdat[, idvar]

  v.names <- setdiff(v.names, intvar)
  nv.names <- setdiff(nv.names, intvar)

  newdat <- newdat[, setdiff(colnames(newdat), intcols)]

  out <- array(NA, c(NROW(newdat), length(nv.names) + length(v.names), NCOL(int)))
  dimnames(out) <- list(newdat[, idvar], c(nv.names, v.names), colnames(int))

  for (i in 1:dim(out)[2]) {
    nam.tmp <- dimnames(out)[[2]][i]
    here <- substr(colnames(newdat), 1, nchar(nam.tmp)) == nam.tmp
    here <- here & (substr(colnames(newdat), nchar(nam.tmp) + 1, nchar(nam.tmp) + 1) == "." | substr(colnames(newdat), nchar(nam.tmp) +
                                                                                                       1, nchar(nam.tmp) + 1) == "")
    tmp <- newdat[, here, drop = FALSE]
    if (sum(here, na.rm = TRUE) > 1) {
      colnames(tmp) <- gsub(paste(dimnames(out)[[2]][i], ".", sep = ""), "", colnames(tmp))
      tmpcols1 <- as.numeric(colnames(tmp))
      tmp <- tmp[, order(tmpcols1, decreasing = FALSE)]
    }
    out[, i, ] <- as.matrix(tmp)
  }

  nv.names <- NULL
  if (dim(out)[3] > 1) {
    for (i in 1:dim(out)[2]) {
      diff.out <- out[, i, ] - rowMeans(out[, i, ])
      diff.out <- rowMeans(diff.out^2)
      if (sum(diff.out, na.rm = TRUE) == 0) {
        nv.names <- c(nv.names, dimnames(out)[[2]][i])
      }
    }
  } else {
    nv.names <- dimnames(out)[[2]]
  }

  v.names <- setdiff(dimnames(out)[[2]], nv.names)

  return(list(bigdat = out, Intervention = int, nv.names = nv.names, v.names = v.names))
}



# Sub-function of microsynth(); create results file (e.g., as CSV or XLSX)
out.results <- function(results, end.pre, end.post = names(results), file = NULL, printFlag = TRUE) {
  use.xlsx <- length(end.post) > 1

  if (substr(file, nchar(file) - 3, nchar(file)) == ".csv") {
    file <- substr(file, 1, nchar(file) - 4)
    if (use.xlsx) {
      if(printFlag){message("WARNING: Cannot return .csv file since length(end.post) > 1.  Returning .xlsx file instead.\n",
                            sep = "", appendLF = FALSE)}
    }
  } else if (substr(file, nchar(file) - 3, nchar(file)) == ".xls") {
    file <- substr(file, 1, nchar(file) - 4)
    use.xlsx <- TRUE
  } else if (substr(file, nchar(file) - 4, nchar(file)) == ".xlsx") {
    file <- substr(file, 1, nchar(file) - 5)
    use.xlsx <- TRUE
  }

  if (use.xlsx) {
    xlsx.loaded <- is.element("xlsx", loadedNamespaces())
    file <- paste(file, ".xlsx", sep = "")
    if (!requireNamespace("xlsx", quietly = TRUE)) {
      stop("The xlsx package is needed when saving output with multiple post-intervention times, or when the file name specifies a .xlsx extension. Please install xlsx, or if you'd like to write to .csv, only select a single post-intervention time and append the filename appropriately.",
           call. = FALSE)
    }
    if (!xlsx.loaded) {
      attachNamespace("xlsx")
    }
    end.post <- as.numeric(end.post)
    wb <- xlsx::createWorkbook()
    cspValColumn <- xlsx::CellStyle(wb, dataFormat = xlsx::DataFormat("0.0000"))
    csOtherColumn <- xlsx::CellStyle(wb, dataFormat = xlsx::DataFormat("0.00"))
    csPercColumn <- xlsx::CellStyle(wb, dataFormat = xlsx::DataFormat("0.0%"))

    for (i in 1:length(end.post)) {
      out <- results[[i]]
      keep <- rowMeans(is.na(out)) < 1
      out <- out[keep, , drop = FALSE]

      is.pct <- grepl("pct", tolower(colnames(out))) | grepl("lower", tolower(colnames(out))) | grepl("upper", tolower(colnames(out)))
      is.pval <- grepl("pval", tolower(colnames(out)))
      is.oth <- !is.pct & !is.pval

      pct.list <- rep(list(csPercColumn), sum(is.pct))
      names(pct.list) <- as.character(which(is.pct))
      pval.list <- rep(list(cspValColumn), sum(is.pval))
      names(pval.list) <- as.character(which(is.pval))
      oth.list <- rep(list(csOtherColumn), sum(is.oth))
      names(oth.list) <- as.character(which(is.oth))
      dum1 <- which(rownames(out) == "Omnibus")
      rm.omni <- c("Naive.Norm", "Pct.Chng")
      dum2 <- which(is.element(colnames(out), rm.omni))
      if (length(dum1) > 0 & length(dum2) > 0) {
        out[dum1, dum2] <- NA
      }
      p <- NCOL(out)
      n <- NROW(out)

      nam <- paste("Max time = ", end.post[i], sep = "")
      sheet <- xlsx::createSheet(wb, sheetName = nam)
      xlsx::setColumnWidth(sheet, colIndex = 1:(p + 1), colWidth = 13)
      xlsx::addDataFrame(out, sheet, colStyle = c(pct.list, pval.list, oth.list))

      bord1 <- xlsx::Border(color = "black", position = "BOTTOM", pen = "BORDER_THIN")
      cb <- xlsx::CellBlock(sheet, startRow = 1, startColumn = 1, noRows = (n + 1), noColumns = (p + 1), create = FALSE)
      xlsx::CB.setBorder(cb, border = bord1, rowIndex = 1, colIndex = 1:(p + 1))

      font <- xlsx::Font(wb, isBold = TRUE)
      xlsx::CB.setFont(cb, font = font, rowIndex = 1, colIndex = 2:(p + 1))

    }

    xlsx::saveWorkbook(wb, file = file)
    rm(wb)
  } else {
    file <- paste(file, ".csv", sep = "")
    utils::write.csv(results[[1]], file, na = "")
  }

}



# Sub-function of microsynth(); calculate confidence intervals
make.ci <- function(means, ses, alpha = 0.05) {
  z.score <- stats::qnorm(1 - alpha/2)

  lower <- means - ses * z.score
  upper <- means + ses * z.score

  out <- cbind(exp(lower) - 1, exp(upper) - 1)
  rownames(out) <- names(means)
  colnames(out) <- c("Lower", "Upper")
  return(out)
}



# Sub-function of microsynth(); calculate confidence intervals
make.ci2 <- function(stats, delta.out, alpha = 0.05) {
  means <- stats[1, ]
  sds <- sqrt(delta.out[1, ])

  z.scores <- stats/sqrt(delta.out)
  quants <- apply(z.scores, 2, stats::quantile, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)

  lower <- means - sds * quants[2, ]
  upper <- means - sds * quants[1, ]

  out <- cbind(exp(lower) - 1, exp(upper) - 1)
  rownames(out) <- colnames(stats)
  colnames(out) <- c("Lower", "Upper")
  return(out)
}



# Sub-function of microsynth(); allow parellalization
msCluster <- function(n) {

  requireNamespace("parallel", quietly = TRUE)

  if (is.logical(n)) {
    if(n) {
      n1 <- parallel::detectCores()
      # n1 <- future::availableCores()
    } else {
      n1 <- 1
    }
  } else if (length(n) == 0) {
    n1 <- parallel::detectCores() - 1
    # n1 <- future::availableCores() - 1
  } else {
    n1 <- n
  }

  # use 2 cores if CRAN, otherwise, use number detected
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  cores <- ifelse(chk==TRUE, 2, n1)

  return(cores)
}



# Sub-function of microsynth()
remove.vars <- function(vars, nams, objnam = "result.var", printFlag = TRUE) {
  if (is.list(vars)) {
    vars1 <- names(vars)
  } else {
    vars1 <- vars
  }

  rm.vars <- setdiff(vars1, nams)

  if (length(rm.vars) > 0) {
    rm.here <- which(is.element(vars1, rm.vars))
    vars <- vars[-rm.here]
    if(printFlag){message("WARNING: The following variables will be removed from ", objnam,
                          " since they are not in the dataset: \n",
                          sep = "", appendLF = FALSE)}
    if(printFlag){message(paste(rm.vars, collapse = ", ", sep = ""), "\n\n",
                          sep = "", appendLF = FALSE)}
  }

  return(vars)
}



# Sub-function of microsynth()
check.matchout <- function (match.out, end.pre) {
  bad <- FALSE
  for (i in 1:length(match.out)) {
    tmp <- match.out[[i]]
    if (end.pre < sum(tmp)) {
      bad <- TRUE
      cum.tmp <- cumsum(tmp)
      cum.tmp <- cum.tmp <= end.pre
      match.out[[i]] <- c(tmp[cum.tmp])
      if (end.pre - sum(tmp[cum.tmp]) > 0) {
        match.out[[i]] <- c(match.out[[i]], end.pre - sum(tmp[cum.tmp]))
      }
    }
  }
  return(list(match.out, bad))
}

