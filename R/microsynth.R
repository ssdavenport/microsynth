#' Synthetic control methods for micro- and meso-level data.
#'
#' Implements the synthetic control method for micro-level data as outlined in
#' Robbins, Saunders, and Kilmer (2017).  \code{microsynth} is a generalization
#' of \code{Synth} (see Abadie and Gardeazabal (2003) and Abadie, Diamond,
#' Hainmueller (2010, 2011, 2014)) that is designed for data at a more granular
#' level (e.g., micro-level). \code{microsynth} may also be used to calculate
#' propensity score-type weights. For more details see the help vignette:
#' \code{vignette('microsynth', package = 'microsynth')}
#'
#' \code{microsynth} develops a synthetic control group by searching for weights
#' that exactly match the treatment group to the synthetic control group across
#' a number of variables while also minimizing the discrepancy between the
#' synthetic control group and the treatment group across a set second set of
#' variables.  \code{microsynth} works in three primary steps: 1) calculation of
#' weights, 2) plotting of treatment vs. synthetic control for pertinent
#' outcomes, and 3) calculation of results.
#' Variables are categorized as outcomes (which are time-variant) and covariates
#' (which are time-invariant).  Using the respective inputs \code{match.covar}
#' and \code{match.out}, the user specifies across which covariates and outcomes
#' (and which pre-intervention time points of the outcomes) treatment is to be
#' exactly matched to synthetic control.  The inputs \code{match.covar.min} and
#' \code{match.out.min} are similar but instead specify variables across which
#' treatment is to be matched to synthetic control as closely as possible.  If
#' there are no variables specified in \code{match.covar.min} and
#' \code{match.out.min}, the function \code{calibrate()} from the \code{survey}
#' package is used to calculate weights otherwise, the function
#' \code{LowRankQP()} from the package of the same name is used.  In the event
#' that the model specified by \code{match.covar} and \code{match.out} is not
#' feasible (i.e., weights do not exist that exactly match treatment and
#' synthetic control subject to the given constraints), a less restrictive
#' backup model is used. \code{microsynth} has the capability to perform
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
#' The software provides the user the option to create time series plots of
#' outcomes for both the treatment and outcome groups during the pre- and
#' post-intervention time periods and to output overall findings in an Excel
#' file.  For each outcome variable, the results list the estimated treatment
#' effect, as well as confidence intervals of the effect and p-values of a
#' hypothesis test that assesses whether the effect is zero.   Such results are
#' produced as needed for each of the three methods of statistical inference
#' noted above.  \code{microsynth} can also apply an omnibus test that examines
#' the presence of a treatment effect jointly across several outcomes.
#'
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
#' @param data A longitudinal data frame.  The data must be entered in tall
#'   format (e.g., at the case/time-level with one row for each time period for
#'   each case).  Missingness is not allowed.  All individuals must have non-NA
#'   values of all variables at all time points.
#'
#' @param idvar A character string that gives that gives the variable in
#'   \code{data} that identifies multiple records from the same case.
#'
#' @param intvar A character string that gives the variable in \code{data} that
#'   corresponds to the intervention variable.  The intervention variable
#'   indicates which cases (and times) have received the intervention.  The
#'   variable should be binary, with a 1 indicating treated.  If \code{int.time}
#'   is specified, a case is considered treated if there is 1 or more non-zero
#'   entries in the column indicated by \code{intvar} for that case (at any time
#'   point).  If \code{int.time} is not specified, an attempt will be made to
#'   use\code{intvar} to determine which time the intervention occurred (i.e.,
#'   the intervention is considered to have occurred at the time of the first
#'   non-zero entry in \code{intvar}).
#'
#' @param timevar A character string that gives the variable in
#'   \code{data} that differentiate multiple records from the same case.  Can be
#'   set to \code{NULL} only when used with cross-sectional data (i.e., with one
#'   observation per entry in \code{idvar}).
#'
#' @param w A list of the form as returned by a prior application of
#'   \code{microsynth}.  If \code{w = NULL}, weights are calculated from
#'   scratch.  Entering a \code{non-NULL} value affords the user the ability to
#'   use previously calculated  weights.
#'
#' @param match.out Either A) logical, B) a vector of variable names that
#'   indicates across which time-varying variables treatment is to be exactly matched
#'   to synthetic control pre-intervention, or, to allow more flexibility, or C) a
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
#'   time point from \code{start.time} to \code{int.time}. Otherwise, to allow
#'   more flexibility, \code{match.out} may also be a list that gives
#'   an outcome-based model outlining more specific constraints that are to be
#'   exactly satisfied within calibration weighting.  In this case, each entry
#'   of \code{match.out} is a vector of integers, and the names of entries of
#'   \code{match.out} are the outcome variables to which the vectors correspond.
#'   Each element of the vectors gives a number of time points that are to be
#'   aggregated for the respective outcome, with the first element indicating
#'   time points immediately prior the time of the intervention.  The sum of
#'   the elements in each vector should not exceed the number of
#'   pre-intervention time periods in the data.
#'
#'   The following examples show the proper formatting of \code{match.out} as a
#'   list.  Assume that there are two outcomes, Y1 and Y2 (across which
#'   treatment is to be matched to synthetic control), and that the intervention
#'   occurs at time 10.   Let \code{match.out = list('Y1' = c(1, 3, 3), 'Y2'=
#'   c(2,5,1))}.  According to this specification, treatment is to be matched to
#'   synthetic control across: a) The value of Y1 at time 10; b) the sum of Y1
#'   across times 7, 8 and 9; c) the sum of Y1 across times 4, 5 and 6; e) The
#'   sum of Y2 across times time 9 and 10; e) the sum of Y2 across times 4, 5,
#'   6, 7, and 8; f) the value of Y2 at time 3.  Likwise, if \code{match.out =
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
#' @param int.time A logical, or an integer that gives the final time point at
#'   which treatment and synthetic control will be matched to one another.  If an
#'   an integer, \code{int.time} indicates the start of the evaluation period.

#'   If the intervention is expected to have an
#'   instantaneous effect, \code{int.time} should be set to when the intervention
#'   occurred. Otherwise, \code{int.time} should be set to the time immediately
#'   prior to the intervention. If \code{intvar} is formatted such that
#'   all treatment cases receive 0 pre-intervention and 1 post-intervention, and
#'   all control cases receive 0 at all time points, then \code{int.time} may
#'   be logical. Setting \code{int.time = TRUE} indicates that effects are
#'   expected to be instantaneous, i.e., \code{int.time} is set to the time of
#'   the last 0 observed among the treatment column. Setting
#'   \code{int.time = FALSE} or \code{int.time = NULL} will set \code{int.time}
#'   as the time of the first 1 in the treatment column.
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
#'   for constraints that are to be exactly satisfied.  If max.mse is not
#'   satisfied by these constraints, and either \code{check.feas = TRUE} or
#'   \code{use.backup = TRUE}, then back-up models are used.
#'
#' @param result.var A vector of variable names giving the outcome
#'   variables for which results will be reported.  Time-varying covariates
#'   should be excluded from \code{result.var}.  If \code{result.var = TRUE}
#'   (the default), \code{result.var} is set as being equal to all time-varying
#'   variables that appear in \code{data}.  If \code{result.var = NULL} or
#'   \code{result.var = FALSE}, results are not tabulated.
#' @param omnibus.var A vector of variable names that indicates the outcome
#'   variables that are to be used within the calculation of the omnibus
#'   statistic.  Can also be a logical indicator.  When \code{omnibus.var =
#'   TRUE}, it is reset as being equal to \code{result.var}.  When
#'   \code{omnibus.var = NULL} or \code{omnibus = FALSE}, no omnibus statistic
#'   is calculated.
#' @param max.time An integer that gives the maximum post-intervention time that
#'   is taken into when compiling results.  That is, the treatment and synthetic
#'   control groups are compared across the outcomes listed in \code{result.var}
#'   from the first time following the intervention up to \code{max.time}.  Can
#'   be a vector (ordered, increasing) giving multiple values of
#'   \code{max.time}.  In this case, the results will be compiled for each entry
#'   in \code{max.time}.  When \code{max.time = NULL} (the default), it is reset
#'   to the maximum time that appears in the column given by \code{timevar}.
#' @param period An integer that gives the granularity of the data that will be
#'   used for plotting and compiling results; if \code{match.out} and
#'   \code{match.out.min} are provided a vector of variable names, it will also
#'   affect the calculation of weights used for matching. In this case, matching
#'   of treatment and synthetic control is performed at a temporal granularity
#'   defined by \code{period}. For instance, if monthly data are provided and
#'   \code{period = 3}, data are aggregated to quarters for plots and results
#'   (and weighting unless otherwise specified). If \code{match.out} and
#'   \code{match.out.min} are provided a list, \code{period} only affects plots
#'   and how results are displayed.
#' @param cut.mse The maximum error (given as mean-squared error) permissible
#'   for permutation groups.  Permutation groups with a larger than permissible
#'   error are dropped when calculating results.  The mean-squared error is only
#'   calculated over constraints that are to be exactly satisfied.
#' @param test The type of hypothesis test (one-sided lower, one-sided upper, or
#'   two-sided) that is used when calculating p-values.  Entries of
#'   \code{'lower'}, \code{'upper'}, and \code{'twosided'} are recognized.
#' @param result.file A character string giving the name of a file that will be
#'   created in the home directory containing results.  If \code{result.file =
#'   NULL} (the default), no file is created.  If \code{max.time} has length 1,
#'   a \code{.csv} file is created.  If \code{max.time} has length greater than
#'   one, a formatted \code{.xlsx} file is created with one tab for each element
#'   of \code{max.time}.  If \code{result.file} has a \code{.xlsx} (or
#'   \code{.xls}) extension (e.g., the last five characters of result.file are
#'   '.xlsx'), an \code{.xlsx} file is created regardless of the length of
#'   \code{max.time}.
#' @param use.survey If \code{use.survey = TRUE}, Taylor series linearization is
#'   applied to the estimated treatment effect within each permutation group.
#'   Setting \code{use.survey = TRUE} makes for better inference but increases
#'   computation time substantially.  Confidence intervals for permutation
#'   groups are calculated only when \code{use.survey = TRUE}.
#'
#' @param confidence The level of confidence for confidence intervals.
#'
#' @param plot.var A vector of variable names giving the outcome variables that
#'   are shown in plots.  If \code{plot.var = NULL} or \code{plot.var = FALSE}, no
#'   plots are generated.  If \code{plot.var = TRUE}, it is reset as equaling all
#'   time variant variables in \code{data}.
#' @param plot.file A character string giving the name of file that will be
#'   created in the home directory containing plots (if \code{plot.var} is
#'   non-\code{NULL}).  The name should have a \code{.pdf} extension.
#'
#' @param start.time An integer indicating the earliest time shown in the plots;
#'   if match.out and match.out.min are provided a vector of variable names, it
#'   will also determine the beginning of the pre-intervention period used for
#'   matching. When \code{start.time = NULL} (default), it is reset to the
#'   minimum time appearing in the column given by \code{timevar}.
#'
#' @param sep If \code{sep = TRUE}, separate plots will be generated for each
#'   outcome.  Applicable only if plotting variables are specified (
#'   \code{plot.var} is \code{non-NULL}) and plots are saved to file (
#'   \code{plot.file} is \code{non-NULL}). To change display of plots produced
#'   as output, use \code{\link[graphics]{par}}.
#'
#' @param legend.spot The location of the legend in the plots.
#'
#' @param scale.var  A variable name.  When comparing the treatment group to all
#'   cases, the latter is scaled to the size of the former with respect to the
#'   variable indicated by \code{scale.var}.  Defaults to the number of units
#'   receiving treatment (i.e., the intercept).
#'
#' @param maxit The maximum number of iterations used within the calibration
#'   routine (\code{survey::calibrate()} from the \code{survey} package) for
#'   calculating weights.
#' @param cal.epsilon The tolerance used within the calibration routine
#'   (\code{survey::calibrate()} from the \code{survey} package) for calculating
#'   weights.
#' @param calfun The calibration function used within the calibration routine
#'   (\code{survey::calibrate()} from the \code{survey} package) for calculating
#'   weights.
#' @param bounds Bounds for calibration weighting (fed into the
#'   \code{survey::calibrate()} from the \code{survey} package).
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
#' @return \code{microsynth} returns a list with up to four elements: a)
#'   \code{w}, b) \code{Results}, c) \code{svyglm.stats},
#'   and c) \code{Plot.Stats}.  The fourth element is returned only if
#'   \code{plot.var} is not \code{NULL} or \code{FALSE}.
#'
#'   \code{w} is a list with five elements: a) \code{Weights}, b)
#'   \code{Intervention},
#'   c) \code{MSE}, d) \code{Model}, and e) \code{Summary}.  Assume there are
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
#'   weights.  \code{w$Summary} is a three-column matrix that (for treatment,
#'   synthetic control, and the full dataset), shows aggregate values
#'   of the variables across which treatment and synthetic control are matched.
#'   The summary, which is tabulated only for the primary weights, is also
#'   printed by microsynth while weights are being calculated.
#'
#'   Further, \code{Results} is a list where each element gives the final
#'   results for each value of \code{max.time}.  Each element of \code{Results}
#'   is itself a matrix with each row corresponding to an outcome variable (and
#'   a row for the omnibus test, if used) and each column denotes estimates of
#'   the intervention effects and p-values, upper, and lower bounds of
#'   confidence intervals as found using Taylor series linearization (Linear),
#'   jackknife (jack), and permutation (perm) methods where needed.
#'
#'   In addition, \code{svyglm.stats} is a list where each element is a
#'   matrix that includes the output from the regression models run using the
#'   \code{svyglm()} function to estimate the treatment effect.  The list has one
#'   element for each value of \code{max.time}, and the matrices each have
#'   one row per variable in \code{result.var}.
#'
#'   Lastly, \code{Plot.Stats} contains the data that are displayed in the
#'   plots.  \code{Plot.Stats} is a list with four elements (Treatment, Control,
#'   All, Difference).  The first three elements are matrices with one row per
#'   outcome variable and one column per time point.  The last element (which
#'   gives the treatment minus control values) is an array that contains data
#'   for each permutation group in addition to the true treatment area.
#'   Specifically, \code{Plot.Stats$Difference[,,1]} contains the time series of
#'   treatment minus control for the true intervention group;
#'   \code{Plot.Stats$Difference[,,i+1]} contains the time series of treatment
#'   minus control for the i^th permutation group.
#'
#'   A summary of weighted matching variables and of results can be viewed
#'   using \code{\link{summary}}
#'
#' @references Abadie A, Diamond A, Hainmueller J (2010). “Synthetic control
#'   methods for comparative case studies: Estimating the effect of California’s
#'   tobacco control program.” \emph{Journal of the American Statistical
#'   Association}, 105(490), 493-505.
#'
#'   Abadie A, Diamond A, Hainmueller J (2011). “Synth: An R Package for
#'   Synthetic Control Methods in Comparative Case Studies.” \emph{Journal
#'   of Statistical Software}, 42(13), 1-17.
#'
#'   Abadie A, Diamond A, Hainmueller J (2015). “Comparative politics and the
#'   synthetic control method. \emph{American Journal of Political Science},
#'   59(2), 495-510.
#'
#'   Abadie A, Gardeazabal J (2003). “The economic costs of conflict: A case
#'   study of the Basque Country.” \emph{American Economic Review}, pp. 113-132.
#'
#'   Hainmueller, J. (2012), “Entropy Balancing for Causal Effects: A
#'   Multivariate Reweighting Method to Produce Balanced Samples in
#'   Observational Studies,” \emph{Political Analysis}, 20, 25–46.
#'
#'   Robbins MW, Saunders J, Kilmer B (2017). “A framework for synthetic control
#'   methods withhigh-dimensional, micro-level data: Evaluating a neighborhood-
#'   specific crime intervention,” \emph{Journal of the American Statistical
#'   Association}, 112(517), 109-126.
#'
#' @examples
#' set.seed(99199) # for reproducibility
#'
#' # Use seattledmi, block-level panel data, to evaluate a crime intervention.
#'
#' # Declare time-variant (outcome) and time-invariant variables for matching
#' cov.var <- c('TotalPop', 'BLACK', 'HISPANIC', 'Males_1521',
#'        'HOUSEHOLDS', 'FAMILYHOUS', 'FEMALE_HOU', 'RENTER_HOU', 'VACANT_HOU')
#' match.out <- c('i_felony', 'i_misdemea', 'i_drugs', 'any_crime')

#' # Perform matching and estimation, without permutations or jackknife
#' sea1 <- microsynth(seattledmi, idvar='ID', timevar='time',
#'        intvar='Intervention', start.time=1, int.time=NULL, max.time=16,
#'        match.out=match.out, match.covar=cov.var, result.var=match.out,
#'        omnibus.var=match.out, plot.var=match.out, test='lower')
#'
#' # View results
#' summary(sea1)
#'
#' # Repeat matching and estimation, with permutations and jackknife
#' # Set permutations and jack-knife to very few groups (2) for
#' # quick demonstration only.
#' sea2 <- microsynth(seattledmi, idvar='ID', timevar='time',
#'         intvar='Intervention', start.time=1, int.time=12, max.time=c(14, 16),
#'         match.out=match.out, match.covar=cov.var, result.var=match.out,
#'         omnibus.var=match.out, plot.var=match.out, test='lower', perm=250,
#'         jack=TRUE, plot.file=NULL, sep = TRUE, result.file='ExResults2.xlsx')
#'
#' # View results
#' summary(sea2)
#'
#' # Specify additional outcome variables for matching, which makes
#' # matching harder.
#' match.out <- c('i_robbery','i_aggassau','i_burglary','i_larceny',
#'        'i_felony','i_misdemea','i_drugsale','i_drugposs','any_crime')
#'
#' # Perform matching, setting check.feas = T and use.backup = T
#' # to ensure model feasibility
#' sea3 <- microsynth(seattledmi, idvar='ID', timevar='time',
#'         intvar='Intervention', match.out=match.out, match.covar=cov.var,
#'         result.var=match.out, plot.var=match.out, perm=250, jack=TRUE,
#'         test='lower', check.feas=TRUE, use.backup = TRUE,
#'         plot.file=NULL, result.file='ExResults3.xlsx')
#'
#' # Aggregate outcome variables before matching, to boost model feasibility
#' match.out <- list( 'i_robbery'=rep(2, 6), 'i_aggassau'=rep(2, 6),
#'          'i_burglary'=rep(1, 12), 'i_larceny'=rep(1, 12),
#'          'i_felony'=rep(2, 6), 'i_misdemea'=rep(2, 6),
#'          'i_drugsale'=rep(4, 3), 'i_drugposs'=rep(4, 3),
#'          'any_crime'=rep(1, 12))
#'
#' # After aggregation, use.backup and cheack.feas no longer needed
#' sea4 <- microsynth(seattledmi, idvar='ID', timevar='time',
#'          intvar='Intervention', match.out=match.out, match.covar=cov.var,
#'          result.var=names(match.out), omnibus.var=names(match.out),
#'          plot.var=names(match.out), perm=250, jack = 0, test='lower',
#'          plot.file='ExPlots4.pdf', result.file='ExResults4.xlsx')
#'
#' # Generate weights only (for four variables)
#' match.out <- c('i_felony', 'i_misdemea', 'i_drugs', 'any_crime')
#' sea5 <- microsynth(seattledmi,  idvar='ID', timevar='time',
#'          intvar='Intervention', match.out=match.out, match.covar=cov.var,
#'          result.var=FALSE, plot.var=FALSE, perm=250, jack=TRUE)
#'
#' # View weights
#' summary(sea5)
#'
#' # Generate plots only using previous weights
#' sea6 <- microsynth(seattledmi,  idvar='ID', timevar='time',
#'           intvar='Intervention', result.var=FALSE, plot.var=match.out[1:2],
#'           w=sea5$w)
#'
#' # Generate results only
#' sea7 <- microsynth(seattledmi, idvar='ID', timevar='time',
#'           intvar='Intervention', max.time=c(14, 16),
#'           result.var=match.out, plot.var=FALSE, test='lower',
#'           w=sea5$w, result.file='ExResults7.xlsx')
#'
#' # View results (including previously-found weights)
#' summary(sea7)
#'
#' # Apply microsynth in the traditional setting of Synth
#' # Create macro-level (small n) data, with 1 treatment unit
#' set.seed(66872)
#' ids.t <- names(table(seattledmi$ID[seattledmi$Intervention==1]))
#' ids.c <- names(table(seattledmi$ID[seattledmi$Intervention==0]))
#' ids.synth <- c(sample(ids.t, 1), sample(ids.c, 100))
#' seattledmi.one <- seattledmi[is.element(seattledmi$ID,
#'            as.numeric(ids.synth)), ]
#'
#' # Apply microsynth to the new macro-level data
#' sea8 <- microsynth(seattledmi.one, idvar='ID', timevar='time',
#'            intvar='Intervention', match.out=match.out[4],
#'            match.covar=cov.var, result.var=match.out[4],
#'            plot.var=match.out[4], test='lower', perm=250, jack=FALSE,
#'            check.feas=TRUE, use.backup=TRUE)
#'
#' # Use microsynth to calculate propensity score-type weights
#' # Prepare cross-sectional data at time of intervention
#' seattledmi.cross <- seattledmi[seattledmi$time==12, ]
#' seattledmi.cross <- seattledmi.cross[, colnames(seattledmi)!='time']
#'
#' # Apply microsynth to find propensity score-type weights
#' sea9 <- microsynth(seattledmi.cross, idvar='ID', intvar='Intervention',
#'              match.out=FALSE, match.covar=cov.var, result.var=match.out,
#'              test='lower', perm=250, jack=TRUE)
#'
#' # View results
#' summary(sea9)
#'
#'
#' @export

microsynth <- function (data, idvar, intvar, timevar = NULL, start.time = NULL,
                        int.time = NULL, max.time = NULL, match.out = TRUE, match.covar = TRUE,
                        match.out.min = NULL, match.covar.min = NULL, result.var = TRUE,
                        omnibus.var = result.var, plot.var = TRUE, period = 1, scale.var = "Intercept",
                        confidence = 0.9, test = "twosided", perm = 0, jack = 0,
                        use.survey = TRUE, cut.mse = Inf, check.feas = FALSE, use.backup = FALSE,
                        w = NULL, max.mse = 0.01, maxit = 250, cal.epsilon = 1e-04,
                        calfun = "linear", bounds = c(0, Inf), result.file = NULL,
                        plot.file = NULL, sep = FALSE, legend.spot = "bottomleft")
{
  all.tmp <- proc.time()
  if (length(timevar) == 0) {
    if (length(table(data[, idvar])) < NROW(data)) {
      stop("Data are not cross-sectional.  Please specify timevar.")
    }
    else {
      data$Time <- 1
      timevar <- "Time"
      int.time <- 1
    }
  }
  time.tmp <- data[,timevar]
  time.names <- names(table(time.tmp))
  data[,timevar] <- match(as.character(time.tmp), time.names)
  if (length(start.time) > 0 & !is.logical(start.time)) {
    start.time <- match(as.character(start.time), time.names)
  }
  if (length(int.time) > 0 & !is.logical(int.time)) {
    int.time <- match(as.character(int.time), time.names)
  }
  if (length(max.time) > 0 & !is.logical(max.time)) {
    max.time <- match(as.character(max.time), time.names)
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
  if (is.logical(plot.var)) {
    if (!plot.var) {
      reset.plot.var <- TRUE
    }
    else {
      reset.plot.var <- FALSE
    }
    plot.var <- NULL
  }
  else if (length(plot.var) == 0) {
    reset.plot.var <- TRUE
  }
  else {
    reset.plot.var <- FALSE
  }
  match.out <- remove.vars(match.out, dimnames(data)[[2]],
                           "match.out")
  match.out.min <- remove.vars(match.out.min, dimnames(data)[[2]],
                               "match.out.min")
  match.covar <- remove.vars(match.covar, dimnames(data)[[2]],
                             "match.covar")
  match.covar.min <- remove.vars(match.covar.min, dimnames(data)[[2]],
                                 "match.covar.min")
  result.var <- remove.vars(result.var, dimnames(data)[[2]],
                            "out.covar")
  if (!is.logical(omnibus.var)) {
    omnibus.var <- remove.vars(omnibus.var, dimnames(data)[[2]],
                               "omnibus.var")
  }
  if (!is.logical(plot.var)) {
    plot.var <- remove.vars(plot.var, dimnames(data)[[2]],
                            "plot.var")
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
  data <- newreshape(data, nv.names = nv.names, v.names = v.names,
                     timevar = timevar, idvar = idvar, intvar = intvar)
  if (length(result.var) == 0) {
    result.var <- data[[4]]
    if (!reset.result.var) {
      message("result.var = TRUE.  Resetting: \n", appendLF = FALSE)
      message("result.var = c(\"", paste(result.var, collapse = "\",\"",
                                         sep = ""), "\")\n\n", sep = "", appendLF = FALSE)
    }
  }
  if (length(plot.var) == 0) {
    plot.var <- data[[4]]
    if (!reset.plot.var) {
      message("plot.var = TRUE.  Resetting: \n", appendLF = FALSE)
      message("plot.var = c(\"", paste(plot.var, collapse = "\",\"",
                                       sep = ""), "\")\n\n", sep = "", appendLF = FALSE)
    }
  }
  if (length(match.covar) == 0) {
    match.covar <- data[[3]]
    if (!reset.match.covar) {
      message("match.covar = TRUE.  Resetting: \n", appendLF = FALSE)
      message("match.covar = c(\"", paste(match.covar,
                                          collapse = "\",\"", sep = ""), "\")\n\n", sep = "",
              appendLF = FALSE)
    }
  }
  Intervention <- data[[2]]
  data <- data[[1]]
  times <- as.numeric(colnames(Intervention))
  if (length(int.time) == 0) {
    int.time <- FALSE
  }
  if (is.logical(int.time)) {
    int.time1 <- int.time
    int.time <- which(colSums(Intervention) != 0)
    if (length(int.time) == 0) {
      stop("There are no intervention cases.\n")
    }
    else {
      if (int.time1) {
        int.time <- min(times[int.time]) - 1
      }
      else {
        int.time <- min(times[int.time])
      }
      if (int.time == 0) {
        stop("There are no pre-intervention time points.")
      }
      message("Setting int.time = ", time.names[int.time], ".\n\n", sep = "",
              appendLF = FALSE)
    }
  } else {
    Intervention[] <- Intervention[, NCOL(Intervention)]
  }
  if (length(start.time) == 0) {
    start.time <- min(times)
  }
  if (length(max.time) == 0) {
    max.time <- max(times)
  }
  if (length(match.out) == 0) {
    match.out <- result.var
    if (!reset.match.out) {
      message("match.out = TRUE.  Resetting: \n", appendLF = FALSE)
      message("match.out = c(\"", paste(match.out, collapse = "\",\"",
                                        sep = ""), "\")\n\n", sep = "", appendLF = FALSE)
    }
  }
  if (!is.list(match.out) & length(match.out) > 0) {
    match.out.tmp <- match.out
    match.out <- list()
    for (i in 1:length(match.out.tmp)) {
      match.out[[i]] <- rep(period, (int.time - start.time +
                                       1)%/%period)
    }
    names(match.out) <- match.out.tmp
  } else if (length(match.out) > 0) {
    match.out <- check.matchout(match.out, int.time - min(times) + 1)
    if (match.out[[2]]) {
      message("WARNING: match.out calls on time periods that are beyond the data range.\n",
              sep ="", appendLF = FALSE)
      message("match.out is being reset accordingly.\n\n", sep="", appendLF = FALSE)
    }
    match.out <- match.out[[1]]
  }
  if (!is.list(match.out.min) & length(match.out.min) > 0) {
    match.out.tmp1 <- match.out.min
    match.out.min <- list()
    for (i in 1:length(match.out.tmp1)) {
      match.out.min[[i]] <- rep(period, (int.time - start.time +
                                           1)%/%period)
    }
    names(match.out.min) <- match.out.tmp1
  } else if (length(match.out.min) > 0) {
    match.out.min <- check.matchout(match.out.min, int.time - min(times) + 1)
    if (match.out.min[[2]]) {
      message("WARNING: match.out.min calls on time periods that are beyond the data range. \n",
              sep ="", appendLF = FALSE)
      message("match.out.min is being reset accordingly.\n\n", sep="", appendLF = FALSE)
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
        message("omnibus.var = TRUE.  Resetting: \n",
                appendLF = FALSE)
        message("omnibus.var = c(\"", paste(omnibus.var,
                                            collapse = "\",\"", sep = ""), "\")\n\n", sep = "",
                appendLF = FALSE)
      }
    }
    else {
      omnibus.var <- NULL
    }
  }
  int.num <- 1
  dum <- max(colSums(Intervention == 1))
  dum <- min(NROW(Intervention) - dum, dum)
  dum1 <- ((max(max.time) - (max(max.time) - int.time)%%period -
              int.time)/period)
  if (dum <= dum1 + 1) {
    message("WARNING: There is a low number (", dum, ") of cases in the treatment or intervention group.\n",
            sep = "", appendLF = FALSE)
    if (jack > 0) {
      jack <- 0
      message("setting jack = 0.\n", appendLF = FALSE)
    }
    if (use.survey) {
      use.survey <- FALSE
      message("Setting use.survey = FALSE.\n", appendLF = FALSE)
    }
    message("Be cautious of results involving linearization or confidence intervals.\n\n",
            appendLF = FALSE)
  }
  if (reset.match.covar) {
    match.covar <- NULL
  }
  if (reset.match.out) {
    match.out <- NULL
  }
  if (reset.plot.var) {
    plot.var <- NULL
  }
  if (length(rm.col) > 0) {
    for (i in 1:length(rm.col)) {
      if (is.element(rm.col[i], result.var)) {
        message("WARNING: ", rm.col[i], " is a non-numeric variable.  It will be removed from result.var. \n",
                sep = "", appendLF = FALSE)
        result.var <- setdiff(result.var, rm.col[i])
      }
      if (is.element(rm.col[i], omnibus.var)) {
        message("WARNING: ", rm.col[i], " is a non-numeric variable.  It will be removed from omnibus.var. \n",
                sep = "", appendLF = FALSE)
        omnibus.var <- setdiff(omnibus.var, rm.col[i])
      }
      if (is.element(rm.col[i], match.covar)) {
        message("WARNING: ", rm.col[i], " is a non-numeric variable.  It will be removed from match.covar. \n",
                sep = "", appendLF = FALSE)
        match.covar <- setdiff(match.covar, rm.col[i])
      }
      if (is.element(rm.col[i], match.covar.min)) {
        message("WARNING: ", rm.col[i], " is a non-numeric variable.  It will be removed from match.covar.min. \n",
                sep = "", appendLF = FALSE)
        match.covar.min <- setdiff(match.covar.min, rm.col[i])
      }
      if (is.element(rm.col[i], names(match.out))) {
        message("WARNING: ", rm.col[i], " is a non-numeric variable.  It will be removed from match.out. \n",
                sep = "", appendLF = FALSE)
        rm.li <- which(is.element(names(match.out), rm.col[i]))
        match.out <- match.out[-rm.li]
      }
      if (is.element(rm.col[i], names(match.out.min))) {
        message("WARNING: ", rm.col[i], " is a non-numeric variable.  It will be removed from match.out.min. \n",
                sep = "", appendLF = FALSE)
        rm.li1 <- which(is.element(names(match.out.min),
                                   rm.col[i]))
        match.out.min <- match.out.min[-rm.li1]
      }
    }
  }
  if (length(w) == 0) {
    tmp <- proc.time()
    message("Calculating weights...", "\n", appendLF = FALSE)
    w <- get.w(data, match.covar, match.covar.min, match.out,
               match.out.min, boot = perm, jack = jack, Int = Intervention[,
                                                                           as.character(int.time)], int.val = int.num, trim = NULL,
               int.time = int.time, cal.epsilon = cal.epsilon, maxit = maxit,
               bounds = bounds, calfun = calfun, check.feas = check.feas,
               scale.var = scale.var, cut.mse = max.mse, use.backup = use.backup, time.names = time.names)
    tmp <- proc.time() - tmp
    message("Calculation of weights complete: Total time = ",
            round(tmp[3], 2), "\n\n", sep = "", appendLF = FALSE)
  }
  else {
    message("Weights have been provided.  Will not calculate weights.\n",
            appendLF = FALSE)
    is.correct.w <- is.list(w)
    if (is.correct.w) {
      is.correct.w <- is.correct.w & sum(names(w) != c("Weights",
                                                       "Intervention", "MSE", "Model", "Summary")) ==
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
    message("Setting jack = ", jack, "\n", sep = "", appendLF = FALSE)
    message("Setting perm = ", perm, "\n\n", sep = "", appendLF = FALSE)
  }
  not.jack <- !grepl("Jack", colnames(w$Weights))
  if (is.logical(plot.var)) {
    if (plot.var) {
      plot.var <- result.var
    }
    else {
      plot.var <- NULL
    }
  }
  else {
    plot.var <- intersect(plot.var, dimnames(data)[[2]])
  }
  max.time <- max.time - (max.time - int.time)%%period
  if (length(plot.var) > 0 & length(times) == 1) {
    plot.var <- NULL
    message("There is only one time point in the data.\n            Will not generate plots.\n\n",
            appendLF = FALSE)
  }
  if (!reset.result.var | !reset.plot.var) {
    stats <- list()
    stats1 <- list()
    stats2 <- list()
    delta.out <- list()
    dof <- list()
    out.coefs <- list()
    results <- list()
    for (i in 1:length(max.time)) {
      tmp <- proc.time()
      is.graph <- is.graph1 <- ""
      if (!reset.plot.var & !reset.result.var) {
        is.graph <- "Making graphs and calculating basic statistics"
        is.graph1 <- "Completed graphs and calculation of basic statistics"
      }
      else if (reset.plot.var & !reset.result.var) {
        is.graph <- "Calculating basic statistics"
        is.graph1 <- "Completed calculation of basic statistics"
      }
      else if (!reset.plot.var & reset.result.var) {
        is.graph <- "Making graphs"
        is.graph1 <- "Completed graphs"
      }
      message(is.graph, " for max.time = ", time.names[max.time[i]],
              "...", "\n", sep = "", appendLF = FALSE)
      stats[[i]] <- get.stats(data, w$Weights, w$Intervention,
                              w$MSE[1, ], result.var, int.time = int.time,
                              period = period, plot.it = plot.var, max.time = max.time[i],
                              file = plot.file, omnibus.var = omnibus.var,
                              sep = sep, start.time = start.time, legend.spot = legend.spot,
                              cut.mse = cut.mse, twosided = twosided, time.names = time.names)
      if (i == which.max(max.time)) {
        plot.stats <- stats[[i]][[5]]
      }
      stats[[i]] <- stats[[i]][-5]
      tmp <- proc.time() - tmp
      message(is.graph1, " for max.time = ", time.names[max.time[i]],
              ".  Time = ", round(tmp[3], 2), "\n\n", sep = "",
              appendLF = FALSE)
      if (!reset.result.var) {
        w.tmp <- w$Weights
        Inter.tmp <- w$Intervention
        mse.tmp <- w$MSE[1, ]
        if (!use.survey) {
          keep.surv <- !grepl("Perm", colnames(w.tmp))
          w.tmp <- w.tmp[, keep.surv, drop = FALSE]
          Inter.tmp <- Inter.tmp[, keep.surv, drop = FALSE]
          mse.tmp <- mse.tmp[keep.surv]
        }
        tmp <- proc.time()
        message("Calculating survey statistics for max.time = ",
                time.names[max.time[i]], "...", "\n", sep = "", appendLF = FALSE)
        stats.tmp <- get.stats1(data, w.tmp, Inter.tmp,
                                mse.tmp, result.var, int.time = int.time, period = period,
                                max.time = max.time[i], omnibus.var = omnibus.var,
                                cut.mse = cut.mse, twosided = twosided)
        stats1[[i]] <- stats.tmp[[1]]
        stats2[[i]] <- stats.tmp[[2]]
        delta.out[[i]] <- stats.tmp[[3]]
        dof[[i]] <- stats.tmp[[4]]
        out.coefs[[i]] <- stats.tmp[[5]]
        tmp <- proc.time() - tmp
        message("Completed calculation of survey statistics for max.time = ",
                time.names[max.time[i]], ".  Time = ", round(tmp[3], 2),
                "\n\n", sep = "", appendLF = FALSE)
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
        names(results)[i] <- names(out.coefs)[i] <- time.names[max.time[i]]
      }
    }
  }
  if (reset.plot.var) {
    message("No plotting variables specified (e.g., plot.var = NULL).\n",
            appendLF = FALSE)
    message("Plots will not be created.\n", appendLF = FALSE)
  }
  if (reset.result.var) {
    message("No outcome variables specified (e.g., result.var = NULL).\n",
            appendLF = FALSE)
    message("Results will not be tablulated.\n", appendLF = FALSE)
  }
  if (reset.result.var & reset.plot.var) {
    message("Returning weights only.\n", appendLF = FALSE)
  }
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
  if (!reset.plot.var) {
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
    out.results(results, int.time, max.time = names(results), result.file)
  }
  all.tmp <- proc.time() - all.tmp
  message("microsynth complete: Overall time = ", round(all.tmp[3],
                                                        2), "\n\n", sep = "", appendLF = FALSE)
  out <- makemicrosynth(out)
  return(out)
}


my.delta <- function(mu, Sigma) {
  inds1 <- which(rownames(Sigma) == "treat")
  inds2 <- which(rownames(Sigma) != "treat")
  inds <- c(inds1, inds2)
  mu <- mu[inds]
  Sigma <- Sigma[inds, inds]
  mu1 <- mu[1]
  mu2 <- mu[2:length(mu)]
  m <- length(mu2)
  MU2 <- mean(mu2)
  D <- rep(NA, length(mu))
  D[1] <- 1/(mu1 + MU2)
  D[2:length(D)] <- (D[1] - (1/MU2))/m
  D <- as.matrix(D)
  out <- t(D) %*% Sigma %*% D
  return(c(out))
}


make.ci <- function(means, ses, alpha = 0.05) {
  z.score <- stats::qnorm(1 - alpha/2)

  lower <- means - ses * z.score
  upper <- means + ses * z.score

  out <- cbind(exp(lower) - 1, exp(upper) - 1)
  rownames(out) <- names(means)
  colnames(out) <- c("Lower", "Upper")
  return(out)
}


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


get.w <- function (bigdat, covar.var, covar.var1 = NULL, dum, dum1 = NULL,
                   boot = 0, jack = 0, Int, int.val = 1, trim = NULL, maxit = 500,
                   cal.epsilon = 1e-04, int.time, bounds = c(-Inf, Inf), calfun = "raking",
                   qpmeth = "LowRankQP", check.feas = FALSE, use.backup = TRUE,
                   scale.var = "Intercept", cut.mse = 1, time.names = NULL)
{
  n <- dim(bigdat)[1]
  n.int <- sum(Int == int.val)
  back.state <- back.state1 <- back.state2 <- back.state3 <- ""
  fin.boots <- FALSE
  if (boot > 0) {
    n.choose <- choose(n, n.int)
    if (boot > n.choose - 1) {
      boot <- n.choose - 1
      message("Resetting perm = ", boot, "\n", sep = "",
              appendLF = FALSE)
    }
    if (n.choose - 1 <= max(1e+06, boot)) {
      boots <- utils::combn(1:n, n.int)
      check.combn <- function(x, y) {
        return(sum(!is.element(x, y)))
      }
      is.trt.area <- which(apply(boots, 2, check.combn,
                                 x = which(Int == int.val)) == 0)
      boots <- boots[, sample((1:NCOL(boots))[-is.trt.area],
                              boot), drop = FALSE]
      fin.boots <- TRUE
    }
  }
  use.model <- 1
  if (check.feas & use.backup) {
    tmp <- proc.time()
    message("Checking feasibility of first model...", "\n",
            sep = "", appendLF = FALSE)
    is.sol <- is.feasible(bigdat, covar.var, dum, Int = Int,
                          int.val = int.val, int.time = int.time, eps = 1e-04)
    tmp <- proc.time() - tmp
    if (!is.sol) {
      use.model <- 2
      message("First model is infeasible: Time = ", round(tmp[3],
                                                          2), "\n", sep = "", appendLF = FALSE)
      dum.tmp <- merge.dums(dum, dum1)
      tmp <- proc.time()
      message("Checking feasibility of second model...",
              "\n", sep = "", appendLF = FALSE)
      is.sol <- is.feasible(bigdat, covar.var, dum.tmp[[1]],
                            Int = Int, int.val = int.val, int.time = int.time,
                            eps = 1e-04)
      tmp <- proc.time() - tmp
      if (!is.sol) {
        use.model <- 3
        message("Second model is infeasible: Time = ",
                round(tmp[3], 2), "\n", sep = "", appendLF = FALSE)
        message("Will use third model.", "\n", sep = "",
                appendLF = FALSE)
      }
      else {
        message("Second model is feasible: Time = ",
                round(tmp[3], 2), "\n", sep = "", appendLF = FALSE)
      }
    }
    else {
      message("First model is feasible: Time = ", round(tmp[3],
                                                        2), "\n", sep = "", appendLF = FALSE)
    }
  }
  newdat <- get.newdat(bigdat, dum = dum, dum1 = dum1, covar.var = covar.var,
                       covar.var1 = covar.var1, int.time = int.time, time.names = time.names)
  newdat1 <- newdat[[2]]
  newdat <- newdat[[1]]
  duma <- merge.dums(dum, dum1)
  newdata <- get.newdat(bigdat, dum = duma[[1]], dum1 = duma[[2]], covar.var = covar.var,
                        covar.var1 = covar.var1, int.time = int.time, time.names = time.names)
  newdat1a <- newdata[[2]]
  newdata <- newdata[[1]]
  dumb <- merge.dums(duma[[1]], duma[[2]])
  covar.var1b <- union(covar.var, covar.var1)
  covar.varb <- NULL
  newdatb <- get.newdat(bigdat, dum = dumb[[1]], dum1 = dumb[[2]], covar.var = covar.varb,
                        covar.var1 = covar.var1b, int.time = int.time, time.names = time.names)
  newdat1b <- newdatb[[2]]
  newdatb <- newdatb[[1]]
  colnam <- "Main"
  if (is.logical(jack)) {
    if (jack) {
      jack <- min(table(Int == int.val))
    }
    else {
      jack <- 0
    }
  }
  if (jack > min(table(Int == int.val))) {
    jack <- min(table(Int == int.val))
    message("Resetting jack = ", jack, "\n", sep = "", appendLF = FALSE)
  }
  if (jack > 0) {
    rep.G <- assign.groups(Int == int.val, G = jack)
    colnam <- c(colnam, paste("Jack", 1:jack, sep = ""))
  }
  if (boot > 0) {
    colnam <- c(colnam, paste("Perm", 1:boot, sep = ""))
  }
  wghts <- Inter <- matrix(NA, n, boot + jack + 1)
  mse <- matrix(NA, 6, boot + jack + 1)
  mod <- rep(NA, boot + jack + 1)
  rownames(wghts) <- rownames(Inter) <- dimnames(bigdat)[[1]]
  rownames(mse) <- c("First model: Primary variables", "First model: Second variables",
                     "Secondary model: Primary variables", "Second model: Secondary variables",
                     "Third model: Primary variables", "Third model: Secondary variables")
  names(mod) <- colnames(wghts) <- colnames(Inter) <- colnames(mse) <- colnam
  inds <- 1:(boot + jack + 1)
  jack.lower <- 2
  jack.upper <- 1 + jack
  boot.lower <- 2 + jack
  boot.upper <- 1 + jack + boot
  for (i in inds) {
    use.model.i <- use.model
    tmp <- proc.time()
    if (i == 1) {
      samp <- Int == int.val
      use <- rep(TRUE, n)
    }
    else if (grepl("Jack", colnam[i])) {
      samp <- Int == int.val
      g <- as.numeric(gsub("Jack", "", colnam[i]))
      use <- rep.G != g
      rep.meth <- "jackknife"
    }
    else if (grepl("Perm", colnam[i])) {
      g <- as.numeric(gsub("Perm", "", colnam[i]))
      if (!fin.boots) {
        samp <- sample(1:n, sum(Int == int.val), replace = FALSE,
                       prob = NULL)
        samp <- is.element(1:n, samp)
      }
      else {
        samp <- is.element(1:n, boots[, g])
      }
      use <- rep(TRUE, n)
      rep.meth <- "permutation"
    }
    Inter[samp & use, i] <- TRUE
    Inter[!samp & use, i] <- FALSE
    if (use.model.i == 1) {
      mod[i] <- "First"
      ws <- get.w.sub(newdat = newdat, newdat1 = newdat1,
                      int.time = int.time, samp = samp, use = use,
                      n = NROW(newdat), maxit = maxit, calfun = calfun,
                      bounds = bounds, epsilon = cal.epsilon, trim = trim,
                      qpmeth = qpmeth, scale.var = scale.var)
      if (ws$mse > cut.mse | is.na(ws$mse)) {
        if (use.backup) {
          cat.back <- ".  Using second model."
          cat.back1 <- ".  Used second model."
        }
        else {
          cat.back <- cat.back1 <- ".  Consider setting use.backup = TRUE."
        }
        if (i > 1) {
          back.state1 <- paste("First model was infeasible for ",
                               rep.meth, " group ", g, cat.back1, "\n",
                               sep = "")
        }
        else {
          back.state1 <- paste("First model is infeasible for primary weights",
                               cat.back, "\n", sep = "")
        }
        if (use.backup) {
          use.model.i <- 2
        }
      }
    }
    if (use.model.i == 2 & use.backup) {
      mod[i] <- "Second"
      ws <- get.w.sub(newdat = newdata, newdat1 = newdat1a,
                      int.time = int.time, samp = samp, use = use,
                      n = NROW(newdat), maxit = maxit, calfun = calfun,
                      bounds = bounds, epsilon = cal.epsilon, trim = trim,
                      qpmeth = qpmeth, scale.var = scale.var)
      if (ws$mse > cut.mse | is.na(ws$mse)) {
        if (i > 1) {
          back.state2 <- paste("Second model was infeasible for ",
                               rep.meth, " group ", g, ".  Used third model.",
                               "\n", sep = "")
          back.state3 <- paste("First and second model were infeasible for ",
                               rep.meth, " group ", g, ".  Used third model.",
                               "\n", sep = "")
        }
        else {
          back.state2 <- paste("Second model is infeasible for primary weights.  Will use third model.",
                               "\n", sep = "")
          back.state3 <- paste("First and second model are infeasible for primary weights.  Will use third model.",
                               "\n", sep = "")
        }
        use.model.i <- 3
      }
    }
    if (use.model.i == 3 & use.backup) {
      mod[i] <- "Third"
      ws <- get.w.sub(newdat = newdatb, newdat1 = newdat1b,
                      int.time = int.time, samp = samp, use = use,
                      n = NROW(newdat), maxit = maxit, calfun = calfun,
                      bounds = bounds, epsilon = cal.epsilon, trim = trim,
                      qpmeth = qpmeth, scale.var = scale.var)
    }
    if (back.state1 != "" & back.state2 != "") {
      back.state.final <- back.state3
    }
    else if (back.state1 != "" & back.state2 == "") {
      back.state.final <- back.state1
    }
    else if (back.state1 == "" & back.state2 != "") {
      back.state.final <- back.state2
    }
    else {
      back.state.final <- ""
    }
    back.state <- paste(back.state, back.state.final, sep = "")
    back.state1 <- back.state2 <- back.state3 <- ""
    mses <- get.mse(newdat, newdat1, samp, use, ws$wghts,
                    ws$wghts.init, ws$scale.by)
    msesa <- get.mse(newdata, newdat1a, samp, use, ws$wghts,
                     ws$wghts.init, ws$scale.by)
    msesb <- get.mse(newdatb, newdat1b, samp, use, ws$wghts,
                     ws$wghts.init, ws$scale.by)
    mse[1, i] <- mses$mse
    mse[2, i] <- mses$mse1
    mse[3, i] <- msesa$mse
    mse[4, i] <- msesa$mse1
    mse[5, i] <- msesb$mse
    mse[6, i] <- msesb$mse1
    wghts[, i] <- ws$wghts
    tmp <- proc.time() - tmp
    if (i == 1) {
      if (back.state != "")
        message(back.state, appendLF = FALSE)
      back.state <- ""
      use.model <- use.model.i
      message("Created main weights for synthetic control: Time = ",
              round(tmp[3], 2), "\n\n", sep = "", appendLF = FALSE)
      message("Matching summary for main weights:\n", appendLF = FALSE)
      if (use.model.i == 1) {
        printstuff <- mses$printstuff
        message(paste0(capture.output(round(printstuff,
                                            4)), collapse = "\n"), appendLF = FALSE)
        message("\n", appendLF = FALSE)
      }
      else if (use.model.i == 2) {
        printstuff <- msesa$printstuff
        message(paste0(capture.output(round(printstuff,
                                            4)), collapse = "\n"), appendLF = FALSE)
        message("\n", appendLF = FALSE)
      }
      else if (use.model.i == 3) {
        printstuff <- msesb$printstuff
        message(paste0(capture.output(round(printstuff,
                                            4)), collapse = "/n"), appendLF = FALSE)
        message("\n", appendLF = FALSE)
      }
      if (jack > 0) {
        message("Calculating weights for jackknife replication groups...\n",
                appendLF = FALSE)
        tmp.jack <- proc.time()
      }
      else if (boot > 0) {
        message("Calculating weights for permutation groups...\n",
                appendLF = FALSE)
        tmp.boot <- proc.time()
      }
    }
    else if (i >= jack.lower & i <= jack.upper) {
      if (i == jack.lower) {
        message("Completed weights for jackknife group:\n",
                i - 1, sep = "", appendLF = FALSE)
      }
      else if ((i - 1)%%20 != 1 & i != jack.upper) {
        message(", ", i - 1, sep = "", appendLF = FALSE)
      }
      else if ((i - 1)%%20 == 1 & i != jack.upper) {
        message(", \n", i - 1, sep = "", appendLF = FALSE)
      }
      else if ((i - 1)%%20 != 1 & i == jack.upper) {
        message(", ", i - 1, "\n", sep = "", appendLF = FALSE)
      }
      else {
        message(", \n", i - 1, "\n", sep = "", appendLF = FALSE)
      }
      if (i == jack.upper) {
        if (i == jack.lower) {
        }
        if (back.state != "")
          message(back.state, appendLF = FALSE)
        back.state <- ""
        tmp.jack <- proc.time() - tmp.jack
        message("Completed weights for all jackknife replication groups: Time = ",
                round(tmp.jack[3], 2), "\n\n", sep = "", appendLF = FALSE)
        if (boot > 0) {
          message("Calculating weights for permutation groups...\n",
                  appendLF = FALSE)
          tmp.boot <- proc.time()
        }
      }
    }
    else if (i >= boot.lower & i <= boot.upper) {
      if (i == boot.lower) {
        message("Completed weights for permutation group:\n",
                i - jack - 1, sep = "", appendLF = FALSE)
      }
      else if ((i - jack - 1)%%20 != 1 & i != boot.upper) {
        message(", ", i - jack - 1, sep = "", appendLF = FALSE)
      }
      else if ((i - jack - 1)%%20 == 1 & i != boot.upper) {
        message(", \n", i - jack - 1, sep = "", appendLF = FALSE)
      }
      else if ((i - jack - 1)%%20 != 1 & i == boot.upper) {
        message(", ", i - jack - 1, "\n", sep = "", appendLF = FALSE)
      }
      else {
        message(", \n", i - jack - 1, "\n", sep = "",
                appendLF = FALSE)
      }
      if (i == boot.upper) {
        if (i == boot.lower) {
          message("\n", appendLF = FALSE)
        }
        tmp.boot <- proc.time() - tmp.boot
        message("Completed weights for all permutation groups: Time = ",
                round(tmp.boot[3], 2), "\n\n", sep = "", appendLF = FALSE)
      }
    }
  }
  out <- list(Weights = wghts, Intervention = Inter, MSE = mse,
              Model = mod, Summary = printstuff)
  return(out)
}


get.w.sub <- function(newdat = NULL, newdat1 = NULL, bigdat = NULL, dum = NULL, dum1 = NULL, covar.var = NULL, covar.var1 = NULL,
                      int.time, samp, use, n = NROW(newdat), maxit = 500, calfun = "raking", bounds = c(-Inf, Inf), epsilon = 1e-04, trim = NULL, qpmeth = "LowRankQP",
                      scale.var = "Intercept") {
  tmp <- proc.time()

  if (length(newdat) == 0) {
    newdat <- get.newdat(bigdat = bigdat, dum = dum, dum1 = dum1, covar.var = covar.var, covar.var1 = covar.var1, int.time = int.time)
    newdat1 <- newdat[[2]]
    newdat <- newdat[[1]]
  }

  mult <- sum(samp)/sum(use & samp)
  # mult <- 1

  if (length(newdat1) > 0) {
    newdat1 <- as.matrix(newdat1)
  }

  intdat <- newdat[samp & use, , drop = FALSE]
  scale.by <- sum(intdat[, scale.var])
  condat <- newdat[!samp & use, , drop = FALSE]

  alldat <- newdat[, , drop = FALSE]
  scale.by <- scale.by/sum(alldat[, scale.var])

  targets <- colSums(intdat)
  targets <- mult * targets
  # targets <- colSums(newdat[samp, , drop = FALSE])
  init <- rep(mult * NROW(intdat)/NROW(condat), NROW(condat))

  if (length(newdat1) > 0) {
    intdat1 <- newdat1[samp & use, , drop = FALSE]
    condat1 <- newdat1[!samp & use, , drop = FALSE]
    alldat1 <- newdat1[, , drop = FALSE]
    targets1 <- colSums(intdat1)
    targets1 <- mult * targets1
    # targets1 <- colSums(newdat1[samp, , drop = FALSE])
  }

  usevars <- colnames(condat)

  rem <- find.sing(t(as.matrix(condat)) %*% as.matrix(condat))
  keep <- !is.element(1:NCOL(condat), rem)

  form <- paste("~", paste(usevars[keep], collapse = "+", sep = ""), "-1", sep = "")
  form <- stats::formula(form)

  ws.init <- ws <- init

  if (NCOL(newdat) > 1) {
    tryCatch({
      caldesign <- survey::svydesign(ids = ~0, data = condat[, keep], weights = init)
      cali2 <- survey::calibrate(design = caldesign, maxit = maxit, formula = form, population = targets[keep], data = condat[,
                                                                                                                              keep], calfun = calfun, bounds = bounds, force = TRUE, epsilon = epsilon)
      ws.init <- ws <- stats::weights(cali2)
    }, error = function(e) {
      message("ERROR :", conditionMessage(e), "\n", appendLF=FALSE)
    })
  }

  if (length(newdat1) > 0) {
    tmp1 <- proc.time() - tmp
    tmp2 <- proc.time()

    condat2 <- condat1
    targets2 <- targets1

    ws <- my.qp(b.init = ws.init, X = t(condat2), a = targets2, Y = t(condat[, keep, drop = FALSE]), c = targets[keep], qpmeth = qpmeth)
    tmp2 <- proc.time() - tmp2
  }

  if (length(trim) > 0) {
    if (length(trim) == 1) {
      trim = c(trim, 1 - trim)
    }
    cali2 <- survey::trimWeights(cali2, upper = stats::quantile(ws, max(trim)), lower = stats::quantile(ws, min(trim)))
    ws <- stats::weights(cali2)
  }

  mse1 <- mean((colSums(ws * condat[, keep, drop = FALSE]) - targets[keep])^2)
  wghts1 <- rep(NA, n)
  wghts1[!samp & use] <- ws
  wghts1[samp & use] <- mult
  wghts.init1 <- rep(NA, n)
  wghts.init1[!samp & use] <- ws.init
  wghts.init1[samp & use] <- mult

  out <- list(wghts = wghts1, wghts.init = wghts.init1, mse = mse1, scale.by = scale.by)

  return(out)
}


get.newdat <- function(bigdat, dum = NULL, dum1 = NULL, covar.var = NULL, covar.var1 = NULL, int.time, time.names = NULL) {
  n <- dim(bigdat)[1]

  if (length(time.names) == 0) {
    time.names <- as.character(1:dim(bigdat)[3])
  }

  result.var <- names(dum)
  newdat <- list()

  if (length(result.var) > 0) {
    for (j in 1:length(result.var)) {
      dum.tmp <- dum[[j]]
      high <- int.time
      low <- int.time
      i <- 0
      while (low > 0 & i < length(dum.tmp)) {
        i <- i + 1
        low <- high - dum.tmp[i] + 1
        if (i == 1) {
          low.i <- low
          high.i <- high
        }

        if (i == 1) {
          newdat[[j]] <- integer(0)
        }
        newdat[[j]] <- cbind(newdat[[j]], rowSums(as.matrix(bigdat[, result.var[j], low:high])))
        if (low == high) {
          colnames(newdat[[j]])[NCOL(newdat[[j]])] <- paste(result.var[j], ".", time.names[low], sep = "")
        } else {
          colnames(newdat[[j]])[NCOL(newdat[[j]])] <- paste(result.var[j], ".", time.names[low], ".", time.names[high], sep = "")
        }
        high <- low - 1
      }
    }
  }

  covar.dat <- NULL
  if (length(covar.var) > 0) {
    covar.dat <- bigdat[, covar.var, 1]
    covar.dat <- as.matrix(covar.dat)
    colnames(covar.dat) <- covar.var
  }

  if (length(newdat) == 0) {
    newdat <- NULL
  } else {
    newdat <- data.frame(newdat, check.names = FALSE)
  }

  if (length(covar.dat) > 0 & length(newdat) > 0) {
    newdat <- data.frame(covar.dat, newdat, check.names = FALSE)
  } else if (length(covar.dat) > 0 & length(newdat) == 0) {
    newdat <- data.frame(covar.dat, check.names = FALSE)
  }

  if (length(newdat) > 0) {
    newdat <- data.frame(Intercept = rep(1, n), newdat, check.names = FALSE)
  } else {
    newdat <- data.frame(Intercept = rep(1, n))
  }

  newdat1 <- list()
  if (length(dum1) > 0) {
    result.var1 <- names(dum1)
    for (j in 1:length(result.var1)) {
      dum.tmp <- dum1[[j]]
      high <- int.time
      low <- int.time
      i <- 0
      while (low > 0 & i < length(dum.tmp)) {
        i <- i + 1
        low <- high - dum.tmp[i] + 1
        if (i == 1) {
          low.i <- low
          high.i <- high
        }

        if (i == 1) {
          newdat1[[j]] <- integer(0)
        }
        newdat1[[j]] <- cbind(newdat1[[j]], rowSums(as.matrix(bigdat[, result.var1[j], low:high])))
        if (low == high) {
          colnames(newdat1[[j]])[NCOL(newdat1[[j]])] <- paste(result.var1[j], ".", time.names[low], sep = "")
        } else {
          colnames(newdat1[[j]])[NCOL(newdat1[[j]])] <- paste(result.var1[j], ".", time.names[low], ".", time.names[high], sep = "")
        }
        high <- low - 1
      }
    }
  }

  covar.dat1 <- NULL
  if (length(covar.var1) > 0) {
    covar.dat1 <- bigdat[, covar.var1, 1]
    covar.dat1 <- as.matrix(covar.dat1)
    colnames(covar.dat1) <- covar.var1
  }

  if (length(newdat1) == 0) {
    newdat1 <- NULL
  } else {
    newdat1 <- data.frame(newdat1, check.names = FALSE)
  }

  if (length(covar.dat1) > 0 & length(newdat1) > 0) {
    newdat1 <- data.frame(covar.dat1, newdat1, check.names = FALSE)
  } else if (length(covar.dat1) > 0 & length(newdat1) == 0) {
    newdat1 <- data.frame(covar.dat1, check.names = FALSE)
  }

  if (length(newdat1) > 0) {
    newdat1 <- data.frame(newdat1, check.names = FALSE)
  } else {
    newdat1 <- NULL
  }

  return(list(newdat = newdat, newdat1 = newdat1))

}


get.stats <- function (bigdat, w, inter, mse, result.var = dimnames(bigdat)[[2]],
                       int.time, period = 1, plot.it = result.var, max.time = 80,
                       plot.first = 100, file = NULL, sep = TRUE, start.time = 25,
                       legend.spot = "bottomleft", omnibus.var = result.var, cut.mse = 1,
                       scale.var = "Intercept", twosided = FALSE, time.names = NULL)
{

  if (length(time.names) == 0) {
    time.names <- as.character(1:dim(bigdat)[3])
  }

  use.omnibus <- length(omnibus.var) > 0
  all.var <- union(result.var, plot.it)
  all.var <- union(all.var, omnibus.var)
  stat5 <- stat4 <- stat2 <- stat1 <- mu <- matrix(NA, NCOL(w),
                                                   length(result.var) + sum(use.omnibus))
  rownames(stat5) <- rownames(stat4) <- rownames(stat2) <- rownames(stat1) <- rownames(mu) <- colnames(w)
  if (use.omnibus) {
    colnames(stat5) <- colnames(stat4) <- colnames(stat2) <- colnames(stat1) <- colnames(mu) <- c(result.var,
                                                                                                  "Omnibus")
  }
  else {
    colnames(stat5) <- colnames(stat4) <- colnames(stat2) <- colnames(stat1) <- colnames(mu) <- c(result.var)
  }
  bigdat1 <- make.quarter3(bigdat, period = period, int.time = int.time)
  keep <- mse < cut.mse & !is.na(mse)
  keep[1] <- TRUE
  synth <- list()
  inter <- inter[, !grepl("Jack", colnames(w)), drop = FALSE]
  mse <- mse[!grepl("Jack", colnames(w))]
  keep <- keep[!grepl("Jack", colnames(w))]
  w <- w[, !grepl("Jack", colnames(w)), drop = FALSE]
  for (i in 1:NCOL(w)) {
    Inter <- inter[, i]
    use.con <- !is.na(Inter) & Inter == FALSE
    use.tre <- !is.na(Inter) & Inter == TRUE
    w.tmp <- w[use.con, i]
    condat1 <- bigdat1[use.con, all.var, , drop = FALSE]
    intdat1 <- bigdat1[use.tre, all.var, , drop = FALSE]
    test2 <- apply(w.tmp * condat1, c(2, 3), sum)
    test1 <- apply(intdat1, c(2, 3), sum)
    if (i == 1) {
      if (scale.var == "Intercept") {
        scale.by <- sum(use.tre)/sum(use.tre | use.con)
      }
      else {
        scale.by <- sum(bigdat1[use.tre, scale.var, 1])/sum(bigdat1[use.tre |
                                                                      use.con, scale.var, 1])
      }
      alldat1 <- bigdat1[use.tre | use.con, all.var, ,
                         drop = FALSE]
      test3 <- apply(alldat1, c(2, 3), sum)
    }
    if (i == 1) {
      xnams <- as.numeric(colnames(test1))
      tuse <- xnams <= max.time & xnams >= start.time
      use <- xnams <= int.time
      nuse <- xnams >= int.time
      if (max.time > int.time) {
        fuse <- !use & xnams <= max.time
      }
      else if (max.time == int.time) {
        fuse <- xnams >= int.time & xnams <= max.time
      }
      else {
        stop("max.time is less than int.time")
      }
      if (length(plot.it) > 0) {
        no.jack <- which(!grepl("Jack", colnames(w)))
        plot.first <- min(plot.first, NCOL(w) - 1)
        keep1 <- keep[no.jack]
        plotdat.d <- array(NA, c(length(plot.it), length(no.jack),
                                 length(xnams)))
        dimnames(plotdat.d) <- list(plot.it, colnames(w)[no.jack],
                                    time.names[xnams])
        plotdat.a <- plotdat.t <- plotdat.c <- array(NA,
                                                     c(length(plot.it), length(xnams)))
        dimnames(plotdat.a) <- dimnames(plotdat.t) <- dimnames(plotdat.c) <- list(plot.it,
                                                                                  time.names[xnams])
      }
      else {
        plotdat.t <- plotdat.c <- plotdat.a <- plotdat.d <- NULL
      }
    }
    if (use.omnibus) {
      use.cols <- 1:(NCOL(stat1) - 1)
    }
    else {
      use.cols <- 1:NCOL(stat1)
    }
    mu[i, use.cols] <- sum(fuse) * rowMeans(test1[result.var,
                                                  tuse, drop = FALSE])
    stat4[i, use.cols] <- rowSums(test1[result.var, fuse,
                                        drop = FALSE])
    stat5[i, use.cols] <- rowSums(test2[result.var, fuse,
                                        drop = FALSE])
    stat1[i, use.cols] <- stat4[i, use.cols, drop = FALSE] -
      stat5[i, use.cols, drop = FALSE]
    stat2[i, use.cols] <- stat1[i, use.cols, drop = FALSE]/rowSums(test2[result.var,
                                                                         fuse, drop = FALSE])
    if (length(plot.it) > 0) {
      if (i == 1) {
        plotdat.t[plot.it, ] <- test1[plot.it, ]
        plotdat.c[plot.it, ] <- test2[plot.it, ]
        plotdat.a[plot.it, ] <- test3[plot.it, ]
        i1 <- i
      }
      if (is.element(i, no.jack)) {
        plotdat.d[plot.it, i1, ] <- test1[plot.it, ] -
          test2[plot.it, ]
        i1 <- i1 + 1
      }
    }
    if (use.omnibus) {
      if (!twosided) {
        stat1[i, NCOL(stat1)] <- sum(stat1[i, omnibus.var])
        stat2[i, NCOL(stat2)] <- sum(stat1[i, omnibus.var])/sum(mu[i,
                                                                   omnibus.var])
      }
      else {
        stat1[i, NCOL(stat1)] <- sum((stat1[i, omnibus.var])^2)
        stat2[i, NCOL(stat2)] <- sum((stat1[i, omnibus.var]/mu[i,
                                                               omnibus.var])^2)
      }
    }
  }
  if (length(plot.it) > 0) {
    plotdat.d <- plotdat.d[, keep1, , drop = FALSE]
    mu <- mu[no.jack, , drop = FALSE]
    mu <- mu[keep1, , drop = FALSE]
    if (length(file) == 0) {
      graphics::par(mfrow = c(1, 2), ask = FALSE)
    }
    else {
      if (substr(file, nchar(file) - 3, nchar(file)) ==
          ".pdf") {
        file <- substr(file, 1, nchar(file) - 4)
      }
      if (!sep) {
        file <- paste(file, ".pdf", sep = "")
        grDevices::pdf(file = file, width = 8, height = 10.5)
        graphics::par(mfrow = c(3, 2), ask = FALSE)
      }
    }
    for (j in 1:length(plot.it)) {
      ylab1 <- main1 <- main <- plot.it[j]
      ylab2 <- "Treatment - Control"
      use.mu <- which(mu[, plot.it[j]] > 0)
      for (i in 1:dim(plotdat.d)[2]) {
        tmp <- plotdat.d[plot.it[j], i, ]
        if (i == 1) {
          if (sep & length(file) > 0) {
            grDevices::pdf(file = paste(file, "_", main1,
                                        "_TC.pdf", sep = ""), width = 5, height = 5)
          }
          tmp1 <- plotdat.t[plot.it[j], ]
          tmp2 <- plotdat.c[plot.it[j], ]
          tmp3 <- scale.by * plotdat.a[plot.it[j], ]
          lty <- c(1, 2, 4)
          col <- c(2, 1, 3)
          lwd <- c(2, 2, 2)
          ylim <- c(min(tmp1, tmp2), max(tmp1, tmp2))
          ylim <- c(min(ylim, tmp3), max(ylim, tmp3))
          ylim[2] <- 1.2 * ylim[2]
          xxnams1 <- as.numeric(as.character(time.names[xnams[tuse]]))
          xxnams2 <- as.numeric(as.character(time.names[xnams[tuse & use]]))
          xxnams3 <- as.numeric(as.character(time.names[xnams[tuse & nuse]]))
          iint.time <- as.numeric(as.character(time.names[int.time]))
          if (sum(is.na(xxnams1)) > 0) {
            xxnams1 <- xnams[tuse]
            xxnams2 <- xnams[tuse & use]
            xxnams3 <- xnams[tuse & nuse]
            iint.time <- int.time
          }
          xlim <- c(min(xxnams1), max(xxnams1))
          graphics::plot(xxnams1, tmp1[tuse], type = "l",
                         lty = lty[1], col = col[1], lwd = lwd[1],
                         xlim = xlim, xlab = "", ylab = ylab1, main = main,
                         ylim = ylim)
          graphics::abline(v = iint.time, lty = 2, col = 2)
          graphics::lines(xxnams1, tmp2[tuse], type = "l",
                          lty = lty[2], col = col[2], lwd = lwd[2])
          graphics::lines(xxnams1, tmp3[tuse], type = "l",
                          lty = lty[3], col = col[3], lwd = lwd[3])
          leg <- c("Treatment", "Synthetic Control",
                   "All blocks (scaled)")
          graphics::legend(legend.spot, legend = leg,
                           col = col, lty = lty, lwd = lwd, cex = 0.8,
                           bty = "n")
          if (sep & length(file) > 0) {
            grDevices::dev.off()
          }
          if (is.element(i, use.mu)) {
            if (sep & length(file) > 0) {
              grDevices::pdf(file = paste(file, "_",
                                          main1, "_Diff.pdf", sep = ""), width = 5,
                             height = 5)
            }
            ylim1 <- c(min(tmp), max(tmp))
            bigtmp <- plotdat.d[plot.it[j], use.mu, ]
            ylim2 <- 2 * c(stats::quantile(bigtmp, 0.05,
                                           na.rm = TRUE), stats::quantile(bigtmp,
                                                                          0.95, na.rm = TRUE))
            ylim <- c(min(ylim1[1], ylim2[1], na.rm = TRUE),
                      max(ylim1[2], ylim2[2], na.rm = TRUE))
            xlim <- c(min(xxnams1), max(xxnams1))
            graphics::plot(xxnams2, tmp[tuse &
                                          use], type = "l", ylim = ylim, xlim = xlim,
                           col = 2, lty = 2, main = main, xlab = "",
                           ylab = ylab2)
          }
          else {
            if (!sep & length(file) > 0) {
              graphics::plot(1, 1, main = plot.it[j])
            }
          }
        }
        else {
          if (i <= min((dim(plotdat.d)[2]), plot.first)) {
            if (is.element(i, use.mu) & is.element(1,
                                                   use.mu)) {
              graphics::lines(xxnams1, tmp[tuse],
                              col = "azure3")
            }
          }
        }
        if (i == dim(plotdat.d)[2] & is.element(1, use.mu)) {
          tmp <- plotdat.d[plot.it[j], 1, ]
          graphics::lines(xxnams2, tmp[tuse &
                                         use], lty = 1, col = 1, lwd = 2)
          graphics::lines(xxnams3, tmp[tuse &
                                         nuse], col = 2, lwd = 2)
          graphics::abline(v = iint.time, col = 2, lwd = 1,
                           lty = 2)
          graphics::abline(h = 0, lwd = 1, lty = 2)
          if (sep & length(file) > 0) {
            grDevices::dev.off()
          }
        }
      }
    }
    if (!sep & length(file) > 0) {
      grDevices::dev.off()
    }
  }
  stats <- list(stat1[keep, , drop = FALSE],
                stat2[keep, , drop = FALSE],
                stat4[keep, , drop = FALSE],
                stat5[keep, , drop = FALSE],
                list(Treatment = plotdat.t, Control = plotdat.c, All = plotdat.a, Difference = plotdat.d))
  return(stats)
}


get.stats1 <- function(bigdat, w, inter, mse, all.var, int.time, period = 1, max.time = 80, mfrow = c(1, 3), plot.first = 100, omnibus.var = NULL,
                       cut.mse = 1, G = 25, twosided = FALSE) {
  use.omnibus <- length(omnibus.var) > 0
  dof <- NA

  jack <- sum(grepl("Jack", colnames(w)))
  use.jack <- as.numeric(jack > 0)
  mse <- mse[!grepl("Jack", colnames(w))]
  w.jack <- w[, grepl("Jack", colnames(w))]
  w <- w[, !grepl("Jack", colnames(w)), drop = FALSE]
  inter.jack <- inter[, grepl("Jack", colnames(inter))]
  inter <- inter[, !grepl("Jack", colnames(inter)), drop = FALSE]

  if (length(G) == 0) {
    G <- min(table(inter[, 1]))
  }

  keep <- mse < cut.mse & !is.na(mse)
  keep[1] <- TRUE

  if (sum(!keep) > 0) {
    w <- w[, keep]
    inter <- inter[, keep]
    mse <- mse[keep]
  }

  boot <- sum(grepl("Perm", colnames(w)))
  boot.nams <- colnames(w)[grepl("Perm", colnames(w))]

  if (use.jack == 1) {
    all.nams <- c(colnames(w)[1], "Jackknife", boot.nams)
  } else {
    all.nams <- colnames(w)
  }

  tot.col <- 1 + use.jack + boot
  boot.lower <- 2 + use.jack
  boot.upper <- 1 + use.jack + boot

  delta.out <- stat2 <- stat1 <- matrix(NA, tot.col, length(all.var) + sum(use.omnibus))
  rownames(delta.out) <- rownames(stat2) <- rownames(stat1) <- all.nams
  if (use.omnibus) {
    colnames(delta.out) <- colnames(stat2) <- colnames(stat1) <- c(all.var, "Omnibus")
  } else {
    colnames(delta.out) <- colnames(stat2) <- colnames(stat1) <- c(all.var)
  }

  for (i in 1:tot.col) {
    G.tmp <- G
    if (i == 1) {
      i.tmp <- i
    } else {
      i.tmp <- i - use.jack
    }
    is.jack <- use.jack == 1 & i == 2
    tmp <- proc.time()
    Inter <- inter[, i.tmp]
    use <- !is.na(Inter)
    w.tmp <- w[use, i.tmp]
    is.tre <- !is.na(Inter) & Inter == TRUE
    is.tre <- is.tre[use]

    if (i == 1) {
      test <- make.quarter2(bigdat[use, , , drop = FALSE], tre = is.tre, w = w.tmp, period = period, int.time = int.time)
      if (max.time > int.time) {
        use.test <- test[, 1] > int.time & test[, 1] <= max.time
      } else if (max.time == int.time) {
        use.test <- test[, 1] >= int.time & test[, 1] <= max.time
      } else {
        stop("max.time is less than int.time")
      }
      time <- test[use.test, 1, drop = FALSE]
      time.tmp <- length(table(time))
      time.tmp1 <- time.tmp
      treat <- test[use.test, 2, drop = FALSE]
      w.tmp <- test[use.test, 3, drop = FALSE]
      treat1 <- treat
      w.tmp1 <- w.tmp
      test <- test[use.test, -(1:3), drop = FALSE]
    } else if (is.jack) {
      w.tmp <- w.tmp1
      treat <- treat1
      remain <- int.time%%period
      newcol <- (dim(bigdat)[3] - remain)%/%period
      time.tmp <- time.tmp1
      w.jack.tmp <- do.call("rbind", rep(list(w.jack), newcol))
      w.jack.tmp <- w.jack.tmp[use.test, , drop = FALSE]
      w.jack.tmp[is.na(w.jack.tmp)] <- 0
      G.tmp <- NCOL(w.jack.tmp)
    } else {
      remain <- int.time%%period
      newcol <- (dim(bigdat)[3] - remain)%/%period
      time.tmp <- time.tmp1
      w.tmp <- rep(w.tmp, newcol)
      w.tmp <- w.tmp[use.test]
      treat <- rep(as.numeric(is.tre), newcol)
      treat <- treat[use.test]
    }

    for (j in 1:length(all.var)) {
      test.tmp <- data.frame(y = test[, all.var[j]], treat = as.numeric(treat), time = factor(time))
      time.tmp <- 1
      if (!is.jack) {
        design <- survey::svydesign(ids = ~0, data = test.tmp, weights = w.tmp)
      } else {
        sup.out <- suppressWarnings({
          design <- survey::svrepdesign(data = test.tmp, repweights = w.jack.tmp, weights = w.tmp, combined.weights = TRUE,
                                        type = "JK1", mse = TRUE)
        })
      }
      usevars <- colnames(test.tmp)
      if (length(levels(factor(time))) == 1) {
        usevars <- setdiff(usevars, "time")
      }
      for (k in 2:length(usevars)) {
        if (k == 2) {
          form = paste("y~", usevars[k], sep = "")
        } else {
          form <- paste(form, "+", usevars[k], sep = "")
        }
      }
      if (!length(levels(factor(time))) == 1) {
        form <- paste(form, "-1", sep = "")
      } else {
        time.tmp <- 2
      }
      if (j == 1) {
        form1 <- substr(form, 2, 10000)
      }
      form <- stats::formula(form)
      sup.out <- suppressWarnings({
        mod <- survey::svyglm(form, design, family = stats::gaussian())
      })
      coefs <- summary(mod)$coefficients
      if(i == 1) {
        if(j == 1) {
          out.coefs <- matrix(NA, length(all.var), NCOL(coefs))
          rownames(out.coefs) <- all.var
          colnames(out.coefs) <- colnames(coefs)
        }
        out.coefs[j, ] <- coefs[time.tmp, ]
      }
      if (sum(coefs[, "Std. Error"]) == Inf | sum(coefs[, "Std. Error"]) == -Inf) {
        is.inf <- TRUE
      } else {
        is.inf <- FALSE
      }
      # stat1[i, j] <- coefs[time.tmp, 't value']
      stat1[i, j] <- coefs[time.tmp, "Estimate"]/coefs[time.tmp, "Std. Error"]
      yC <- mean(coefs[-time.tmp, "Estimate"])
      yT <- yC + coefs[time.tmp, "Estimate"]
      if (yT/yC > 0) {
        stat2[i, j] <- log(yT/yC)
      } else {
        stat2[i, j] <- NA
      }
      delta.out[i, j] <- my.delta(mu = coefs[, "Estimate"], Sigma = stats::vcov(mod))
    }
    if (use.omnibus) {
      if (i == 1) {
        dum.tmp <- stat1[i, omnibus.var]
        keep <- !is.na(dum.tmp)
        if (sum(!keep) > 0) {
          message("\nThe following variables yield survey statistics with value NA. \nThese will be removed from the omnibus statistic: \n", appendLF=FALSE)
          message(paste0(capture.output(omnibus.var[!keep]), collapse='\n'))
          message("\n", appendLF=FALSE)
        }
        omnibus.var <- omnibus.var[keep]
      }
      test.tmp <- data.frame(test[, omnibus.var, drop = FALSE], treat = as.numeric(treat), time = factor(time))
      form2 <- paste(omnibus.var, sep = "", collapse = ",")
      form2 <- paste("cbind(", form2, ")", sep = "")
      form1 <- paste(form2, form1, sep = "")
      form1 <- stats::formula(form1)
      if (i == 1) {
        reps <- assign.groups(as.numeric(treat), G = G.tmp)
      }
      coefs <- matrix(NA, length(omnibus.var), G.tmp)
      rownames(coefs) <- omnibus.var
      colnames(coefs) <- 1:G.tmp
      Sigma <- matrix(0, length(omnibus.var), length(omnibus.var))
      rownames(Sigma) <- colnames(Sigma) <- omnibus.var
      if (!is.jack) {
        for (g in 1:G.tmp) {
          allmod <- stats::lm(form1, data = test.tmp, weights = w.tmp, subset = which(reps != g))
          coefs[, g] <- as.matrix(stats::coef(allmod))["treat", ]
        }
      } else {
        for (g in 1:G.tmp) {
          test.tmp.tmp <- test.tmp[!is.na(inter.jack[, g]), ]
          w.jack.tmp.tmp <- w.jack.tmp[!is.na(inter.jack[, g]), g]
          allmod <- stats::lm(form1, data = test.tmp.tmp, weights = w.jack.tmp.tmp)
          coefs[, g] <- as.matrix(stats::coef(allmod))["treat", ]
        }
      }
      thetas <- rowMeans(coefs)
      dof <- length(thetas)
      coefs <- coefs - thetas
      thetas <- as.matrix(thetas)
      for (g in 1:G.tmp) {
        Sigma <- Sigma + coefs[, g, drop = FALSE] %*% t(coefs[, g, drop = FALSE])
      }
      Sigma <- (G.tmp - 1)/G.tmp * Sigma
      if (i == 1) {
        rm.var <- find.sing(Sigma)
        keep.var <- 1:NCOL(Sigma)
        if (length(rm.var) > 0) {
          keep.var <- keep.var[-rm.var]
        }
      }
      thetas <- thetas[keep.var, , drop = FALSE]
      Sigma <- Sigma[keep.var, keep.var, drop = FALSE]
      if (!twosided) {
        a <- as.matrix(diag(Sigma)^-0.5)
        stat1[i, NCOL(stat1)] <- t(a) %*% thetas/sqrt(t(a) %*% Sigma %*% a)
      } else {
        stat1[i, NCOL(stat1)] <- t(thetas) %*% solve(Sigma) %*% thetas
      }
    }
    tmp <- proc.time() - tmp
    if (i == 1) {
      message("Completed survey statistics for main weights: Time = ", round(tmp[3], 2), "\n",
              sep = "", appendLF=FALSE)
      if (jack == 0 & boot > 0) {
        message("Calculating survey statistics for permutation groups...\n", appendLF=FALSE)
        tmp.boot <- proc.time()
      }
      if (is.inf) {
        message("WARNING: Infinite standard errors yielded by main weights.\n", appendLF=FALSE)
      }
    } else if (is.jack) {
      message("Completed survey statistics for jackknife: Time = ", round(tmp[3], 2), "\n",
              sep = "", appendLF=FALSE)
      if (boot > 0) {
        message("Calculating survey statistics for permutation groups...\n", appendLF=FALSE)
        tmp.boot <- proc.time()
      }
      if (is.inf) {
        message("WARNING: Infinite standard errors yielded by jackknife weights.\n", appendLF=FALSE)
      }
    } else {
      if (i == boot.lower) {
        message("Completed survey statistics for permutation group:\n", i - use.jack - 1,
                sep = "", appendLF=FALSE)
      } else if ((i - use.jack - 1)%%20 != 1 & i != boot.upper) {
        message(", ", i - use.jack - 1, sep = "", appendLF=FALSE)
      } else if ((i - use.jack - 1)%%20 == 1 & i != boot.upper) {
        message(", \n", i - use.jack - 1, sep = "", appendLF=FALSE)
      } else if ((i - use.jack - 1)%%20 != 1 & i == boot.upper) {
        message(", ", i - use.jack - 1, "\n", sep = "", appendLF=FALSE)
      } else {
        message(", \n", i - use.jack - 1, "\n", sep = "", appendLF=FALSE)
      }
      if (i == boot.upper) {
        if (i == boot.lower) {
          message("\n", appendLF=FALSE)
        }
        tmp.boot <- proc.time() - tmp.boot
        message("Completed survey statistics for permutation groups: Time = ", round(tmp.boot[3], 2), "\n",
                sep = "", appendLF=FALSE)
      }
      if (is.inf) {
        message("WARNING: Infinite standard errors yielded by permutation group ", i, ".\n",
                sep = "", appendLF=FALSE)
      }
    }
  }
  stat2[stat2 == -Inf | stat2 == Inf] <- NA
  return(list(stat1, stat2, delta.out, dof, out.coefs))
}


my.qp <- function(b.init, X, Y, a, c, M = 10000, qpmeth = "LowRankQP", maxit = 1000) {
  q <- NCOL(X)
  n <- NROW(Y)

  a1 <- iso.nam(a, sep = ".")
  c1 <- iso.nam(c, sep = ".")

  a1[a1 == 0] <- 1
  c1[c1 == 0] <- 1

  a1 <- sqrt(a1)
  c1 <- sqrt(c1)

  X <- X/a1
  Y <- Y/c1
  a <- a/a1
  c <- c/c1

  if (qpmeth == "nleqslv") {
    requireNamespace("nleqslv", quietly = TRUE)
    XtX <- t(X) %*% X

    M <- 1/M
    jacobian <- function(all) {
      b <- all[1:q]
      lambda <- all[(q + 1):(q + n)]
      lowerright <- matrix(0, n, n)
      upperleft <- diag(1/b) + XtX
      upperright <- t(Y)
      lowerleft <- Y
      out <- rbind(cbind(upperleft, upperright), cbind(lowerleft, lowerright))
      return((out))
    }

    f1 <- function(all) {
      b <- all[1:q]
      lambda <- all[(q + 1):(q + n)]
      out1 <- (M * log(b)) + t(X) %*% (X %*% b - a) - t(Y) %*% lambda
      out2 <- Y %*% b - c
      return(c(out1, out2))
    }
    b.init[b.init < 1e-06] <- 1e-06
    lambda.init <- solve(Y %*% t(Y)) %*% Y %*% ((M * log(b.init)) + t(X) %*% (X %*% b.init - a))
    all.init <- c(b.init, lambda.init)

    all.root <- nleqslv::nleqslv(all.init, fn = f1, jac = NULL, method = c("Broyden"), global = c("gline"), xscalm = c("fixed"),
                                 control = list(maxit = maxit, xtol = 1e-11, ftol = 1e-11, btol = 1e-06, cndtol = 1e-13))

    message("Number of Jacobian evaluations = ", all.root$njcnt, ". \n", sep = "", appendLF=FALSE)
    message("Number of function evaluations = ", all.root$nfcnt, ". \n", sep = "", appendLF=FALSE)
    message("Number of iterations = ", all.root$iter, ". \n", sep = "", appendLF=FALSE)
    all.root <- all.root$x

    b <- all.root[1:q]
    b[b < 0] <- 0
  } else if (qpmeth == "LowRankQP") {
    requireNamespace("LowRankQP", quietly = TRUE)

    rem <- find.sing(Y %*% t(Y))
    leave <- setdiff(1:NROW(Y), rem)
    sup.out <- utils::capture.output(all.root <- LowRankQP::LowRankQP(Vmat = t(X), dvec = -t(X) %*% a, Amat = Y[leave, , drop = FALSE],
                                                                      bvec = c[leave], uvec = rep(M, q), method = "SMW", verbose = FALSE, niter = maxit))
    b <- all.root$alpha
  } else if (qpmeth == "ipop") {
    requireNamespace("kernlab", quietly = TRUE)

    all.root <- kernlab::ipop(c = -t(X) %*% a, H = t(X) %*% X, b = c, A = Y, l = rep(0, q), r = rep(0, n), u = rep(M, q))
    b <- all.root$alpha
  } else {
    stop("qpmeth not recognized.")
  }

  X <- X * a1
  Y <- Y * c1
  a <- a * a1
  c <- c * c1

  return(b)
}


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
    if (sum(here) > 1) {
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
      if (sum(diff.out) == 0) {
        nv.names <- c(nv.names, dimnames(out)[[2]][i])
      }
    }
  } else {
    nv.names <- dimnames(out)[[2]]
  }

  v.names <- setdiff(dimnames(out)[[2]], nv.names)

  return(list(bigdat = out, Intervention = int, nv.names = nv.names, v.names = v.names))
}


make.quarter2 <- function(dat, tre, w, period = 1, int.time) {
  n <- dim(dat)[1]
  p <- dim(dat)[2]
  q <- dim(dat)[3]
  add.back <- min(as.numeric(dimnames(dat)[[3]]))
  remain <- int.time%%period
  newcol <- (q - remain)%/%period
  times <- period * (1:newcol) + remain
  times <- times + add.back - 1

  out <- matrix(NA, n * length(times), p + 3)
  colnames(out) <- c("Time", "Treatment", "w", dimnames(dat)[[2]])
  start <- remain + 1
  for (i in 1:length(times)) {
    stop <- start + period - 1
    tmp <- dat[, , start:stop, drop = FALSE]
    here <- ((i - 1) * n + 1):(i * n)
    out[here, 1] <- times[i]
    out[here, 2] <- as.numeric(tre)
    out[here, 3] <- w
    out[here, 4:NCOL(out)] <- apply(tmp, c(1, 2), sum)
    start <- stop + 1
  }

  return(out)
}


make.quarter3 <- function(dat, period = 1, int.time) {
  n <- dim(dat)[1]
  p <- dim(dat)[2]
  q <- dim(dat)[3]
  add.back <- min(as.numeric(dimnames(dat)[[3]]))

  remain <- int.time%%period
  newcol <- (q - remain)%/%period
  times <- period * (1:newcol) + remain
  times <- times + add.back - 1

  out <- array(NA, c(n, p, newcol))
  dimnames(out) <- list(dimnames(dat)[[1]], dimnames(dat)[[2]], times)
  start <- remain + 1
  for (i in 1:newcol) {
    stop <- start + period - 1
    out[, , i] <- apply(dat[, , start:stop, drop = FALSE], c(1, 2), sum)
    start <- stop + 1
  }
  return(out)
}


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


iso.nam <- function(a, sep = "-") {
  nams <- names(a)
  locs <- regexpr(sep, nams, fixed = TRUE)
  outs <- nams
  outs[locs > 1] <- substr(outs[locs > 1], 1, locs[locs > 1] - 1)

  out <- rep(NA, length(outs))
  for (i in names(table(outs))) {
    out[outs == i] <- mean(a[outs == i])
  }
  return(out)
}


out.results <- function(results, int.time, max.time = names(results), file = NULL) {
  use.xlsx <- length(max.time) > 1

  if (substr(file, nchar(file) - 3, nchar(file)) == ".csv") {
    file <- substr(file, 1, nchar(file) - 4)
    if (use.xlsx) {
      message("WARNING: Cannot return .csv file since length(max.time) > 1.  Returning .xlsx file instead.\n",
              sep = "", appendLF=FALSE)
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
    max.time <- as.numeric(max.time)
    wb <- xlsx::createWorkbook()
    cspValColumn <- xlsx::CellStyle(wb, dataFormat = xlsx::DataFormat("0.0000"))
    csOtherColumn <- xlsx::CellStyle(wb, dataFormat = xlsx::DataFormat("0.00"))
    csPercColumn <- xlsx::CellStyle(wb, dataFormat = xlsx::DataFormat("0.0%"))

    for (i in 1:length(max.time)) {
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

      nam <- paste("Max time = ", max.time[i], sep = "")
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


merge.dums <- function(dum, dum1) {
  nam1 <- names(dum)
  nam2 <- names(dum1)

  if (length(nam1) > 0) {
    nam <- union(nam1, nam2)

    dum.out <- list()
    dum.out[[length(nam)]] <- NA
    dum.new <- list()
    dum.new[[length(nam1)]] <- NA
    names(dum.new) <- nam1

    dum[[length(nam) + 1]] <- NA
    dum1[[length(nam) + 1]] <- NA

    max.l.dum <- 0

    for (i in 1:length(nam)) {
      here1 <- which(nam1 == nam[i])
      here2 <- which(nam2 == nam[i])

      if (length(here1) == 0 & length(here2) == 0) {
        dum.out[[i]] <- NULL
      } else if (length(here1) > 0 & length(here2) == 0) {
        max.l.dum <- max(max.l.dum, length(dum[[here1]]))
        dum.new[[here1]] <- sum(dum[[here1]])
        dum.out[[i]] <- dum[[here1]]
      } else if (length(here1) == 0 & length(here2) > 0) {
        dum.out[[i]] <- dum1[[here2]]
      } else {
        max.l.dum <- max(max.l.dum, length(dum[[here1]]))
        dum.new[[here1]] <- sum(dum[[here1]])

        dum1c <- cumsum(dum[[here1]])
        dum2c <- cumsum(dum1[[here2]])

        dumc <- union(dum1c, dum2c)
        dumc <- dumc[order(dumc, decreasing = FALSE)]

        if (length(dumc) > 1) {
          dumc <- c(dumc[1], dumc[2:length(dumc)] - dumc[1:(length(dumc) - 1)])
        }

        dum.out[[i]] <- dumc
      }
    }

    if (max.l.dum <= 1) {
      dum.new <- list()
    }

    names(dum.out) <- nam
  } else {
    dum.new <- dum
    dum.out <- dum1
  }
  return(list(dum = dum.new, dum1 = dum.out))
}


is.feasible <- function(bigdat, covar.var, dum, Int, int.val = 1, int.time, eps = 0.001) {
  n <- dim(bigdat)[1]

  newdat <- get.newdat(bigdat, dum = dum, covar.var = covar.var, int.time = int.time)[[1]]

  intdat <- newdat[Int == int.val, ]
  condat <- newdat[Int != int.val, ]
  targets <- colSums(intdat)
  is.sol <- check.feasible2(t(condat), targets, eps = eps)
  return(is.sol)
}


find.sing <- function(X) {
  Z <- rref(X)

  rem <- NULL
  cont <- TRUE
  while (cont) {
    Z.tmp <- diag(Z)
    if (sum(Z.tmp == 0) == 0) {
      cont <- FALSE
    }
    if (cont) {
      Z.tmp <- min(which(Z.tmp == 0))
      Z <- Z[-NROW(Z), -Z.tmp]
      rem <- c(rem, Z.tmp + length(rem))
    }
  }
  return(rem)
}


check.feasible2 <- function(A, b, eps = 1e-07, M = 10000, meth = "LowRankQP") {
  rem <- find.sing(A %*% t(A))
  leave <- setdiff(1:NROW(A), rem)
  A <- A[leave, ]
  b <- b[leave]
  A3 <- cbind(A, diag(NROW(A)), (-1) * diag(NROW(A)))
  a <- c(rep(0, NCOL(A)), rep(1, 2 * NROW(A)))

  if (meth == "LowRankQP") {
    requireNamespace("LowRankQP", quietly = TRUE)
    Vmat <- matrix(0, NCOL(A3), 1)
    uvec <- rep(M, NCOL(A3))
    sup.out <- utils::capture.output(all.root <- LowRankQP::LowRankQP(Vmat = Vmat, dvec = a, Amat = A3, bvec = b, uvec = uvec,
                                                                      method = "SMW", verbose = FALSE, niter = maxit))
    sol <- all.root$alpha
  }
  if (meth == "simplex") {
    requireNamespace("boot", quietly = TRUE)
    simp <- boot::simplex(a = a, A3 = A3, b3 = b)
    sol <- simp$soln
  }
  # w <- sol[1:NCOL(A)]
  sol <- sol[-(1:NCOL(A))]
  sol <- sqrt(mean(sol * sol))
  if (sol >= eps) {
    out <- FALSE
  } else {
    out <- TRUE
  }
  return(out)
}


assign.groups <- function(strata = NULL, n = length(strata), G = min(table(strata))) {

  states <- names(table(strata))

  Gs <- sample(1:G, G)
  Gs1 <- Gs

  rep.G <- rep(NA, n)

  for (i in 1:length(states)) {
    here <- which(strata == states[i])

    n <- length(here)
    samp <- sample(1:n, n)

    J <- floor(n/G)

    rep.G1 <- rep(NA, n)

    for (g in 1:G) {
      here1 <- samp[1:J]
      samp <- samp[-(1:J)]
      rep.G1[here1] <- g
    }

    if (length(is.na(rep.G1)) > 0) {
      here2 <- which(is.na(rep.G1))
      if (length(here2) > length(Gs1)) {
        Gs1 <- c(Gs1, Gs)
      }
      rep.G1[here2] <- Gs1[1:length(here2)]
      Gs1 <- Gs1[-(1:length(here2))]
    }

    rep.G[here] <- rep.G1
  }

  return(rep.G)
}


get.mse <- function(newdat, newdat1 = NULL, samp, use, ws, ws.init, scale.by) {

  mult <- sum(samp)/sum(use & samp)

  intdat <- newdat[samp & use, , drop = FALSE]
  condat <- newdat[!samp & use, , drop = FALSE]
  alldat <- newdat[, , drop = FALSE]
  targets <- colSums(intdat)
  targets <- mult * targets
  ws <- ws[!samp & use]
  mse <- mean((colSums(ws * condat) - targets)^2)

  mse1 <- NA
  if (length(newdat1) != 0) {
    intdat1 <- newdat1[samp & use, , drop = FALSE]
    condat1 <- newdat1[!samp & use, , drop = FALSE]
    alldat1 <- newdat1[, , drop = FALSE]
    targets1 <- colSums(intdat1)
    targets1 <- mult * targets1
    ws.init <- ws.init[!samp & use]
    mse1 <- mean((colSums(ws * condat1) - targets1)^2)
  }

  if (length(newdat1) == 0) {
    printstuff <- cbind(Targets = targets, Weighted.Control = colSums((ws) * (condat)), All.scaled = (scale.by) * colSums(alldat))
  } else {
    # printstuff <- cbind(Targets = c(targets, targets1), Initial.Weighted.Control = colSums(ws.init * cbind(condat, condat1)),
    # Final.Weighted.Control = colSums(ws * cbind(condat, condat1)), All.scaled = scale.by * colSums(cbind(alldat, alldat1)))
    printstuff <- cbind(Targets = c(targets, targets1), Final.Weighted.Control = colSums(ws * cbind(condat, condat1)), All.scaled = scale.by *
                          colSums(cbind(alldat, alldat1)))
  }

  return(list(mse = mse, mse1 = mse1, printstuff = printstuff))
}


rref <- function(A, tol = sqrt(.Machine$double.eps), verbose = FALSE, fractions = FALSE) {
  ## A: coefficient matrix tol: tolerance for checking for 0 pivot verbose: if TRUE, print intermediate steps fractions: try to
  ## express nonintegers as rational numbers Written by John Fox Modified by Geoffrey Brent 2014-12-17 to fix a bug
  if (fractions) {
    mass <- requireNamespace("MASS", quietly = TRUE)
    if (!mass)
      stop("fractions=TRUE needs MASS package")
  }
  if ((!is.matrix(A)) || (!is.numeric(A)))
    stop("argument must be a numeric matrix")
  n <- nrow(A)
  m <- ncol(A)
  x.position <- 1
  y.position <- 1
  # change loop:
  while ((x.position <= m) & (y.position <= n)) {
    col <- A[, x.position]
    col[1:n < y.position] <- 0
    # find maximum pivot in current column at or below current row
    which <- which.max(abs(col))
    pivot <- col[which]
    if (abs(pivot) <= tol)
      x.position <- x.position + 1  # check for 0 pivot
    else {
      if (which > y.position)
      {
        A[c(y.position, which), ] <- A[c(which, y.position), ]
      }  # exchange rows
      A[y.position, ] <- A[y.position, ]/pivot  # pivot
      row <- A[y.position, ]
      A <- A - outer(A[, x.position], row)  # sweep
      A[y.position, ] <- row  # restore current row
      if (verbose)
        if (fractions) {
          message(paste0(capture.output(fractions(A)), collapse="\n"), appendLF=FALSE)
        } else {
          message(paste0(capture.output(round(A, round(abs(log(tol, 10))))), collapse="\n"), appendLF=FALSE)
        }
      x.position <- x.position + 1
      y.position <- y.position + 1
    }
  }
  for (i in 1:n) if (max(abs(A[i, 1:m])) <= tol)
    A[c(i, n), ] <- A[c(n, i), ]  # 0 rows to bottom
  if (fractions)
    fractions(A) else round(A, round(abs(log(tol, 10))))
}


remove.vars <- function(vars, nams, objnam = "result.var") {
  if (is.list(vars)) {
    vars1 <- names(vars)
  } else {
    vars1 <- vars
  }

  rm.vars <- setdiff(vars1, nams)

  if (length(rm.vars) > 0) {
    rm.here <- which(is.element(vars1, rm.vars))
    vars <- vars[-rm.here]
    message("WARNING: The following variables will be removed from ", objnam,
            " since they are not in the dataset: \n",
            sep = "", appendLF=FALSE)
    message(paste(rm.vars, collapse = ", ", sep = ""), "\n\n",
            sep = "", appendLF=FALSE)
  }

  return(vars)
}

check.matchout <- function (match.out, int.time) {
  bad <- FALSE
  for (i in 1:length(match.out)) {
    tmp <- match.out[[i]]
    if (int.time < sum(tmp)) {
      bad <- TRUE
      cum.tmp <- cumsum(tmp)
      cum.tmp <- cum.tmp <= int.time
      match.out[[i]] <- c(tmp[cum.tmp])
      if (int.time - sum(tmp[cum.tmp]) > 0) {
        match.out[[i]] <- c(match.out[[i]], int.time - sum(tmp[cum.tmp]))
      }
    }
  }
  return(list(match.out, bad))
}


