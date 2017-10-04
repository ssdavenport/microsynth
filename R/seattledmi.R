#' Data for a crime intervention in Seattle, Washington
#'
#' The dataset contains information used to evaluate a Drug Market Intervention
#' (DMI) occurring in parts of Seattle, Washington in 2013. The data include
#' 2010 block-level Census data and counts of crime reported by the Seattle
#' Police, by crime type. Crime data are available for one year prior to the
#' intervention and two years after. DMIs are an intervention intended to
#' disrupt drug markets by targeting enforcement priorities at specific market
#' participants. The intervention was applied to 39 blocks in Seattle's
#' International District.
#'
#' @format A data frame with 154,272 rows and 22 columns, consisting of 9,642
#'   unique blocks with 16 (quarterly) observations each. It contains the
#'   following variables:
#'
#' \describe{
#'   \item{ID}{unique Census block ID}
#'
#'   \item{time}{time unit (in quarters)}
#'
#'   \item{Intervention}{time-variant binary indicator; all treated units
#'   receive 0 pre-intervention and 1 from the start of the intervention onward,
#'   while untreated cases receive 0s throughout}
#'
#'   \item{i_robbery}{number of robberies reported in that block-quarter
#'   (time-variant)}
#'
#'   \item{i_aggassau}{number of aggravated assaults reported}
#'
#'   \item{i_burglary}{number of burglaries reported}
#'
#'   \item{i_larceny}{number of larcenies reported}
#'
#'   \item{i_felony}{number of felony crimes reported}
#'
#'   \item{i_misdemea}{number of misdemeanor crimes reported}
#'
#'   \item{i_drugsale}{number of drug sales reported}
#'
#'   \item{i_drugposs}{number of drug possession incidents reported}
#'
#'   \item{i_drugs}{number of drug sale or possession incidents reported}
#'
#'   \item{any_crime}{number of all crimes reported}
#'
#'   \item{TotalPop}{number of residents}
#'
#'   \item{BLACK}{number of African American residents}
#'
#'   \item{HISPANIC}{number of Hispanic residents}
#'
#'   \item{Males_1521}{number of male residents aged 15-21}
#'
#'   \item{HOUSEHOLDS}{number of households}
#'
#'   \item{FAMILYHOUS}{number of family households}
#'
#'   \item{FEMALE_HOU}{number of female-headed households}
#'
#'   \item{RENTER_HOU}{number of households occupied by renters}
#'
#'   \item{VACANT_HOU}{number of vacant housing units}
#' }
#' @source Demographic data obtained from the 2010 Census, and administrative
#'   crime data from the Seattle Police Department.
"seattledmi"
