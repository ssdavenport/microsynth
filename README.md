
# microsynth

The goal of microsynth is to implement the synthetic control method for
micro-level data as outlined in Robbins, Saunders, and Kilmer (2017). 
is designed for use in assessment of the effect of an intervention using
longitudinal data. However, it may also be used to calculate propensity
score-type weights in cross-sectional data.  is a generalization of 
(see Abadie and Gardeazabal (2003) and Abadie, Diamond, Hainmueller
(2010, 2011, 2014)) that is designed for data at a more granular level
(e.g., micro-level). For more details see the help vignette: .

## Installation

You can install microsynth from github with:

``` r
# install.packages("devtools")
devtools::install_github("ssdavenport/microsynth")
```

## Example

This is a basic example which shows you how to run a simple microsynth
model, using block-level panel data to evaluate a crime intervention:

``` r
library("microsynth")
# Use seattledmi, block-level panel data, to evaluate a crime intervention.

# Declare time-variant (outcome) and time-invariant variables for matching
cov.var <- c('TotalPop', 'BLACK', 'HISPANIC', 'Males_1521',
       'HOUSEHOLDS', 'FAMILYHOUS', 'FEMALE_HOU', 'RENTER_HOU', 'VACANT_HOU')


match.out <- c('i_felony', 'i_misdemea', 'i_drugs', 'any_crime')

set.seed(99199) # for reproducibility

# Perform matching and estimation, without permutations or jackknife
# runtime: < 1 min
sea1 <- microsynth(seattledmi,
                  idvar="ID", timevar="time", intvar="Intervention",
                  start.pre=1, end.pre=12, end.post=16,
                  match.out=match.out, match.covar=cov.var,
                  result.var=match.out, omnibus.var=match.out,
                  test="lower")
#> Calculating weights...
#> Created main weights for synthetic control: Time = 1.68
#> Matching summary for main weights:
#>               Targets Weighted.Control All.scaled
#> Intercept          39          39.0002    39.0000
#> TotalPop         2994        2994.0519  2384.7477
#> BLACK             173         173.0010   190.5224
#> HISPANIC          149         149.0026   159.2682
#> Males_1521         49          49.0000    97.3746
#> HOUSEHOLDS       1968        1968.0340  1113.5588
#> FAMILYHOUS        519         519.0108   475.1876
#> FEMALE_HOU        101         101.0010    81.1549
#> RENTER_HOU       1868        1868.0203   581.9340
#> VACANT_HOU        160         160.0115    98.4222
#> i_felony.12        14          14.0000     4.9023
#> i_felony.11        11          11.0002     4.6313
#> i_felony.10         9           9.0000     3.0741
#> i_felony.9          5           5.0000     3.2642
#> i_felony.8         20          20.0000     4.4331
#> i_felony.7          8           8.0000     3.7617
#> i_felony.6         13          13.0000     3.0012
#> i_felony.5         20          20.0007     3.1549
#> i_felony.4         10          10.0000     4.0246
#> i_felony.3          7           7.0000     3.3693
#> i_felony.2         13          13.0002     3.2803
#> i_felony.1         12          12.0000     3.4381
#> i_misdemea.12      15          15.0002     4.2470
#> i_misdemea.11      12          12.0000     4.6070
#> i_misdemea.10      12          12.0000     4.0772
#> i_misdemea.9       14          14.0000     3.7414
#> i_misdemea.8       12          12.0000     3.9680
#> i_misdemea.7       20          20.0000     4.2551
#> i_misdemea.6       16          16.0005     3.5594
#> i_misdemea.5       24          24.0000     3.5635
#> i_misdemea.4       21          21.0002     4.3360
#> i_misdemea.3       21          21.0000     4.3846
#> i_misdemea.2       14          14.0000     3.5352
#> i_misdemea.1       16          16.0000     4.1540
#> i_drugs.12         13          13.0000     1.6543
#> i_drugs.11          8           8.0000     1.5128
#> i_drugs.10          3           3.0000     1.3227
#> i_drugs.9           4           4.0000     0.9788
#> i_drugs.8           4           4.0000     1.1123
#> i_drugs.7          10          10.0000     1.0516
#> i_drugs.6           4           4.0000     1.2377
#> i_drugs.5           2           2.0000     1.2296
#> i_drugs.4           1           1.0000     1.1245
#> i_drugs.3           5           5.0000     1.3550
#> i_drugs.2          12          12.0000     1.1366
#> i_drugs.1           8           8.0002     1.3591
#> any_crime.12      272         272.0012    65.3398
#> any_crime.11      227         227.0017    64.2396
#> any_crime.10      183         183.0010    55.6929
#> any_crime.9       176         176.0005    53.2377
#> any_crime.8       228         228.0005    55.8143
#> any_crime.7       246         246.0024    55.8062
#> any_crime.6       200         200.0010    52.8292
#> any_crime.5       270         270.0014    50.6531
#> any_crime.4       250         250.0010    57.2946
#> any_crime.3       236         236.0010    58.8681
#> any_crime.2       250         250.0012    51.5429
#> any_crime.1       242         242.0010    55.1145
#> Calculation of weights complete: Total time = 2.04
#> 
#> Calculating basic statistics for end.post = 16...
#> Completed calculation of basic statistics for end.post = 16.  Time = 3.85
#> 
#> Calculating survey statistics for end.post = 16...
#> Completed survey statistics for main weights: Time = 8.6
#> Completed calculation of survey statistics for end.post = 16.  Time = 8.6
#> 
#> microsynth complete: Overall time = 17.17

sea1 # View description of microsynth object and results
#>  microsynth object
#> 
#> Scope:
#>  Units:          Total: 9642 Treated: 39 Untreated: 9603
#>  Study Period(s):    Pre-period: 1 - 12  Post-period: 13 - 16
#>  Constraints:        Exact Match: 58     Minimized Distance: 0
#> Time-variant outcomes:
#>  Exact Match: i_felony, i_misdemea, i_drugs, any_crime (4)
#>  Minimized Distance: (0)
#> Time-invariant covariates:
#>  Exact Match: TotalPop, BLACK, HISPANIC, Males_1521, HOUSEHOLDS, FAMILYHOUS, FEMALE_HOU, RENTER_HOU, VACANT_HOU (9)
#>  Minimized Distance: (0)
#> 
#> Results:
#> end.post = 16
#>            Trt    Con Pct.Chng Linear.pVal Linear.Lower Linear.Upper
#> i_felony    46  68.22   -32.6%      0.0109       -50.3%        -8.4%
#> i_misdemea  45  71.80   -37.3%      0.0019       -52.8%       -16.7%
#> i_drugs     20  23.76   -15.8%      0.2559       -46.4%        32.1%
#> any_crime  788 986.44   -20.1%      0.0146       -32.9%        -4.9%
#> Omnibus     --     --       --      0.0011           --           --
plot_microsynth(sea1) # Plot results
```

![](README-example-1.png)<!-- -->![](README-example-2.png)<!-- -->![](README-example-3.png)<!-- -->![](README-example-4.png)<!-- -->
