/* Se, PPV and FDR */

clear
cd "/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/06. Data/06. Results/"
use "/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/06. Data/06. Results/kraken2.dta"

drop Matches_found Total_PSVs Percentage Total_PSVs_filtered
gen method = 0 

rename pcv N_PSV 
replace enrichment = 1 if enrichment == .3
replace enrichment = 2 if enrichment == .5
replace enrichment = 3 if enrichment == .7
replace enrichment = 4 if enrichment == .9


expand 4

bysort niter: generate kraken2_status=_n
generate kraken2_weight=truepos if kraken2_status==1   /* True positive */
replace kraken2_weight=falsepos if kraken2_status==2  /* False positive weight */
replace kraken2_weight=falseneg if kraken2_status==3  /* False negative */
replace kraken2_weight=trueneg if kraken2_status==4  /* True negative weight */

save kraken2_temp.dta, replace

clear
use "/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/06. Data/06. Results/msweep.dta"

drop Filtered_matches Total_PSV Percentage_in_filtered Num_filtered_PSVs
gen method = 1

rename pcv N_PSV 
replace enrichment = 1 if enrichment == .3
replace enrichment = 2 if enrichment == .5
replace enrichment = 3 if enrichment == .7
replace enrichment = 4 if enrichment == .9


expand 4

bysort niter: generate msweep_status=_n
generate msweep_weight=truepos if msweep_status==1   /* True positive */
replace msweep_weight=falsepos if msweep_status==2  /* False positive weight */
replace msweep_weight=falseneg if msweep_status==3  /* False negative */
replace msweep_weight=trueneg if msweep_status==4  /* True negative weight */

append using kraken2_temp.dta
sort niter method




/** MODEL Se **/


* Define the outcome for sensitivity (1 for TP and FP, 0 for FN and TN)
*msweep 
generate sim_outcome=1 if msweep_status==1 | msweep_status==2 & method == 1 /* TP and FP */
replace sim_outcome=0 if msweep_status==3 | msweep_status==4 & method == 1 /* FN and TN */
*kraken2 
replace sim_outcome=1 if kraken2_status==1 | kraken2_status==2 & method == 0 /* TP and FP */
replace sim_outcome=0 if kraken2_status==3 | kraken2_status==4 & method == 0 /* FN and TN */

* Define true_status as 1 for actual positives (TP and FN)
*msweep 
generate true_status=1 if msweep_status==1 | msweep_status==3 & method == 1 /* TP and FN */
replace true_status=0 if msweep_status==2 | msweep_status==4 & method == 1 /*FP and TN */
*kraken2
replace true_status=1 if kraken2_status==1 | kraken2_status==3 & method == 0/* TP and FN */
replace true_status=0 if kraken2_status==2 | kraken2_status==4 & method == 0/*FP and TN */

gen weight = msweep_weight if method == 1
replace weight = kraken2_weight if method == 0

* Model
logistic sim_outcome i.enrichment##i.N_PSV##i.method if true_status == 1 [fweight=weight] /* not converging */
logistic sim_outcome i.enrichment##i.method i.N_PSV##i.method if true_status == 1 [fweight=weight]
testparm i.enrichment#i.method /* p = 0.92 */
testparm i.N_PSV#i.method /* p < 0.001 */
logistic sim_outcome i.enrichment i.N_PSV##i.method if true_status == 1 [fweight=weight]
testparm i.enrichment

* Margins table
margins i.N_PSV#i.method
marginsplot

* Pairwise comparisons
margins i.N_PSV#i.method, pwcompare(effects)
*margins i.N_PSV, at(method = (0 1)) pwcompare
*margins i.N_PSV if method == 0, pwcompare mcompare(sidak)
*margins i.N_PSV if method == 1, pwcompare mcompare(sidak)
*margins i.N_PSV if method ==1

contrast N_PSV@method


/** MODEL PPV **/
logistic true_status i.enrichment##i.N_PSV##i.method if sim_outcome == 1 [fweight=weight] 
testparm i.enrichment#i.N_PSV#i.method /* p = 0.91 */

logistic true_status i.enrichment##i.method i.N_PSV##i.method if sim_outcome == 1 [fweight=weight]
testparm i.enrichment#i.method /* p = 0.98 */
testparm i.N_PSV#i.method /* p < 0.001 */
logistic true_status i.enrichment i.N_PSV##i.method if sim_outcome == 1 [fweight=weight]
testparm i.enrichment

* Margins table
margins i.N_PSV#i.method
marginsplot 

* Pairwise comparisons
margins i.N_PSV#i.method, pwcompare(effects)


/** MODEL FDR **/
sum FDR
gen fdr_status = 1 - true_status
logistic fdr_status i.enrichment##i.N_PSV##i.method if sim_outcome == 1 [fweight=weight]
testparm i.enrichment#i.N_PSV#i.method /* p = 0.81 */

logistic fdr_status i.enrichment##i.method i.N_PSV##i.method if sim_outcome == 1 [fweight=weight]
testparm i.enrichment#i.method /* p = 0.98 */
testparm i.N_PSV#i.method /* p < 0.001 */
logistic fdr_status i.enrichment i.N_PSV##i.method if sim_outcome == 1 [fweight=weight]

* Margins table
margins i.N_PSV#i.method
marginsplot 

* Pairwise comparisons
margins i.N_PSV#i.method, pwcompare(effects)

bysort N_PSV: summarize FDR







/* Sp and NPV */

clear
cd "/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/06. Data/06. Results/"
use "/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/06. Data/06. Results/kraken2.dta"

drop Matches_found Total_PSVs Percentage Total_PSVs_filtered
gen method = 0 

rename pcv N_PSV 
replace enrichment = 1 if enrichment == .3
replace enrichment = 2 if enrichment == .5
replace enrichment = 3 if enrichment == .7
replace enrichment = 4 if enrichment == .9


expand 4

bysort niter: generate kraken2_status=_n
generate kraken2_weight=truepos if kraken2_status==1   /* True positive */
replace kraken2_weight=falsepos if kraken2_status==2  /* False positive weight */
replace kraken2_weight=falseneg if kraken2_status==3  /* False negative */
replace kraken2_weight=trueneg if kraken2_status==4  /* True negative weight */

save kraken2_temp.dta, replace

clear
use "/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/06. Data/06. Results/msweep.dta"

drop Filtered_matches Total_PSV Percentage_in_filtered Num_filtered_PSVs
gen method = 1

rename pcv N_PSV 
replace enrichment = 1 if enrichment == .3
replace enrichment = 2 if enrichment == .5
replace enrichment = 3 if enrichment == .7
replace enrichment = 4 if enrichment == .9


expand 4

bysort niter: generate msweep_status=_n
generate msweep_weight=truepos if msweep_status==1   /* True positive */
replace msweep_weight=falsepos if msweep_status==2  /* False positive weight */
replace msweep_weight=falseneg if msweep_status==3  /* False negative */
replace msweep_weight=trueneg if msweep_status==4  /* True negative weight */

append using kraken2_temp.dta
sort niter method


/** MODEL Sp **/

* Define the outcome for specificity (1 for TN and FP, 0 for TP and FN)
*msweep
generate sim_outcome=1 if msweep_status==3 | msweep_status==4 & method == 1 /* TP and FP */
replace sim_outcome=0 if msweep_status==1 | msweep_status==2 & method == 1 /* FN and TN */
*kraken2
replace sim_outcome=1 if kraken2_status==3 | kraken2_status==4 & method == 0 /* TP and FP */
replace sim_outcome=0 if kraken2_status==1 | kraken2_status==2 & method == 0 /* TP and FP */

* Define true_status as 1 for actual negatives (TP and FN)
*msweep 
generate true_status=1 if msweep_status==4 | msweep_status==2 & method == 1 /* True negatives (TN) or false positives (FP) */
replace true_status=0 if msweep_status==1 | msweep_status==3 & method == 1 /* True positives (TP) or false negatives (FN) */
*kraken2 
replace true_status=1 if kraken2_status==4 | kraken2_status==2 & method == 0 /* True negatives (TN) or false positives (FP) */
replace true_status=0 if kraken2_status==1 | kraken2_status==3 & method == 0 /* True positives (TP) or false negatives (FN) */

gen weight = msweep_weight if method == 1
replace weight = kraken2_weight if method == 0

* Model
logistic sim_outcome i.enrichment##i.N_PSV##i.method if true_status == 1 [fweight=weight]
testparm i.enrichment#i.N_PSV#i.method /*p=0.9458*/


logistic sim_outcome i.enrichment##i.method i.N_PSV##i.method if true_status == 1 [fweight=weight]
testparm i.enrichment#i.method /* p = 0.7354 */
testparm i.N_PSV#i.method /* p < 0.0001 */

logistic sim_outcome i.enrichment i.N_PSV##i.method if true_status == 1 [fweight=weight]
testparm i.enrichment /*0.7935*/
testparm i.N_PSV#i.method /* p < 0.0001 */

* Margins table
margins i.N_PSV#i.method
marginsplot

* Pairwise comparisons
margins i.N_PSV#i.method, pwcompare(effects)



/** MODEL NPV **/
logistic true_status i.enrichment##i.N_PSV##i.method if sim_outcome == 1 [fweight=weight] 
testparm i.enrichment#i.N_PSV#i.method /*p=0.0473, but several level of the triple interaction dropped due to collinearity*/

logistic true_status i.enrichment##i.method i.N_PSV##i.method if sim_outcome == 1 [fweight=weight]
testparm i.enrichment#i.method /* p = 0.8184 */
testparm i.N_PSV#i.method /* p < 0.0001 */

logistic true_status i.enrichment i.N_PSV##i.method if sim_outcome == 1 [fweight=weight] /*perfect prediction in two levels of the interaction N_PSV#method*/
testparm i.enrichment /*0.9939*/
testparm i.N_PSV#i.method /* p < 0.0001 */

* Margins table
margins i.N_PSV#i.method
marginsplot

* Pairwise comparisons
margins i.N_PSV#i.method, pwcompare(effects)

