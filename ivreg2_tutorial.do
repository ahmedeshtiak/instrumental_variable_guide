/*===========================================================================
  COMPLETE ivreg2 TUTORIAL
  Using Simulated Large Household Survey Data from an RCT
  
  Context: A poverty graduation programme (like BRAC TUP/DIUPG) where:
    - Instrument (Z)     = Randomized treatment assignment (intent-to-treat)
    - Endogenous var (D) = Actual programme participation (take-up)
    - Outcome (Y)        = Monthly household income (BDT)
  
  Topics covered:
    0.  Setup & data generation
    1.  OLS baseline (naive estimator, biased)
    2.  First-stage regression & instrument relevance
    3.  Basic ivreg2 (2SLS)
    4.  Weak instrument tests (F-stat, Kleibergen-Paap, Stock-Yogo)
    5.  Overidentification tests (Sargan / Hansen J)
    6.  Endogeneity test (Durbin-Wu-Hausman)
    7.  Heteroskedasticity-robust & cluster-robust SE
    8.  LIML (Limited Information Maximum Likelihood)
    9.  CUE (Continuously Updating Estimator)
   10.  GMM (Generalized Method of Moments)
   11.  Multiple endogenous variables
   12.  Partial R² and concentration parameter
   13.  ivreg2 with absorbing fixed effects (ivreg2h / absorb option)
   14.  Anderson-Rubin confidence sets (weak-instrument-robust inference)
   15.  Storing & presenting results with esttab
  
  Requirements:
    - Stata 14+
    - ssc install ivreg2
    - ssc install ranktest
    - ssc install estout
    - ssc install weakiv   (for Section 14)
  
  Author:  Tutorial prepared for BIGD Research Associates
  Date:    June 2026
===========================================================================*/


*===========================================================================
* SECTION 0: SETUP & SIMULATED DATA GENERATION
*===========================================================================

clear all
set more off
capture log close
log using "ivreg2_tutorial_log.txt", replace text

* ── Install required packages (comment out if already installed) ──────────
// ssc install ivreg2,   replace
// ssc install ranktest, replace
// ssc install estout,   replace
// ssc install weakiv,   replace

* ── Seed for full reproducibility ────────────────────────────────────────
set seed 20260606          // YYYYMMDD of today – always document your seed!

* ── Data dimensions ─────────────────────────────────────────────────────
* 3 000 households, 50 villages (clusters), 5 districts
local N        = 3000
local n_vill   = 50
local n_dist   = 5

set obs `N'

* ─────────────────────────────────────────────────────────────────────────
* 0.1  Geography
* ─────────────────────────────────────────────────────────────────────────
gen village_id  = ceil(_n / (`N'/`n_vill'))   // 60 HHs per village
gen district_id = ceil(village_id / 10)        // 10 villages per district

label var village_id  "Village identifier (1–50)"
label var district_id "District identifier (1–5)"

* ─────────────────────────────────────────────────────────────────────────
* 0.2  Household-level characteristics (pre-treatment baseline controls)
* ─────────────────────────────────────────────────────────────────────────
* Head's years of education (0–12, right-skewed for ultra-poor)
gen edu_head    = max(0, round(rnormal(3, 2.5)))
replace edu_head = min(edu_head, 12)

* Household size (2–10 members)
gen hh_size     = 2 + int(8 * rbeta(2, 4))

* Land owned (decimals; ~60 % landless, rest up to 100)
gen land_owned  = (runiform() > 0.60) * int(runiform()*100)

* Female-headed household
gen female_head = (runiform() < 0.35)

* Baseline (pre-programme) income in BDT/month – log-normal
gen ln_base_inc = rnormal(8.5, 0.6)           // mean ≈ BDT 4 900
gen base_income = exp(ln_base_inc)

* Unobserved ability/motivation (correlated with take-up – source of endogeneity)
gen u_ability   = rnormal(0, 1)

label var edu_head    "Head's education (years)"
label var hh_size     "Household size (# members)"
label var land_owned  "Land owned (decimals)"
label var female_head "Female-headed HH (1=yes)"
label var base_income "Baseline monthly income (BDT)"
label var u_ability   "Unobserved motivation (latent)"

* ─────────────────────────────────────────────────────────────────────────
* 0.3  INSTRUMENT: Randomized treatment assignment (Z)
*       Village-level randomization – 25 treatment villages
* ─────────────────────────────────────────────────────────────────────────
* Create a village-level random draw, then assign top-25 to treatment
gen vill_rand = .
forvalues v = 1/`n_vill' {
    quietly replace vill_rand = runiform() if village_id == `v' & _n == ///
        (select(_n, village_id == `v'))[1]
    // simpler: one draw per village then broadcast
}

* Cleaner approach: village-level random score, broadcast to HH level
bysort village_id: gen vill_draw = runiform() if _n == 1
bysort village_id: replace vill_draw = vill_draw[1]

* Top-half of villages assigned to treatment arm
xtile vill_rank = vill_draw, nq(2)
gen Z = (vill_rank == 2)            // Binary instrument: assignment to treatment
label var Z "Treatment assignment (village-level RCT, Z=1 if assigned)"

* ─────────────────────────────────────────────────────────────────────────
* 0.4  ENDOGENOUS VARIABLE: Actual programme participation (D)
*       Take-up is imperfect → D ≠ Z (classic non-compliance / partial compliance)
*       D is a function of Z AND unobserved motivation u_ability
* ─────────────────────────────────────────────────────────────────────────
* Latent propensity to participate
gen d_latent = -0.5 + 1.2*Z + 0.4*u_ability + 0.05*female_head ///
               - 0.03*edu_head + rnormal(0, 0.5)
gen D = (d_latent > 0)              // Binary: actual programme participation
label var D "Actual programme participation (1=yes)"

* Take-up rates for sanity check
tab Z D, row nofreq

* ─────────────────────────────────────────────────────────────────────────
* 0.5  SECOND INSTRUMENT (for overidentification tests later)
*       Household-level distance to nearest programme centre (km)
*       Inversely related to D, conditionally independent of outcome
* ─────────────────────────────────────────────────────────────────────────
gen dist_centre  = max(0.5, rnormal(5, 2) - 1.5*Z + rnormal(0, 0.3))
label var dist_centre "Distance to programme centre (km)"

* For the IV we use an indicator: more than 5 km away (reduces take-up)
gen far_from_centre = (dist_centre > 5)
label var far_from_centre "Far from centre >5 km (2nd instrument)"

* ─────────────────────────────────────────────────────────────────────────
* 0.6  OUTCOME: Monthly household income at endline (BDT)
*       True causal effect of participation: β_D = 2 000 BDT
*       Endogeneity: motivated HHs (high u_ability) also earn more → OLS biased UP
* ─────────────────────────────────────────────────────────────────────────
gen ln_income = 8.4 + 2000/exp(8.4)*D + 0.08*edu_head ///
                + 0.05*ln(base_income) - 0.02*hh_size  ///
                + 0.3*u_ability                       ///   ← endogeneity channel
                + rnormal(0, 0.4)

* Rescale to make income in BDT for intuitive interpretation
gen income = exp(ln_income) * 100   // endline monthly income, BDT

label var income    "Endline monthly income (BDT)"
label var ln_income "Log endline monthly income"

* Quick descriptives
sum income D Z edu_head hh_size land_owned female_head base_income

* Save dataset
save "rct_household_survey.dta", replace


*===========================================================================
* SECTION 1: OLS BASELINE (NAIVE – BIASED ESTIMATOR)
*===========================================================================
/*
  OLS regresses Y on D directly. Because D is correlated with u_ability
  (motivated people both participate AND earn more), OLS overestimates
  the causal effect of D. This is the classic upward ability bias.
*/

di _newline(2)
di as txt "========== SECTION 1: OLS (BIASED) =========="

reg income D edu_head hh_size land_owned female_head base_income i.district_id, ///
    vce(cluster village_id)

estimates store OLS
/*
  → Coefficient on D is upward biased because D is endogenous.
  → We need an instrument that shifts D without directly affecting income.
  → Z (randomized assignment) is our candidate instrument.
*/


*===========================================================================
* SECTION 2: FIRST-STAGE REGRESSION & INSTRUMENT RELEVANCE
*===========================================================================
/*
  The first stage regresses D on Z (+ controls).
  Key diagnostic: F-statistic on excluded instrument(s).
    - Rule of thumb: F > 10 (Stock & Yogo 2005 threshold for 10% maximal IV bias)
    - Kleibergen-Paap rk Wald F used when SE are clustered.
*/

di _newline(2)
di as txt "========== SECTION 2: FIRST STAGE =========="

reg D Z edu_head hh_size land_owned female_head base_income i.district_id, ///
    vce(cluster village_id)

estimates store first_stage

* Test joint significance of excluded instrument
test Z
di "First-stage F-stat on Z: " r(F)

* Partial R² of the instrument
gen D_hat = .
quietly reg D Z edu_head hh_size land_owned female_head base_income i.district_id
quietly predict D_hat_full
quietly reg D   edu_head hh_size land_owned female_head base_income i.district_id
quietly predict D_hat_noZ
quietly gen resid_full  = D - D_hat_full
quietly gen resid_noZ   = D - D_hat_noZ
quietly corr D_hat_full D_hat_noZ
drop D_hat D_hat_full D_hat_noZ resid_full resid_noZ

/*
  Interpretation:
  - A large, significant coefficient on Z confirms Z is a strong predictor of D.
  - High partial F ensures the instrument is not weak.
*/


*===========================================================================
* SECTION 3: BASIC ivreg2 — TWO-STAGE LEAST SQUARES (2SLS)
*===========================================================================
/*
  Syntax:
    ivreg2 Y X1 X2 ... (D = Z1 Z2 ...) [, options]
  
  The parentheses ( ) specify:
    Left  of =  →  endogenous regressors
    Right of =  →  excluded instruments (IVs)
  
  Variables OUTSIDE the parentheses are exogenous controls (included instruments).
*/

di _newline(2)
di as txt "========== SECTION 3: BASIC 2SLS =========="

ivreg2 income (D = Z) edu_head hh_size land_owned female_head base_income ///
       i.district_id, first

/*
  Output breakdown:
  ─────────────────────────────────────────────────────────────────────
  • Coefficient on D:  LATE (Local Average Treatment Effect)
      = causal effect for COMPLIERS (those induced to participate by Z)
  • first option:  also prints the first-stage regression
  • Underidentification test (Kleibergen-Paap rk LM statistic):
      H0: equation is underidentified → reject H0 for valid instrument
  • Weak identification test (Kleibergen-Paap rk Wald F statistic):
      Compare with Stock-Yogo critical values (printed in output)
  • Overidentification test: not applicable with exactly 1 instrument (just-identified)
  ─────────────────────────────────────────────────────────────────────
*/

estimates store iv2sls_basic


*===========================================================================
* SECTION 4: WEAK INSTRUMENT DIAGNOSTICS IN DETAIL
*===========================================================================
/*
  Weak instruments inflate IV standard errors and cause finite-sample bias
  toward OLS. Key statistics:
  
  (a) Stock-Yogo F-critical values (iid errors):
      Default size distortion thresholds at 5% nominal level.
  
  (b) Kleibergen-Paap rk Wald F (with robust/clustered SE):
      Recommended when using vce(cluster ...) or vce(robust).
  
  (c) Cragg-Donald F (iid errors, classical): legacy, avoid with clustered SE.
*/

di _newline(2)
di as txt "========== SECTION 4: WEAK INSTRUMENT TESTS =========="

* Run ivreg2 with cluster-robust SE and request all weak-ID statistics
ivreg2 income (D = Z) edu_head hh_size land_owned female_head base_income ///
       i.district_id,                ///
       robust                        ///   HC1-robust standard errors
       first                         ///   show first stage
       ffirst                        ///   first-stage F for each endogenous var
       savefirst                     //    save first-stage results

* Retrieve saved first-stage estimates
estimates restore _ivreg2_D_first_stage
di "Kleibergen-Paap F (robust): " e(widstat)

/*
  Reading the output:
  ─────────────────────────────────────────────────────────────────────
  Kleibergen-Paap rk Wald F:
    > 16.38  →  <5% size distortion (10% maximal relative bias)
    > 22.30  →  <5% size distortion (5% maximal relative bias)
  Stock-Yogo critical values are printed directly by ivreg2.
  
  If F < 10: instrument is WEAK → consider LIML or Anderson-Rubin CI
  ─────────────────────────────────────────────────────────────────────
*/


*===========================================================================
* SECTION 5: OVERIDENTIFICATION TEST (SARGAN / HANSEN J)
*===========================================================================
/*
  When we have MORE instruments than endogenous variables (overidentified model),
  we can test whether the instruments are exogenous (uncorrelated with ε).
  
  Sargan statistic: assumes homoskedastic errors
  Hansen J statistic: robust to heteroskedasticity (preferred)
  
  H0: all instruments are exogenous (jointly valid)
  Reject H0 → at least one instrument is invalid
  
  NOTE: This test CANNOT tell you WHICH instrument is invalid.
        With one instrument, the model is just-identified → no test possible.
        We need ≥ 2 instruments here.
*/

di _newline(2)
di as txt "========== SECTION 5: OVERIDENTIFICATION TEST =========="

* Now use BOTH instruments: Z and far_from_centre
ivreg2 income (D = Z far_from_centre) ///
       edu_head hh_size land_owned female_head base_income ///
       i.district_id, ///
       robust         ///
       first

/*
  Output now includes:
  • Hansen J statistic + p-value
    - p > 0.05: fail to reject H0 → both instruments appear exogenous (good)
    - p < 0.05: at least one instrument is suspect → re-examine exclusion restriction
  • C-statistic (difference-in-Hansen): tests exogeneity of subset of instruments
*/

* C-statistic: test exogeneity of far_from_centre while maintaining Z as clean
ivreg2 income (D = Z far_from_centre) ///
       edu_head hh_size land_owned female_head base_income ///
       i.district_id, ///
       robust         ///
       orthog(far_from_centre)   // C-stat: is far_from_centre exogenous?

estimates store iv2sls_overid


*===========================================================================
* SECTION 6: ENDOGENEITY TEST (DURBIN-WU-HAUSMAN)
*===========================================================================
/*
  Tests whether D is actually endogenous. If D were exogenous, OLS would
  be consistent AND more efficient than IV → no reason to use IV.
  
  H0: D is exogenous (OLS is consistent)
  
  DWH statistic: Hausman-type, Sargan-Hansen formulation in ivreg2
  
  If p > 0.05: fail to reject → OLS may be preferred (D is exogenous)
  If p < 0.05: reject H0    → D is endogenous → use IV
*/

di _newline(2)
di as txt "========== SECTION 6: ENDOGENEITY (DWH) TEST =========="

ivreg2 income (D = Z) edu_head hh_size land_owned female_head base_income ///
       i.district_id, ///
       robust         ///
       endog(D)       //  request Durbin-Wu-Hausman test for D

/*
  Interpretation:
  • Robust score χ² (Baum-Schaffer-Stillman) or Hausman F-stat
  • p < 0.05: D is endogenous → IV is necessary
  • p > 0.05: OLS is consistent; IV loses precision without gaining consistency
*/


*===========================================================================
* SECTION 7: CLUSTER-ROBUST STANDARD ERRORS
*===========================================================================
/*
  In household surveys, errors are correlated WITHIN clusters (villages).
  Ignoring clustering understates SEs, inflates t-stats, leads to over-rejection.
  
  Rule: cluster at the level of treatment assignment.
  Here: treatment is assigned at village level → cluster by village_id.
  
  ivreg2 supports:
    vce(robust)           HC1-robust (heteroskedasticity only)
    vce(cluster varname)  cluster-robust (heteroskedasticity + within-cluster correlation)
    vce(hc2) / vce(hc3)   finite-sample corrected HC
*/

di _newline(2)
di as txt "========== SECTION 7: CLUSTER-ROBUST SE =========="

* (a) Homoskedastic (default – for comparison)
ivreg2 income (D = Z) edu_head hh_size land_owned female_head base_income ///
       i.district_id
estimates store iv_homo

* (b) HC1-robust
ivreg2 income (D = Z) edu_head hh_size land_owned female_head base_income ///
       i.district_id, robust
estimates store iv_robust

* (c) Cluster-robust (recommended for RCT with village-level treatment)
ivreg2 income (D = Z) edu_head hh_size land_owned female_head base_income ///
       i.district_id, cluster(village_id)
estimates store iv_cluster

* Compare SEs side by side
esttab iv_homo iv_robust iv_cluster, ///
       keep(D) se star(* 0.10 ** 0.05 *** 0.01) ///
       mtitles("Homoskedastic" "HC-Robust" "Cluster-Robust") ///
       title("Table: Effect of Participation on Income — SE comparison")

/*
  Key insight: cluster-robust SEs are typically LARGER (more conservative).
  Under RCTs with few clusters, consider wild cluster bootstrap (boottest package).
*/


*===========================================================================
* SECTION 8: LIML — LIMITED INFORMATION MAXIMUM LIKELIHOOD
*===========================================================================
/*
  LIML is the maximum-likelihood equivalent of 2SLS for a single equation.
  
  Advantages over 2SLS:
  • Less finite-sample bias when instruments are weak
  • Median-unbiased under weak-instrument asymptotics
  
  Disadvantage:
  • Larger variance than 2SLS when instruments are strong
  
  When to use: F < 10 (weak instruments) → prefer LIML or Fuller's k-class
  
  Syntax: add liml option to ivreg2
*/

di _newline(2)
di as txt "========== SECTION 8: LIML =========="

ivreg2 income (D = Z far_from_centre) ///
       edu_head hh_size land_owned female_head base_income ///
       i.district_id, ///
       cluster(village_id) ///
       liml                //  LIML estimator

estimates store iv_liml

* Fuller's LIML (reduces finite-sample bias further, k parameter)
ivreg2 income (D = Z far_from_centre) ///
       edu_head hh_size land_owned female_head base_income ///
       i.district_id, ///
       cluster(village_id) ///
       fuller(4)           //  Fuller k-class, k=4 minimises MSE under normality

estimates store iv_fuller

/*
  In the output:
  • LIML eigenvalue (κ): close to 1 → 2SLS ≈ LIML (strong instruments)
                          κ >> 1   → instruments are weak, LIML preferred
*/


*===========================================================================
* SECTION 9: CUE — CONTINUOUSLY UPDATING ESTIMATOR
*===========================================================================
/*
  CUE (Hansen, Heaton & Yaron 1996) iterates jointly over β and the weight
  matrix, unlike GMM which uses a two-step weight matrix.
  
  Advantage:  better finite-sample properties than 2-step GMM
  Disadvantage: can have multiple local optima; slower convergence
  
  Syntax: add cue option
*/

di _newline(2)
di as txt "========== SECTION 9: CUE =========="

ivreg2 income (D = Z far_from_centre) ///
       edu_head hh_size land_owned female_head base_income ///
       i.district_id, ///
       cluster(village_id) ///
       cue                 //  CUE estimator

estimates store iv_cue


*===========================================================================
* SECTION 10: GMM — GENERALIZED METHOD OF MOMENTS
*===========================================================================
/*
  GMM exploits moment conditions: E[Z'ε] = 0
  
  With overidentified model (more IVs than endogenous vars), GMM is more
  efficient than 2SLS by optimally weighting the moment conditions.
  
  gmm2s = two-step efficient GMM (robust to heteroskedasticity)
  
  Syntax: add gmm2s option
  
  Note: with homoskedastic errors, 2SLS = GMM (optimal weighting is the same)
*/

di _newline(2)
di as txt "========== SECTION 10: GMM =========="

* (a) One-step GMM (uses identity weighting = same as 2SLS with robust)
ivreg2 income (D = Z far_from_centre) ///
       edu_head hh_size land_owned female_head base_income ///
       i.district_id, ///
       robust

estimates store iv_gmm1s

* (b) Two-step efficient GMM (optimal weighting matrix)
ivreg2 income (D = Z far_from_centre) ///
       edu_head hh_size land_owned female_head base_income ///
       i.district_id, ///
       robust gmm2s        //  efficient GMM with robust weight matrix

estimates store iv_gmm2s

* (c) Cluster-robust two-step GMM
ivreg2 income (D = Z far_from_centre) ///
       edu_head hh_size land_owned female_head base_income ///
       i.district_id, ///
       cluster(village_id) gmm2s

estimates store iv_gmm_cl

/*
  GMM output adds:
  • J-statistic (Hansen J test of overidentification)
  • Weight matrix used
  Two-step GMM is generally preferred over 2SLS in the overidentified case
  with heteroskedasticity.
*/


*===========================================================================
* SECTION 11: MULTIPLE ENDOGENOUS VARIABLES
*===========================================================================
/*
  Suppose both programme participation (D) AND asset transfer receipt (A)
  are endogenous. We need at least as many instruments as endogenous variables.
  
  Here we construct:
    D = programme participation     (endogenous)
    A = asset transfer value (BDT)  (endogenous – correlated with D and u_ability)
  
  Instruments:
    Z               = village-level assignment
    far_from_centre = distance-based instrument
*/

di _newline(2)
di as txt "========== SECTION 11: MULTIPLE ENDOGENOUS VARS =========="

* Generate second endogenous variable: asset value received
gen asset_value = D * (5000 + 2000*u_ability + rnormal(0, 500))
replace asset_value = 0 if D == 0
label var asset_value "Asset transfer value (BDT)"

* Regenerate income with both D and asset_value as causal inputs
gen income2 = exp(8.4 + 0.4*D + 0.0001*asset_value + 0.08*edu_head ///
                  + 0.05*ln(base_income) - 0.02*hh_size             ///
                  + 0.3*u_ability + rnormal(0, 0.4)) * 100

label var income2 "Endline income (2-endogenous model)"

* 2SLS with two endogenous variables and two instruments
ivreg2 income2 (D asset_value = Z far_from_centre) ///
       edu_head hh_size land_owned female_head base_income ///
       i.district_id, ///
       cluster(village_id) ///
       first ffirst        //  show both first stages

/*
  Interpretation:
  • Two first-stage regressions are shown (one for D, one for asset_value)
  • Kleibergen-Paap rk statistic now tests JOINT rank condition for both IVs
  • Model is exactly identified (2 IVs, 2 endogenous) → no overidentification test
  
  IMPORTANT: Each endogenous variable needs its own "relevant" instrument.
  Cross-equation partial F-stats (ffirst) tell you how well identified each is.
*/

estimates store iv_multi


*===========================================================================
* SECTION 12: PARTIAL R² AND CONCENTRATION PARAMETER
*===========================================================================
/*
  Partial R²: share of variation in D explained by Z AFTER partialling out X.
  
  Concentration parameter (μ²): measure of instrument strength in finite samples.
  Rule of thumb: μ² > 10 for reliable inference.
  
  ivreg2 reports these under the first stage with the 'first' option.
  We can also compute manually for pedagogical purposes.
*/

di _newline(2)
di as txt "========== SECTION 12: PARTIAL R² & CONCENTRATION PARAMETER =========="

* Obtain partial R² via ivreg2 with partialr2 option
ivreg2 income (D = Z) edu_head hh_size land_owned female_head base_income ///
       i.district_id, ///
       cluster(village_id) ///
       first partialr2

/*
  Partial R² from first stage output:
  → Proportion of endogenous variation in D explained by Z
  → Higher = stronger instrument
  
  Concentration parameter ≈ N × (partial R²) / (1 - partial R²)
  → Rule of thumb: > 10 for reliable 2SLS inference
*/

* Manual computation for illustration
quietly reg D edu_head hh_size land_owned female_head base_income i.district_id
predict D_resid, resid

quietly reg D Z edu_head hh_size land_owned female_head base_income i.district_id
predict D_hat_with_Z

gen Z_contribution = D_hat_with_Z - D + D_resid   // variation due to Z
quietly corr Z_contribution D
local partial_r2 = r(rho)^2
di as result "Manual Partial R² of Z: " %6.4f `partial_r2'

local conc_param = `N' * `partial_r2' / (1 - `partial_r2')
di as result "Concentration parameter μ²: " %8.2f `conc_param'

drop D_resid D_hat_with_Z Z_contribution


*===========================================================================
* SECTION 13: FIXED EFFECTS ABSORB (LARGE-FE MODELS)
*===========================================================================
/*
  In large datasets with many fixed effects (e.g., 500 villages), including
  them as dummies is slow. The absorb() option demeans by FE without creating
  dummy variables — equivalent to the within estimator.
  
  Syntax: absorb(fe_var)
  
  Note: absorb() requires ivreg2 ≥ 3.1 and ranktest.
        The absorbed FEs are NOT reported in output but are partialled out.
*/

di _newline(2)
di as txt "========== SECTION 13: ABSORBING FIXED EFFECTS =========="

* Without absorb (slow with many FEs):
ivreg2 income (D = Z) edu_head hh_size land_owned female_head base_income ///
       i.district_id, ///
       cluster(village_id)
estimates store iv_fedir

* With absorb (preferred for large FE):
ivreg2 income (D = Z) edu_head hh_size land_owned female_head base_income, ///
       cluster(village_id) ///
       absorb(village_id)   //  absorb village FEs (within-village identification)
estimates store iv_feabs

/*
  Notes on absorb():
  • Partials out ALL within-village variation
  • Treatment must vary WITHIN the absorbed unit (village) for identification
  • Our Z varies at village level → no within-village variation → 
    absorb(village_id) would absorb Z! Use district or subdistrict FEs instead.
  
  Let's demonstrate with district FEs absorbed:
*/
ivreg2 income (D = Z) edu_head hh_size land_owned female_head base_income, ///
       cluster(village_id) ///
       absorb(district_id)
estimates store iv_dist_abs

* Compare district dummies vs absorb
esttab iv_fedir iv_dist_abs, ///
       keep(D) se star(* 0.10 ** 0.05 *** 0.01) ///
       mtitles("i.district dummies" "absorb(district_id)") ///
       title("Table: District FE — dummy vs absorb approach (should match)")


*===========================================================================
* SECTION 14: WEAK-INSTRUMENT-ROBUST INFERENCE (ANDERSON-RUBIN)
*===========================================================================
/*
  Standard 2SLS t-tests are INVALID when instruments are weak.
  The Anderson-Rubin (1949) test and projection-based confidence sets
  remain valid regardless of instrument strength.
  
  Use the weakiv package (Finlay & Magnusson):
    ssc install weakiv
  
  Produces:
  • AR confidence interval (valid CI even with weak IVs)
  • K statistic (Kleibergen 2002 score test)
  • CLR statistic (conditional likelihood ratio)
*/

di _newline(2)
di as txt "========== SECTION 14: WEAK-INSTRUMENT-ROBUST CI (AR / CLR) =========="

* First run ivreg2 and then pass it to weakiv
ivreg2 income (D = Z) edu_head hh_size land_owned female_head base_income ///
       i.district_id, ///
       cluster(village_id)

* Attempt weakiv if installed (comment out if not installed)
capture {
    weakiv ivreg2 income (D = Z) edu_head hh_size land_owned female_head ///
           base_income i.district_id, ///
           cluster(village_id)
}
if _rc == 0 {
    di as result "weakiv ran successfully. Check AR and CLR confidence sets above."
}
else {
    di as txt "Note: weakiv not installed. Run: ssc install weakiv"
    di as txt "The Anderson-Rubin CI is robust to weak instruments."
    di as txt "Strongly recommended when Kleibergen-Paap F < 10."
}

/*
  Anderson-Rubin confidence set:
  • Often wider than 2SLS CI when instruments are weak
  • Has correct coverage regardless of instrument strength
  • May be unbounded (entire real line) if instrument is very weak
  
  Interpretation:
  • If AR CI ≈ 2SLS CI → instruments are strong, 2SLS is reliable
  • If AR CI >> 2SLS CI → instruments are weak, report AR CI instead
*/


*===========================================================================
* SECTION 15: STORING & PRESENTING RESULTS WITH esttab
*===========================================================================
/*
  Professional output table comparing OLS, 2SLS, LIML, and GMM.
*/

di _newline(2)
di as txt "========== SECTION 15: RESULTS TABLE =========="

* Re-run main specifications cleanly for table
quietly reg income D edu_head hh_size land_owned female_head base_income ///
    i.district_id, vce(cluster village_id)
estimates store T_ols

quietly ivreg2 income (D = Z) edu_head hh_size land_owned female_head ///
    base_income i.district_id, cluster(village_id)
estimates store T_2sls

quietly ivreg2 income (D = Z far_from_centre) edu_head hh_size land_owned ///
    female_head base_income i.district_id, cluster(village_id) liml
estimates store T_liml

quietly ivreg2 income (D = Z far_from_centre) edu_head hh_size land_owned ///
    female_head base_income i.district_id, cluster(village_id) gmm2s
estimates store T_gmm

* Print to screen
esttab T_ols T_2sls T_liml T_gmm, ///
    keep(D edu_head hh_size female_head) ///
    b(%9.1f) se(%9.1f) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2 widstat, labels("N" "R²" "KP-F (1st stage)") fmt(%9.0f %9.3f %9.2f)) ///
    mtitles("OLS" "2SLS" "LIML" "GMM") ///
    title("Table 1: Effect of Programme Participation on Monthly Income (BDT)") ///
    note("Cluster-robust SE in parentheses (clustered at village level)." ///
         "OLS is biased due to selection on unobserved motivation." ///
         "2SLS uses randomized village assignment (Z) as instrument." ///
         "LIML and GMM use Z and distance-to-centre as instruments." ///
         "KP-F = Kleibergen-Paap rk Wald F for weak-ID test.")

* Export to Excel
esttab T_ols T_2sls T_liml T_gmm using "ivreg2_results.csv", ///
    keep(D edu_head hh_size female_head) ///
    b(%9.1f) se(%9.1f) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2 widstat, labels("N" "R²" "KP-F") fmt(%9.0f %9.3f %9.2f)) ///
    mtitles("OLS" "2SLS" "LIML" "GMM") ///
    replace


*===========================================================================
* SECTION 16: DIAGNOSTICS CHECKLIST SUMMARY
*===========================================================================

di _newline(2)
di as txt "=========================================================="
di as txt " ivreg2 DIAGNOSTICS CHECKLIST"
di as txt "=========================================================="
di as txt ""
di as txt " □ 1.  Instrument relevance:"
di as txt "       Kleibergen-Paap F > 10 (ideally > 16 for 5% size distortion)"
di as txt ""
di as txt " □ 2.  Underidentification test:"
di as txt "       KP rk LM p-value < 0.05 (reject H0: underidentified)"
di as txt ""
di as txt " □ 3.  Overidentification (just-id: not testable):"
di as txt "       Hansen J p-value > 0.10 (fail to reject: IVs appear exogenous)"
di as txt ""
di as txt " □ 4.  Endogeneity test (Durbin-Wu-Hausman):"
di as txt "       p < 0.05 → D is endogenous → use IV"
di as txt "       p > 0.05 → D is exogenous → OLS is consistent"
di as txt ""
di as txt " □ 5.  Weak-IV-robust CI (Anderson-Rubin):"
di as txt "       Always report if KP-F is borderline (<16)"
di as txt ""
di as txt " □ 6.  Cluster SE at level of treatment assignment"
di as txt ""
di as txt " □ 7.  Exclusion restriction (UNTESTABLE by data alone):"
di as txt "       Z must not affect income EXCEPT through D"
di as txt "       Argue theoretically; use partial tests (C-stat, placebo)"
di as txt ""
di as txt " □ 8.  Monotonicity (for LATE interpretation):"
di as txt "       No defiers: Z=1 never reduces participation vs Z=0"
di as txt "=========================================================="


*===========================================================================
* SECTION 17: KEY EQUATIONS REFERENCE
*===========================================================================
/*
  ────────────────────────────────────────────────────────────────────
  2SLS in Two Stages:
    Stage 1:  D̂ᵢ = π₀ + π₁Zᵢ + Xᵢ'γ + νᵢ
    Stage 2:  Yᵢ = β₀ + β₁D̂ᵢ + Xᵢ'δ + εᵢ
  
  β₁ (2SLS) = Cov(Zᵢ, Yᵢ) / Cov(Zᵢ, Dᵢ)   (just-identified)
            = ITT estimate / First-stage compliance rate
            = LATE (under monotonicity)
  
  GMM moment condition: E[Zᵢ'εᵢ] = 0
  
  LIML: β_LIML = argmin  ε'PZε / ε'MXε   (eigenvector solution)
  
  Concentration parameter: μ² = π₁² × Var(Z) / Var(ν) × N
  ────────────────────────────────────────────────────────────────────
  
  ESTIMATOR CHOICE GUIDE:
  ──────────────────────────────────────────────────────────────────────
  Condition                       │ Recommended estimator
  ────────────────────────────────┼─────────────────────────────────────
  Strong IVs, homoskedastic       │ 2SLS
  Strong IVs, heteroskedastic     │ 2SLS + robust SE
  Strong IVs, clustered           │ 2SLS + cluster SE  ← MOST COMMON
  Overidentified, hetero.         │ Efficient GMM (gmm2s)
  Weak IVs (F < 10)               │ LIML or Fuller(4)
  Weak IVs, inference             │ Anderson-Rubin CI (weakiv)
  Many weak IVs (JIVE context)    │ HLIM or JIVE (see ivreg2 help)
  ──────────────────────────────────────────────────────────────────────
*/


log close

di as result _newline "Tutorial complete. Dataset saved: rct_household_survey.dta"
di as result "Results exported: ivreg2_results.csv"
di as result "Log file: ivreg2_tutorial_log.txt"
