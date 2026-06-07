/*=============================================================================
  ivreg2 TUTORIAL FOR STATA
  Using Household Survey Dataset of Randomized Control Trial Intervention 

  Author  : Ahmed Eshtiak | BIGD, BRAC University
  Dataset : rct_hh_survey.dta  (5,000 HHs · 90 villages · 3 districts)
            

  IV SETUP IN THIS DATASET
  ────────────────────────────────────────────────────────────────────────────
  Suppose programme TAKE-UP is imperfect and endogenous:
    · treat     = village-level random assignment (our INSTRUMENT, Z)
    · takeup    = actual participation (ENDOGENOUS variable, D)  ← we create
    · Outcome Y = income_el, pce_el, days_work_el, etc.

  Endogeneity story: motivated households both self-select into the programme
  AND earn more regardless → OLS(takeup → income) is upward-biased.
  Randomised village assignment (treat) is a valid instrument because:
    Relevance:  treat strongly predicts takeup (first stage F >> 10)
    Exclusion:  treat affects income ONLY through programme participation
    Exogeneity: treat is randomised → uncorrelated with unobservables

  Second instrument: vill_market_dist (distance to market) negatively
  predicts takeup (farther = harder to attend field sessions), and is
  plausibly exogenous to endline income conditional on controls.

  SECTIONS
  ────────────────────────────────────────────────────────────────────────────
  Part 0  — Setup, globals, and endogenous take-up variable
  Part 1  — OLS (naive, biased baseline)
  Part 2  — First-stage regression and instrument relevance
  Part 3  — Basic ivreg2 (just-identified 2SLS)
  Part 4  — Weak instrument diagnostics (KP-F, Stock-Yogo, Cragg-Donald)
  Part 5  — Overidentification tests (Hansen J, C-statistic)
  Part 6  — Endogeneity test (Durbin-Wu-Hausman)
  Part 7  — Standard error options (homoskedastic, robust, cluster-robust)
  Part 8  — LIML (Limited Information Maximum Likelihood)
  Part 9  — Fuller k-class estimator
  Part 10 — Two-step efficient GMM
  Part 11 — CUE (Continuously Updating Estimator)
  Part 12 — Multiple endogenous variables
  Part 13 — Absorbing fixed effects (absorb option)
  Part 14 — Weak-instrument-robust inference (Anderson-Rubin via weakiv)
  Part 15 — Multiple outcomes table with outreg2
  Part 16 — Diagnostics checklist and estimator choice guide
=============================================================================*/


*********
**# SETUP
*********

clear all
set more off

**# Directory Setup
global PATH "D:\Ahmed Eshtiak\Local Disk D\STATA Github\instrumental_variable_guide"

global RAW "$PATH\0_raw"
global DO "$PATH\1_do"
global RESULT "$PATH\2_result"
global CLEAN "$PATH\3_clean"

**# Required packages 
//ssc install ivreg2, replace
//ssc install ranktest, replace
//ssc install estout, replace
//ssc install weakiv, replace


**# Load the dataset 
use "$RAW/rct_hh_survey.dta", clear



**# CREATE ENDOGENOUS TAKE-UP VARIABLE
//80 % take-up rate in treated villages
// 10 % "leakage" in control villages (e.g., NGO spillovers)

set seed 20250601    

* Latent unobserved motivation (correlated with both takeup and income)
gen u_motivation = rnormal(0, 1)
label var u_motivation "Unobserved motivation (latent, for simulation)"

* Actual programme take-up (endogenous)
gen takeup_latent = -0.5 + 2.0*treat - 0.04*vill_market_dist + 0.3*u_motivation + rnormal(0, 0.4)
gen takeup = (takeup_latent > 0)
label var takeup "Actual programme take-up (1=Yes) [ENDOGENOUS]"

* Verify plausible take-up rates
tab treat takeup, row nofreq
quietly sum income_el
replace income_el = income_el + 200*u_motivation

* Convenience: standardise market distance for interaction models later
egen vill_market_std = std(vill_market_dist)
label var vill_market_std "Market distance (standardised)"

* Second binary instrument: far from market (>6 km)
gen far_market = (vill_market_dist > 6)
label var far_market "Far from market: >6 km (binary instrument)"

save "$CLEAN/rct_hh_survey_iv.dta", replace


**********************************************************
**# OLS — BIASED ESTIMATOR (Intention to Treatment Effect)
**********************************************************
regress income_el treat hh_size hh_female_hh hh_age_head hh_edu_head hh_land_owned asset_index_bl i.district, vce(cluster village_id)

estimates store ITT

//Coefficient on treat > 463, this is downward bias



**#  FIRST-STAGE REGRESSION
/*
  The first stage regresses the ENDOGENOUS variable (takeup) on the IV (treat) and all exogenous controls.

  What we look for:
  ─────────────────────────────────────────────────────────────────────────
  1. Coefficient on treat: large, positive, and significant
  2. F-statistic on excluded instrument(s): > 10 (rule of thumb)
     Kleibergen-Paap rk Wald F when using clustered SE (preferred)
  3. Partial R²: proportion of variation in takeup explained by treat
     after partialling out controls
  ─────────────────────────────────────────────────────────────────────────
*/

regress takeup treat hh_size hh_female_hh hh_age_head hh_edu_head hh_land_owned asset_index_bl i.district, vce(cluster village_id)
estimates store first_stage

* Test excluded instrument only
test treat //First-stage F on treat: 17533.14| F >> 10: treat is a strong instrument → 2SLS is reliable
//F < 10: instrument is weak → use LIML or Fuller 


*************************************
**# Instrumental Variable  Regression
*************************************
ivreg2 income_el (takeup = treat) hh_size hh_female_hh hh_age_head hh_edu_head hh_land_owned asset_index_bl i.district, cluster(village_id) first // also print first stage

estimates store LATE
// β₂SLS = LATE (Local Average Treatment Effect for compliers)

esttab ITT LATE, b(%9.1f) se(%9.1f) star(* 0.10 ** 0.05 *** 0.01)  mtitles("OLS (biased)" "2SLS (LATE)") title("OLS vs 2SLS: Effect of Take-up on Monthly Income") // Compare OLS vs 2SLS


**#Tests whether takeup is actually endogenous (DURBIN-WU-HAUSMAN)
//If takeup were exogenous, OLS would be consistent AND more efficient than IV → no reason to use IV.

ivreg2 income_el (takeup = treat)  hh_size hh_female_hh hh_age_head hh_edu_head hh_land_owned asset_index_bl i.district, cluster(village_id)  endog(takeup) 

ivreg2 income_el (takeup = treat) hh_size hh_female_hh hh_age_head hh_edu_head hh_land_owned asset_index_bl, cluster(village_id)


**# LIML — LIMITED INFORMATION MAXIMUM LIKELIHOOD
// When instruments are WEAK (F < 10), prefer LIML

ivreg2 income_el (takeup = treat far_market) hh_size hh_female_hh hh_age_head hh_edu_head hh_land_owned asset_index_bl i.district, cluster(village_id) liml


**# FULLER k-CLASS ESTIMATOR
//when instruments are borderline weak 

ivreg2 income_el (takeup = treat far_market) hh_size hh_female_hh hh_age_head hh_edu_head hh_land_owned asset_index_bl i.district,cluster(village_id) fuller(4)


**************************
**# TWO-STEP EFFICIENT GMM
**************************
// With homoskedastic errors: 2SLS = efficient GMM (same weight matrix)
// With heteroskedasticity:  GMM > 2SLS in efficiency (smaller SE)

* One-step GMM 
ivreg2 income_el (takeup = treat far_market) hh_size hh_female_hh hh_age_head hh_edu_head hh_land_owned asset_index_bl i.district,robust
estimates store GMM_1step

* Two-step efficient GMM with cluster-robust weighting
ivreg2 income_el (takeup = treat far_market) hh_size hh_female_hh hh_age_head hh_edu_head hh_land_owned asset_index_bl i.district, cluster(village_id) gmm2s
estimates store GMM_2step

