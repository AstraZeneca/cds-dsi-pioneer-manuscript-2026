# PIONEER: Bayesian Joint Modelling of Mechanistic Tumour Growth and Time-to-Event Endpoints for Dynamic Prediction of Ongoing Oncology Trials

**Authors.** Karim Naguib, Roger Berché, Lu Li, Antonia Bevan, and Paul Metcalfe.
**Affiliation.** AstraZeneca, Oncology R&D, Data Science & AI.
**Year.** 2026.

**Working title (long form).** *PIONEER: Bayesian Joint Modelling of Mechanistic Tumour Growth and Time-to-Event Endpoints for Dynamic Prediction of Ongoing Oncology Trials*.

**Working title (short).** *PIONEER: Joint mechanistic-multistate Bayesian forecasting for oncology trials.*

**Target.** Methodology / biostatistics journal with strong oncology readership; primary candidates: *Statistics in Medicine*, *Biostatistics*, *CPT: Pharmacometrics & Systems Pharmacology*, *Journal of the Royal Statistical Society Series C*. Decide after first complete draft.

---

## 0. Status

This file is the **living outline** of the PIONEER methods paper. Each section gets fleshed out in a later commit; for now we have:

- §1 (Introduction + Literature Review) — first draft below.
- §§2–7 — subsection bullets, anchor points and figures-to-include only.
- §8 (To-Dos) — open questions and missing material.

---

## 1. Introduction

### 1.1 Why we wrote this paper

Modern oncology trials, particularly those evaluating immuno-oncology and targeted agents, are run under intense pressure to deliver early, defensible read-outs. Sponsors must commit to dose, schedule, expansion cohorts and indication choices on the basis of immature data — typically a Phase II or interim Phase III data cut-off (DCO) where the longitudinal tumour size record is rich but progression-free survival (PFS) and overall survival (OS) are heavily right-censored. Three needs follow from this state of affairs.

**(i) Forecasting from immature follow-up.** Decision makers want a credible projection of what PFS and OS curves will look like at the protocol-specified mature DCO, not just an estimate of the current Kaplan–Meier (KM) curve.

**(ii) Coherent uncertainty.** A point estimate of "median PFS at month 24" without a calibrated prediction interval is of no use for go/no-go logic; what is needed is a posterior predictive distribution over the full survival curve, jointly with all other endpoints.

**(iii) Mechanistic transportability.** Patient-level dynamics observed in one trial, indication or treatment line should be re-usable as priors or hierarchical components in another. Hazards-only models tied to baseline covariates do not deliver this — biological parameters do.

PIONEER is the framework AstraZeneca has built to satisfy these three needs simultaneously. The model couples a *mechanistic* two-compartment state-space submodel of longitudinal tumour size (sum of longest diameters, SLD) to a *multistate proportional-hazard* submodel for the clinical events of interest (progression, death without progression, post-progression death, dropout, off-trial death) inside a single Bayesian hierarchical framework. All five clinical endpoints — target lesion progression, multistate PFS, OS, objective response rate, and off-treatment-time (OTT)-style summaries — are derived in the same generated-quantities pass from the joint posterior, so every reported quantity automatically inherits the model's full uncertainty. The mechanistic submodel feeds the multistate submodel as a set of *time-varying covariates* (standardised log SLD, decrease and growth rates, optionally SLD velocity); the multistate submodel feeds back through the joint likelihood, so neither half is fitted in isolation.

### 1.2 What is novel about the approach

The literature on PFS/OS prediction from tumour kinetics is by now substantial (see §1.3), and it is therefore worth being precise about what is and is not new in PIONEER. The novelty has four facets, in decreasing order of how distinctive each one is.

**(N1) Mechanistic dynamics as structured time-varying covariates of a *multistate* hazard.** Existing tumour-dynamics-to-OS frameworks (Stein–Bruno–Claret lineage; Tardivon 2019; Chan 2021; Kerioui 2020) almost universally feed a one-dimensional tumour metric — typically the model-predicted SLD or a tumour growth rate constant — into a *single* survival hazard, usually a Cox or parametric proportional-hazards model for OS. PIONEER replaces that single hazard with the full illness-death (and dropout) multistate, and *each* transition is allowed its own time-varying tumour-derived covariate vector with transition-specific coefficients. This means the model can express, for example, that growth rate predicts progression hazard but only modestly affects post-progression death hazard, which depends instead on the patient's tumour state at progression. Standard joint models cannot make that distinction because they have only one event channel.

**(N2) Two-compartment decrease/growth dynamics inside the joint model.** The decrease/growth (D+G) decomposition that descends from Stein 2008 / Claret 2009 is well-known as a *summary* model fit to SLD trajectories, and joint papers usually link a fitted D+G to survival in a *two-stage* fashion (estimate dynamics first, plug into survival second). In PIONEER the two-compartment decomposition is the *latent state* of a state-space model that is fitted *jointly* with the multistate, in log-baseline-normalised space, with a Student-t observation likelihood, optional AR(1) process noise at population and patient level, and explicit handling of left-censoring at the SLD limit of detection. We are not aware of a published joint-model implementation that combines the two-compartment latent-state structure, multistate survival, hierarchical borrowing across trials, and full Bayesian inference in one self-contained framework.

**(N3) Generic multi-level hierarchy with QR-rotated covariate machinery.** The hierarchy is *level-agnostic*: the same code accommodates "population → trial → patient", "population → trial → arm → patient", or "population → indication → trial → patient", with each module — tumour regression rate (`tr`), growth fraction (`frac`), initial state (`init`) and each multistate transition — independently choosing whether random intercepts and random slopes are active at each level. Borrowing across trials (e.g. learning growth-rate priors from HISTORICAL and applying them to SCLC-01) drops out of this design without needing a separate meta-analytic step. Population covariate effects are entered via QR-decomposed design matrices, and elicited priors are mapped through the QR rotation, which preserves interpretability while avoiding the ill-conditioning that bites when 20+ correlated baseline covariates enter the same linear predictor.

**(N4) End-to-end uncertainty propagation from priors to clinical readouts.** Because all endpoints are derived in `generated quantities` from the joint posterior, the posterior predictive distribution of *any* clinical summary — median PFS, twelve-month OS rate, objective response rate, hazard ratio between two arms at a fixed DCO — is obtained by reading off the corresponding draws and summarising. There are no point estimates "fixed in" and no separate post-hoc sensitivity analyses for upstream uncertainty. The model also produces unconditional ("super-population", `spop_*`) and observed-data-conditioned (`sample_*`) predictive distributions side-by-side, which lets us answer *both* "what would happen to a future patient like ours?" and "what will happen to *these* patients past the data cut-off?" within the same fit.

We argue (§5–§6) that these four ingredients together change the *use case* the model can serve: instead of a one-shot OS extrapolation at end of Phase II, PIONEER is fit repeatedly over the life of a trial, providing dynamically updated forecasts of mature PFS/OS as more data accrues and as the model is borrowed across trials in the same indication.

### 1.3 Literature review

We organise the prior literature into six arcs, each of which feeds either the methods (§3, §4) or the discussion (§7).

#### A. Tumour growth-inhibition (TGI) and the decrease/growth decomposition

The starting point for any longitudinal-tumour-to-survival framework is the empirical observation that tumour size kinetics correlate strongly with patient survival, and that this correlation is captured well by a two-exponential "decrease + growth" decomposition. Stein and colleagues introduced the formulation in advanced prostate cancer [@stein2008; Stein 2011], showing that the growth-rate constant in particular is a strong correlate of OS. Wang 2009 extended the framework to non-small-cell lung cancer (NSCLC). Claret and co-authors [Claret 2009; Claret 2012] turned the framework into a forecasting tool, demonstrating that Phase II tumour dynamics can be combined with a parametric OS submodel to predict Phase III OS, with calibration successes in colorectal and NSCLC contexts. Bruno 2014 consolidated this body of work into a methodological review, and Bruno 2020 argued explicitly for tumour-dynamic modelling as the standard supporting tool for go/no-go and dose-decision use-cases. The framework is by now the dominant TGI/OS paradigm in pharmacometrics. PIONEER's `tr` (total rate) and `frac` (decrease/growth allocation) modules implement the same biological idea, but as the *latent state* of a state-space model fitted jointly with multistate survival rather than as a deterministic regression target.

#### B. Population mixed-effects modelling of tumour size

Mixed-effects (population-PK/PD-style) modelling of tumour size — typically NONMEM-style nonlinear mixed-effects models — is the methodological foundation that the TGI literature rests on. The canonical methodological review is Ribba 2014 (CPT:PSP), which catalogues the family of models (Simeoni, Claret, Stein, Wang) and the population-level inference machinery. Holford 2013 and the broader pop-PK lineage justify the multi-level structure used in §3. More recent contributions include Tardivon 2019, who applied a population NLME tumour-size model jointly with an OS model in atezolizumab-treated urothelial carcinoma; and Desmée 2017, who introduced HMC-based Bayesian inference for nonlinear joint models of tumour kinetics and survival in metastatic prostate cancer. PIONEER's multi-level hierarchy generalises the trial-and-patient split of these works to an arbitrary number of levels with module-specific enable flags.

#### C. Joint longitudinal-survival modelling

Since the PFS/OS endpoints depend on a longitudinal biomarker that is itself patient-specific and noisily observed, naïve insertion of the biomarker as a time-varying covariate in a Cox model is biased: the biomarker is endogenous because it shares patient-level frailty with the survival outcome. The textbook treatment is Rizopoulos 2012, with Hickey 2016 reviewing more recent extensions to multiple longitudinal outcomes and competing risks. Kerioui 2020 is the closest published methodological cousin to PIONEER: a Bayesian HMC joint model in an immuno-oncology setting, with a tumour-size submodel feeding an OS hazard. Kerioui 2022 gives a non-technical introduction to the family. PIONEER differs from these in three ways: (a) it uses a *multistate* hazard rather than a single OS hazard; (b) the longitudinal submodel is *two-compartment state-space* rather than nonlinear regression; (c) the longitudinal-to-hazard linkage uses *multiple* time-varying tumour features (log SLD, log decrease rate, log growth rate, optional velocity) with transition-specific coefficients.

#### D. Immuno-oncology-specific tumour kinetics

Immuno-oncology (IO) regimens exhibit dynamics — delayed responses, pseudo-progression, slow regrowth after a deep response — that monomerically-summarised (single-time-point or single-rate-constant) models can fail to capture. Mistry 2019 examined resistance kinetics in NSCLC under EGFR inhibition. Yin 2019 fitted population PK and tumour growth dynamics for pembrolizumab. Chan 2021 developed a TGI-to-OS framework specifically for atezolizumab across solid tumours, which is the closest published analogue to the use case PIONEER serves for durvalumab (SCLC-01) and other AstraZeneca IO assets. Bruno 2023 used the framework explicitly for decision support on a Phase Ib/II combination, with retrospective re-sampling of IMpower150 to validate. PIONEER inherits this lineage; the novelty over Chan 2021 is the multistate decomposition of OS into competing pathways, which we argue is more clinically faithful for IO trials in which a non-trivial fraction of OS events occur without prior PD.

#### E. State-space / mechanistic modelling rationale

Phenomenological state-space models of tumour size have a smaller but consistent literature in oncology pharmacometrics, and Bayesian state-space approaches are relatively recent. Dean 2010 made the original clinical-decision case for tumour growth rate as a routine summary statistic. The case for *state-space* (as opposed to nonlinear regression) is that the latent state can absorb measurement noise and limit-of-detection censoring coherently; PIONEER's Student-t observation density, log-baseline normalisation, and AR(1) process-noise option (§3) sit in this tradition.

#### F. Multistate survival models in oncology

Multistate (illness-death) extensions of competing-risks survival have been used extensively in chronic-disease epidemiology and increasingly in oncology — particularly when post-progression survival is itself of interest, when dropout is informative, or when *PFS-as-surrogate-for-OS* arguments must be examined explicitly (Prentice 1989; Buyse et al. 2000; recent overviews in Hougaard 1999). Our model implements an illness-death structure with five transitions ($0\to1$ progression, $0\to2$ direct death, $1\to2$ post-progression death, $0\to3$ dropout, $3\to2$ off-trial death), each with its own Gaussian-process baseline hazard. Allowing arm or treatment effects on the post-progression hazards $1\to2$ and direct-death hazard $0\to2$ — versus pooling them, which structurally enforces PFS surrogacy — is a configurable model choice, and we will use it in §6 to discuss surrogacy explicitly. Compared with prior IO joint models (which collapse OS to a single hazard), this formulation is what permits us to *test* surrogacy as a model-fit comparison rather than assume it.

#### G. RECIST and the regulatory anchor

Since SLD and RECIST PD are the regulatory observables, no longitudinal-tumour modelling paper can avoid Eisenhauer 2009 (RECIST 1.1). PIONEER computes RECIST classifications inside the Stan generated-quantities block from the latent state trajectory, which means PFS-by-RECIST is itself a posterior quantity. The FDA's MIDD (model-informed drug development) program documents [FDA MIDD references — TBD] explicitly endorses tumour-dynamic modelling as supporting evidence for end-of-Phase-II decisions, which is the primary use case PIONEER targets.

### 1.4 How novel is PIONEER, candidly?

There is no single ingredient in PIONEER that has not appeared somewhere in the literature: two-compartment dynamics (Stein/Claret/Bruno), Bayesian HMC joint modelling (Desmée/Kerioui), multistate models in oncology (Hougaard; multistate-survival textbooks), GP baseline hazards (Riihimäki & Vehtari and the general GP-survival literature), QR-rotated regression in Stan (Stan User's Guide), TGI→OS in IO (Chan; Bruno 2023). What is novel is the combination — and specifically, three combinatorial claims that we will defend in §5:

1. The two-compartment latent state and the multistate hazard are fitted *jointly*, not in stages, so endogeneity bias is correctly handled.
2. Each multistate transition has its own time-varying-covariate coefficients on the tumour-derived features, so the model can express clinically distinct mechanisms — e.g. that growth rate matters for $0\to1$ but not for $1\to2$.
3. All endpoints (PFS, OS, ORR, OTT, hazard ratios at any DCO) are read off the joint posterior, with two predictive flavours (super-population and observed-conditioned) available without re-fitting.

We will avoid two claims that would be unsafe under reviewer scrutiny:

- We will *not* claim the two-compartment decomposition recovers two real clonal populations from SLD data. It does not. The argument for the decomposition is forecasting performance, parameter transportability, and clinical interpretability — not biological mechanism. (See `docs/lu_learn/why_mechanistic.md` §1.)
- We will *not* claim the model improves *in-sample* fit relative to a flexible empirical (spline / GAM) baseline. It almost certainly does not. The argument is *out-of-sample forecasting and parameter transport*, which we evidence with leave-future-out cross-validation (LFO; §6) and across-trial transfer (§6).

### 1.5 Roadmap

§2 fixes notation. §3 specifies the mechanistic submodel. §4 specifies the multistate submodel and the time-varying covariate bridge. §5 describes priors and the multi-level hierarchy. §6 reports applications: SCLC-01 (durvalumab combination, NSCLC), with secondary illustrations from LUNG and Pioneer; LFO performance; cross-trial borrowing. §7 discusses limitations, the surrogacy question, and outlook for routine deployment. §8 concludes.

---

## 2. Notation and Setup

*Outline only — to be written after §3, §4 stabilise.*

- Patients $i=1,\dots,N$. Trial / arm / level indices.
- Two clocks: study time $t$ (week 0 = first treatment) and calendar time.
- Visits: ragged; visit set $\mathcal{V}_i$, observed SLDs $y_{ij}$ at week $t_{ij}$.
- Observed events: $T^{01}_i, T^{02}_i, T^{12}_i, T^{03}_i, T^{32}_i$, censoring indicators, final state $S_i \in \{0,1,2,3\}$.
- Covariates: design matrix $X$, QR-decomposed $Q$, $R$.
- Hierarchy: $n_\text{levels}$, $n_\text{groups\_per\_level}$, `patient_level_groups`.

---

## 3. Mechanistic Submodel: Two-Compartment State-Space SLD

*Outline.*

- Latent log-states $D_i(t)$, $G_i(t)$.
- Total SLD reconstruction: $\log(\text{SLD}_i(t)/\text{SLD}_i(0)) = \log(\exp D_i + \exp G_i)$.
- Baseline normalisation; why log-baseline-normalised space.
- Dynamics: $D_{t+1} = D_t - r^\text{dec}_i(t)$, $G_{t+1} = G_t + r^\text{grow}_i(t)$.
- `tr` × `frac` decomposition: $\log r^\text{dec} = \mu^\text{tr} + \log\sigma(\eta^\text{frac})$ etc.
- Initial state via `init_logit_loc_patient`.
- Optional AR(1) process noise (population and patient).
- Observation density: Student-t on log-normalised SLD, left-censored at LoD.
- Outputs that flow to multistate: `states_full_grid`, `patient_log_decrease_rate`, `patient_log_growth_rate`, derived $\log\text{SLD}$ and SLD velocity.

**Figures.**
- Fig. 3.1 — model-architecture TikZ (already in `quarto/website/images/model-diagram-tikz.svg`).
- Fig. 3.2 — example of D + G trajectories with observed SLD overlaid; one IO patient with delayed response and one fast progressor.

---

## 4. Multistate Submodel and the Tumour→Hazard Bridge

*Outline.*

- Illness-death + dropout state diagram (already in `quarto/website/images/multistate-diagram.svg`).
- Five transitions $0\to1, 0\to2, 1\to2, 0\to3, 3\to2$; each independently switchable.
- Per-transition hazard:
  $$\log\lambda_{jk,i}(t) = \log\lambda^\text{pop}_{jk}(t) + \sum_\ell \log\lambda^{\text{level }\ell}_{jk}(t) + \beta^\text{tv}_{jk} W_i(t) + (Q\beta^\text{ti(QR)}_{jk})_i + \text{level slopes}.$$
- Gaussian-process baseline on coarse knot grid (`ms_gp_grid_step`), pushed to weekly resolution.
- Time-scale options for $1\to2$: Markov, semi-Markov, or extended (sojourn + clock).
- Time-varying tumour covariates $W_i(t)$: standardised log SLD, log decrease rate, log growth rate, optional SLD velocity.
- `ms_prog_deterministic`: how we avoid double-counting RECIST-PD events between mechanistic and hazard halves.
- Likelihood with interval censoring on $0\to1$ (`ms_ic_gap_01`).

**Figures.**
- Fig. 4.1 — multistate state diagram (existing).
- Fig. 4.2 — illustrative posterior baseline hazards $\lambda_{01}, \lambda_{02}, \lambda_{12}$ from SCLC-01.
- Fig. 4.3 — coefficient forest plot for time-varying tumour covariates by transition.

---

## 5. Priors, QR Machinery, and the Multi-Level Hierarchy

*Outline.*

- Hierarchy as level-agnostic: see §1 of `docs/multi_level_hierarchy_design.md`.
- Module-level enable flags (`enable_level_intercept_*`, `enable_level_cov_*`).
- Non-centred parameterisation everywhere; raw $\sim N(0,1)$.
- QR-rotated covariates: elicited prior $(\mu, \Sigma)$ rotated by $R$.
- Prior choices: `tr_loc_pop`, `frac_logit_loc_pop`, `init_logit_loc_pop` (with explanation of why the original wide priors were tightened).
- GP hyper-priors: `inv_gamma` on length-scale and marginal SD.
- Process-noise priors when AR(1) is enabled.
- Cross-trial borrowing example: HISTORICAL-derived priors transferred to SCLC-01.

**Figures.**
- Fig. 5.1 — prior vs posterior for each population-level intercept (`tr_loc_pop` etc.).
- Fig. 5.2 — schematic of the multi-level hierarchy with module enable matrix.

---

## 6. Applications

*Outline.*

### 6.1 SCLC-01 (primary case study)
- Trial description; immature DCO at the time of the prediction.
- Mature DCO comparison (held-out): predicted vs observed PFS and OS KM.
- Predicted ORR, durable response, twelve-month OS rate.

### 6.2 Borrowing across trials (HISTORICAL → SCLC-01)
- Hierarchical pooling at the trial level.
- Quantitative impact on posterior tightness vs no-borrow alternative.

### 6.3 Leave-future-out cross-validation
- Definition (existing `quarto/.../leave-future-out-cross-validation.qmd`).
- ELPD comparison vs pure baseline-covariate hazard model and vs single-hazard joint model.
- Argument: out-of-sample is where the mechanistic + multistate combination earns its keep.

### 6.4 Surrogacy and treatment-effect decomposition (illustrative, not regulatory)
- Configuration with arm-level random effects on $1\to2$ and $0\to2$.
- Posterior decomposition of $\Delta\text{OS}$ into mediated-via-PFS vs direct components.
- Compare with Prentice-style criteria.

### 6.5 Pioneer (PSA + OS, brief)
- Same machinery applied to prostate cancer with PSA replacing SLD.
- One- or two-page demonstration that the framework generalises.

**Figures.**
- Fig. 6.1 — predicted-vs-observed PFS / OS KM at mature DCO (SCLC-01).
- Fig. 6.2 — borrowing impact (posterior on $\mu^\text{tr}$ with vs without HISTORICAL).
- Fig. 6.3 — ELPD-LFO bar chart, three model variants.
- Fig. 6.4 — surrogacy decomposition.

---

## 7. Discussion

*Outline (bullets).*

- What PIONEER buys you and what it does not.
- Limitations: phenomenological D+G; multistate identifiability under sparse trial-level data; computational cost.
- When to prefer a single-hazard joint model (small trial, target endpoint = OS only, no dropout signal).
- Future work: PSA / circulating-tumour-DNA dynamics; multi-biomarker latent states; integration with model-based meta-analysis.
- Regulatory positioning vs FDA MIDD.

---

## 8. Open To-Dos for the Manuscript

- [ ] Pull exact abstracts/quotes from Kerioui 2020 and Bruno 2020.
- [ ] Add Hougaard / multistate-survival textbook reference; check Buyse 2000 for surrogacy.
- [ ] Fill in the FDA MIDD reference placeholder (currently `[FDA MIDD references — TBD]`).
- [ ] Decide citation style (Vancouver vs APA) and convert references in `bib/pioneer.bib`.
- [ ] Add a model-comparison subsection (§6.3) showing mechanistic out-forecasts a baseline-covariate hazard model — already partially evidenced in `quarto/website/analysis/model-validation.qmd`.
- [ ] Confirm author affiliations (single-affiliation footnote vs multiple).
- [ ] Decide whether Pioneer (PSA / OS) gets its own §6.5 or moves to a separate paper.
- [ ] Write a short pre-registration paragraph for the SCLC-01 forecast: which DCO is the held-out target, what is the locked model spec.
- [ ] Add reproducibility statement: GitHub orgs (`azu-oncology-rd`), Stan model files, frozen renv.lock.
- [ ] Decide on inclusion of a sensitivity analysis to the AR(1) process-noise toggle.
