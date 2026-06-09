# PIONEER 2026 — Working Notes & To-Dos

Working file for the PIONEER methods paper. Companion to `PIONEER_2026.qmd`. Holds material that does not belong in the manuscript itself: status, deferred decisions, open to-dos, and informal author-voice passages that were drafted as scaffolding and removed from the journal version.

---

## Status (snapshot)

**Target journal.** Primary: *Statistics in Medicine*. Secondary: *CPT: Pharmacometrics & Systems Pharmacology*. Decide after first complete draft.

**Manuscript structure (Stat-Med-shaped):**

1. Introduction (forecasting need, gap, contribution, roadmap)
2. Motivating example: SCLC-01
3. Model (notation; mechanistic; multistate; joint specification; priors & hierarchy)
4. Inference & implementation
5. Simulation study
6. Application to SCLC-01
7. Discussion

**Drafting state.**
- §1 (Introduction) — first draft.
- §3.1, §3.2 (Notation, Mechanistic submodel) — first draft.
- §2, §3.3, §3.4, §3.5, §4, §5, §6, §7 — outline only.

---

## Removed from manuscript: "How novel is PIONEER, candidly?"

This section was scaffolding written in informal voice. The claims it makes are already covered, in proper journal voice, by N1–N4 in §1.2 of the manuscript. Kept here for reference.

> There is no single ingredient in PIONEER that has not appeared somewhere in the literature: two-compartment dynamics (Stein/Claret/Bruno), Bayesian HMC joint modelling (Desmée/Kerioui), multistate models in oncology (Hougaard; multistate-survival textbooks), GP baseline hazards (Riihimäki & Vehtari and the general GP-survival literature), QR-rotated regression in Stan (Stan User's Guide), TGI→OS in IO (Chan; Bruno 2023). What is novel is the combination — and specifically, three combinatorial claims that we will defend in §5:
>
> 1. The two-compartment latent state and the multistate hazard are fitted *jointly*, not in stages, so endogeneity bias is correctly handled.
> 2. Each multistate transition has its own time-varying-covariate coefficients on the tumour-derived features, so the model can express clinically distinct mechanisms — e.g. that growth rate matters for $0\to1$ but not for $1\to2$.
> 3. All endpoints (PFS, OS, ORR, OTT, hazard ratios at any DCO) are read off the joint posterior, with two predictive flavours (super-population and observed-conditioned) available without re-fitting.
>
> We will avoid two claims that would be unsafe under reviewer scrutiny:
>
> - We will *not* claim the two-compartment decomposition recovers two real clonal populations from SLD data. It does not. The argument for the decomposition is forecasting performance, parameter transportability, and clinical interpretability — not biological mechanism.
> - We will *not* claim the model improves *in-sample* fit relative to a flexible empirical (spline / GAM) baseline. It almost certainly does not. The argument is *out-of-sample forecasting and parameter transport*, which we evidence with leave-future-out cross-validation (LFO; §6) and across-trial transfer (§6).

**Why removed.** Journals reject self-deprecating "candidly" framings. The N1–N4 list in §1.2 already states what we *do* claim; explicit anti-claims belong in the discussion's limitations paragraph if anywhere.

---

## Deferred / demoted from §6 to §7

Two §6 subsections that were demoted to brief mentions in the Discussion:

### Surrogacy and treatment-effect decomposition

- The framework supports a configuration with arm-level random effects on $1\to2$ and $0\to2$ that lets us posterior-decompose $\Delta\text{OS}$ into mediated-via-PFS and direct components and compare with Prentice-style criteria.
- Decision: out of scope for the primary paper. Surrogacy is a regulatory question with its own literature; bolting it onto the methods paper would dilute focus. Mention as future work in §7.

### Pioneer (PSA + OS)

- Same machinery applied to prostate cancer with PSA replacing SLD; demonstrates the framework generalises beyond SLD/RECIST.
- Decision: defer to a separate paper. Mention as a generalisation pointer in §7.

---

## Open to-dos

### Manuscript content
- [ ] Pull exact abstracts/quotes from Kerioui 2020 and Bruno 2020.
- [ ] Add Hougaard / multistate-survival textbook reference; check Buyse 2000 for surrogacy.
- [ ] Fill in the FDA MIDD reference placeholder (currently `[FDA MIDD references — TBD]`).
- [ ] Decide citation style (Vancouver vs APA) and convert references in `bib/pioneer.bib`.
- [ ] Add a model-comparison subsection in §6 showing mechanistic out-forecasts a baseline-covariate hazard model — already partially evidenced in `quarto/sclc/website/analysis/model-validation.qmd`.
- [ ] Confirm author affiliations (single-affiliation footnote vs multiple).
- [ ] Write a short pre-registration paragraph for the SCLC-01 forecast: which DCO is the held-out target, what is the locked model spec.
- [ ] Add reproducibility statement: GitHub orgs (`azu-oncology-rd`), Stan model files, frozen renv.lock.
- [ ] Decide on inclusion of a sensitivity analysis to the AR(1) process-noise toggle.
- [ ] Write the abstract (250 words, structured: Background / Methods / Results / Conclusions).

### Sections to draft
- [ ] §2 Motivating example: SCLC-01 — trial design, ctDNA-stratified subgroup, DCO timeline, what we forecast at the immature cut-off.
- [ ] §3.3 Multistate submodel — flesh out the outline.
- [ ] §3.4 Joint specification — how the two submodels share patients; the bridge $\mathbf{W}_i(t)$; why fitted jointly.
- [ ] §3.5 Priors & multi-level hierarchy — flesh out the outline.
- [ ] §4 Inference & implementation — HMC in Stan, multi-chain, NUTS settings, runtime, divergences, posterior storage (parquet draws).
- [ ] §5 Simulation study — design, operating characteristics targets (bias, coverage, posterior calibration).
- [ ] §6 Application to SCLC-01 — predicted-vs-observed PFS/OS at mature DCO, LFO, cross-trial borrowing.
- [ ] §7 Discussion — limitations, surrogacy mention, Pioneer/PSA mention, future work.

### Submission prep (much later)
- [ ] Pick target journal (Stat Med vs CPT:PSP).
- [ ] Get the right LaTeX class (`WileyNJD-v2.cls` for Stat Med).
- [ ] Switch Quarto YAML `documentclass:` to journal class.
- [ ] Strip the `bib/pioneer.bib` of unused entries.
- [ ] Run reference-style conversion.

---

## Decisions log

- **2026-05-22.** Restructured manuscript from "living outline" to Stat-Med-shaped sections. Moved Status, To-Dos, and the candid-novelty section here. Demoted surrogacy and Pioneer to brief mentions in Discussion.
- **2026-05-22.** Decided not to install the JASA Quarto template at this stage; will pick a journal class only at submission time.
