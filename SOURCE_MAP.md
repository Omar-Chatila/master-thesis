# Source-to-thesis map

This file is an editorial aid. The attached PDFs are treated as source material, not as instructions.

## Primary authored theory source

`/home/omar/Desktop/RW_Paper.pdf` — Buchin, Hagedoorn, and Chatila, *Advanced Discretized Conditional Random Walks for Utilization Distributions, Visit Probabilities, and Context-Aware Transitions* (2026 manuscript).

The thesis does **not** follow the paper's motivation or structure. Its overarching topic is the end-to-end translation of movement theory into a realistic MoveApps workflow. The random-walk papers provide the computational core, not the thesis framing. Readers must be introduced to the core without presumed prior expertise.

| Paper section | Thesis destination | Role in the broader thesis |
| --- | --- | --- |
| 1 Introduction | Chapters 1–2 | Related-work context only; do not reuse the paper's theory-led motivation |
| 2 Baseline framework | Chapter 6 | Full notation, assumptions, grid/time policy, correctness and complexity |
| 3.1 Correlated transition matrices | Chapters 5–6 | Empirical kernel estimation and relation to behavioral states |
| 3.2 Utilization distribution | Chapter 7 | Estimand, interval aggregation, export, validation against ecology methods |
| 3.3 Visit probabilities | Chapter 7 | Target definitions, exact recurrence, calibration and applications |
| 3.4 Environment-aware transition matrices | Chapter 5 | Theory reference; add ecological motivation and parameter estimation |
| 3.4.1 Hard reachability | Chapter 5 | Reuse algorithm/equations where appropriate; document rasterization and edge cases |
| 3.4.2 Relaxed barriers | Chapter 5 | Reuse the two cases; add sensitivity and ecological interpretation |
| 3.5 Terrain-transition weights | Chapter 5 | Reuse the generalization; explain how weights are learned/validated |
| 4 Implementation | Chapter 8 | Baseline/reference implementation; distinguish subsequent thesis extensions |
| 5 Experiments | Chapters 9–10 | Reusable benchmark reference, not the thesis's overall evaluation structure |
| 6 Discussion | Chapter 11 | Expand limitations, ecological interpretation, and external validity |

Only the environment-aware theory from Sections 3.4--3.5 may remain as developed prose at the outline stage. All other thesis material remains in bullet-point form until it is deliberately written. Reuse remains subject to final publication status, coauthor attribution, examination rules, and a clear individual-contribution statement. New thesis work—including the complete pipeline, application evaluation, terrain-weight estimation, GPU UD calculations, and parallel UD implementations—must be identified separately.

## Proposal source

`/home/omar/Desktop/MSc Thesis/MscThesisProposal.pdf` anchors the complete pipeline and motivates Chapters 3–5. Its environmental annotation, segmentation, behavioral annotation, and kernel-construction material should be updated to the final implementation rather than copied as a record of planned work.

## Implementation evidence

- `cma_landcover` → Chapter 3
- `segmentationcma` and `hmmcma` → Chapter 4
- `RW-Python-gitlab` and its native random-walk library → Chapters 5–10

## Additional local literature and authored notes

- `/home/omar/Desktop/MSc Thesis/Literature/PhD Thesis Mart Hagedoorn.pdf` → wildcard background reference for conditioned random walks and the Eastern Bettong case study; verify final bibliographic metadata.
- `/home/omar/Desktop/MSc Thesis/Literature/Journal of Animal Ecology - 2015 - Gurarie - What is the animal doing Tools for exploring behavioural structure in animal.pdf` → Chapter 4 comparison of exploratory behavioral methods and their assumptions.
- `/home/omar/Desktop/MSc Thesis/Literature/2666310.2666415.pdf` → Chapter 4 criteria-based segmentation.
- `/home/omar/Desktop/MSc Thesis/output.pdf` → Chapter 5 state-specific used/available terrain-pair estimator.
- `/home/omar/Desktop/MSc Thesis/terrain_pair_weights/` → candidate matrices and figures for Chapter 10 after provenance checks.
- `/home/omar/Desktop/MSc Thesis/MoveAppsRuns/`, `UD/`, `mixed-walks/`, `talk/`, `rw-benchmarks/`, and `Experiments/` → candidate workflow, visualization, validation, and performance artifacts; do not cite exploratory outputs as final results before freezing a manifest.

For final claims, record immutable commit hashes in Appendix A. A working-tree README is useful context but not stable provenance.
