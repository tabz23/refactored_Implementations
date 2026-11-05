# Original runs of the paper's results

These are the runs from which the submitted paper was made, including the runtimes it
reports. They were made between 2026-10-24 and 2026-11-03 (each folder's `provenance.txt`
and logs record the commit and date):

| Folder | Paper results |
|---|---|
| `results_hscc_2027_newbellmans_12workers_current_code/` | Figs. 2–5 and their runtimes (solver commit `8751e29`) |
| `results_quantitative_verification/` | §9.2.5 (solver commit `d8463b7`, verification scripts of 2026-11-03) |
| `results_bound_verification/` | Appendix H (Python only) |

**They give the same results as the current folders in `../`.** On 2026-11-04, all three were
regenerated from scratch with the current code (`reproduce_all.sh --updated`, solver commit
`9f6e912`; `di14` and `evasion12` were finished at `1572bab`, whose `cpp/` is the same). The
comparison logs are in this folder: `reproduce_updated.log` (first run, stopped when the disk
filled during `di14`) and `reproduce_updated_part2.log` (the rest). The result:

* Figs. 2–5: every phase of all eight cases has the same classification and the same
  `phases.csv` counts, and all 79 PNGs are pixel-identical.
* §9.2.5 and Appendix H: every output is identical (`verification_ra`, `verification_avoid`,
  `verification_ra_violation`, the error bounds and the unclassified-cell check of Theorem 3).

The regenerated folders replaced these ones in `../`. What changed is only the code:

* **Made more readable** (2026-11-04): the stopping tolerances are passed by their paper names,
  `--delta-lower` and `--delta-upper` (δ̲\* and δ̄\*), instead of `--delta-min`; `--gamma` is
  given only for `ra_discount`, since the other modes use γ = 1; `--tau` defaults to `--dt`;
  code the paper does not use was removed (the "conservative" stopping mode and
  `--delta-max`, the separate Algorithm-1 driver, `van_der_pol_avoid`; the Python port moved
  to `../../old_results/`); run folders have shorter names, e.g.
  `dubins_res20_epsmin1e-06_dlower1e-08_dupper1e-08_gamma0.96_dt0.3` instead of
  `dubins_res20_eps0.05_dmin1e-08_dmax1e-08_gamma0.96_dt0.3_tau0.3_std`.
* **The fixes we made** (2026-10-26 to 2026-11-04): Algorithm 2's threshold ε_min is given
  directly (`--eps-min`); the §9.2.5 checks use the policy and the floor time bounds of
  Proposition 4 as stated in the paper, simulate every sampled cell (example paths included)
  from its own random state, and check the second (violation) claim of Proposition 4;
  `reproduce_all.sh` compares every phase and PNG and can resume an interrupted run.

None of these changes alters a result, as the comparison above shows.

**Runtimes.** The paper keeps the runtimes of these original runs. The regenerated runs of
`vdp18` (417 s vs 292 s), `evasion12` (1,682 s vs 1,301 s), `di14` (4,225 s vs 1,910 s) and
`bicycle8` (79 s vs 64 s) were slower because the machine was swapping during the regeneration
(up to 30 GB of swap); the other four cases took about as long as here.

The bicycle out-of-memory record (`bicycle13_out_of_memory/`) moved from this `bicycle8` run
folder into the current one.
