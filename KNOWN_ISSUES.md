# Known issues

Issues found in the review of the migration to the sample x subunit / three-column distance scheme.
No correctness problems were found in the indices themselves; every item below is minor.

- **[decision] MAD with a representative absent from the sample.** The mean distance of the cluster then leaves out the representative (its zero self-distance), while a present representative is included, as in the legacy implementation. Check against the definition in Finn (2024). 
- **[Docs] Vignette explains the effect of `sig` backwards.** A smaller `sig` does increase distance-based multiplicity, but the text says this means "less diversity lost to clustering", while `multiplicity.distance` documents higher values as more diversity lost. [vignettes/divermeta.Rmd:227](vignettes/divermeta.Rmd#L227)
- **[Check] `R CMD check` NOTEs from the visualisation functions.** Undefined global variables `value` and `mag`, and a help link to `ade4`, which is not a dependency. [R/visualize_dist.R:30](R/visualize_dist.R#L30)
- **[Check] ggplot2 is used without checking it is installed.** It is only in Suggests. The `visualize_*` examples are now guarded with `@examplesIf`, but the vignette still calls them unconditionally, so it fails to build where ggplot2 is missing.
- **[Site] The committed pkgdown site (`docs/`) is stale.** It has pages for removed functions (`cluster_distance_matrix`, `*.by_blocks`, `dist_quadratic_form`, `convert_to_dist_indices`) and none for the new ones.
- **[Site] The pkgdown site publishes internal notes.** [docs/INSTRUCTIONS.md](docs/INSTRUCTIONS.md) (an old internal prompt) and `docs/TODO.md` are published. Make sure the current `INSTRUCTIONS.md` and this file are not published on the next build.
- **[Release] Breaking API change without a changelog.** Samples are now in rows, `clusters` became `clust` and distances are a three-column table, but there is no `NEWS.md` and the version is still 0.0.3.
- **[Repo] The migration is not committed.** The new files in `R/`, `man/` and `tests/` are untracked and `R/support_functions.R` is only staged for deletion.
