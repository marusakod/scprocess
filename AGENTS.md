# AGENTS.md

Instructions for coding agents working in this repository.

## Versioning, task scope and analysis provenance

### Required approval before main or release publication

Commits and pushes to `dev` are allowed within the authorized task. They do not
authorize promotion to `main` or a release. Before pushing/merging to `main`,
creating or pushing a release tag, or publishing any GitHub/Conda release,
obtain explicit approval from **both Will Macnair and Marusa** for the concrete
candidate commit, version and release contents. If either approval is absent or
unclear, leave the reviewed work on `dev` and ask; never infer Marusa's approval
from Will's approval, a generic request to fix/commit work, or passing tests.
Preparing and incrementing VERSION on `dev` is not release approval. Review
results must be available before seeking the two approvals. Material changes
after approval require renewed approval of the changed candidate.

Automatic publication on tag pushes is disabled. Do not re-enable publication
automation without an approval mechanism that requires both maintainers. The
user's explicit instruction to withdraw the premature releases and restore main
is authorization for that rollback only, not a new release.

Follow [the versioning and development policy](resources/mkdocs/docs/development/versioning.md).
At the start of a new task, inspect the checkout and active jobs, and establish
with the user whether software work is a patch, small additional functionality,
or major rewriting/breaking changes. Reuse a classification already agreed in
the conversation. For an analysis using unchanged software, identify the pinned
software version instead of inventing a release bump.

Keep implementation, meaningful tests and relevant docs in neat, compact commits.
Squash unpublished fixup/WIP commits when useful; never rewrite published,
shared, tagged or analysis-cited history without explicit authorization. Preserve
unrelated user edits. Prefer clean, committed checkouts for processing, training
and report runs; document any explicitly chosen dirty-run exception with a base
SHA, complete patch and relevant untracked source snapshots before execution.

`main` is the tested release branch; `dev` integrates development. For 0.x,
compatible fixes use patch releases, while features, breaks and intentional
scientific-default changes use minor releases. After 1.0, use normal semantic
versioning. VERSION, annotated release tag and GitHub release must agree. Record
breaking changes and changes affecting results separately. Models have immutable
independent IDs/checksums. Record analysis provenance per stage; a later commit
does not retroactively identify earlier dirty runs.


## Project-specific Instructions for Development

- When you have finished working on some functionality, you MUST check that the docs (in resources/mkdocs) are up-to-date.

Inspect this branch directly before development: dev may contain functionality
that has not been reviewed or promoted to main.
