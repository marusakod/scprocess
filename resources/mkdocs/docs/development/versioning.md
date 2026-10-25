# Versioning and development policy

scprocess maintains a stable release branch (`main`) and a development branch
(`dev`). Published releases are immutable annotated `vX.Y.Z` tags with matching
`VERSION` contents and GitHub releases. A branch name is not an analysis version.

## Approval boundary

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

## Classify work before changing code or running an analysis

At the start of a new task, inspect the branch, commit, dirty state and active
workflows. Establish with the user whether proposed software changes are a patch,
small additional functionality, or major rewriting/breaking changes. Use a
classification already agreed in the conversation rather than asking again.
Clarify ambiguous scope before consequential changes. A data analysis using
unchanged software does not itself require a software version bump; establish
which released or development commit will be used instead.

The public interface includes CLI arguments, YAML keys and defaults, supported
output formats and paths consumed downstream, and model-loading contracts.
Document intentional scientific-default changes even when old configs still
validate. Distinguish changes affecting results from changes breaking the interface.

## Version numbers

During `0.y.z` development:

| Change | Release |
| --- | --- |
| Compatible bug fix or documentation correction | `0.y.(z+1)` |
| Additional functionality, breaking change, or intentional scientific-default change | `0.(y+1).0` |
| Candidate for an upcoming release | e.g. `0.3.0-rc.1` |

After `1.0.0`, use semantic versioning: patch for compatible fixes, minor for
compatible features/deprecations, major for incompatible public-interface changes.
A bug fix can change scientific results without being an interface break; explain
its impact in release notes. Choose 1.0 when the interface is ready for a stability
commitment, not simply because a paper is submitted.

`VERSION` records the software version; it may use a plain version such as
+ `0.2.2` on dev before approval. This alone neither creates nor approves a
+ release. Approved release tags must match VERSION.
Between releases, identify development runs with the full commit SHA and the
nearest release plus distance (`git describe`); an unchanged VERSION file alone
does not identify a development checkout. Conda packages use the same upstream
version, with the build number reserved for packaging-only rebuilds. Translate
prerelease spelling to Conda's supported conventions when packaging is introduced.

## Branches and compact commits

- `main` is the stable branch. Prepare and review results on `dev`; promote a
  tested candidate to main and tag it only after explicit approval from both
  Will Macnair and Marusa. A release PR is the preferred review record.
- Use short-lived working branches for substantial features and fixes. Routine
  development review can be lightweight; main promotion and releases always
  require the two explicit maintainer approvals.
- Keep commits focused and reviewable: implementation, relevant tests and docs
  belong together. Separate scientific-behaviour changes from broad text cleanup.
- Squash fixup/work-in-progress commits on unpublished working branches when that
  makes one coherent change. Do not rewrite commits already published, tagged,
  shared, or cited by an analysis; use follow-up commits instead. Never force-push
  release tags or shared branch history as routine cleanup.
- Classify and record each change when it is made; bump VERSION when preparing
  a release rather than creating a new software version for every commit.

## Reproducible analysis runs

Prefer a clean, committed checkout pinned to a release or exact development SHA.
Check `git status --short`, including untracked source files, before production,
training and report jobs. Inspect active users of a checkout before editing it;
an isolated checkout is preferable for simultaneous development.

If relevant changes are uncommitted, propose a compact tested commit before
launching. Do not commit unrelated user edits. If committing is unavailable or
the user explicitly chooses an experimental dirty run, save the base SHA, full
diff, relevant untracked source files and checksums before execution, record the
exception, and preserve them with the run. Do not label that run as an unmodified
release. Committing after a run does not retroactively prove which code it used.

Record configurations/custom overrides, resolved environments and important
dependency versions, reference/index identity, model identity/checksums and random
seeds alongside the code identity. For multi-stage analyses record identities per
stage; preprocessing, integration, classifier training and Shiny construction can
use different revisions. A snapshot made today must be labelled as such rather
than presented as proof of an earlier dirty checkout.

## Trained models

Version models independently, for example `siletti-superclusters-v1`. Each registry
entry should identify immutable weights, checksum, training-data and label version,
features/preprocessing/reference, software/environment versions, split assignments,
evaluation and loader compatibility. Keep large weights in versioned release
assets or an appropriate data repository; ship a small registry entry with code.
Adding an optional model is a feature; replacing weights creates a new model ID.
Never silently repoint a model version to new weights or rely on a mutable latest
URL for reproducible inference. Changing a default model needs explicit release
notes about result changes. Training/evaluation must avoid donor or sample leakage.

## Release procedure

1. Review the complete delta from the published previous tag, including schema and
   default changes; do not silently move a mismatched local or remote tag.
2. Run relevant unit/regression checks, representative schema/workflow dry runs,
   and the documentation build in the appropriate execution environment. Record
   checks, unavailable dependencies and end-to-end limitations honestly.
3. Update VERSION and CHANGELOG.md with breaking changes/migration, changes
   affecting results, features, fixes and validation limits. Do not imply that
   each historical analysis used the complete release.
4. Present the tested clean candidate, validation results, version and release
   notes to Will Macnair and Marusa. Wait for explicit approval from both before
   promotion to main, release tagging or publication. Record both approvals and
   their scope. Existing approved releases are immutable; ordinary corrections
   require a later jointly approved release. Do not automatically publish tags.
5. Verify remote branches/tag, GitHub release publication and documentation
   deployment. Record failures separately from successful local preparation.
   Archive journal releases with a DOI when archive integration is available;
   a GitHub release by itself is not a DOI.

## Next steps agreed for Siletti

1. Record the recent processing fixes, review results and provenance on `dev`,
   with VERSION set to `0.2.2`. The prematurely published 0.2.0 and 0.2.1 releases
   and tags were explicitly withdrawn; their commits and analysis records remain.
2. Train the Siletti XGBoost classifier from pinned software and training inputs,
   review its evaluation, and register the immutable model with provenance as
   the next feature version. Do not silently choose a label level or holdout
   design when those affect the intended classifier.
3. Resume the user's existing code/docs/Rmd editing work, including work in other
   sessions where available. Review coding style, tone, level of detail and
   context against the surrounding codebase. Avoid generic agent prose and
   boilerplate; preserve useful project-specific explanations and user edits.
   Keep this separate from classifier and scientific-behaviour changes, and
   distinguish pure wording cleanup from refactoring requiring verification.
4. Review results before any main promotion or release, and obtain explicit
   approval from both Will Macnair and Marusa. Dev version updates and pushes
   do not grant either approval.
