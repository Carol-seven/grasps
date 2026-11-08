# Testing policy

- Batch related changes before verification to save time and tokens.
- During development, run only the targeted tests needed for the
  changes.
- Do not run the full test suite or thorough package checks after each
  edit or commit. Reserve thorough testing, including package build and
  `R CMD check`, for immediately before pushing.
- Reuse valid check results. Repeat checks only when subsequent changes
  or failures make the earlier results insufficient.
- Documentation-only changes generally do not require tests.
- Report what was checked and what was deferred.
