---
name: commit
description: Run the complete gs-waveform-measurements quality gate, review and stage the diff deliberately, update the dated journal, and create a Conventional Commit. Use whenever committing changes in this repository.
---

# Commit gs-waveform-measurements changes

1. Read `AGENTS.md`, the newest journal entry, and `docs/gotchas/INDEX.md`.
2. Inspect `git status`, `git diff`, and any already-staged diff. Check for secrets, caches,
   debug output, and unrelated user changes.
3. Run the standard gate:

   ```bash
   cargo fmt --all -- --check
   cargo clippy --all-targets -- -D warnings
   cargo test --all-targets
   cargo doc --no-deps --lib
   ```

4. Run `pre-commit run --all-files` (if installed). Fix the cause of every failure.
5. Obtain the real time with `date '+%Y-%m-%d %H:%M'` and add a newest-first entry to
   `docs/journal/YYYY-MM.md` describing the request, findings, changes, verification, and
   remaining work. Include it in the same changeset so a journal-only follow-up does not create
   an infinite bookkeeping loop.
6. Review the final diff, then stage explicit paths.
7. Commit with `type(scope): imperative summary`. Allowed types are `feat`, `fix`, `docs`,
   `style`, `refactor`, `perf`, `test`, `build`, `ci`, `chore`, and `revert`.
8. Confirm `git status --short` is empty. Never use `--no-verify` or `SKIP=`.
