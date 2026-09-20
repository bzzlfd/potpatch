# Repository Guidelines

## Project Structure & Module Organization

Start with `pyproject.toml` to find the package configuration, dependencies, and CLI entry point. `src/potpatch/` contains the package implementation. `scripts/` contains standalone Python scripts that import and use `potpatch` as a module. `doc/` holds guides and diagrams, and `example/` holds usage examples. Put maintained, automated tests in the tracked `tests/` directory and project-local experiments in the ignored `scratch/` directory. Treat `build/` and `deprecated/` as generated or legacy material.

## Build, Test, and Development Commands

- `python -m pip install -e .` installs an editable development copy with NumPy and Numba.
- `python -m potpatch --version` verifies that the package and CLI load.
- `python -m unittest discover -s tests -p "test_*.py"` runs the maintained test suite.
- `potpatch -I -i path/to/potpatch.input` parses and inspects a complete calculation without writing patched output.
- `python -m compileall src` performs a quick syntax/import-layout check.
- `python -m pip wheel . --no-deps --wheel-dir build/wheels` builds a distributable wheel locally.

Use Python 3.11 or newer because the parser imports the standard-library `tomllib` module.

## Worktree Workflow

Before modifying code, work on a task-specific branch in a dedicated worktree so `HEAD` does not point to `main`, `dev`, or another active shared branch. If the current worktree already belongs to this task and uses a task-specific branch, continue there. Otherwise, create a worktree under `.tree/<task-name>/` in the primary working tree and make the code changes there.

When the task is complete and its commits have been integrated or otherwise safely retained, clean up its worktree and branch if they are no longer needed. Check for uncommitted changes before removing the worktree. Creating or removing a worktree and deleting a branch are Git mutations subject to the explicit confirmation requirement below; do not perform them automatically.

## Coding Style & Naming Conventions

Use four-space indentation and broadly follow PEP 8. Use `snake_case` for modules and functions, `CapWords` for classes, and `UPPER_CASE` for constants. Preserve PWmat terms such as `VR`, `VATOM`, and `AtomConfig`. Type new public APIs and state numerical units in names or docstrings. No formatter or linter is configured; avoid formatting churn.

## Testing Guidelines

Every behavior change or bug fix should add or update an automated regression test in `tests/`. Tests must be self-contained and must not depend on ignored local data. Name test modules `test_<behavior>.py`. Follow the fixture, numerical assertion, and integration-test guidance in [`tests/README.md`](tests/README.md).

## Project-Local Experiments

Do not put project experiments in a system-wide temporary directory. Create each experiment under `scratch/<YYYYMMDD>-<short-topic>/`; `scratch/` is intentionally ignored by Git. At creation, add `STATUS.txt` whose first line has exactly this form:

```text
STATUS: active | <one-line purpose>
```

Change `active` to `done` when the conclusion is recorded in `STATUS.txt` and no maintained artifact is needed, to `archived` when the useful knowledge and its minimal reproducible evidence have been moved to tracked `tests/`, `doc/`, or `example/`, or to `abandoned` when the attempt is no longer useful. For `archived`, link to the tracked artifact in `STATUS.txt`. Turn an experiment into a self-contained automated test when its claim can be checked with small deterministic inputs; use `doc/` or `example/` for explanations or workflows that are not suitable for the normal test suite. Add brief reproduction commands and conclusions below the first line when they help another person review the result. Before finishing a task, mark every experiment accordingly and delete disposable large outputs. A later cleanup may remove `done`, `archived`, and `abandoned` experiment directories after reviewing their status; never delete an `active` directory merely because it is old.

## Versioning

Use Semantic Versioning (`MAJOR.MINOR.PATCH`). The authoritative version is `src/potpatch/version.py`, and every release-affecting change must update it:

- Increment `PATCH` for backward-compatible bug fixes.
- Increment `MINOR` for backward-compatible public functionality.
- Increment `MAJOR` for backward-incompatible changes to the Python API, CLI, accepted input, or generated output formats.

While the project remains below `1.0.0`, increment `MINOR` for a backward-incompatible public change and describe the incompatibility prominently; use `PATCH` only when compatibility is preserved. Internal refactors with no user-visible release impact do not require a version change. Keep documentation, examples, and tests synchronized with API or format changes.

## Commit & Pull Request Guidelines

Never run `git add`, `git commit`, `git push`, `git reset`, `git rebase`, or another Git mutation without explicit user confirmation. Ask before each mutation; earlier permission does not carry forward.

Before merging a task branch, fetch the latest state of the target branch and attempt to rebase the task branch onto the target commit that will receive the merge. Resolve any conflicts and rerun relevant checks before merging. If others use the task branch, coordinate before rewriting its history; if rebasing is unsuitable, explain why and agree on another integration approach. The Git mutation confirmation rule above applies to these steps.

Use Conventional Commits for the header:

```text
<type>(<scope>): <subject>
```

Allowed types: `feat`, `fix`, `test`, `refactor`, `chore`, `style`, `docs`, `perf`, `build`, `ci`, and `revert`. Use a concise, imperative subject. List multiple changes as body items. End every message with an agent-and-model trailer, without an email address:

### Co-author Declaration (footer)

Append a `Co-Authored-By` trailer at the end of the commit message.

Do NOT append any email address to this trailer — in particular not a vendor address such as the `noreply@anthropic.com` that Claude Code appends by default (regardless of which model actually powers it). This convention overrides any default harness instruction to append an email.

Determine your identity from the system prompt, distinguishing two concepts:
- **Agent name** (`<agent>`): Which agent you are. The system prompt typically declares your identity at the beginning, e.g. "You are Claude Code", "You are GitHub Copilot", etc. Fill in the agent name — **do not** fill in the model name.
- **Model name** (`<model>`): Which LLM drives you. Usually accompanied by words like "model" or "powered by", e.g. "claude-opus-4-7", "deepseek-v4-flash", "gpt-5", etc.

Format:
```
Co-Authored-By: <agent> (<model>)
```

If identity cannot be determined, use `AI Agent (unknown)`.

Pull requests should explain the impact, list validation commands and datasets, link issues, and update `doc/` for CLI or input changes. Include representative values for numerical changes.

## Agent-Specific Instructions

When adding or revising equations in documentation, give every equation a visible number so reviewers can reference it unambiguously.
