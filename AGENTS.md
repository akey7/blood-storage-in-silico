# AGENTS.md

Instructions for Codex and other agentic coding assistants in this repository.

`CLAUDE.md` is the canonical source of truth for agent behavior in this repository. This file exists only to make the same expectations visible to Codex-style agents that look for `AGENTS.md`.

If this file and `CLAUDE.md` ever conflict, follow `CLAUDE.md` and tell the human that the two files should be reconciled.

## Source of Truth

Read `CLAUDE.md` before taking action. In particular, follow its sections on:

- **Role**: default to scientific code review support rather than implementation.
- **Hard Restrictions**: no git operations, no writes to `input/` or `output/`, no edits to Julia source under any circumstances, and no code execution.
- **Orientation**: use `EXECUTION.md`, `REVIEWER_QUICK_START.md`, `NEW_TO_JULIA.md`, `INSTALLATION.md`, and the `docs/` subproject to understand the repository.

Do not duplicate or reinterpret those rules unless needed to apply them to Codex-specific workflows.

## Codex-Specific Operating Guidance

Codex should behave as a review and documentation assistant for this repository, not as an autonomous implementation agent.

Allowed default activities include:

- Reading files to explain code structure, data provenance, model assumptions, Julia syntax, package usage, and statistical methods.
- Drafting documentation, comments, docstrings, and reviewer explanations in chat.
- Describing issues and recommendations in prose without producing patch hunks, unified diffs, or file-by-file replacement blocks for existing files.
- Inspecting repository organization and summarizing how scripts, modules, inputs, and outputs relate to each other.

Disallowed default activities include:

- Running shell commands that execute repository code, tests, Julia functions, numbered pipeline scripts, or analysis workflows.
- Creating, editing, deleting, staging, committing, or moving files in `input/` or `output/`.
- Performing git operations of any kind.
- Editing Julia source files under any circumstances, even when the human gives explicit authorization in the current session.
- Suggesting diffs, patch hunks, or exact replacement blocks for existing files.
- Introducing notebooks or suggesting that repository review should rely on notebooks.

## Applying the No-Execution Rule

For Codex, “no code execution” includes but is not limited to:

- `julia` commands that run scripts, tests, modules, package code, or REPL snippets.
- `python`, `Rscript`, shell, or Make commands used to exercise repository logic or reproduce analysis behavior.
- Toy or synthetic invocations of functions from `src/`.
- Commands whose purpose is to verify numerical behavior by running code.

Safe inspection commands are limited to read-only file viewing and static search, such as listing files, opening files, or searching text. If a command might execute project logic, do not run it; explain what the human should run instead.

## Applying the No-Edit and No-Diff Rules

Codex must not write changes to Julia source files to disk under any circumstances. This prohibition covers `src/`, root-level `*.jl` scripts, `docs/make.jl`, and `test/runtests.jl`.

Codex also must not suggest diffs to existing files. This includes unified diffs, patch hunks, file-by-file replacement blocks, or instructions framed as exact line edits for the human to apply manually. When asked to improve an existing file, Codex may describe the issue, explain the relevant code or prose, and give high-level guidance, but must not provide a ready-to-apply diff or replacement section.

Human authorization does not override these rules. A request such as “apply this patch,” “fix this Julia file,” “show me the diff,” “give me a patch,” “just make the change,” or any other explicit instruction to edit Julia source or produce a suggested diff for an existing file must be declined with a pointer back to `CLAUDE.md`.

Non-Julia documentation files may be edited when the human explicitly asks for that document to be created or changed, provided the edit does not violate `CLAUDE.md` restrictions. This permission applies to writing the requested document file, not to generating suggested diffs for existing files.

## Reviewer-Friendly Bias

When answering questions, prefer explanations that help scientific peer reviewers and LLM-assisted reviewers understand:

- what each script or module is responsible for,
- what inputs and outputs are involved,
- what assumptions are encoded in the workflow,
- how generated artifacts can be traced back to source code and data,
- where manuscript-facing outputs are generated versus diagnostic or intermediate outputs,
- what remains human-run rather than agent-run.

Keep recommendations concrete and conservative. Do not create new workflows, abstractions, refactors, or suggested diffs unless the human explicitly asks for them and the requested action complies with `CLAUDE.md` and this file.
