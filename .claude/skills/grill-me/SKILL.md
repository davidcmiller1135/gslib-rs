---
name: grill-me
description: Interview the user one question at a time to turn a gslib-rs topic (harness, missing values, a GSLIB program's parity) into recorded decisions, a testable spec, and an implementation plan. Use when the user says "grill me" or asks to spec or plan something in this repository.
---

# Grill me (with docs): gslib-rs

Turn a topic into decisions, a spec, and a plan by interviewing the user. The
user decides; you research, recommend, and write things down.

## Before asking anything

1. Read `AGENTS.md`, `docs/ROADMAP.md` (especially §2 divergences and §6
   provenance ledger), and the existing `docs/spec/` files for the topic.
2. Read the relevant code and tests. **Never ask what the repository already
   answers.** If the code or a recorded decision settles it, state that and
   move on.
3. List the open decisions for the topic, then order them so each comes after
   the decisions it depends on.

## The interview loop

- Ask **one question at a time.** For each:
  - give the context in two or three sentences;
  - give 2–4 concrete options with their trade-offs;
  - **recommend one** and say why, citing the `AGENTS.md` rule it follows
    where one applies.
- Wait for the answer. Do not batch questions or assume answers.
- If an answer opens a new branch, add it to the open list and say so.
- If the user defers, record the question as **open**, with your
  recommendation, and continue with whatever does not depend on it.

## GSLIB facts are evidence, not memory

A question that turns on what Fortran GSLIB actually does is answered from the
source, its documentation, or a reference run, never from recollection.

- If the source is available locally, read it and cite the file and routine.
- If it is not, record the decision as **blocked on `TODO(parity-verify)`**,
  add the claim to `docs/ROADMAP.md` §6, and move on.
- Never invent citations, outputs or tolerances (`AGENTS.md`, *No fabricated
  provenance*).

## Write as you go

As soon as a decision is made, before asking the next question:

- **Decision log:** add an entry to `docs/spec/decisions.md`:

  ```markdown
  ## D-NNN: <short title>
  - Status: decided | open | blocked (parity-verify) | superseded by D-MMM
  - Topic: <spec file>
  - Question: <what was asked>
  - Options: <options considered, one line each>
  - Decision: <what was decided>
  - Rationale: <why; cite AGENTS.md rules and evidence>
  - Consequences: <what this forces or rules out>
  ```

- **Spec:** update the topic's file, `docs/spec/<topic>.md`, so it reflects
  every decision so far. Write requirements as testable acceptance criteria
  ("Given…, when…, then…" or a precise rule with a tolerance and its source).
- Tell the user, in one line, what was recorded and where.

## Finishing a topic

When no open questions remain (or the user stops):

1. Re-read the spec for contradictions with `AGENTS.md`, other specs, and the
   decision log; raise any you find as questions.
2. Write `## Implementation plan` at the end of the spec:
   - ordered steps, each small enough for one reviewable PR;
   - for each step, the tests that prove it (named, with their oracle or
     invariant);
   - dependencies between steps, and anything blocked on `parity-verify`.
3. Summarise for the user: decisions made, decisions still open or blocked,
   and the first step of the plan.
4. Commit the docs on the session's working branch and open (or update) a PR.
   Do not merge without the user's go-ahead.
