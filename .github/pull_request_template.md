## What and why

<!-- What changes, and the reason. Link the roadmap item or issue. -->

## GSLIB parity

<!-- Which GSLIB program/behaviour this touches, and whether it matches the
Fortran or diverges on purpose (if so, add it to docs/ROADMAP.md §2). -->

## Testing

<!-- New/changed tests. For numerical changes: analytic cases, properties,
or comparison against Fortran output. -->

## Checklist

- [ ] `cargo fmt --all --check`
- [ ] `cargo clippy --workspace --all-targets -- -D warnings`
- [ ] `cargo test --workspace`
- [ ] `cargo deny check` (if dependencies changed)
- [ ] Docs and `docs/ROADMAP.md` updated
