# Vendored GSLIB 2.0 (Fortran 77)

The original GSLIB 2.0 source, used as the reference ("oracle") for parity
testing (see `docs/spec/harness.md`) and as the authority for questions about
GSLIB's behaviour (`AGENTS.md`, *No fabricated provenance*).

## Provenance

- **Archive:** `gslib77.zip`, provided by the repository owner on 2026-11-27.
- **SHA-256:** `f8a4bd8fe2b90b26f6b3f33ccf795f9d73cb2560c524128c7865b050a5f763d1`
- **Version:** GSLIB 2.0. Programs declare `VERSION=2.000`; the `README`
  reads "DOSEXEC GSLIB 2.0, Last Revised: February 14, 1997".
- Details: finding F-003 in `docs/spec/decisions.md`.

## What is here

`src/` holds every file of the archive **except `BIN/`** (1997 Windows/DOS
executables, which are never used as an oracle; decision D-005):

- one directory per program (`.FOR` source, `.INC` includes, `.PAR`
  parameter files);
- the library routines in `GSLIB/`;
- the book data sets in `DATA/`;
- the original `README`.

## Rules

- **Never modify these files.** They are byte-identical to the archive,
  including DOS CRLF line endings; `MANIFEST.sha256` lists each file's
  checksum, and `.gitattributes` stops Git from converting them. Changes
  needed for the reference build (e.g. widening output formats) are applied
  by the recipe to a copy at build time.
- To verify: `cd vendor/gslib-2.0/src && sha256sum -c ../MANIFEST.sha256`.

## Licence

GSLIB 2.0 is © 1996 The Board of Trustees of the Leland Stanford Junior
University. It is distributed with permission to copy, modify and
redistribute on condition that its notice remains intact; the full notice is
in each program file and in the repository's `NOTICE`. The `GSLIB/` library
routines and the data sets carry no notice of their own and are covered by
the same distribution. This project's own code is MIT-licensed (`LICENSE`).
