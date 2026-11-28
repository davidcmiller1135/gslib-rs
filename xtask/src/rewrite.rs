//! Output-format widening for the reference build (`docs/spec/harness.md`,
//! R1 step 2; `AGENTS.md`, *Reference build*).
//!
//! GSLIB writes results with fixed-point formats such as `f9.3`, which
//! truncate a double-promoted build's precision. These rules rewrite those
//! formats in a *copy* of the source, never the vendored files.
//!
//! # Rules (version [`RULES_VERSION`])
//!
//! Derived from a survey of the vendored GSLIB 2.0 source: every numeric
//! output descriptor is `fW.D` (313 in FORMAT statements, 72 in quoted WRITE
//! formats); there are no `e`/`g`/`d` output descriptors and no scale
//! factors; no READ references a FORMAT label; and the only numeric READ
//! format is an internal read in `vargplt`.
//!
//! 1. `[r]fW.D` becomes `[r]es25.16e3`, keeping any repeat count `r`. That is
//!    17 significant digits (enough to round-trip an `f64`) and a field never
//!    narrower than its content plus one blank, so columns stay separated.
//! 2. Rule 1 applies to FORMAT statements and to quoted formats in WRITE
//!    statements, **except** internal writes (the unit is a `character`
//!    variable declared in the same file), and FORMAT statements whose label
//!    such an internal write uses. Widening those could overflow the string.
//! 3. READ statements and everything else are left alone.
//!
//! Lines may grow past column 72, so the build uses
//! `-ffixed-line-length-none` (no vendored code line passes column 72).

use std::collections::BTreeSet;

/// Version of these rules, recorded in every fixture's provenance.
pub(crate) const RULES_VERSION: u32 = 1;

/// The replacement edit descriptor.
const WIDE: &str = "es25.16e3";

/// What rewriting one file did.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct Report {
    /// Descriptors widened.
    pub widened: usize,
    /// Statements skipped because they write to a character variable.
    pub internal_skipped: usize,
}

/// A fixed-form statement: its line indices and its text (columns 7+ joined).
struct Statement {
    lines: Vec<usize>,
    label: String,
    text: String,
}

fn is_comment(line: &str) -> bool {
    line.trim().is_empty() || line.starts_with(['c', 'C', '*', '!'])
}

fn is_continuation(line: &str) -> bool {
    !is_comment(line) && line.chars().nth(5).is_some_and(|c| c != ' ' && c != '0')
}

/// Columns 7 onwards (or nothing for a short line).
fn body(line: &str) -> &str {
    line.char_indices()
        .nth(6)
        .map_or("", |(i, _)| line.get(i..).unwrap_or(""))
}

fn statements(lines: &[&str]) -> Vec<Statement> {
    let mut out: Vec<Statement> = Vec::new();
    for (i, line) in lines.iter().enumerate() {
        if is_comment(line) {
            continue;
        }
        if is_continuation(line) {
            if let Some(last) = out.last_mut() {
                last.lines.push(i);
                last.text.push_str(body(line));
            }
            continue;
        }
        let label: String = line.chars().take(5).filter(char::is_ascii_digit).collect();
        out.push(Statement {
            lines: vec![i],
            label,
            text: body(line).to_owned(),
        });
    }
    out
}

/// Lower-case text with blanks removed (fixed form ignores blanks).
fn squeeze(text: &str) -> String {
    text.chars()
        .filter(|c| !c.is_whitespace())
        .collect::<String>()
        .to_ascii_lowercase()
}

/// Names declared in `character` statements.
fn character_variables(stmts: &[Statement]) -> BTreeSet<String> {
    let mut names = BTreeSet::new();
    for s in stmts {
        let t = squeeze(&s.text);
        let Some(rest) = t.strip_prefix("character") else {
            continue;
        };
        // Drop a leading length spec: `*40` or `*(40)`.
        let rest = rest.strip_prefix('*').map_or(rest, |r| {
            r.trim_start_matches(|c: char| c.is_ascii_digit() || c == '(' || c == ')')
        });
        let mut depth = 0_i32;
        let mut current = String::new();
        for c in rest.chars().chain(std::iter::once(',')) {
            match c {
                '(' => depth += 1,
                ')' => depth -= 1,
                ',' if depth == 0 => {
                    let name: String = current
                        .chars()
                        .take_while(|ch| ch.is_ascii_alphanumeric() || *ch == '_')
                        .collect();
                    if !name.is_empty() {
                        names.insert(name);
                    }
                    current.clear();
                    continue;
                }
                _ => {}
            }
            if depth == 0 {
                current.push(c);
            }
        }
    }
    names
}

/// For `write(unit, fmt, ...)`: the unit's variable name and the format
/// argument, from squeezed text.
fn write_unit_and_format(squeezed: &str) -> Option<(String, String)> {
    let args = squeezed.strip_prefix("write(")?;
    let mut depth = 0_i32;
    let mut parts = vec![String::new()];
    let mut in_quote = false;
    for c in args.chars() {
        match c {
            '\'' => in_quote = !in_quote,
            '(' if !in_quote => depth += 1,
            ')' if !in_quote && depth == 0 => break,
            ')' if !in_quote => depth -= 1,
            ',' if !in_quote && depth == 0 => {
                parts.push(String::new());
                continue;
            }
            _ => {}
        }
        if let Some(p) = parts.last_mut() {
            p.push(c);
        }
    }
    let unit: String = parts
        .first()?
        .chars()
        .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
        .collect();
    Some((unit, parts.get(1).cloned().unwrap_or_default()))
}

fn is_delimiter(c: Option<char>) -> bool {
    c.is_none_or(|c| matches!(c, '(' | ',' | '/' | ':' | ' '))
}

/// Widen `[r]fW.D` descriptors in `segment`. Returns the new text and count.
fn widen(segment: &str) -> (String, usize) {
    let chars: Vec<char> = segment.chars().collect();
    let mut out = String::with_capacity(segment.len() + 32);
    let mut count = 0;
    let mut i = 0;
    while let Some(&c) = chars.get(i) {
        if c == 'f' || c == 'F' {
            // Look back over a repeat count to a delimiter.
            let start = chars.get(..i).map_or(0, |before| {
                i - before
                    .iter()
                    .rev()
                    .take_while(|d| d.is_ascii_digit())
                    .count()
            });
            let before = start.checked_sub(1).and_then(|k| chars.get(k)).copied();
            // Look ahead for W.D.
            let w = chars
                .get(i + 1..)
                .map_or(0, |a| a.iter().take_while(|d| d.is_ascii_digit()).count());
            let dot = chars.get(i + 1 + w).copied();
            let d = chars
                .get(i + 2 + w..)
                .map_or(0, |a| a.iter().take_while(|d| d.is_ascii_digit()).count());
            let after = chars.get(i + 2 + w + d).copied();
            let ends_cleanly = after.is_none_or(|a| !a.is_ascii_alphanumeric());
            if is_delimiter(before) && w > 0 && dot == Some('.') && d > 0 && ends_cleanly {
                out.push_str(WIDE);
                count += 1;
                i += 2 + w + d;
                continue;
            }
        }
        out.push(c);
        i += 1;
    }
    (out, count)
}

/// Widen descriptors outside quotes (`outside = true`) or inside quoted
/// strings that look like formats (`outside = false`).
fn widen_line(line: &str, outside: bool) -> (String, usize) {
    let Some((head_end, _)) = line.char_indices().nth(6) else {
        return (line.to_owned(), 0);
    };
    let (head, rest) = line.split_at(head_end);
    let mut out = String::from(head);
    let mut count = 0;
    // Split on quotes: even segments are outside, odd are inside.
    for (k, seg) in rest.split('\'').enumerate() {
        if k > 0 {
            out.push('\'');
        }
        let inside = k % 2 == 1;
        let target = if outside {
            !inside
        } else {
            inside && seg.trim_start().starts_with('(')
        };
        if target {
            let (w, n) = widen(seg);
            out.push_str(&w);
            count += n;
        } else {
            out.push_str(seg);
        }
    }
    (out, count)
}

/// Apply the rules to one fixed-form source file.
pub(crate) fn widen_formats(source: &str) -> (String, Report) {
    let crlf = source.contains("\r\n");
    let text = source.replace("\r\n", "\n");
    let mut lines: Vec<String> = text.split('\n').map(str::to_owned).collect();
    let borrowed: Vec<&str> = lines.iter().map(String::as_str).collect();
    let stmts = statements(&borrowed);
    let chars = character_variables(&stmts);

    let mut report = Report::default();
    let mut internal_labels = BTreeSet::new();
    let mut internal_writes = BTreeSet::new();
    for (k, s) in stmts.iter().enumerate() {
        if let Some((unit, fmt)) = write_unit_and_format(&squeeze(&s.text))
            && chars.contains(&unit)
        {
            internal_writes.insert(k);
            report.internal_skipped += 1;
            if fmt.chars().all(|c| c.is_ascii_digit()) && !fmt.is_empty() {
                internal_labels.insert(fmt);
            }
        }
    }

    for (k, s) in stmts.iter().enumerate() {
        let t = squeeze(&s.text);
        let mode = if t.starts_with("format(") {
            (!internal_labels.contains(&s.label)).then_some(true)
        } else if t.starts_with("write(") {
            (!internal_writes.contains(&k)).then_some(false)
        } else {
            None
        };
        let Some(outside) = mode else { continue };
        for &i in &s.lines {
            if let Some(line) = lines.get_mut(i) {
                let (new, n) = widen_line(line, outside);
                *line = new;
                report.widened += n;
            }
        }
    }

    let joined = lines.join("\n");
    (
        if crlf {
            joined.replace('\n', "\r\n")
        } else {
            joined
        },
        report,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    // Fixture snippets are written for these tests; they are not GSLIB code.

    #[test]
    fn widens_format_statements_keeping_repeat_counts() {
        let src = "      write(lout,100) a,b\n  100 format(f9.3,1x,3f12.4,i5,a)\n";
        let (out, r) = widen_formats(src);
        assert_eq!(
            out,
            "      write(lout,100) a,b\n  100 format(es25.16e3,1x,3es25.16e3,i5,a)\n"
        );
        assert_eq!(r.widened, 2);
    }

    #[test]
    fn widens_quoted_write_formats_but_not_data_strings() {
        let src = "      write(lout,'(f9.3,1x,f9.3)') x,y\n      write(lout,*) 'f1.2 is text'\n";
        let (out, r) = widen_formats(src);
        assert!(out.contains("'(es25.16e3,1x,es25.16e3)'"), "{out}");
        assert!(out.contains("'f1.2 is text'"), "{out}");
        assert_eq!(r.widened, 2);
    }

    #[test]
    fn leaves_reads_alone() {
        let src = "      read(title,'(3x,7f10.3)') v\n      read(lin,*) a\n";
        let (out, r) = widen_formats(src);
        assert_eq!(out, src);
        assert_eq!(r.widened, 0);
    }

    #[test]
    fn skips_internal_writes_and_their_format_labels() {
        let src = "      character str*12, title*80\n      write(str,'(f12.2)') v\n      write(title,200) v\n  200 format(f8.3)\n      write(lout,'(f8.3)') v\n";
        let (out, r) = widen_formats(src);
        assert!(out.contains("write(str,'(f12.2)')"), "{out}");
        assert!(out.contains("200 format(f8.3)"), "{out}");
        assert!(out.contains("write(lout,'(es25.16e3)')"), "{out}");
        assert_eq!(r.internal_skipped, 2);
    }

    #[test]
    fn handles_continuation_lines_and_crlf() {
        let src = "  101 format(' X =',f8.2,\r\n     +       ' Y =',f8.2)\r\n";
        let (out, r) = widen_formats(src);
        assert_eq!(r.widened, 2);
        assert!(out.contains("\r\n"), "line endings preserved");
        assert!(out.contains("+       ' Y =',es25.16e3)"), "{out}");
        assert!(out.contains("' X ='"), "quoted text untouched");
    }

    #[test]
    fn does_not_touch_identifiers_or_comments() {
        let src = "c     format(f8.3) in a comment\n      buf12 = 1.5\n      x = af1.2\n";
        let (out, r) = widen_formats(src);
        assert_eq!(out, src);
        assert_eq!(r.widened, 0);
    }

    #[test]
    fn character_declarations_with_lengths() {
        let stmts = statements(&["      character*40 a, b(10)", "      character c*12, d"]);
        let names = character_variables(&stmts);
        assert_eq!(names.into_iter().collect::<Vec<_>>(), ["a", "b", "c", "d"]);
    }
}
