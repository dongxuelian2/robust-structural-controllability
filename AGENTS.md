# Research instructions

## Zero-defensive-prose rule — mandatory

Apply this rule to all repository-authored English prose: README, mathematical research notes, `PROOF.md`, LaTeX exposition, documentation, pull-request descriptions, and code comments. Enforce it when auditing historical material and when writing anything new.

**Prohibited rhetorical constructions (including grammatical variants):**
- Chinese examples to recognize and eliminate in translated source material: `不是……`, `并不是……`, `不构成……`, `还没有做……`, `尚未……`, `而不是……`, `并未……`, `不意味着……`, `不代表……`, `不能据此……`, `这并不证明……`, `我们不声称……`, `不保证……`, `仅仅是……`.
- English equivalents to eliminate from authored explanations: `not a ...`, `is not ...`, `does not ...`, `do not ...`, `has not ...`, `have not ...`, `not yet ...`, `rather than ...`, `instead of ...`, `does not constitute ...`, `does not establish ...`, `does not imply ...`, `should not be interpreted as ...`, `we do not claim ...`, `no claim is made ...`, `merely ...`, `only a preliminary ...`, `beyond the scope ...`, `not intended to ...`.
- Delete disclaimers, self-exculpation, premature qualification, repetitive status caveats, and defensive contrasts. Do not relocate them to a different paragraph or replace them with synonyms serving the same rhetorical purpose. Zero such constructions are acceptable in explanatory prose.

**Required replacement style:** Lead with the actual theorem, deliverable, tested behavior, evidence, applicability conditions, exact scope, and concrete open mathematical question, using affirmative statements. Examples:
- `We prove the classification for n >= N under hypotheses H1 and H2.`
- `Verified cases: n = 2,...,7. General-case proof target: ...`
- `CI checks: 16 passed. Proof: Proposition 4. Research priority: external review.`
- `Algorithm: FPT for disjoint ordinary cardinality quotas. Next target: robust mixed-status source-SCC accessibility.`

**Mathematical truth takes precedence over stylistic cleanup.** Preserve essential logical negation in theorems, inequalities, counterexamples, constraints, source quotations, and executable code where changing it would alter correctness. In ordinary explanation, express boundaries by positive quantified scope and precise next objectives. If a claim is unsupported, correct or remove it; never invent evidence to avoid a negative sentence.

Before opening or updating a PR, review every added or modified paragraph for prohibited rhetoric and its semantic equivalents. Record the direct result and evidence; keep mathematical limitations in compact, exact statements of domain and open problems.

## Mandatory English and proof-quality policy

- English is the sole language for repository-authored prose: all Markdown documentation, proof notes, README, AGENTS, manuscript exposition, source comments and docstrings, and PR descriptions. Translate existing Chinese prose throughout tracked research files, including nested directories and archives. Keep original titles/quotations and notation only where necessary, with English context; do not create parallel Chinese documentation.
- Audit the **entire existing research corpus**, not only changed files. Check every mathematical claim against its exact definitions, quantifiers, hypotheses, boundary cases, proof dependencies, citations, computational evidence, and the stated scope of conclusions. Repair gaps with complete arguments, correct mistaken claims, and add essential missing definitions, lemmas, references, or counterexamples. Do not inflate conclusions or silently turn an incomplete argument into a theorem.
- A translation must preserve equations, proof structure, assumptions, labels, bibliographic facts, and mathematical strength. If a proof is incomplete or ambiguous, resolve it mathematically where possible; otherwise identify the precise unresolved step in clear English and keep the claim at its justified status.
- Submit focused, reviewable PRs per repository. Give exact file locations, substantive corrections, and validation steps; do not mix unrelated research with translation. Preserve cross-references and all functional code, run existing relevant checks, and manually inspect for residual Chinese prose. A character scan alone is insufficient.
- Apply these standards to every future contribution. Write concise, direct English; avoid ceremonial documentation, defensive boilerplate, and unsupported novelty or correctness claims.

This is a public, mathematics-first research repository.

- Work directly on the stated problem. Verify the exact model and literature before claiming novelty.
- Prioritize theorems, proofs, constructive algorithms, counterexamples, and precise complexity bounds.
- Use short, readable notes and source files. Add exact-arithmetic scripts only when they materially verify a claim.
- Keep experimental findings distinct from proven statements. Record counterexamples with enough data to reproduce them.
- Treat intermediate rigorous results as valuable independent outputs.
- Avoid CI scaffolding, hashes, verbose process reports, and unnecessary project-management files.
- Use portable Python and Linux-compatible shell commands when needed.
- State results directly and accurately. Do not claim a new bound until established.
