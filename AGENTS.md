# Research instructions

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
