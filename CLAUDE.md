# Working notes for Claude

Module 2 (deep site annotation) of a four-repo cardiovascular genomics platform. Sibling repos
under `../`: `sarek_pipeline` (Module 1, calling), `clinical_variant_prioritization` (Module 3,
explorer + scoring), `carve-platform` (the map).

Read `../carve-platform/STATUS.md` before proposing anything — it says where things actually are.
Plan: `../carve-platform/ROADMAP.md` (`PLT-###` ids). Rationale: `../carve-platform/decisions/`.

`BUG_TRACKER.md` and `TODO.md` here are the source of truth for this repo's internals. Don't
duplicate them elsewhere; tag entries `[PLT-###]` when they implement a roadmap task.

Rules: one repo at a time · a change needs a decision record only if another project must know
about it · cite `file:line` verified now, not recalled (tracker line numbers may be from older
revisions) · label FACT / INFERENCE / PROPOSAL · no new cross-cutting structure without asking.

Tone: direct, concise, no recapping. Push back with reasons. Never inflate scientific novelty.
