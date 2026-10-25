# CHANGES entry template (spec §0.3.8)

Copy this file to `NNNN_short_slug.md` (next free number). Every deviation from
`SINA_Implementation_Plan_v4_FINAL.md` needs an entry in the same commit as — or
before — the code/doc change that embodies it.

```markdown
# NNNN — <short title>

- **change:** What diverges from the spec (quote the spec section it touches).
- **rationale:** Why the deviation is necessary or better.
- **requester:** Who asked for it (owner / agent / reviewer).
- **impact:** Consequences for downstream phases, gates, budgets, claims.
- **date:** YYYY-MM-DD
- **approved:** owner — YYYY-MM-DD | pending
```

Entries are numbered sequentially, never renumbered, never deleted. A rejected or
superseded change keeps its entry and gains a `superseded-by:` line pointing at the
newer entry.
