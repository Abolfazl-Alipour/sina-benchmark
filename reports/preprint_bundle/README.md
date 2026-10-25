<!-- track_scope: benchmark -->
# Preprint bundle (preview — NOT submitted)

- Source of truth: `reports/SINA_PAPER.md` (edit that file, regenerate this PDF).
- PDF: `SINA_preprint_2026-10-07.pdf` (29 pages, pandoc → HTML → weasyprint from the source at commit `f904822`+ abstract tightening; verify with `pypdf`).
- Status: **final preview for owner review only. Nothing has been uploaded anywhere.** Submission (bioRxiv/ChemRxiv) is the owner's action.
- Regenerate: `pandoc reports/SINA_PAPER.md -f markdown -t html5 --standalone -o /tmp/sina.html && weasyprint /tmp/sina.html reports/preprint_bundle/SINA_preprint_<date>.pdf` (pandoc + weasyprint installed in the `sina` env 2026-10-07).
- Before any submission: owner reads the PDF end-to-end; confirm the claim-ceiling lint is still clean (90 files, 0 violations at bundle time); T1 credential hygiene (token revocation + history purge decision) should precede making anything public.
