# EB-FPCA Background Learning

## Goal

Create a durable first-pass knowledge base for developing EB-FPCA on top of
`EBSmoothr`, starting from the Zotero screenshot in `EBFPCA/related-paper/`.

## Plan

- Identify the papers listed in `EBFPCA/related-paper/list-of-paper.png`.
- Query local Zotero metadata and attachment records to find stored PDFs.
- Copy available Zotero PDFs into `EBFPCA/related-paper/`.
- Record papers that do not currently have local Zotero PDF attachments.
- Extract readable text from the copied PDFs and local Zotero HTML/cache where
  needed.
- Summarize paper-level lessons and translate them into concrete implications
  for the current `EBSmoothr` package.
- Create `EBFPCA/knowledge.md` as an iteratively maintained document.

## Validation Plan

- Confirm the number of copied PDFs.
- Confirm the knowledge document includes all screenshot-listed papers.
- Avoid modifying existing package code or unrelated dirty worktree files.
