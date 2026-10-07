# Source provenance and rights documentation gaps (leishmaniasis child, 2026-10-08 audit record)

This file records what this repository's own documents say about source hashes and terms. It is a documentation record, not a license verdict, not an integrity certificate, and not a clearance of any recorded rights hold. A matching claim below means two recorded hash strings are equal at a named field in pinned repository documents. It does not mean any assay bytes were re-downloaded or verified. No rights record was found for these sources, so reuse and redistribution are NOT cleared by this file.

Audited child commit: `7e9ce0f479b75ede8dd5a6ec594a105c6437a69d`  
Compared shared-parent repo: `mega27-25-biomarkers-underserved-diseases` at `36ff87995f4fe0e1d08f9e2b2cd885d341b88937`

Rights: No source-specific rights check found in inspected disease-child documentation. Hash matches are not legal clearance or byte verification.

## 1. Result-file hash records in this child (6)

Hash values are shown as 12-hex display prefixes only; the full value is in the cited file and field. Classes are kept separate on purpose: an expression/assay source hash, a GEO series-matrix hash (metadata that may include expression), a reference annotation, and a published artifact are different kinds of record.

- `results/leish_p28_result.json` `/sha256/GSE216638_series_matrix.txt.gz` hash `d67a275a4e9b...` (GSE216638); class: GEO_series_matrix_file_includes_metadata_and_may_include_expression; equal to shared-parent field
- `results/leish_p28_result.json` `/sha256/GSE216638_resultsAll_Counts.xlsx` hash `29d3d2b37392...` (GSE216638); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/leish_p35_result.json` `/sha256/normalized` hash `ac826a7f4b4f...` (GSE127831); class: source_hash_record_target_needs_context; equal to shared-parent field
- `results/leish_p35_result.json` `/sha256/design` hash `1544d482c0e7...` (GSE127831); class: sample_or_reference_metadata_not_assay; equal to shared-parent field
- `results/leish_p30_result.json` `/sha256/GSE214397_series_matrix.txt.gz` hash `7d12af7c4465...` (GSE214397); class: GEO_series_matrix_file_includes_metadata_and_may_include_expression; equal to shared-parent field
- `results/leish_p30_result.json` `/sha256/GSE214397_processed_log2cpm.txt.gz` hash `c065bd30189c...` (GSE214397); class: source_hash_record_target_needs_context; equal to shared-parent field

## 4. Metadata crosswalk tables (not assay bytes)

These per-sample tables carry hash columns describing GEO source-response metadata. They are not blanket expression-matrix integrity.

- `projects/leishmaniasis/sources/GSE55664_used_sample_crosswalk.csv`: columns ['source_sha256'], 35 records; file identical to a shared-parent file: True; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity
- `projects/leishmaniasis/sources/GSE80008_used_sample_crosswalk.csv`: columns ['source_sha256'], 30 records; file identical to a shared-parent file: True; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity
- `projects/leishmaniasis/sources/independence/GSE63931_used_sample_crosswalk.csv`: columns ['source_sha256'], 16 records; file identical to a shared-parent file: False; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity

## Coverage boundary

- Source accessions listed in the child manifest: 7.
- Of those, 4 have no mapped result-file payload hash in this audit: GSE125993, GSE55664, GSE63931, GSE80008.
- 0 hash fields inherited from other diseases' records are not counted toward this child.
- Records absent from child manifest can still have result-file hashes. Neither presence nor absence proves assay acquisition by this scout. Other-disease copied hashes never count toward this child.
