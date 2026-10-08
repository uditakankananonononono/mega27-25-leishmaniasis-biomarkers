# Record-unit derivation (2026-10-08)

Documentation note, not a license verdict and not an integrity certificate. It explains the "disease-tagged record units" count in README.md.

## What a record unit is
One row of `results/disease_tagged_accession_manifest.csv` (columns: accession, source, role). A row is one GEO accession record that this project used or fetched for a stated role: an individual GSM sample record or a GSE series record. Rows are not patients, not independent datasets and not independent cohorts. Many GSM rows are nested in a parent GSE and the `role` text says so. No per-patient or per-cell unit is counted.

## Current count
196 rows, 196 unique accessions: {'GSE series': 7, 'GSM sample': 188, 'other': 1}. Non-GSE/GSM row: MONDO_0011989.

## Why README said 179
The README count was written when the manifest had 179 rows and was not updated when later commits added rows:
manifest rows by commit: cf9019f=179, 590a665=196.
- 590a665 (2026-09-28) +17 rows: 8 GSE63931 healthy-skin and 8 lesion GSM expression-probe rows, plus 1 GSE record for the post-selection PON2/CACNA2D2 same-tissue follow-up.
Those commits changed the manifest, not the README sentence, so the README lagged by 17. The earlier 179 equals the manifest at commit cf9019f.

## Method
`git show <commit>:results/disease_tagged_accession_manifest.csv`, row count and (accession, role) set difference between commits; no row was removed between commits.

## Boundary
The recount does not change what the rows mean. The README caveat stays: these are nested record units, not independent datasets or people.
