import csv
import gzip
import hashlib
from pathlib import Path
from run_gse63931 import run

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]


def test_postselection_source_followup():
    result=run()
    assert (result['source_records'],result['cases'],result['healthy'],result['n_fixed_measured'])==(16,8,8,2)
    assert [g['gene'] for g in result['genes']]==['PON2','CACNA2D2']
    assert result['both_directional_pass']
    cross=list(csv.DictReader(open(HERE/'GSE63931_used_sample_crosswalk.csv')))
    assert len({r['gsm'] for r in cross})==16
    for r in cross:
        assert hashlib.sha256(gzip.open(HERE/'GSE63931'/f'{r["gsm"]}.soft.txt.gz','rb').read()).hexdigest()==r['source_sha256']
    for name in ('dataset_manifest.csv','disease_tagged_accession_manifest.csv'):
        manifest=list(csv.DictReader(open(ROOT/'results'/name)))
        ids={r['accession'] for r in manifest}
        assert len(ids)==len(manifest)==196
        assert {'GSE63931'}|{r['gsm'] for r in cross}<=ids
