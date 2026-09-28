"""Post-selection same-tissue follow-up of two already selected leishmaniasis genes."""
import csv
import gzip
import hashlib
import io
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, t

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
from ubiomark.geo import parse_series_matrix
from ubiomark.stats import hedges_g

SOURCE='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63931'
PIN={'matrix':'21f4424b5b8197eb41f53ef394ee70b57f016ce4ff7d4ac2937b6b05f865b003',
     'platform':'3f967e8c265fe56343f5b79e21b522128d987804fae8a4245ed2e69b03b037f1'}
GENES=('PON2','CACNA2D2')

def run():
    src=ROOT/'data/geo/gse63931'
    for key,f in [('matrix',src/'GSE63931_series_matrix.txt.gz'),('platform',src/'GPL17077_table.txt')]:
        if hashlib.sha256(f.read_bytes()).hexdigest()!=PIN[key]:
            raise ValueError(f'input {key} SHA-256 changed')
    mat,ann,_=parse_series_matrix(str(src/'GSE63931_series_matrix.txt.gz'))
    if mat.shape!=(50561,16) or set(ann.Sample_platform_id)!={'GPL17077'}:
        raise ValueError('source series has unexpected probe/sample count or platform')
    series=HERE/'GSE63931.series.soft.txt'
    if hashlib.sha256(series.read_bytes()).hexdigest()!='b4df34eab61c5653605deae371d63457eeb434ed849eee95da4eb75e41841669':
        raise ValueError('series metadata changed')
    lines=[line for line in (src/'GPL17077_table.txt').open(errors='replace') if not line.startswith(('^','!','#'))]
    platform=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',dtype=str,low_memory=False)
    if len(platform)!=50739 or not {'ID','GENE_SYMBOL','CONTROL_TYPE'}<=set(platform.columns):
        raise ValueError('unexpected platform annotation')
    probe_map={}
    for gene in GENES:
        rows=platform.loc[(platform.GENE_SYMBOL==gene)&(platform.CONTROL_TYPE=='FALSE')]
        if len(rows)!=1:
            raise ValueError(f'{gene} has {len(rows)} annotated target probes; no phenotype-selected probe allowed')
        probe=rows.ID.iloc[0]
        if probe not in mat.index:raise ValueError(f'{gene} unique platform probe absent in matrix')
        probe_map[gene]=probe
    cases=[];ctrl=[];rows=[]
    for gsm,a in ann.iterrows():
        title=a.Sample_title
        if title.startswith('Cutaneous lesion') and a.Sample_source_name_ch1.startswith('Cutaneous lesion (L. braziliensis)'):
            group='lesion';cases.append(gsm)
        elif title.startswith('Normal skin') and 'non-infected' in a.Sample_source_name_ch1:
            group='healthy_skin';ctrl.append(gsm)
        else:raise ValueError(f'unknown source phenotype {gsm}: {title}')
        url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=full'
        file=HERE/'GSE63931'/f'{gsm}.soft.txt.gz'
        if not file.exists():
            with urllib.request.urlopen(url,timeout=35) as resp: raw=resp.read()
            if len(raw)<20000:raise ValueError(f'short source record {gsm}')
            with gzip.open(file,'wb',compresslevel=9) as f:f.write(raw)
        else:raw=gzip.open(file,'rb').read()
        text=raw.decode().replace('\r\n','\n')
        for line in [f'^SAMPLE = {gsm}',f'!Sample_geo_accession = {gsm}',
                     f'!Sample_title = {title}',f'!Sample_source_name_ch1 = {a.Sample_source_name_ch1}',
                     '!Sample_organism_ch1 = Homo sapiens','!Sample_platform_id = GPL17077',
                     '!Sample_series_id = GSE63931','!Sample_characteristics_ch1 = tissue: skin']:
            if line not in text.splitlines():raise ValueError(f'{gsm}: unexpected field {line}')
        start=text.index('!sample_table_begin\n')+len('!sample_table_begin\n')
        end=text.index('!sample_table_end',start)
        sample=pd.read_csv(io.StringIO(text[start:end]),sep='\t',index_col=0,low_memory=False)
        if len(sample)!=50561 or set(sample.index)!=set(mat.index):
            raise ValueError(f'{gsm} mismatched probe set')
        # Check all per-sample probe values against the pinned series matrix.
        if not np.allclose(sample.loc[mat.index,'VALUE'].to_numpy(float),
                           mat[gsm].to_numpy(float),atol=1e-5,rtol=1e-6,equal_nan=True):
            raise ValueError(f'{gsm} source values differ')
        rows.append({'gsm':gsm,'source_url':url,'source_sha256':hashlib.sha256(raw).hexdigest(),
                     'title':title,'source_group':group,'probes_matched':len(sample)})
    if len(cases)!=8 or len(ctrl)!=8 or len(rows)!=16:raise ValueError('unexpected case/control counts')
    pd.DataFrame(rows).to_csv(HERE/'GSE63931_used_sample_crosswalk.csv',index=False)
    result=[]
    for gene in GENES:
        probe=probe_map[gene]
        a=mat.loc[probe,cases].to_numpy(dtype=float)
        b=mat.loc[probe,ctrl].to_numpy(dtype=float)
        g,v=hedges_g(a[None,:],b[None,:])
        welch=ttest_ind(a,b,equal_var=False)
        delta=a.mean()-b.mean()
        s1=a.var(ddof=1)/len(a);s0=b.var(ddof=1)/len(b)
        se=np.sqrt(s1+s0);df=(s1+s0)**2/(s1*s1/(len(a)-1)+s0*s0/(len(b)-1))
        ci=[float(delta-t.ppf(.975,df)*se),float(delta+t.ppf(.975,df)*se)]
        # One-sided p for frozen negative direction; two-sided p needed for context.
        pneg=float(t.cdf(delta/se,df))
        result.append({'gene':gene,'probe':probe,'n_lesion':len(a),'n_healthy':len(b),
                       'mean_lesion':float(a.mean()),'mean_healthy':float(b.mean()),
                       'mean_difference':float(delta),'ci95_difference':ci,
                       'hedges_g':float(g[0]),'g_variance':float(v[0]),
                       'welch_two_sided_p':float(welch.pvalue),
                       'one_sided_negative_p':pneg,'bonferroni_directional_pass':bool(delta<0 and pneg<=.025)})
    rec={'source':SOURCE,'sha256':PIN,'source_records':16,'cases':8,'healthy':8,'n_fixed_measured':2,
         'genes':result,'both_directional_pass':all(x['bonferroni_directional_pass'] for x in result),
         'status':'Post-selection follow-up after P30/P35, not untouched validation. Participant independence versus Brazilian cohorts unverified; same-task published benchmark and clinical classifier absent.'}
    (HERE/'GSE63931_two_gene_result.json').write_text(json.dumps(rec,indent=2)+'\n')
    return rec

if __name__=='__main__':print(json.dumps(run(),indent=2))
