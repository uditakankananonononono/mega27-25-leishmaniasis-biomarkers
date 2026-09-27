"""Post-audit sensitivity, not a new untouched validation or a replacement of the frozen run.

GSE125993 is the union of GSE125991/125992. Use its expression once. Exclude HighAb
samples whose source phenotype does not assert active visceral leishmaniasis.
"""
import sys,json
import numpy as np,pandas as pd
sys.path.insert(0,'src')
from ubiomark import geo,stats
from scipy.stats import norm
G='GSE125993';expr,ann,_=geo.parse_series_matrix(geo.download_matrices(G)[0])
lab=pd.read_csv(f'results/series/leishmaniasis__{G}.labels.csv').set_index('gsm').label
state=ann.Sample_characteristics_ch1_0.str.lower()
active=state.eq('disease state: active visceral leishmaniasis case')
healthy=state.eq('disease state: endemic healthy control (serology and quantiferon negative)')
highab=state.eq('disease state: high anti-leishmanial antibody levels by direct agglutination test')
assert (active.sum(),healthy.sum(),highab.sum())==(21,16,8)
assert set(lab[lab=='case'].index)==set(ann.index[active|highab])
assert set(lab[lab=='control'].index)==set(ann.index[healthy])
assert len(set(ann.index))==84 and set(expr.columns)==set(ann.index)
gene=geo.to_gene_level(expr,geo.probe_to_symbol('GPL10558'))
g,v=stats.hedges_g(gene.loc[:,active].to_numpy(float),gene.loc[:,healthy].to_numpy(float))
effect=pd.DataFrame({'g':g,'v':v},index=gene.index).replace([np.inf,-np.inf],np.nan).dropna()
effect.index.name='gene';effect.to_csv('results/leish_GSE125993_strict_active_effects.csv.gz')
disc=pd.read_csv('results/meta_discovery/leishmaniasis.csv.gz',index_col=0);disc=disc[np.isfinite(disc.p)]
top=disc.sort_values('p').head(50)
common=disc.index.intersection(effect.index);selected=top.index.intersection(common)
sg=np.sign(disc.loc[selected,'mu'].to_numpy());z=effect.loc[selected,'g'].to_numpy()/np.sqrt(effect.loc[selected,'v'].to_numpy())
p=2*norm.sf(np.abs(z));agree=np.sign(z)==sg
allz=effect.loc[common,'g'].to_numpy()/np.sqrt(effect.loc[common,'v'].to_numpy());allp=2*norm.sf(np.abs(allz));alld=np.sign(disc.loc[common,'mu'].to_numpy())
rng=np.random.default_rng(20260925);B=10000;n=len(selected);null_sign=np.empty(B);null_hits=np.empty(B)
for i in range(B):
 idx=rng.choice(len(common),n,replace=False);a=np.sign(allz[idx])==alld[idx];null_sign[i]=a.mean();null_hits[i]=np.sum(a&(allp[idx]<.05))
result={'audit_type':'post-outcome phenotype sensitivity, not preregistered confirmation','source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125993','n_source_GSM':84,'n_active':int(active.sum()),'n_endemic_healthy':int(healthy.sum()),'n_excluded_highab_previously_case':int(highab.sum()),'n_top50_measured':n,'sign_agreement':int(agree.sum()),'sign_fraction':float(agree.mean()),'null_sign_mean':float(null_sign.mean()),'empirical_p_sign':float((1+(null_sign>=agree.mean()).sum())/(B+1)),'direction_p05_hits':int(np.sum(agree&(p<.05))),'null_hits_mean':float(null_hits.mean()),'empirical_p_hits':float((1+(null_hits>=np.sum(agree&(p<.05))).sum())/(B+1)),'old_result':'results/replication_results.json; 29 case labels include eight HighAb; do not erase'}
with open('results/leish_GSE125993_strict_active_sensitivity.json','w') as f:json.dump(result,f,indent=2)
print(json.dumps(result,indent=2))
