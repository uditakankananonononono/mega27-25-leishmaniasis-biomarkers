"""P28 pre-registered cross-parasite leishmaniasis lesion sign transport."""
import csv,hashlib,json,re,sys,urllib.request,warnings
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,pandas as pd
from ubiomark import geo,stats
root=Path(__file__).resolve().parents[1];folder=root/'data/geo/p28';folder.mkdir(parents=True,exist_ok=True)
sources=[('GSE216638_series_matrix.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE216nnn/GSE216638/matrix/GSE216638_series_matrix.txt.gz','d67a275a4e9b6b3f097f43d4afed82b82ad6daf11bc9c81bd39dbd1b93c21868'),('GSE216638_resultsAll_Counts.xlsx','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE216nnn/GSE216638/suppl/GSE216638_resultsAll_Counts.xlsx','29d3d2b37392d4f5a9f4461801d8173def9a5d634e08add50433500cd3923944')]
for name,url,sha in sources:
 p=folder/name
 if not p.exists():urllib.request.urlretrieve(url,p)
 assert hashlib.sha256(p.read_bytes()).hexdigest()==sha
_,ann,_=geo.parse_series_matrix(str(folder/sources[0][0]));assert len(ann)==22 and ann.index.is_unique
case=['A','B','C','D','E','F','G','H','I','J','K','L','M','O','P','Q'];ctrl=['H1','H2','H3','H4','H5','H7'];cols=case+ctrl
mapping={r.Sample_title.split(':')[0]:gsm for gsm,r in ann.iterrows()};assert set(mapping)==set(cols) and len(mapping)==22
prior={r['accession'] for r in csv.DictReader((root/'results/dataset_manifest.csv').open())};assert not prior.intersection(ann.index)
for code in cols:
 row=ann.loc[mapping[code]];status=row.Sample_characteristics_ch1_1
 assert status.startswith('disease state: ')
 assert ('healthy' in status)==(code in ctrl)
 assert row.Sample_characteristics_ch1_2==('treatment: healthy' if code in ctrl else 'treatment: NO prior anti-leishmanial treatment')
assert sum('/UCL' in ann.loc[mapping[c],'Sample_title'] for c in case)==7
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message='Unknown extension is not supported')
 x=pd.read_excel(folder/sources[1][0],engine='openpyxl')
assert list(x.columns[:3])==['ENS','symbol','entrez'] and set(cols)==set(x.columns[3:25]) and len(x.columns)==43 and len(x)==44135
assert x.ENS.is_unique and x.ENS.astype(str).str.fullmatch(r'ENSG\d+').all()
values=x[cols].to_numpy(dtype=float);assert np.isfinite(values).all() and (values>=0).all() and np.equal(values,np.floor(values)).all()
notnull=x.symbol.notna();valid=x.symbol.astype(str).str.fullmatch(r'[A-Za-z][A-Za-z0-9_.-]*');notnull=notnull&valid # exclude the slash-delimited ambiguous THRA1/BTR row
known=x.loc[notnull,'symbol'].astype(str);dup=set(known[known.duplicated(keep=False)])
expr=x.loc[notnull&~x.symbol.isin(dup),['symbol']+cols].set_index('symbol').astype(float)
assert expr.index.is_unique
lib=values.sum(axis=0);assert (lib>0).all();cpm=expr.div(pd.Series(lib,index=cols),axis=1)*1e6
filtered=cpm.loc[(cpm>1).mean(axis=1)>=.2];log=np.log2(filtered+1)
g,v=stats.hedges_g(log[case].to_numpy(),log[ctrl].to_numpy());effect=pd.DataFrame({'g':g,'v':v},index=log.index).replace([np.inf,-np.inf],np.nan).dropna()
D=pd.read_csv(root/'results/meta_discovery/leishmaniasis.csv.gz',index_col=0);C=pd.read_csv(root/'results/leish_p28_fixed_panel.csv').set_index('gene');assert len(C)==20 and (C.k>=2).all()
obs=C.index.intersection(effect.index);cp=int((C.loc[obs,'mu']>0).sum());cn=len(obs)-cp
agreement=np.sign(C.loc[obs,'mu'])==np.sign(effect.loc[obs,'g'])
pool=D.index.intersection(effect.index);pool=pool[(D.loc[pool,'k']>=2)&~pool.isin(C.index)];up=np.array(pool[D.loc[pool,'mu']>0]);down=np.array(pool[D.loc[pool,'mu']<0]);assert len(up)>=cp and len(down)>=cn
rng=np.random.default_rng(20260925);null=np.zeros(10000,dtype=int)
for i in range(10000):
 picks=list(rng.choice(up,cp,replace=False))+list(rng.choice(down,cn,replace=False))
 null[i]=int((np.sign(D.loc[picks,'mu'])==np.sign(effect.loc[picks,'g'])).sum())
from scipy.stats import ttest_ind
rawp=ttest_ind(log.loc[obs,case].to_numpy(),log.loc[obs,ctrl].to_numpy(),axis=1,equal_var=False).pvalue
adjusted=stats.bh_fdr(rawp)
with (root/'results/leish_p28_genes.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['gene','discovery_mu','cohort_g','cohort_v','agrees','welch_p','panel_bh_q','status'],lineterminator='\n');w.writeheader()
 for gene in C.index:
  if gene in obs:
   i=obs.get_loc(gene);w.writerow(dict(gene=gene,discovery_mu=float(C.loc[gene,'mu']),cohort_g=float(effect.loc[gene,'g']),cohort_v=float(effect.loc[gene,'v']),agrees=int(agreement.loc[gene]),welch_p=float(rawp[i]),panel_bh_q=float(adjusted[i]),status='measured'))
  else:w.writerow(dict(gene=gene,discovery_mu=float(C.loc[gene,'mu']),cohort_g='',cohort_v='',agrees='',welch_p='',panel_bh_q='',status='not uniquely mapped or filtered'))
with (root/'results/leish_p28_samples.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['gsm','column','title','disease','treatment'],lineterminator='\n');w.writeheader()
 for code in cols:
  gsm=mapping[code];r=ann.loc[gsm];w.writerow(dict(gsm=gsm,column=code,title=r.Sample_title,disease=r.Sample_characteristics_ch1_1,treatment=r.Sample_characteristics_ch1_2))
pd.DataFrame({'matched_random_agreements':null}).to_csv(root/'results/leish_p28_null.csv.gz',index=False)
count=int(agreement.sum());p=float((1+(null>=count).sum())/10001)
rec=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216638',sha256={n:s for n,_,s in sources},n_lesion=16,n_UCL=7,n_NUCL=9,n_healthy=6,n_unique_mapped_genes=len(expr),n_effect_genes=len(effect),duplicate_symbol_count=len(dup),n_fixed_measured=len(obs),missing_fixed=list(C.index.difference(obs)),n_positive=cp,n_negative=cn,n_sign_match=count,null_mean=float(null.mean()),empirical_p=p,registered_descriptive_support=bool(len(obs)>=15 and count>=15 and p<.025),panel_bh_hits=int((adjusted<.05).sum()))
(root/'results/leish_p28_result.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
