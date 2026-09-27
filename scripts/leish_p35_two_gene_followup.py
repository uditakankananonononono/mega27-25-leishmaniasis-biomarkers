"""P35 frozen post-selection lesion two-gene recurrence in GSE127831."""
import csv,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats as ss
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ubiomark import stats
root=Path(__file__).resolve().parents[1];folder=root/'data/geo/p35'
sha={'normalized':'ac826a7f4b4f472285d3644b43aa6e9ce3212c1b70ae585fecb39bb7d34b19d8','design':'1544d482c0e7f4d53a4275a169f6e01c6f217d791e233d3eb7addc86f96deab6'}
files={'normalized':folder/'GSE127831_Amorim_GEO_normalized.txt.gz','design':folder/'GSE127831_studydesign.csv.gz'}
for key,p in files.items():assert hashlib.sha256(p.read_bytes()).hexdigest()==sha[key]
r=list(csv.DictReader((root/'results/leish_p35_samples.csv').open()));assert len(r)==28 and len({z['gsm'] for z in r})==28 and len({z['column'] for z in r})==28
manifest=list(csv.DictReader((root/'results/dataset_manifest.csv').open()));assert all(sum(m['accession']==z['gsm'] for m in manifest)==1 for z in r)
x=pd.read_csv(files['normalized'],sep='\t',index_col=0);d=pd.read_csv(files['design'],encoding='utf-8-sig').set_index('sample');assert set(x.columns)==set(d.index)=={z['column'] for z in r} and len(d)==28
for z in r:
 key=z['column'];assert d.loc[key,'disease']==z['diagnosis'] and str(d.loc[key,'Treat. outcome'])==z['outcome'].split(': ')[-1]
 assert (key.startswith('CL'))==(z['diagnosis']=='cutaneous')
assert np.isfinite(x.to_numpy(dtype=float)).all() and not x.index.isna().any()
case=[z['column'] for z in r if z['diagnosis']=='cutaneous'];ctrl=[z['column'] for z in r if z['diagnosis']=='control'];assert len(case)==21 and len(ctrl)==7
result=[]
for gene in ('PON2','CACNA2D2'):
 assert sum(x.index==gene)==1
 a=x.loc[gene,case].to_numpy(dtype=float);b=x.loc[gene,ctrl].to_numpy(dtype=float)
 g,v=stats.hedges_g(a.reshape(1,-1),b.reshape(1,-1));welch=ss.ttest_ind(a,b,equal_var=False);diff=float(np.mean(a)-np.mean(b));se=np.sqrt(np.var(a,ddof=1)/len(a)+np.var(b,ddof=1)/len(b));df=(np.var(a,ddof=1)/len(a)+np.var(b,ddof=1)/len(b))**2/((np.var(a,ddof=1)/len(a))**2/(len(a)-1)+(np.var(b,ddof=1)/len(b))**2/(len(b)-1));ci=ss.t.interval(.95,df,loc=diff,scale=se)
 result.append(dict(gene=gene,n_lesion=len(a),n_healthy=len(b),lesion_mean=float(np.mean(a)),healthy_mean=float(np.mean(b)),normalized_difference=diff,welch_ci95=list(map(float,ci)),hedges_g=float(g[0]),g_variance=float(v[0]),welch_p=float(welch.pvalue),down=bool(diff<0),bonferroni_pass=bool(diff<0 and welch.pvalue<.025)))
rec=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127831',sha256=sha,n_lesion=21,n_healthy=7,n_symbols=len(x),duplicated_symbols=int(x.index.duplicated(keep=False).sum()),treatment_outcomes={k:int(sum(z['outcome']==f'treatment_outcome: {k}' for z in r)) for k in ['failure','cure']},genes=result,registered_descriptive_support=bool(all(z['bonferroni_pass'] for z in result)),selection_status='Post-P28 and P30 selected; recurrence only, overlap with other Bahia study patients unverified')
(root/'results/leish_p35_result.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
