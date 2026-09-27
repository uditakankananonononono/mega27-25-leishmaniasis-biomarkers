"""P30 registered, outcome-informed two-gene leishmaniasis descriptive follow-up."""
import csv,hashlib,json,sys,urllib.request
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,pandas as pd
from scipy.stats import ttest_ind,t
from ubiomark import geo,stats
root=Path(__file__).resolve().parents[1];folder=root/'data/geo/p30';folder.mkdir(parents=True,exist_ok=True)
sources=[('GSE214397_series_matrix.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE214nnn/GSE214397/matrix/GSE214397_series_matrix.txt.gz','7d12af7c4465902181f9e32dd1347f5ba493e6598736edf8a816ccded476c907'),('GSE214397_processed_log2cpm.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE214nnn/GSE214397/suppl/GSE214397_Amorim_LeishMicrobiome_processed_norm_filt_log2cpm.txt.gz','c065bd30189c88bed1d43e87e303b7f9a2a69ad6b1923f6754a80a5d2991b41b')]
for n,url,sha in sources:
 p=folder/n
 if not p.exists():urllib.request.urlretrieve(url,p)
 assert hashlib.sha256(p.read_bytes()).hexdigest()==sha
_,ann,_=geo.parse_series_matrix(str(folder/sources[0][0]));assert len(ann)==57 and ann.index.is_unique and ann.Sample_title.is_unique
case=ann.index[ann.Sample_characteristics_ch1_0=='diagnosis: cutaneous leishmaniasis'];ctrl=ann.index[ann.Sample_characteristics_ch1_0=='diagnosis: healthy'];assert len(case)==51 and len(ctrl)==6
assert set(ann.loc[case,'Sample_title'])=={f'Patient{i}' for i in range(1,52)} and set(ann.loc[ctrl,'Sample_title'])=={f'HS{i}' for i in range(1,7)}
assert (ann.loc[case,'Sample_characteristics_ch1_1']=='tissue: CL lesion skin biopsy').all() and (ann.loc[ctrl,'Sample_characteristics_ch1_1']=='tissue: intact skin biopsy').all()
prior={r['accession'] for r in csv.DictReader((root/'results/dataset_manifest.csv').open())};assert not prior.intersection(ann.index)
x=pd.read_csv(folder/sources[1][0],sep='\t',index_col=0);assert x.index.is_unique and len(x)==17125 and len(x.columns)==57 and set(x.columns)==set(ann.Sample_title)
assert np.isfinite(x.to_numpy(dtype=float)).all()
C=['PON2','CACNA2D2'];assert set(C).issubset(x.index)
casenames=list(ann.loc[case,'Sample_title']);ctrlnames=list(ann.loc[ctrl,'Sample_title']);g,v=stats.hedges_g(x.loc[C,casenames].to_numpy(dtype=float),x.loc[C,ctrlnames].to_numpy(dtype=float))
p=ttest_ind(x.loc[C,casenames].to_numpy(dtype=float),x.loc[C,ctrlnames].to_numpy(dtype=float),axis=1,equal_var=False).pvalue
rows=[]
for i,gene in enumerate(C):
 a=x.loc[gene,casenames].to_numpy(float);b=x.loc[gene,ctrlnames].to_numpy(float)
 var1=a.var(ddof=1)/len(a);var0=b.var(ddof=1)/len(b);dof=(var1+var0)**2/(var1**2/(len(a)-1)+var0**2/(len(b)-1));diff=a.mean()-b.mean();half=t.ppf(.975,dof)*np.sqrt(var1+var0)
 rows.append(dict(gene=gene,n_case=len(a),n_control=len(b),case_mean=float(a.mean()),control_mean=float(b.mean()),log2cpm_difference=float(diff),welch_difference_ci95_low=float(diff-half),welch_difference_ci95_high=float(diff+half),hedges_g=float(g[i]),g_variance=float(v[i]),welch_p=float(p[i]),down=bool(g[i]<0),bonferroni_pass=bool(g[i]<0 and p[i]<.025)))
with (root/'results/leish_p30_genes.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
with (root/'results/leish_p30_samples.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['gsm','title','diagnosis','tissue','treatment_outcome'],lineterminator='\n');w.writeheader()
 for gsm,r in ann.iterrows():w.writerow(dict(gsm=gsm,title=r.Sample_title,diagnosis=r.Sample_characteristics_ch1_0,tissue=r.Sample_characteristics_ch1_1,treatment_outcome=r.Sample_characteristics_ch1_7))
result=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214397',sha256={n:s for n,_,s in sources},n_lesion=51,n_healthy=6,unique_symbols=len(x),genes=rows,registered_descriptive_support=all(r['bonferroni_pass'] for r in rows),selection_status='post-P28 selected; not untouched validation or novel discovery')
(root/'results/leish_p30_result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
