"""Prepare already-downloaded authoritative raw sources; no network at application runtime.
Usage: python scripts/prepare_full_sources.py HENKE_TAR FFAST_CHECKOUT
Downloads: https://henke.lbl.gov/optical_constants/sf.tar.gz and
https://github.com/usnistgov/FFAST.jl (pin recorded commit before rebuilding).
"""
import hashlib,json,sys,tarfile,subprocess,shutil
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import xraydb

root=Path(__file__).resolve().parents[1]
tarpath,ffast=Path(sys.argv[1]),Path(sys.argv[2])
out=root/'data/sources/full';out.mkdir(parents=True,exist_ok=True)
commit=subprocess.check_output(['git','-C',str(ffast),'rev-parse','HEAD'],text=True).strip()
stamp=datetime.now(timezone.utc).isoformat()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest={'schema_version':1,'retrieved_utc':stamp,'NIST_commit':commit,'Henke_archive_sha256':sha(tarpath),'sources':{}}
shells=pd.read_csv(ffast/'data/shelldata.csv')
for source in ('Henke','NIST'):
 (out/source).mkdir(exist_ok=True)
 manifest['sources'][source]={}
with tarfile.open(tarpath,'r:gz') as archive:
 members={Path(m.name).stem.capitalize():m for m in archive.getmembers() if m.isfile() and m.name.endswith('.nff')}
 for symbol,member in members.items():
  (out/'Henke'/f'{symbol}.nff').write_bytes(archive.extractfile(member).read())
shutil.copy2(ffast/'data/shelldata.csv',out/'NIST/shelldata.csv')
shutil.copy2(ffast/'LICENSE',out/'NIST/LICENSE')
manifest['NIST_shelldata_sha256']=sha(out/'NIST/shelldata.csv')
for source in ('Henke','NIST'):
 for z in range(1,93):
  element=xraydb.atomic_symbol(z)
  if source=='Henke':
   path=out/source/f'{element}.nff'
   raw=np.loadtxt(path,skiprows=1)[:,:3]
   raw=raw[np.argsort(raw[:,0],kind='stable')]
   _,inverse,counts=np.unique(raw[:,0],return_inverse=True,return_counts=True)
   # Conflicting duplicate-energy observations have no unique tabulated factor.
   # Remove all such knots and forbid interpolation across the resulting gap.
   valid=np.isfinite(raw).all(axis=1)&(raw[:,1]>-9999)&(raw[:,2]>=0)&(counts[inverse]==1)
   # Undefined f1 sentinels must not be crossed: boundary intervals remain unavailable.
   indices=np.flatnonzero(valid);a=raw[valid]
   gaps=[(float(a[i,0]),float(a[i+1,0]),'undefined-factor gap') for i in range(len(a)-1) if indices[i+1]!=indices[i]+1]
   # Henke explicitly places narrow pairs about absorption discontinuities.
   d=np.diff(a[:,0])
   for i in range(1,len(d)-1):
    if d[i]<0.2*min(d[i-1],d[i+1]): gaps.append((float(a[i,0]),float(a[i+1,0]),'Henke narrow edge bracket'))
   # Check every independently tabulated atomic shell. If a close Henke bracket
   # is unavailable, conservatively forbid its enclosing raw interval.
   for label,edge in xraydb.xray_edges(element).items():
    e=edge.energy
    if not a[0,0]<e<a[-1,0]: continue
    candidates=[g for g in gaps if abs((g[0]+g[1])/2-e)<=max(2.,.03*e)]
    if not candidates:
     j=int(np.searchsorted(a[:,0],e));gaps.append((float(a[j-1,0]),float(a[j,0]),f'{label}: XrayDB shell locator (conservative)'))
   url=f'https://henke.lbl.gov/optical_constants/sf/{element.lower()}.nff'
   version='Henke/Gullikson/Davis; current official sf.tar.gz snapshot'
   edges={'bracket_'+str(i):[lo,hi] for i,(lo,hi,_) in enumerate(gaps)}
  else:
   path=out/source/f'{element}.csv';shutil.copy2(ffast/f'data/mac[{z}].csv',path)
   raw=pd.read_csv(path)[['E','f1','f2']].to_numpy();raw[:,0]*=1000
   valid=np.isfinite(raw).all(axis=1)&(raw[:,2]>=0);a=raw[valid]
   gaps=[];edges={}
   for label in shells.columns[7:]:
    e=shells.iloc[z-1][label]
    if not np.isfinite(e) or not a[0,0]<e<a[-1,0]: continue
    edges[label]=float(e)
    j=int(np.searchsorted(a[:,0],e))
    # If an exact edge knot exists its value is usable, but both adjoining
    # intervals are conservatively excluded from interpolation.
    end=j+1 if a[j,0]==e and j+1<len(a) else j
    gaps.append((float(a[j-1,0]),float(a[end,0]),f'{label}: archived shelldata edge'))
   url=f'https://raw.githubusercontent.com/usnistgov/FFAST.jl/{commit}/data/mac%5B{z}%5D.csv'
   version=f'usnistgov/FFAST.jl {commit}'
  if np.any(np.diff(a[:,0])<=0) or not np.isfinite(a).all() or len(a)<4: raise ValueError(f'{source}/{element} invalid raw knots')
  # log(f2) is undefined at zero; no interpolation in intervals changing
  # between zero and positive. Exact zero knots and all-zero intervals are valid.
  for i in np.flatnonzero((a[:-1,2]==0)!=(a[1:,2]==0)):
   gaps.append((float(a[i,0]),float(a[i+1,0]),'zero/positive f2 boundary'))
  manifest['sources'][source][element]={'Z':z,'file':f'{source}/{path.name}','url':url,'version':version,'sha256':sha(path),'normalization':'stable energy sort; all duplicate-energy knots excluded; gaps blocked (raw file unchanged)','raw_rows':len(raw),'valid_rows':len(a),'domain_eV':[float(a[0,0]),float(a[-1,0])],'forbidden_intervals_eV':[list(g) for g in sorted(set(gaps))],'edges_eV':edges,'interpolation':'linear full f1(E) minus Z; log-log positive f2; zero on zero-zero intervals; no crossing forbidden intervals; exact knots retained','f0':'common native XrayDB/Waasmaier–Kirfel neutral atom at q/(4π)','factor_scope':'both selected elements use this source anomalous factors'}
(out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
target=root/'src/atrxl/data_sources/full'
shutil.copytree(out,target,dirs_exist_ok=True)
print(json.dumps({'counts':{s:len(v) for s,v in manifest['sources'].items()},'nist_commit':commit,'henke_archive_sha256':sha(tarpath)},indent=2))
