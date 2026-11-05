"""Descriptive scientific figure; preserve per-pocket and fixed-N comparisons."""
from common import *
import os
os.environ.setdefault('MPLCONFIGDIR',str(OUT/'matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def main():
 d=load(OUT/'summary.json');ss=d['summary'];ns=load(OUT/'PER_N_DIVERSITY.json');labels=['Top-1','Categorical','T=1.5','Top-k 50','Top-p .95','Prior + iter.','Alternating'];x=np.arange(len(ss));colors=['#c84b45','#377eb8','#4daf4a','#986bb8','#db953d','#3f9e95','#8e6b58'];fig,ax=plt.subplots(1,3,figsize=(16,4.7),layout='constrained')
 ax[0].bar(x-.18,[s['identity_unique_initial'] for s in ss],.36,color='#9cc5e3',label='Initial');ax[0].bar(x+.18,[s['identity_unique_terminal'] for s in ss],.36,color='#285980',label='Terminal');ax[0].set_ylabel('Fragment classes observed');ax[0].set_title('Identity coverage (all attempts)');ax[0].legend(frameon=False)
 mean=[];err=[]
 for s in ss:
  vv=[r['native_final_valid_unique_per_attempt'] for r in ns if r['strategy']==s['strategy'] and r['receptor_id']!='ALL' and r['N_stratum']==6];m=float(np.mean(vv));mean.append(m);err.append([m-min(vv),max(vv)-m])
 ax[1].bar(x,mean,color=colors);ax[1].errorbar(x,mean,yerr=np.array(err).T,fmt='none',ecolor='#333333',capsize=3);ax[1].set_ylim(0,1.05);ax[1].set_title('Fixed N=6, within each pocket');ax[1].set_ylabel('Unique 2D / attempts (mean and range)')
 for shift,stage,color,label in [(-.24,'native_final','#a54b4b','Neural endpoint'),(0,'rigid','#6f8bb3','Rigid repair'),(.24,'local_flex','#4c9a79','Local-flex repair')]:
  ax[2].bar(x+shift,[s[stage+'_high_quality_count']/s['attempts'] for s in ss],.24,color=color,label=label)
 ax[2].set_ylim(0,1.05);ax[2].set_ylabel('HQ_no_water / all attempts');ax[2].set_title('Geometry quality, full denominator');ax[2].legend(frameon=False,fontsize=8)
 for a in ax:
  a.set_xticks(x,labels,rotation=35,ha='right');a.spines[['top','right']].set_visible(False);a.grid(axis='y',alpha=.2);a.set_axisbelow(True)
 fig.savefig(OUT/'DIVERSITY_QUALITY.png',dpi=200);fig.savefig(OUT/'DIVERSITY_QUALITY.svg');plt.close(fig)
if __name__=='__main__':main()

