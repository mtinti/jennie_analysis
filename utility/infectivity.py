
from scipy import stats
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
#import statsmodels.api as sm
fig,ax=plt.subplots()
sns.set(style="ticks")
df = pd.DataFrame.from_csv('anna_data/infectivity_dataset.txt',sep='\t')
df['counts'] = df['count']#*10+1

#df['counts'] = np.log10(df['counts'] )
del df['count']
print df.head()

df['condition'] = ['minus_dox' if n == 'minus_Dox' else n for n in df['condition']]

#sns.stripplot(x="condition", y="count", data=df, palette="PRGn")
sns.boxplot(x="condition", y="counts", data=df,ax=ax)#, palette="PRGn")
sns.swarmplot(x="condition", y="counts", data=df, size=10,color=".3",linewidth=0,ax=ax)
sns.despine(offset=10, trim=True,ax=ax)
ax.set( ylabel='10^7 parasites/ml')
ax.set( xlabel='Condition')
ax.yaxis.label.set_size(20)
ax.xaxis.label.set_size(20)
#ax.set( yscale="log")
#ax.set_ylim(0,1500)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.savefig('Fig2.svg')
'''
x=df[df['condition']=='WT']['count'].values
y=df[df['condition']=='dKO4']['count'].values
res=stats.mannwhitneyu(x=x,y=y)
x1, x2 = 0, 1
y, h, col = x.max() + 0.1, 0.1, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y-0.8], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, 'pvalue='+str(round(res[-1],4)), ha='center', va='bottom', color=col)

x=df[df['condition']=='WT']['count'].values
y=df[df['condition']=='dKO2']['count'].values
res=stats.mannwhitneyu(x=x,y=y)
x1, x2 = 0, 2
y, h, col = x.max() + 0.5, 0.1, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y-1.7], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, 'pvalue='+str(round(res[-1],4)), ha='center', va='bottom', color=col)


x=df[df['condition']=='WT']['count'].values
y=df[df['condition']=='dKO3']['count'].values
res=stats.mannwhitneyu(x=x,y=y)
x1, x2 = 0, 3
y, h, col = x.max() + 1, 0.1, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y-1.7], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, 'pvalue='+str(round(res[-1],4)), ha='center', va='bottom', color=col)


x=df[df['condition']=='plus_dox']['count'].values
y=df[df['condition']=='minus_dox']['count'].values
res=stats.mannwhitneyu(x=x,y=y)
x1, x2 = 4, 5
y, h, col = x.max() + 0.1, 0.1, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y-1], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, 'pvalue='+str(round(res[-1],4)), ha='center', va='bottom', color=col)
'''


