
from scipy import stats
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
#import statsmodels.api as sm
sns.set(color_codes=True)
##gammas = sns.load_dataset("gammas")
##print gammas



df = pd.DataFrame.from_csv('deltaCt.txt',sep='\t')
fig,ax = plt.subplots()
df['remove']=[1 if n%2 == 0 else 0 for n in df['class']]
df = df[df['remove']==0]
#df['value'] = np.log10(df['value'])
sns.boxplot(x="class", y="dCt", data=df, ax=ax)#, palette="PRGn")
sns.swarmplot(x="class", y="dCt", data=df, size=8,color=".3",linewidth=0,ax=ax)
sns.despine(offset=10, trim=True, ax=ax)
ax.set_xticklabels(['+tet', '24h','48h','72h'])
ax.set(xlabel='Condition')
ax.yaxis.label.set_size(20)
ax.xaxis.label.set_size(20)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.savefig('fig_1A.eps')
plt.show()


sns.boxplot(x="condition", y="dCt", data=df)#, palette="PRGn")
sns.swarmplot(x="condition", y="dCt", data=df, size=5, color=".3",linewidth=0)
sns.despine(offset=10, trim=True)
plt.show()

'''
data = [float(n) for  n in open('wb.txt').read().split('\n')]
fig,ax=plt.subplots()      
#ax = sns.swarmplot(data=data, size=10, color=".3",linewidth=0, ax=ax)
ax.boxplot(data)
for y in data:
    ax.plot(1,y,'b.',markersize=15)
plt.show()
'''

fig,ax=plt.subplots()      
df = pd.DataFrame.from_csv('wb2.txt', sep='\t')
df['condition'] = [n[0:-1]+str(int(n[-1])-1) if n != 'SM' else 'WT' for n in df['condition']]
df['condition'] = [n if n != 'dKO1' else 'null-1' for n in df['condition']] 
df['condition'] = [n if n != 'dKO2' else 'null-2' for n in df['condition']]
df['condition'] = [n if n != 'dKO3' else 'null-3' for n in df['condition']]
print df.head()
sns.barplot(x='condition', y='value', data=df,ax=ax)
ax.yaxis.label.set_size(20)
ax.xaxis.label.set_size(20)
ax.set_ylabel('Normalized Intensity')
ax.set_xlabel('Condition')
plt.tick_params(axis='both', which='major', labelsize=14)
plt.savefig('fig_3B.eps')
plt.show()


data = pd.DataFrame.from_csv('anna_data/cell_growth_cKO_2.txt',sep='\t')
#data['count']=np.log10(data['count'])
#data['count']=np.log10(data['count'])
data.columns = ['Days','condition','subject','Counts']
print data.head()
ax = sns.tsplot(time="Days", value="Counts",
                unit="subject", condition="condition",
                 data=data,err_style="unit_points")
ax.set_xlim(1.5,10.5)
ax.yaxis.label.set_size(20)
ax.xaxis.label.set_size(20)
#ax.set( yscale="log")
ax.set_ylabel('Log10( 10^5 parasites/ml )')
plt.tick_params(axis='both', which='major', labelsize=14)
plt.legend(prop={'size': 16},loc=2)
plt.savefig('fig_1B.eps')
plt.show()


data = pd.DataFrame.from_csv('anna_data/cell_growth_nKO_2.txt',sep='\t')
#data['count']=np.log10(data['count'])
#data['count']=np.log10(data['count'])
print data.head()
data['condition'] = [n if n == 'WT' else 'dKO' for n in data['condition']]

data.columns = ['condition','Days','subject','Counts']
ax = sns.tsplot(time="Days", value="Counts",
                unit="subject", condition="condition",
                 data=data,err_style="unit_points")
ax.set_xlim(1.5,10.5)
ax.set_xlim(1.5,10.5)
ax.yaxis.label.set_size(20)
ax.xaxis.label.set_size(20)
#ax.set( yscale="log")
ax.set_ylabel('Log10( 10^4 parasites/ml )')
plt.tick_params(axis='both', which='major', labelsize=14)
plt.legend(prop={'size': 16},loc=2)
plt.savefig('fig_1C.eps')
plt.show()

'''
data = pd.DataFrame.from_csv('anna_data/cell_growth_nKO.txt',sep='\t')
data['condition'] = [n if n == 'WT' else 'dKO' for n in data['condition']]
data.columns = ['condition','Days','subject','Counts']
data['Counts']=data['Counts']*10000
data['Counts']=np.log10(data['Counts'])
print data.head()
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

x = 10 ** np.arange(1, 10)
y = x * 2
#data = pd.DataFrame(data={'x': x, 'y': y})

f, ax = plt.subplots(figsize=(7, 7))
#ax.set(yscale="log")
sns.regplot("Days", "Counts", data[data['condition']=='WT'], ax=ax,  scatter_kws={"s": 100})
'''