from sklearn.decomposition import PCA
from scipy import stats
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
#import statsmodels.api as sm
sns.set(color_codes=True)
from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import TruncatedSVD, PCA, NMF, FactorAnalysis, FastICA, IncrementalPCA, KernelPCA

#df_norm = pd.DataFrame.from_csv('norm_count.csv')
#df_norm.columns=cols
#print df_raw.head()


def make_scatter_matrix(df_raw):
    sns.set(font_scale = 1.5)
    #sns.set(style="white")
    def corrfunc(x, y, **kws):
        r, _ = stats.pearsonr(x, y)
        ax = plt.gca()
        ax.annotate("r = {:.2f}".format(r),
                    xy=(.1, .9), xycoords=ax.transAxes)
    
    selection =  df_raw.iloc[:,0:4]
    selection = selection+1
    selection = selection[selection.sum(axis=1)>5]
    selection = np.log10(selection.dropna())
    print selection.head()
    #selection = selection
    g = sns.PairGrid(selection, palette=["red"]) 
    g.map_upper(plt.scatter, s=10)
    g.map_diag(sns.distplot, kde=False) 
    g.map_lower(sns.kdeplot, cmap="Blues_d") 
    g.map_lower(corrfunc)
    
    g.savefig('q1.svg')
    g.savefig('q1.png')
    g.savefig('q1.pdf')
   
    plt.show()
    
    selection =  df_raw.iloc[:,4:]
    selection.columns =['null-1','null-2','null-3','null-4']
    selection = selection+1
    selection = selection[selection.sum(axis=1)>5]
    selection = np.log10(selection.dropna())
    print selection.head()
    #selection = selection
    g = sns.PairGrid(selection, palette=["red"]) 
    g.map_upper(plt.scatter, s=10)
    g.map_diag(sns.distplot, kde=False) 
    g.map_lower(sns.kdeplot, cmap="Blues_d") 
    g.map_lower(corrfunc)
    g.savefig('q2.svg')
    g.savefig('q2.png')
    g.savefig('q2.pdf')
    plt.show()
    
'''
#make_scatter_matrix()

def make_pca(in_df):
    pca = PCA(n_components=2)
    #in_df = df_raw[df_raw.sum(axis=1)>10]
    #in_df = df_raw.replace([np.inf, -np.inf], np.nan)
    #in_df = in_df.dropna()
    #in_df = in_df[in_df.sum(axis=1)>1]
    #print in_df.head()
    #in_df=in_df+1
    #in_df=in_df/in_df.median()
    #in_df = in_df[in_df.sum(axis=1)>0.00001]
    #print in_df.head()
    #in_df = np.log2(in_df)       
    #in_df = (in_df-in_df.mean())/in_df.std()
    
    
    #in_df = in_df/in_df.median()
    #in_df.to_csv('test.csv')
    #print in_df.head()
    pca.fit(in_df)
    temp_df = pd.DataFrame()
    temp_df['pc_1']=pca.components_[0]
    temp_df['pc_2']=pca.components_[1]
    temp_df.index = cols
    #print temp_df.head()
    print(pca.explained_variance_ratio_)
    fig,ax=plt.subplots()
    temp_df.iloc[0:4,:].plot(kind='scatter',x='pc_1',y='pc_2',s=50, c='b', ax=ax,legend = 'wt')
    temp_df.iloc[4:,:].plot(kind='scatter',x='pc_1',y='pc_2',s=50, c='r', ax=ax, legend='ko')
    for i, txt in enumerate(cols):
        ax.annotate(txt, (temp_df['pc_1'].values[i]+0.0003,temp_df['pc_2'].values[i]))

    #ax.set_xlabel('PC1({})'.format(round(pca.explained_variance_ratio_[0],2)))
    #ax.set_ylabel('PC2({})'.format(round(pca.explained_variance_ratio_[1],2)))
    #ax.set_xlim(0.345,0.360)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    
    ax.yaxis.label.set_size(20)
    ax.xaxis.label.set_size(20)
    plt.tick_params(axis='both', which='major', labelsize=14)
    fig.savefig('pca.svg')
    fig.savefig('pca.png')
    fig.savefig('pca.pdf')
    plt.legend()
    plt.show()
    
def make_mdf(df):
    df=df-df.mean()
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna()
    #X_true=df.values 
    print df.head()
    '''
    similarities = euclidean_distances(X_true)
    mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=9,
                   dissimilarity="precomputed", n_jobs=1)
    pos = mds.fit(similarities).embedding_
    clf = PCA(n_components=2)
    X_true = clf.fit_transform(X_true)
    pos *= np.sqrt((X_true ** 2).sum()) / np.sqrt((pos ** 2).sum())
    pos = clf.fit_transform(pos)
    
    #fig,ax = plt.figure(1)
    #ax = plt.axes([0., 0., 1., 1.])
    plt.scatter(X_true[:, 0], X_true[:, 1], color='navy',  lw=0,
            label='True Position')
    '''
def make_ma(df):
    
    df['log10_basemean']=np.log10(df['baseMean'])
    df['p value']=-np.log10(df['pvalue'])
    print df.head()
    fig,ax=plt.subplots()
    df.plot(kind='scatter',x='log10_basemean', y='log2FoldChange', 
            s=10, c='p value',cmap='Blues',ax=ax)#,alpha=0.3)
    x= df.loc['Tb927.10.4500']['log10_basemean']
    y= df.loc['Tb927.10.4500']['log2FoldChange']
    print x,y
    ax.annotate('TbCMT1', xy=(x, y), xycoords='data', 
                xytext=(3, -1.9), textcoords='data',
                arrowprops=dict(#arrowstyle="->",
                            connectionstyle="arc3",
                            width=1,headwidth=5) ) 

    x= df.loc['Tb10.v4.0040']['log10_basemean']
    y= df.loc['Tb10.v4.0040']['log2FoldChange']
    print x,y
    ax.annotate(' ', xy=(x, y), xycoords='data', 
                xytext=(3, -1.8), textcoords='data',
                arrowprops=dict(#arrowstyle="->",
                            connectionstyle="arc3",
                            width=1,headwidth=5) ) 

    
    x= df.loc['Tb927.6.970']['log10_basemean']
    y= df.loc['Tb927.6.970']['log2FoldChange']
    ax.annotate('Cysteine Peptidase', xy=(x, y), xycoords='data', 
                xytext=(4.5, -1.7), textcoords='data',
                arrowprops=dict(#arrowstyle="->",
                            connectionstyle="arc3",
                            width=1,headwidth=5) )  

    x= df.loc['Tb927.6.1000']['log10_basemean']
    y= df.loc['Tb927.6.1000']['log2FoldChange']
    ax.annotate('', xy=(x, y), xycoords='data', 
                xytext=(4.5, -1.5), textcoords='data',
                arrowprops=dict(#arrowstyle="->",
                            connectionstyle="arc3",
                            width=1,headwidth=5) )
    
    x= df.loc['Tb927.6.1020']['log10_basemean']
    y= df.loc['Tb927.6.1020']['log2FoldChange']
    ax.annotate('', xy=(x, y), xycoords='data', 
                xytext=(4.5, -1.5), textcoords='data',
                arrowprops=dict(#arrowstyle="->",
                            connectionstyle="arc3",
                            width=1,headwidth=5) )
    
    
    
    x= df.loc['Tb927.7.6500']['log10_basemean']
    y= df.loc['Tb927.7.6500']['log2FoldChange']
    ax.annotate('VSGs', xy=(x, y), xycoords='data', 
                xytext=(3, 2), textcoords='data',
                arrowprops=dict(#arrowstyle="->",
                            connectionstyle="arc3",
                            width=1,headwidth=5) )
    
    x= df.loc['Tb427VSG-1336']['log10_basemean']
    y= df.loc['Tb427VSG-1336']['log2FoldChange']    
    ax.annotate(' ', xy=(x, y), xycoords='data', 
                xytext=(3, 2), textcoords='data',
                arrowprops=dict(#arrowstyle="->",
                            connectionstyle="arc3",
                            width=1,headwidth=5) )  


    x= df.loc['Tb427VSG-1583']['log10_basemean']
    y= df.loc['Tb427VSG-1583']['log2FoldChange']    
    ax.annotate(' ', xy=(x, y), xycoords='data', 
                xytext=(3, 2), textcoords='data',
                arrowprops=dict(#arrowstyle="->",
                            connectionstyle="arc3",
                            width=1,headwidth=5) )       

    x= df.loc['Tb427VSG-1381']['log10_basemean']
    y= df.loc['Tb427VSG-1381']['log2FoldChange']    
    ax.annotate(' ', xy=(x, y), xycoords='data', 
                xytext=(3, 2), textcoords='data',
                arrowprops=dict(#arrowstyle="->",
                            connectionstyle="arc3",
                            width=1,headwidth=5) )           
    #g = sns.jointplot("log10_basemean", "log2FoldChange", data=df, color="r",kind='scatter')
                  #xlim=(0, 60), ylim=(0, 12), color="-log10_pvalue", size=7)
    x= df.loc['Tb427VSG-702']['log10_basemean']
    y= df.loc['Tb427VSG-702']['log2FoldChange']    
    ax.annotate(' ', xy=(x, y), xycoords='data', 
                xytext=(3, 2), textcoords='data',
                arrowprops=dict(#arrowstyle="->",
                            connectionstyle="arc3",
                            width=1,headwidth=5) )     
    
    
    
    x= df.loc['Tb427VSG-7565']['log10_basemean']
    y= df.loc['Tb427VSG-7565']['log2FoldChange']
    print x,y
    ax.annotate('VSG', xy=(x, y), xycoords='data', 
                xytext=(-.3, -1.5), textcoords='data',
                arrowprops=dict(#arrowstyle="->",
                            connectionstyle="arc3",
                            width=1,headwidth=5) )     
    '''
    x= df.loc['Tb927.9.350']['log10_basemean']
    y= df.loc['Tb927.9.350']['log2FoldChange']
    print x,y
    ax.annotate('', xy=(x, y), xycoords='data', 
                xytext=(.5, -1.5), textcoords='data',
                arrowprops=dict(#arrowstyle="->",
                            connectionstyle="arc3",
                            width=1,headwidth=5) ) 
    '''
    
    x= df.loc['Tb427VSG-2']['log10_basemean']
    y= df.loc['Tb427VSG-2']['log2FoldChange']
    print x,y
    ax.annotate('MITat 1.2', xy=(x, y), xycoords='data', 
                xytext=(4.5, 1), textcoords='data',
                arrowprops=dict(#arrowstyle="->",
                            connectionstyle="arc3",
                            width=1,headwidth=5) )  
    
    
    ax.set_xlabel('Base Mean')
    ax.set_ylabel('Fold Change')
    
    
    ax.yaxis.label.set_size(20)
    ax.xaxis.label.set_size(20)
    plt.tick_params(axis='both', which='major', labelsize=14)
    #log10 mean normalized counts
    fig.savefig('ma.svg')
    fig.savefig('ma.png')
    fig.savefig('ma.pdf')    
    plt.show()


if __name__ == '__main__':
    cols = ["WT1", "WT2", "WT3", "WT4","null-1", "null-2", "null-3", "null-4"]
    #df_raw = pd.DataFrame.from_csv('voomed_data.csv')
    df_raw = pd.DataFrame.from_csv('raw_count_2.csv')    
    df_raw.columns = cols
    print df_raw.head()
    df_norm = pd.DataFrame.from_csv('norm_count_2.csv')
    df_norm.columns=cols
    
    #make_scatter_matrix(df_raw)
    
    
    make_scatter_matrix(df_raw)

    #df_norm['fc']=(df_norm.iloc[:,0:4].mean(axis=1))/(df_norm.iloc[:,4:].mean(axis=1))
    #df_norm = df_norm[(df_norm['fc']<-1.2)|(df_norm['fc']>1.2)]
    #del df_norm['fc']
    print df_norm.shape  
              
    make_pca(np.log10(df_norm+1).dropna())
    
    df = pd.DataFrame.from_csv('ordered_results_deseq_2.csv')
    #make_ma(df)  
    
    
    #make_mdf(df)
   
#df = pd.DataFrame.from_csv('ordered_results_deseq2.csv')   
#fig,ax=plt.subplots()
#stats.probplot(-np.log10(df['pvalue'].dropna()), dist="norm", plot=plt)
#probplot = sm.ProbPlot(df['pvalue'].dropna(),fit=True)
#probplot.qqplot(line='45')
#plt.show()
#make_ma()
#qqplots
#sns.pairplot(np.log10(df_raw.dropna()+1))
'''

