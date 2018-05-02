from scipy import stats
from scipy.spatial import distance
import pandas as pd
import numpy as np
import re
from sklearn import cluster
from sklearn.metrics import  silhouette_score
import matplotlib
import matplotlib.pyplot as plt
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
import statsmodels.stats.multitest as smm
from sklearn.decomposition import PCA
import seaborn as sns
from statsmodels.sandbox.stats.multicomp import multipletests
#matplotlib.style.use('ggplot')
#cmap = matplotlib.cm.get_cmap('Dark2')

def quantileNormalize(df_input, keep_na=True):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        norm = [rank[i] for i in t]
        if keep_na == True:
            norm = [np.nan if np.isnan(a) else b for a,b in zip(df[col],norm)]
        df[col] =  norm             
    return df

def norm(X):
    return X/X.mean()

def norm_last(X):
    return X/X[-1]

def constandNormalize(df='', n_channels= 9 ):
    n_channels= n_channels
    
    def compute(X):
        X = X / (n_channels * X.mean())
        return X
        
    def step(df):
        df = df.apply(compute, axis=1)
        error_row = abs(df.mean(axis=0) - (float(1)/n_channels )  ).sum()/2
        df = df.apply(compute, axis=0)
        error_col = abs(df.mean(axis=1) - (float(1)/n_channels )  ).sum()/2
        return df, error_row, error_col
        
    def normalize(df):
        error_col = 1
        error_row = 1
        a=0
        while 1:
            a+=1
            df, error_row, error_col = step(df)
            print (a,error_row, error_col)
            if a >50:
                print ('normalized in',a,'steps')
                return df
            if (error_row < 1e-5) and (error_col < 1e-5):
            #if ((error_row + error_col)/2 )< 1e-5: 
            #if (error_row < 1e-5) or (error_col < 1e-5):             
                print ('normalized in',a,'steps')
                return df 
    return normalize(df)



def make_pca(in_df):
    cols = in_df.columns
    pca = PCA(n_components=2)
    pca.fit(in_df)
    temp_df = pd.DataFrame()
    temp_df['pc_1']=pca.components_[0]
    temp_df['pc_2']=pca.components_[1]
    temp_df.index = cols
    print(pca.explained_variance_ratio_)
    fig,ax=plt.subplots()
    temp_df.iloc[0:3,:].plot(kind='scatter',x='pc_1',y='pc_2',s=50, c='b', ax=ax,legend = 'ND')
    temp_df.iloc[3:6,:].plot(kind='scatter',x='pc_1',y='pc_2',s=50, c='r', ax=ax, legend='HF')    
    temp_df.iloc[6:9,:].plot(kind='scatter',x='pc_1',y='pc_2',s=50, c='r', ax=ax, legend='KO')
    for i, txt in enumerate(cols):
        ax.annotate(txt, (temp_df['pc_1'].values[i]+0.0003,temp_df['pc_2'].values[i]))
    #ax.set_xlabel('PC1({})'.format(round(pca.explained_variance_ratio_[0],2)))
    #ax.set_ylabel('PC2({})'.format(round(pca.explained_variance_ratio_[1],2)))
    #ax.set_xlim(0.345,0.360)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    
    ax.yaxis.label.set_size(20)
    ax.xaxis.label.set_size(20)
    #plt.tick_params(axis='both', which='major', labelsize=16, rotation=70, horizontalalignment='right' )
    #plt.setp( axs[1].xaxis.get_majorticklabels(), )
    fig.savefig('pca.svg')
    fig.savefig('pca.png')
    fig.savefig('pca.pdf')
    plt.legend()
    plt.show()
    
def make_scatter_matrix(in_df, fig_tag):
    sns.set(font_scale = 1.5)
    #sns.set(style="white")
    def corrfunc(x, y, **kws):
        r, _ = stats.pearsonr(x, y)
        ax = plt.gca()
        ax.annotate("r = {:.2f}".format(r),
                    xy=(.1, .9), xycoords=ax.transAxes)
    
    g = sns.PairGrid(in_df, palette=["red"]) 
    g.map_upper(plt.scatter, s=10)
    g.map_diag(sns.distplot, kde=False) 
    g.map_lower(sns.kdeplot, cmap="Blues_d") 
    g.map_lower(corrfunc)
    g.savefig(fig_tag+'_sm.svg')
    g.savefig(fig_tag+'_sm.png')
    g.savefig(fig_tag+'_sm.pdf')
    plt.show()

def make_analysis(norm_df, onesample_results, tag):
    norm_df['T_'+tag]=onesample_results[0]
    norm_df['PVAL_'+tag]=onesample_results[1]
    #p-value distributions
    norm_df['PVAL_'+tag].plot(kind='hist',bins=100)
    plt.title('p-val distributions')
    plt.show()

    #correct for multy-hypotesis
    p_adjusted = multipletests(norm_df['PVAL_'+tag], method='bonferroni')
    norm_df['PVAL_'+tag+'_ADJ']=p_adjusted[1]


    norm_df.plot(kind='scatter',x='exp10_'+tag,y='log2_fc_'+tag)
    plt.title('MA plot')
    plt.show()

    norm_df['-log10_PVAL_'+tag] = -np.log10(norm_df['PVAL_'+tag])
    norm_df.plot(kind='scatter',x='log2_fc_'+tag,y='-log10_PVAL_'+tag)
    plt.title('vulcano plot')
    
    
if __name__ == '__main__':
    pass
'''
def make_anova(in_df):
    def interpolate(in_series):
        in_series = in_series.replace(0.0, np.nan)
        in_series = in_series.interpolate() 
        return in_series
    
    res_1=[]
    res_2=[]
    cols_1 = ['E5014_Reporter intensity ' +str(n) for n in range(1,10)]
    cols_2 = ['E5015_Reporter intensity ' +str(n) for n in range(1,10)]
    cols_3 = ['E5016_Reporter intensity ' +str(n) for n in range(1, 10)]
              
    #print in_df.head()
    for prot in in_df.index.values:
        #prot='Tb927.2.5810'
        #print prot
        temp_df=pd.DataFrame()

        val_1 = interpolate(in_df.loc[prot][cols_1])
        val_2 = interpolate(in_df.loc[prot][cols_2])
        val_3 = interpolate(in_df.loc[prot][cols_3]) 
        
        to_use = []
        
        for n in [val_1,val_2,val_3]:
            if n.isnull().values.any():
                pass
            else:
                to_use.append(n)
                
        
        temp_df = pd.DataFrame()
        temp_df['values'] = pd.concat(to_use).values
        replica = []
        a=1
        for n in to_use:
            replica+=[a for n in cols_1]
            a+=1
        temp_df['replica']=replica
        tp = []
        for n in to_use:
            tp+=[n for n in range(1,10)]  
        #print temp_df
        model = ols(data=temp_df, formula='values ~  C(tp)').fit()
        aov_table = anova_lm(model, typ=1)
        res_1.append( aov_table['F'][0] )
        res_2.append( aov_table['PR(>F)'][0]) 
        #print aov_table
        #break
    rej, pval_corr = smm.multipletests(res_2, alpha=0.05, method='bonferroni', returnsorted=False)[:2]
    return res_1,res_2,pval_corr
    

def make_anova2(in_df):
    def interpolate(in_series):
        in_series = in_series.replace(0.0, np.nan)
        in_series = in_series.interpolate() 
        return in_series
    
    res_1=[]
    res_2=[]
    cols_1 = ['E5014_Reporter intensity ' +str(n) for n in range(1,10)]
    cols_2 = ['E5015_Reporter intensity ' +str(n) for n in range(1,10)]
    cols_3 = ['E5016_Reporter intensity ' +str(n) for n in range(1,10)]
              
    #print in_df.head()
    for prot in in_df.index.values:
        #prot='Tb927.2.5810'
        #print prot
        temp_df=pd.DataFrame()

        val_1 = interpolate(in_df.loc[prot][cols_1])
        val_2 = interpolate(in_df.loc[prot][cols_2])
        val_3 = interpolate(in_df.loc[prot][cols_3]) 
        
        to_use = []
        
        for n in [val_1,val_2,val_3]:
            if n.isnull().values.any():
                pass
            else:
                to_use.append(n.values)
        
        temp_df = pd.DataFrame()
        for index,item in enumerate(to_use):
            temp_df['S'+str(index)]=item
        
        #print temp_df
        #temp_df = pd.concat(to_use, axis=1).values
        #print temp_df.iloc[8,:]
        F, p = stats.f_oneway(temp_df.iloc[0,:].values,temp_df.iloc[1,:].values,temp_df.iloc[2,:].values,
                              temp_df.iloc[3,:].values,temp_df.iloc[4,:].values,temp_df.iloc[5,:].values,
                                temp_df.iloc[6,:].values,temp_df.iloc[7,:].values,temp_df.iloc[8,:].values)
        res_1.append(F)
        res_2.append(p)
        #print p
    rej, pval_corr = smm.multipletests(res_2, alpha=0.05, method='bonferroni', returnsorted=False)[:2]
    return res_1,res_2,pval_corr                 

    

def fold_change(X):
    return X.max()/X.min()


#this function perform the quantile normalization of a dataset
#adapted from https://github.com/ShawnLYU/Quantile_Normalize   
def quantileNormalize(df_input, keep_na=True):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        norm = [rank[i] for i in t]
        if keep_na == True:
            norm = [np.nan if np.isnan(a) else b for a,b in zip(df[col],norm)]
        df[col] =  norm             
    return df

#this function remove from the protein txt output file of MaxQuant the 
#raws that are not useful for the analysis, filter out the maxquant rubbish
#columns plus the protein with less than X unique peptides
def clean(df, cols = [], unique_peptide_limits=0):    
    start = df.shape
    print 'starting:', df.shape[0],'ids'
    for col in cols:
        temp = df.shape[0]
        df['filter']=[1 if n == '+' else 0 for n in df[col]]
        df = df[df['filter'] == 0]
        print 'removed ', temp - df.shape[0], col        
    print 'tot ', start[0]-df.shape[0] , ' entries removed'
    
    df['unique_filter'] = [int(n.split(';')[0]) for n 
                           in df['Peptide counts (unique)']]
    temp = df.shape[0]
    df = df[df['unique_filter'] >= unique_peptide_limits]   
    print 'removed uni-pep', temp - df.shape[0]
    return df

#this function remove any text between parenthesis in a text
#used to extract the sequence of the triptic peptide from the 
#probility columns of the MaxQuant output 
def remove_between(text='', start='', end=''):
    new = re.sub(r'\{start}[^{end}]*\)'.format(start=start, end=end), '', text)
    return new

#constand normalization as described in:
#http://www.mcponline.org/content/15/8/2779 
def constandNormalize(df='', n_channels= 9 ):
    n_channels= n_channels
    
    def compute(X):
        X = X / (n_channels * X.mean())
        return X
        
    def step(df):
        df = df.apply(compute, axis=1)
        error_row = abs(df.mean(axis=0) - (float(1)/n_channels )  ).sum()/2
        df = df.apply(compute, axis=0)
        error_col = abs(df.mean(axis=1) - (float(1)/n_channels )  ).sum()/2
        return df, error_row, error_col
        
    def normalize(df):
        error_col = 1
        error_row = 1
        a=0
        while 1:
            a+=1
            df, error_row, error_col = step(df)
            print a,error_row, error_col
            if a >50:
                print 'normalized in',a,'steps'
                return df
            if (error_row < 1e-5) and (error_col < 1e-5):
            #if ((error_row + error_col)/2 )< 1e-5: 
            #if (error_row < 1e-5) or (error_col < 1e-5):             
                print 'normalized in',a,'steps'
                return df 
    return normalize(df)

    
def correlate_pair(df_1 = '', df_2 = '', 
                   ax='', coor_type='',
                   label = ''):
    common = set(df_1.index.values) & set(df_2.index.values)
    values = []
    if coor_type == 'spearman':
        for item in common:
            values_a = df_1.loc[item].values
            values_b = df_2.loc[item].values
            values.append(stats.spearmanr(values_a, values_b)[0])
            
    if coor_type == 'pearson':
        for item in common:
            values_a = df_1.loc[item].values
            values_b = df_2.loc[item].values
            values.append(stats.pearsonr(values_a, values_b)[0])
    
    temp=pd.Series(values)
    label = label+ ' n>0.7='+str(temp[temp>0.7].shape[0])
    temp.plot(kind='kde',ax=ax,label=label)
    return temp


    
def correlate_pair_all(df_1 = '', df_2 = '',  df_3 = '',
                       ax='', cor_type='',
                   label = '', strategy = ''):
    
    if strategy == 'stringent':
        #print 1
        common = set(df_1.index.values) & set(df_2.index.values)  & set(df_3.index.values)
        values = []
        if cor_type == 'spearman':
            for item in common:
                values_a = df_1.loc[item].values
                values_b = df_2.loc[item].values
                values_c = df_3.loc[item].values
                temp_a = stats.spearmanr(values_a, values_b)[0]
                temp_b = stats.spearmanr(values_a, values_c)[0]
                temp_c = stats.spearmanr(values_b, values_c)[0]
                temp = np.array([temp_a,temp_b,temp_c])
                values.append(temp.mean())
                
        if cor_type == 'pearson':
            for item in common:
                #print item
                values_a = df_1.loc[item].values
                values_b = df_2.loc[item].values
                values_c = df_3.loc[item].values
                temp_a = stats.pearsonr(values_a, values_b)[0]
                temp_b = stats.pearsonr(values_a, values_c)[0]   
                temp_c = stats.pearsonr(values_b, values_c)[0]
                temp = np.array([temp_a,temp_b,temp_c])
                values.append(temp.mean())
        temp=pd.Series(values)
        label = label+ ' n>0.7='+str(temp[temp>0.7].shape[0])
        temp.plot(kind='kde',ax=ax,label=label)
        return temp    
    
    
    
def correlate_columns(df_1 = '', df_2 = '', 
                   ax='', coor_type='',
                   label = ''):
    common = set(df_1.index.values) & set(df_2.index.values)
    df_1 = df_1.loc[common]
    df_2 = df_2.loc[common]   
    values = []
    if coor_type == 'spearman':
        for n_1,n_2 in zip(df_1.columns.values,df_2.columns.values):
            values_a = df_1[n_1].values
            values_b = df_2[n_2].values
            #print values_a.shape
            #print values_b.shape
            #print stats.spearmanr(values_a, values_b)
            values.append(stats.spearmanr(values_a, values_b)[0])
            
    if coor_type == 'pearson':
        for n_1,n_2 in zip(df_1.columns.values,df_2.columns.values):
            values_a = df_1[n_1].values
            values_b = df_2[n_2].values
            values.append(stats.spearmanr(values_a, values_b)[0])
    
    temp=pd.Series(values)
    ax.plot(np.arange(1,len(values)+1), temp.values,label=label)
    
    return temp
    
def get_gene_description(in_fasta):
    res = {}
    for l in open(in_fasta):
        if l.startswith('>'):
            a=0
            item_list = l.split(' | ')
            temp_id = strip(item_list[0])[1:]
            for item in item_list:
                if item.startswith('product='):
                    desc = strip(item.split('=')[1])
                    a=1
            if a == 1:
                res[temp_id]=desc
            else:
                res[temp_id]='no_desc'

    return res


    
def plot_bic(df, n_of_cluster=20): 
    
    def compute_bic(kmeans,X):
        centers = [kmeans.cluster_centers_]
        labels  = kmeans.labels_
        m = kmeans.n_clusters
        n = np.bincount(labels)
        N, d = X.shape
        #compute variance for all clusters beforehand
        cl_var = (1.0 / (N - m) / d) * sum([sum(distance.cdist(X[np.where(labels == i)], [centers[0][i]], 'euclidean')**2) for i in range(m)])
        const_term = 0.5 * m * np.log(N) * (d+1)
        BIC = np.sum([n[i] * np.log(n[i]) -
                   n[i] * np.log(N) -
                 ((n[i] * d) / 2) * np.log(2*np.pi*cl_var) -
                 ((n[i] - 1) * d/ 2) for i in range(m)]) - const_term
        return(-BIC)
    
    
    X = df.values
    range_n_clusters = np.arange(1,n_of_cluster,1)
    scores = []
    std_scores = []
    for n in range_n_clusters:
        temp = []
        print n,
        for fold in range(5):
            KMeans = cluster.KMeans(n_clusters = n, init="k-means++").fit(X) 
            BIC =compute_bic(KMeans,X)
            temp.append(BIC)
        scores.append(np.average(temp))
        std_scores.append(np.std(temp))        
    print 'done'        
    scores = np.array(scores)
    std_scores = np.array(std_scores)
    plt.plot(range_n_clusters, scores)
    plt.plot(range_n_clusters,scores+std_scores,c='r',alpha=0.5)
    plt.plot(range_n_clusters,scores-std_scores,c='r',alpha=0.5)
    plt.xticks(range_n_clusters)
    plt.xlabel("# clusters")
    plt.ylabel("# BIC")
    plt.show()
    
    
def siluette_scoring(df,n_of_cluster=20):
    X=df.values 
    range_n_clusters = np.arange(2,n_of_cluster,1)
    scores = []
    std_scores = []
    for n in range_n_clusters:
        temp = []
        print n,
        for item in range(100):
            clusterer = cluster.KMeans(n_clusters=n)
            cluster_labels = clusterer.fit_predict(X)
            silhouette_avg = silhouette_score(X, cluster_labels)
            temp.append(silhouette_avg)
        #print "For n_clusters =", n, "The average silhouette_score is :", np.average(temp)
        scores.append(np.average(temp))
        std_scores.append(np.std(temp))
    print 'done'    
    scores = np.array(scores)
    std_scores = np.array(std_scores)
    temp_res = {}
    temp_res['scores']=scores
    temp_res['std_scores']=std_scores
    temp_res['clusters']=np.arange(2,n_of_cluster,1)
    temp_res = pd.DataFrame.from_dict(temp_res)
    plt.plot(temp_res['clusters'], temp_res['scores'],c='b')
    plt.plot(temp_res['clusters'],temp_res['scores']+temp_res['std_scores'],c='r',alpha=0.5)
    plt.plot(temp_res['clusters'],temp_res['scores']-temp_res['std_scores'],c='r',alpha=0.5)
    plt.xticks(range_n_clusters)
    #plt.xlim((2,7))
    #plt.savefig('silhouette_score.png')
    #plt.savefig('silhouette_score.svg')
    plt.show()
    return temp_res
    #temp_res.to_csv('silhouette_score.txt',sep='\t')


def fuzzy_scoring(df,n_of_cluster=20): 
    fig,ax=plt.subplots(ncols=1,nrows=1)
    range_n_clusters = np.arange(2,n_of_cluster,1)
    X=df.values
    fpcs = []
    fpcs_std = []
    for n in range_n_clusters:
        print n,
        temp = []
        for item in range(100):
            cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(X, n, 2, error=0.005, maxiter=1000, init=None)
            temp.append(fpc)
        fpcs.append(np.average(temp))
        fpcs_std.append(np.std(temp))
    
    print 'done'  
    scores = np.array(fpcs)
    std_scores = np.array(fpcs_std)
    ax.plot(range_n_clusters, scores)
    ax.plot(range_n_clusters,scores+std_scores,c='r',alpha=0.5)
    ax.plot(range_n_clusters,scores-std_scores,c='r',alpha=0.5)
    ax.set_xlabel("Number of centers")
    ax.set_ylabel("Fuzzy partition coefficient")
    #ax2.set_ylim(0.65,1.1)
    ax.set_xticks(range_n_clusters)
    return fig
    
    
    
def get_nofupep(df):
    proteins_id = df.index.values
    n_of_pep = df['Peptide counts (unique)'].values
    res = dict(zip(proteins_id,n_of_pep))
    return res
    
    
def make_K_mean_clustering(df_train, df_test, n_of_cluster=9):
    clusterer = cluster.KMeans(n_clusters=n_of_cluster)
    clusterer.fit(df_train.values )
    cluster_labels = clusterer.predict(df_test.values)
    return cluster_labels
    

def plot_clusters(df,nrows,ncols):
    n_clusters = len(set(df['cluster'].values))
    a=0
    fig,axs = plt.subplots(nrows = nrows, ncols=ncols, figsize=(12,12) )
    for row in np.arange(0,nrows,1):
        for col in np.arange(0,ncols,1):
            if a >=n_clusters:
                continue
            axs[row,col].plot()
            axs[row,col].set_title('cluster: '+str(a))
            temp_df = df[df['cluster']==a].iloc[:,0:9]
            
            #print temp_df.head()
            #del temp_df['cluster']
            #del temp_df['std_']
            #temp_df = temp_df.apply(normalize,1)
            #temp_df.plot(kind = 'box', ax=axs[row,col])#,alpha=0.1,c='r')
            temp_df.columns = np.arange(1,len(temp_df.columns)+1, 1)
            #print temp_df.head()
            color = float(a)/n_clusters
            temp_df.T.plot( ax=axs[row,col], alpha=0.3, c=cmap(color))
            #locs, labels = plt.xticks() 
            temp_df.plot(kind = 'box', ax=axs[row,col])#,positions=x)
            #plt.xticks(locs)
            axs[row,col].legend().set_visible(False)
            axs[row,col].set_ylim(0,1.1)
            axs[row,col].set_xlim(0,10)
            axs[row,col].set_title(str(a))
            #axs[row,col].set_xtick(np.arange(-1,11,1))
            a+=1  
            #res_ref[a]=temp_df.median()
            #print temp_df.median()
    plt.tight_layout()
    #plt.savefig('predict_clusters.png')
    #plt.savefig('predict_clusters.svg')
    return fig
    
def attach_external(df, reference):
    res = []
    for n in df.index.values:
        temp_res = 0
        prots = n.split(';')
        for prot in prots:
            print prot, n
            if prot in reference:
                temp_res=1
                break
        res.append(temp_res)
    return res
'''    
    
    
    
    
    

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        