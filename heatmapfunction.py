#data preprocessing libraries 
import numpy as np
import matplotlib.pyplot as mtp
import pandas as pd
import seaborn as sns
import scipy.stats as st


def heatMap(z):
    gene_id=z.iloc[:,0]
    z_num=z._get_numeric_data()
    col=[]
    
    for i in z_num.columns:
        col.append(i)
        def f_inner():
            try:
                cpm=[]
                for i in range(0,len(z_num.columns)):
                        colsum=z_num.iloc[:,i].sum()

                        c=(z_num.iloc[:,i]/colsum)*10**6
                        cpm.append(c)

            except:
                print("ignore error")
            cpm=np.array(cpm)
            logv=np.log2(cpm+1)
            return logv
    logval=f_inner()
    zscr=st.zscore(logval)
    
    zscr=pd.DataFrame(zscr)
    new=zscr.T
    #new['gene']=id1
    zscr=new.fillna(0)
    #zscr['gene']=id1
    zscr=zscr.values

    mtp.figure(figsize=(20,20))
    ax=sns.heatmap(zscr[0:10,0:10],xticklabels=col[0:10], yticklabels=gene_id[0:10],annot=False,cmap='gnuplot')
    ax.set_xticklabels(ax.get_xticklabels(), fontsize = 7)
    ax.set_yticklabels(ax.get_yticklabels(), rotation = 0, fontsize = 7)
    mtp.show()
    #mtp.savefig("heatmap1.pdf")
   
data=pd.read_csv("GSE205432_NormalizedCounts.csv")
heatMap(data)
 
