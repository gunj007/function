#data preprocessing libraries 
import numpy as np
import matplotlib.pyplot as mtp
import pandas as pd
import seaborn as sns
import scipy.stats as st


def normalization(z):
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
    zscr1=zscr.values
    return zscr1

data=pd.read_csv("GSE205432_NormalizedCounts.csv")
print(normalization(data))
 