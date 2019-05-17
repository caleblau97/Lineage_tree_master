#!/usr/bin/env python
# coding: utf-8

# In[12]:


#get_ipython().run_line_magic('matplotlib', 'inline')
from importlib import reload
import LinVarMap
import pandas as pd
import numpy as np
import sys
#np.set_printoptions(threshold=sys.maxsize)
#import itertools


# In[13]:

### TO TEST GIT HUB
def common_anc(node0,node1):
    z=[node0,node1]
    com=''
    for i in range(len(min(z,key=len))):
        if node0[i] is node1[i]:
            com=com+node0[i]
        else:
            break
    return com

def cov_pattern(X):
    cbins=X.columns
    #if X is not isinstance(X, pd.DataFrame):
    cov_pattern = pd.DataFrame(index=cbins, columns=cbins)
    for i in cbins:
        for j in cbins:
            mcra = common_anc(i,j) # name of mrca of i and j
            i_gen,j_gen,m_gen = len(i),len(j),len(mcra) # generation of i,j and generation mcra of i and j
            cov_pattern.loc[i,j] = str(i_gen)+str(j_gen)+str(m_gen)
            #print(cov_pattern.loc[i,j])
    return cov_pattern,cbins

def cov(Xsim): #pandas dataframe
    X_val=Xsim.values
    n,p=X_val.shape
    node=len(Xsim.columns)
    patt,cbins=cov_pattern(Xsim)
    #creating empty arrays to put calculated sums
    xx_cum=np.zeros([node,node])
    x_cum=np.zeros([1,node])
    
    #average cell size over samples
    d_avg=np.average(X_val,axis=0)

    #calculating sample summations in covariance equation
    for i in range(n):
        xx_cum=xx_cum+np.outer(X_val[i,:],X_val[i,:])
        x_cum=x_cum+(X_val[i,:])
    #for three generations mu_T,X and X_T,mu are the same
    unstr_cov=(1/(n-1))*(xx_cum-2*np.outer(d_avg,x_cum)+n*np.outer(d_avg,d_avg)) #where does n next to avg outer product come from
    pd_unstr_cov=pd.DataFrame(unstr_cov,index=cbins,columns=cbins)
    #calculating averaged covariances
    uniq=np.unique(patt)
    avg_cov=np.zeros([len(uniq)])
    
    for i in range(len(uniq)):
        uniq_pos=np.where(patt==uniq[i]) #to check print uniq[i] and unstr_cov[uniq_pos] in loop
        avg_cov[i]=np.average(unstr_cov[uniq_pos])
    
    avg_cov= pd.DataFrame(avg_cov,index=uniq)
    
    #calculating average averages
    uniq_diag=np.diag(patt)
    avg_mu=np.zeros([len(uniq_diag)])
    
    for i in range(len(uniq_diag)):
        uniq_pos=np.where(uniq_diag==uniq_diag[i]) #to check print uniq[i] and unstr_cov[uniq_pos] in loop
        avg_mu[i]=np.average(d_avg[uniq_pos])

    return pd_unstr_cov,avg_cov,avg_mu


# ### T cell data

# In[14]:


#gmax=3
#X = LinVarMap.X_from_excel('Tcells/Tcell-Area-Gen1to8.xlsx',gmax=gmax)
#print(X)
#print(X)
#crow=X.index
#print(crows)
#print(cbins)
#a= X.columns.shape
#print(a)
#zeros_X=list(range(1,16))
#print(zeros_X)
#c1='1011'
#c2='11'
#print(X)
      
#cov_pattern=cov_pattern(X)
#print(common_anc(c1,c2))

#print(cov_pattern(X.columns))
#cov_pattern=cov_pattern(X.columns)
#cov_pattern = pd.DataFrame(index=X.columns, columns=X.columns) #index=rows, columns = columns the data frame creates an data array table thing 
#print(cov_pattern)
#print(cov_pattern.loc[10,10])
#Xpatt = LinVarMap.make_cov_pattern(X.columns)
#print(Xpatt)
#mvn = LinVarMap.MVN(em_tol=1e-2, em_maxIter=5000, em_verbose=True)
#fit=mvn.fit(X,pattern=Xpatt, MLEalgo='ortho_select',Morder=1)


# ### Simulated data

# In[15]:


sim = LinVarMap.simData()
sim.makeX(n_gens=4, n_samples=5, prob=0.8, missing=0, seed=1)
Xsim = sim.X.copy()
pXsim = sim.cov_pattern.copy()
cov_true = sim.cov_true.copy()
mvn = LinVarMap.MVN(em_tol=1e-3, em_maxIter=10000, em_verbose=True)
mvn.fit(Xsim,pattern=pXsim,MLEalgo='ortho_select',Morder=1)


# In[16]:


unstr_cov,avg_cov,avg_mu=cov(Xsim)
#print(Xsim)
#print(cov_pattern(Xsim))
print(avg_cov)
print(unstr_cov)
print(avg_mu)
print('//')
print(Xsim.cov()-unstr_cov)


# ### Worm data

# In[225]:


#worm = LinVarMap.wormData() #Create wormData object
#worm.load_worm_data(imarker=1) #load (unbalanced) data from txt files
#worm.balance_pedigree(gmax=8)
#mvn = LinVarMap.MVN(em_tol=1e-3, em_maxIter=10000, em_verbose=True) 
## Fit mvn to worm data
#mvn.fit(worm.X, pattern=worm.cov_pattern, MLEalgo='ortho_select',Morder=1)
#
#
## ### Visualize results
#
## In[33]:
#
#
##mvn.graphplot(data='pcorr')
#plt.show()
#
## In[34]:
#
#
#mvn.specgraph(data='Q_Linv')
#
#
## In[35]:
#
#
#mvn.explvariance(output=None)
#
#
## In[ ]:




