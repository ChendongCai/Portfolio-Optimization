#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 11:35:12 2018

@author: chendongcai
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 15:33:08 2017

@author: chendongcai
"""

import pandas as pd
import numpy as np
import random as rd
import numba
import time
import os
import matplotlib.pyplot as plt
from sklearn.utils import shuffle

f=pd.read_table('/Files/INFORMS Competition/Updated/Timeseries_data_SP500.txt')
result_template=pd.read_csv('/Files/INFORMS Competition/Updated/result template/results_template.csv')

rootdir='/Files/INFORMS Competition/Updated/Riskmodels/'
allfile=os.listdir(rootdir)
allfile=sorted(allfile)
cov_file=[]
for i in range(len(allfile)):
    path=os.path.join(rootdir,allfile[i])
    locals()['cov'+str(i)]=pd.read_csv(path)
    cov_file.append(locals()['cov'+str(i)])
  
by_sector=f.groupby(f.SECTOR).sum()
by_date=f.groupby(f.DATE).sum()
by_date.index=pd.to_datetime(by_date.index)
by_date=by_date.sort_index()
date=f.DATE.unique()
date_file=[]
for i in range(len(date)):
    locals()['d'+str(i)]=f[f.DATE==date[i]]
    date_file.append(locals()['d'+str(i)])
##seperate into different date
#initiate gene
#def gen():
    
pc=0.6
#交叉概率
pm=0.01
#变异概率
#%%
def createcovm(date,cov):    
    df=cov
    df_temp = df.pivot_table(index='ROW_INDEX', columns='COLUMN_INDEX', values='VALUE', fill_value = 0)
    df_original = df_temp.copy()
    np.fill_diagonal(df_temp.values, 0)
    final_df = (df_original.T + df_temp).reindex(index = date.SEDOL.unique(), columns = date.SEDOL.unique())    
    if date is d124:
        final_df.loc["BYT3MK1", "BYT3MK1"] = 0.5
        final_df.fillna(0, inplace = True)    
    covm=final_df
    return covm

covm_file=[0]*len(cov_file)
for i in range(len(cov_file)):
    covm_file[i]=createcovm(date_file[i],cov_file[i])   
#%%
def initiate1(date,n,lim):
    length=len(date.SEDOL.unique())
    pool=np.array([i for i in range(length)])#[date.BETA<=1]
    pool=list(pool)
    mcapq=np.array(date.groupby(by='MCAP_Q')['MCAP_Q'].count())
#    pool=[i for i in range(date.index[0],date.index[-1]+1)]
#    origin=[i for i in range(0,length)]
    order={}
    gene_len=len(pool)
    for i in range(n):
        locals()['o'+str(i)]=rd.sample(pool,gene_len)#100--length
        order[i]=locals()['o'+str(i)]
    weight={}
    for i in range(n):
        locals()['w'+str(i)]=[rd.uniform(0,lim) for i in range(gene_len)]#100--length
        weight[i]=locals()['w'+str(i)]
    return (order,weight,pool)

#initiate stragety 
def initiate(date,n,lim):
    length=len(date.SEDOL.unique())
    pool=np.array([i for i in range(length)])[date.BETA<=1.0]
    pool=list(pool)
    mcapq=np.array(date.groupby(by='MCAP_Q')['MCAP_Q'].count())
#    pool=[i for i in range(date.index[0],date.index[-1]+1)]
#    origin=[i for i in range(0,length)]
    order={}
    gene_len=len(pool)
    for i in range(n):
        locals()['o'+str(i)]=rd.sample(pool,gene_len)#100--length
        order[i]=locals()['o'+str(i)]
    weight={}
    for i in range(n):
        locals()['w'+str(i)]=[]#100--length
        weight[i]=locals()['w'+str(i)]
        for j in order[i]:
            weight[i].append(rd.uniform(date.BENCH_WEIGHT.iloc[j]*0.9,date.BENCH_WEIGHT.iloc[j]*1.1))#0.0005 adjusted
    return (order,weight,pool)
#n=20-100
##lim=0.03-0.035
#for i in range(5):
#    locals()['p'+str(i)]=gen(d1,1000)

#def cardinality(floor,ceiling,order,weight):
#    length=len(order)
#    for i in range(length):
#        count=0
#        summ=0
#        for j in range(len(order[i])):
#            summ+=weight[i][j]
#            count+=1
#            if summ>=1:
#                break
#        while count>ceiling:
#            weight[i]=np.array(weight[i])+0.0001
#            for j in range(len(order[i])):
#                summ+=weight[i][j]
#                count+=1
#                if summ>=1:
#                    break
#        while count<floor:
#            weight[i]=np.array(weight[i])-0.0001
#            for j in range(len(order[i])):
#                summ+=weight[i][j]
#                count+=1
#                if summ>=1:
#                    break
#    return weight


def trans1(order,weight):
    length=len(order)
    gene={}
    size={}
    gene_weight={}
    gene_order={}
    pf_len=rd.randint(55,70)
    for i in range(length):
        count=0
        summ=0
        for j in range(len(order[i])):
            summ+=weight[i][j]
            count+=1
            if summ>=1:
                break
        locals()['g_o'+str(i)]=order[i][:count]    
        gene_order[i]=locals()['g_o'+str(i)]
        temp1=weight[i][:count]
        temp=[]
        for k in range(len(temp1)):
            temp.append(temp1[k]/sum(temp1))
        locals()['g_w'+str(i)]=temp    
        gene_weight[i]=locals()['g_w'+str(i)]
        size[i]=count
    return (gene_order,gene_weight,size)

def trans(order,weight,date):
    length=len(order)
    gene={}
    size={}
    gene_weight={}
    gene_order={}
    for i in range(length):
        count=0
        summ=0
        pf_len_start=round(len(order[i])/3)
        pf_len=rd.randint(50,70)
        locals()['g_o'+str(i)]=order[i][pf_len_start:pf_len_start+pf_len]    
        gene_order[i]=locals()['g_o'+str(i)]
        temp1=weight[i][pf_len_start:pf_len_start+pf_len]
        temp=[]
        for k in range(len(temp1)):
            temp.append(temp1[k]/sum(temp1))
        locals()['g_w'+str(i)]=temp    
        gene_weight[i]=locals()['g_w'+str(i)]
        size[i]=pf_len
    for n in gene_weight:
        for i in range(len(gene_weight[n])):
            if gene_weight[n][i]-date.BENCH_WEIGHT.iloc[gene_order[n][i]]>0.05:
                to_be_changed=gene_weight[n][i]
                gene_weight[n][i]=date.BENCH_WEIGHT.iloc[gene_order[n][i]]+0.03
                for k in range(len(gene_weight[n])):
                    if k!=i:
                        gene_weight[n][k]=gene_weight[n][k]+(to_be_changed-date.BENCH_WEIGHT.iloc[gene_order[n][i]]-0.03)/(size[n]-1)
    return (gene_order,gene_weight,size)
    
def createcovm1(date,cov):
    length=len(date)
    covm=np.zeros(shape=(length,length))
    count=-1
    
    for i in range(length):
        for j in np.arange(i,length):
            count+=1
            covm[i][j]=cov.VALUE[count]
    covm+=covm.T
    for i in range(length):
            covm[i][i]=covm[i][i]/2
    covm=pd.DataFrame(covm,index=cov.COLUMN_INDEX[:length],columns=cov.COLUMN_INDEX[:length])
    return covm


#%%   
def turnto01(gene_order,gene_weight,covm,date):
    mtrx_index=list(covm.index)
    gene_order_index={}
    gene_weight_index={}
    gene_alpha_index={}
    for i in gene_order:        
        index=list(date.SEDOL.iloc[np.array(gene_order[i])])
        alpha=np.array(date.ALPHA_SCORE.iloc[gene_order[i]])
        wb=np.array(date.BENCH_WEIGHT.iloc[gene_order[i]])
        d=list(np.array(gene_weight[i])-wb)
#        gene_index=gene_index+[' ']*(len(mtrx_index)-len(gene_index))
        order_index=[]
        weight_index=[]
        alpha_index=[]
        for j in mtrx_index:
            if j in index:
                order_index.append(1)
                weight_index.append(d[index.index(j)])
                alpha_index.append(alpha[index.index(j)])
            else:
                order_index.append(0)
                weight_index.append(0)
                alpha_index.append(0)
        gene_order_index[i]=np.array(order_index)
        gene_weight_index[i]=np.array(weight_index)
        gene_alpha_index[i]=np.array(alpha_index)
    return gene_order_index,gene_weight_index,gene_alpha_index

def violate(gene_order,gene_weight,date):
    sector=date.SECTOR.unique()
    mcapq=date.MCAP_Q.unique()
    c_sector={}
    c_mcapq={}
    for i in gene_order:
        c_mcapq[i]={}
        c_sector[i]={}
        for k in sector:
            c_sector[i][k]=[] 
        for l in mcapq:
            c_mcapq[i][l]=[]
        for j in gene_order[i]:
            c_sector[i][date.SECTOR.iloc[j]].append(j)
            c_mcapq[i][date.MCAP_Q.iloc[j]].append(j)
    cstr_mcapq={}
    for i in c_mcapq:
        cstr_mcapq[i]={}
        for j in c_mcapq[i]:
            con_sum=0
            for k in c_mcapq[i][j]:
                con_sum+=gene_weight[i][gene_order[i].index(k)]-date.BENCH_WEIGHT.iloc[k]
            cstr_mcapq[i][j]=con_sum
    cstr_sector={}
    for i in c_sector:
        cstr_sector[i]={}
        for j in c_sector[i]:
            con_sum=0
            for k in c_sector[i][j]:
                con_sum+=gene_weight[i][gene_order[i].index(k)]-date.BENCH_WEIGHT.iloc[k]
            cstr_sector[i][j]=con_sum   
    return cstr_mcapq,cstr_sector        



def turnover(result_template,gene_order,gene_weight_index,date_i,date_j,last_sol):
    n_mr={}
    n_gwi={}
    n_ls={}
    monthly_return=result_template[result_template.DATE==date.DATE.iloc[0]]
    diff_stock=list(set(date_file[date_i-1].SEDOL.unique()) & set(date_file[date_i].SEDOL.unique()))
    index=date.index-date.index[0]
    index2=date_file[date_i-1].index-date_file[date_i-1].index[0]
    for n in gene_order:
        new_mr[n]=[0]*len(diff_stock)
        new_gwi[n]=[0]*len(diff_stock)
        new_ls[n]=[0]*len(diff_stock)
        for i in range(len(diff_stock)):
            n_mr[n][i]=monthly_return.FOUR_WEEKLY_RETURN[monthly_return.SEDOL==diff_stock[i]].values[0]
            n_gwi[n][i]=gene_weight_index[n][index[date.SEDOL==diff_stock[i]][0]]
            n_ls[n][i]=last_sol[index2[date_file[date_i-1].SEDOL==diff_stock[i]][0]] 
    return n_mr,n_gwi,n_ls

last_sol=np.array(result_template.WEIGHTS[result_template.DATE=='01/03/2007'])
#Error:
def evaluate(last_sol,date_i,date_j,gene_order,gene_weight,gene_weight_index,gene_alpha_index,lam,result_template):
##create a new dataframe to contain covariance---covm
    date=date_file[date_i]
    covm=covm_file[date_j]
    evaluation=[]
    prob={}
    cum_prob={}
    c_m,c_s=violate(gene_order,gene_weight,date)
    monthly_return=result_template[result_template.DATE==date_file[date_i-1].DATE.iloc[0]]
    diff_stock1=list(set(date_file[date_i-1].SEDOL.unique()) - set(date_file[date_i].SEDOL.unique()))
    diff_stock2=list(set(date_file[date_i].SEDOL.unique()) - set(date_file[date_i-1].SEDOL.unique()))
    previous_index=date_file[date_i-1].index-date_file[date_i-1].index[0]
    now_index=date_file[date_i].index-date_file[date_i].index[0]
    new_mr=monthly_return
    new_ls=last_sol
    for k in diff_stock1:
        new_mr=new_mr.FOUR_WEEKLY_RETURN[new_mr.SEDOL!=k]
        new_ls=np.delete(new_ls,previous_index[date_file[date_i-1].SEDOL==k].values[0])
    for n in gene_order:
        new_gwi=np.array(gene_weight_index[n])
        for l in diff_stock2:
            new_gwi=np.delete(new_gwi,now_index[date_file[date_i].SEDOL==l].values[0])
            
        turnover=(abs(new_gwi-new_ls*(1+np.array(new_mr)))).sum()
        risks=np.dot(np.dot(gene_weight_index[n].T,covm),gene_weight_index[n])
        returns=lam*np.dot(gene_weight_index[n].T,gene_alpha_index[n]) 
        c11=np.sqrt(risks)
        penalty=0
        penalty2=0
        for k in c_s[n]:
            penalty+=c_s[n][k]-0.1 if c_s[n][k]>=0.1 else 0
        for m in c_m[n]:
            penalty+=c_m[n][m]-0.1 if c_m[n][m]>=0.1 else 0 
        if c11<0.1 and c11>0.05:
            pass
        elif c11<0.05:
            penalty2=0.05-c11
        elif c11>0.1:
            penalty2=c11-0.1
        fitness=risks-returns+100000*penalty+1000*penalty2#+turnover
        evaluation.append(fitness)
    evaluation=np.array(evaluation)
    sum_eva=np.sum(1/evaluation)
    for i in gene_order:
        prob[i]=(1/evaluation[i])/sum_eva
    prob=sorted(prob.items(),key=lambda d:d[1])   
    p,ind=[],[]
    for i in prob:
        p.append(i[1])
        ind.append(i[0])
    p=np.array(p)
    for i in range(len(p)):
        cum_prob[ind[i]]=sum(p[:i+1])
    cum_prob=sorted(cum_prob.items(),key=lambda d:d[1])
    return prob,cum_prob,evaluation


#def evaluate1(date,covm,gene_order,gene_weight,gene_weight_index,lam):
###create a new dataframe to contain covariance---covm
#    evaluation=[]
#    prob={}
#    cum_prob={}
#    for n in gene_order:
#        alpha=np.array(date.ALPHA_SCORE.iloc[gene_order[n]])
#        wb=np.array(date.BENCH_WEIGHT.iloc[gene_order[n]])
#        sub_covm=np.zeros(shape=(len(gene_order[n]),len(gene_order[n])))
#        length=len(gene_order[n])
#        for i in range(length):
#            for j in range(length):
#                row=date.SEDOL.iloc[gene_order[n][i]]
#                col=date.SEDOL.iloc[gene_order[n][j]]
#                sub_covm[i][j]=covm.loc[row,col] 
#        d=np.array(gene_weight[n])-wb
#        profit=np.dot(np.dot(d.T,sub_covm),d)
#        risk=lam* np.dot(d.T,alpha)   
#        fitness=profit-risk
#        evaluation.append(fitness)
#    evaluation=np.array(evaluation)
#    sum_eva=np.sum(1/evaluation)
#    for i in gene_order:
#        prob[i]=(1/evaluation[i])/sum_eva
#    prob=sorted(prob.items(),key=lambda d:d[1])   
#    p,ind=[],[]
#    for i in prob:
#        p.append(i[1])
#        ind.append(i[0])
#    p=np.array(p)
#    for i in range(len(p)):
#        cum_prob[ind[i]]=sum(p[:i+1])
#    cum_prob=sorted(cum_prob.items(),key=lambda d:d[1])
#    return prob,cum_prob,evaluation

#轮盘选择
def select(prob,cum_prob):
    threshold=rd.random()
    global idx
    for i in range(len(cum_prob)):
        if threshold<=cum_prob[i][1] and threshold>=prob[i][1]:
            idx=prob[i][0]
            break
        else:
            continue
    return idx
#error without global

#def selectbest(generation,order,weight,size):
#    #
#    #
#    return (index1,index2)

#def crossover1(order,weight,index1,index2,size): 
#    point=rd.randint(1,round((size[index1]+size[index2])/2))
#    p1,p2=order[index1],order[index2]
#    p11,p12=p1[:point],p1[point:]
#    p21,p22=p2[:point],p2[point:]
#    w1,w2=weight[index1],weight[index2]
#    w11,w12=w1[:point],w1[point:]
#    w21,w22=w2[:point],w2[point:]
#    l1,l2={},{}
#    for i in p2:
#        if i not in p11:
#            p11.append(i)
#            w11.append(w2[p2.index(i)])
#            
#    for j in p1:
#        if j not in p21:
#            p21.append(j) 
#            w21.append(w1[p1.index(j)])
#    return p11,p21,w11,w21

def crossover(order,weight,index1,index2,size): 
    point=rd.randint(round(len(order[index1])/3),round(len(order[index1])*2/3))
    p1,p2=order[index1],order[index2]
    w1,w2=weight[index1],weight[index2]
    p11,p12=p1[:point],p1[point:]
    p21,p22=p2[:point],p2[point:]
    w11,w12=w1[:point],w1[point:]
    w21,w22=w2[:point],w2[point:]
    l1,l2={},{}
    r1,r2={},{}
#   pool=shuffle([i for i in range(1000,3000)])
#   pool1,pool2=pool[:1000],pool[1000:]
    for i in p12:
        if i in p21:
            l1[p21.index(i)]=i
            r1[p21.index(i)]=w12[p12.index(i)]
        else:
            l1[p12.index(i)+1000]=i
            r1[p12.index(i)+1000]=w12[p12.index(i)]
    for j in p22:
        if j in p11:
            l2[p11.index(j)]=j
            r2[p11.index(j)]=w22[p22.index(j)]
        else:
            l2[p22.index(j)+1000]=j
            r2[p22.index(j)+1000]=w22[p22.index(j)]
    s1,s2=sorted(l1.keys()),sorted(l2.keys())
    t1,t2=sorted(r1.keys()),sorted(r2.keys())
    for m in s1:
        p11.append(l1[m])
    for n in s2:
        p21.append(l2[n])
    for m in t1:
        w11.append(r1[m])
    for n in t2:
        w21.append(r2[n])
#    order[index1],order[index2]=p11,p21
    return p11,p21,w11,w21
#weight是否需要crossover

def mutation(po):
    point1=rd.randint(1,50)#size(index)
    point2=rd.randint(70,len(po)-1)#len(po)
    po[point1],po[point2]=po[point2],po[point1]  
    return po

#def mutation1(po,pw,index,size):
#    point1=rd.randint(1,50)#size(index)
#    point2=rd.randint(51,len(po))#len(po)
#    temp1=po[point1]
#    po[point1]=po[point2]
#    po[point2]=temp1
#    temp2=pw[point1]
#    pw[point1]=pw[point2]
#    pw[point2]=temp2
#    return po,pw
    
def sector(date):
    sectors=date.SECTOR.unique()
    sector_index={}
    for i in sectors:
        sector_index[i]=date.index[date.SECTOR==i]
    return sector_index
     
def MCAPQ(date):
    MCAPQ=date.MCAP_Q.unique()
    MCAPQ_index={}
    for i in MCAPQ:
        MCAPQ_index[i]=date.index[date.MCAP_Q==i]
    return MCAPQ_index 

def GA(iteration,popsize,pc,pm,date_i,date_j,lam):
    #For the first iteration
    #date is d0
    #cov is cov0
#parameter setting:iteration is set to 60; lam is set to 1
#initiate
    date=date_file[date_i]
    cov=cov_file[date_j]
    pop_order,pop_weight,pool=initiate(date,popsize,0.034)
    iter_eva=[]
    covm=createcovm(date,cov)
    for it in range(iteration):        
        new_pop_order,new_pop_weight={},{}
        pop_gene_order,pop_gene_weight,size=trans(pop_order,pop_weight,date)
        goi,gwi,gai=turnto01(pop_gene_order,pop_gene_weight,covm,date)
        prob,cum_prob,eva=evaluate(last_sol,date_i,date_j,pop_gene_order,pop_gene_weight,gwi,gai,lam,result_template)
        best_gene_index=[prob[-i][0] for i in range(1,13)]
        best_order_gene=[]
        best_weight_gene=[]
        best_of=eva.min()
        for index in best_gene_index:
            best_order_gene.append(pop_order[index])
            best_weight_gene.append(pop_weight[index])
        #determine the propotion
        iter_eva.append(np.average(eva))
        #generate offsprings
        for k in range(int(popsize/2)):
            #selection
            index=[]
            for i in range(2):
                idx=select(prob,cum_prob)
                index.append(idx)
            #selection
            # crossover opreration
            if rd.random()<pc:
                p11,p21,w11,w21=crossover(pop_order,pop_weight,index[0],index[1],size)
                new_pop_order[2*k],new_pop_order[2*k+1]=p11,p21
                new_pop_weight[2*k],new_pop_weight[2*k+1]=w11,w21
            else:
                new_pop_order[2*k],new_pop_order[2*k+1]=pop_order[index[0]],pop_order[index[1]]
                new_pop_weight[2*k],new_pop_weight[2*k+1]=pop_weight[index[0]],pop_weight[index[1]]
            # crossover opreration
            # mutation operation
            for i in range(2):
                if rd.random()<pm:
                    po=mutation(new_pop_order[2*k+i]) 
                    new_pop_order[2*k+i]=po
            # mutation operation
# Apply Elitism            
        pop_order=new_pop_order
        pop_weight=new_pop_weight
        pop_gene_order,pop_gene_weight,size=trans(pop_order,pop_weight,date)
        goi,gwi,gai=turnto01(pop_gene_order,pop_gene_weight,covm,date)
        prob,cum_prob,eva=evaluate(last_sol,date_i,date_j,pop_gene_order,pop_gene_weight,gwi,gai,lam,result_template)
        worst_gene_order=[prob[i][0] for i in range(12)]
        best_idx=0
        for index in worst_gene_order:
            pop_order[index]=best_order_gene[best_idx]
            pop_weight[index]=best_weight_gene[best_idx]
            best_idx+=1
# Elitism ends here       
    gene_order,gene_weight,size=trans(pop_order,pop_weight,date)
    goi,gwi,gai=turnto01(gene_order,gene_weight,covm,date)
    prob,cum_prob,final_eva=evaluate(last_sol,date_i,date_j,gene_order,gene_weight,gwi,gai,lam,result_template)
    iter_eva.append(np.average(final_eva))
#0.033 is limit of weights    
    return gene_order,gene_weight,iter_eva,final_eva

def present(date,cov,gene_order,gene_weight,iter_eva,final_eva):
    covm=createcovm(date,cov)
    index=np.where(final_eva==final_eva.min())
    fig=plt.figure(figsize=(9,9))
    plt.plot([i for i in range(1,len(iter_eva)+1)],iter_eva)
    result={}
#   count=0
    for i in gene_order[index[0][0]]:
        result[i]=date.SEDOL.iloc[i]
    result=pd.DataFrame(pd.Series(result)).reset_index(drop=True)
    result.rename(columns=lambda x:str(x).replace('0','STOCK'),inplace=True)
#2.用 Dict 暴力修改
    result['WEIGHT']=np.array(gene_weight[index[0][0]])
    BENCH_WEIGHT=[]
    BETA=[]
    SECTOR=[]
    MCAPQ=[]
    for i in result.STOCK:
        BENCH_WEIGHT.append(date.BENCH_WEIGHT[date.SEDOL==i].values)
        BETA.append(date.BETA[date.SEDOL==i].values)
        SECTOR.append(date.SECTOR[date.SEDOL==i].values)
        MCAPQ.append(date.MCAP_Q[date.SEDOL==i].values)
    result['BENCH_WEIGHT']=np.array(BENCH_WEIGHT)
    result['BETA']=np.array(BETA)
    result['SECTOR']=np.array(SECTOR)
    result['MCAP_Q']=np.array(MCAPQ)
    result['DIFFERENCE']=result['WEIGHT']-result['BENCH_WEIGHT']
    minimum=[]
    for i in result.index:
        minimum.append(min(result.WEIGHT[i],result.BENCH_WEIGHT[i]))
    result['MINIMUM']=np.array(minimum)
    length=len(result.DIFFERENCE)
    sub_covm=np.zeros(shape=(length,length))
    for i in range(length):
        for j in range(length):
            row=result.STOCK[i]
            col=result.STOCK[j]
            sub_covm[i][j]=covm.loc[row,col] 
    d=np.array(result.DIFFERENCE)
    constraint_11=np.sqrt(np.dot(np.dot(d,sub_covm),d.T))
    constraint_10=1-sum(result.MINIMUM)
    constraint_8=(result.MINIMUM*result.BETA).sum()
    constraint_7,constraint_6=violate(gene_order,gene_weight,date) 
    constraint_9=len(result)
    
    return result,constraint_6,constraint_7,constraint_8,constraint_9,constraint_10,constraint_11

#将Dict转化为DataFrame，如果一个key只对应一个value：
#pd.DataFrame.from_dict(my_dict,orient='index').T
#%%

#Testing Running Time

start=time.time()
go,gw,iter_eva,final_eva=GA(1000,40,0.7,0.02,1,1,1)
end=time.time()
run_time=end-start
print(run_time)    
rs,c6,c7,c8,c9,c10,c11=present(d1,cov1,go,gw,iter_eva,final_eva)
# Testing section

#%%
def run(date_file,cov_file,iteration,i):
        go,gw,iter_eva,final_eva=GA(iteration,30,0.8,0.01,date_i,date_j,1)
        rs,c8,c10,c11=present(date_file[date_i],cov_file[date_j],go,gw,iter_eva,final_eva)
                


#2016-07-06 
#BYT3MK1  




def to_result(rs,date,result_template):
    monthly_return=result_template.FOUR_WEEKLY_RETURN[result_template.DATE==date.DATE.iloc[0]]
    sol=[0]*len(date)
    index=date.index-date.index[0]    
    for i in rs.STOCK:
        sol[index[date.SEDOL==i][0]]=rs.WEIGHT[rs.STOCK==i].values[0]        
    return sol
        
    
    
    
    
    
    
    
    
