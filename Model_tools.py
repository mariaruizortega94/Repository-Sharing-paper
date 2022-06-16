#!/usr/bin/env python
# coding: utf-8

import numpy as np
import mpmath
mpmath.dps = 20

def filter_sequences(df):
    
    default_headers = {'Nucleotides': 'cdr3nt', 'cdr3_nt': 'cdr3nt', 'JUNCTION': 'cdr3nt', 'nt':'cdr3nt', 'cdr3nt':'cdr3nt',
                       'Aminoacids': 'cdr3aa', 'cdr3_aa': 'cdr3aa', 'Amino_acids': 'cdr3aa', 'aa': 'cdr3aa', 'cdr3aa':'cdr3aa',
                       'V_CALL': 'v_call', 'v':'v_call', 'v_gene': 'v_call', 'V': 'v_call', 'v_call':'v_call',
                       'J_CALL': 'j_call', 'j':'j_call', 'j_gene': 'j_call', 'J': 'j_call', 'j_call':'j_call'}
    
    df.dropna(inplace=True)
    df = df.rename(columns = default_headers)
    df.drop_duplicates(subset='cdr3aa', keep='first', inplace=True)
    df['junction_length'] = df.cdr3nt.str.len()
    df = df[df['junction_length'] > 6] # Use only longth CDR3
    df = df[df['junction_length'] % 3 == 0] # Get in-frame CDR3s
    df = df[~df['cdr3aa'].str.contains('*', regex=False)] # Get productive CDR3s (aa seqs of the form C...W with no stop codons)
    df = df[df.cdr3aa.str.startswith('C')]
    df = df[df.cdr3aa.str.endswith('W')]
    
    return(df)

def get_shared(df, x_list):
    
    columns = list(df.columns)
    shared = []
    for x in x_list:
        is_shared = False
        if len(x)>1:
            is_shared=True
        
        shared.append(is_shared)
        
    df['shared'] = shared
    df_shared = df.loc[df['shared']==True, columns]
    x_shared = [x for x in x_list if len(x)>1]
    
    return(df_shared, x_shared)


def function_for_p(xs,y,N,f2=1,nb_of_people=10):
    sum0 = 0
    N_sum = sum(N)
    for x in range(1,nb_of_people+1):
        prob = np.exp(-(y/f2) * N[x-1])
        if x in xs:
            sum0 += (N[x-1]/(1-prob))
    return(sum0-N_sum)

def solve_for_p(N,i,x,a=float(1e-55),b=float(1e-6)):
        
    def onevarfunc(y):
        return function_for_p(x,y,N)

    if onevarfunc(a)*onevarfunc(b)<0:
        return (i,brentq(onevarfunc, a, b))
    else:
        p_approx = len(x)/sum(N)
        return (i,p_approx)

def remove_zeros(df, ppost, x):
    
    # Update lists removing those indexes where ppost = 0 
        
    columns = list(df.columns)
    df['ppost'] = ppost
    df = df.loc[df['ppost']!=0.0, columns]

    return(df, ppost, x)

def func(x,q):
    return(x-np.log10(q))

def Likelihood(bin_min, bin_max, x, n):
    prod_x=[]
    bins = np.logspace(bin_min-1, bin_max, base=10, num=20)
    for p in bins:
        prod0 = 1
        for j in range(1,len(n)+1):
            if j in x:
                prod0 = p * n[j-1]*prod0 if p<1e-21 else (1 - float(mpmath.exp(-p * n[j-1])))*prod0
            else:
                prod0 = (1 - (p * n[j-1]))*prod0 if p<1e-21 else float((mpmath.exp(-p* n[j-1])))*prod0
        prod_x.append(prod0)
    return(prod_x, bins)