import hashlib
from functools import lru_cache
from math import comb

import numpy as np
from scipy.stats import norm


@lru_cache(maxsize=128)
def exact_lookup(n,m):
    n,m=min(n,m),max(n,m)
    limit=n*m//2
    coefficients=np.zeros(limit+1)
    coefficients[0]=1
    for step in range(1,n+1):
        for u in range(step,limit+1):
            coefficients[u]+=coefficients[u-step]
    for shift in range(m+1,m+n+1):
        if shift<=limit:
            coefficients[shift:]-=coefficients[:-shift].copy()
    if (coefficients < -.01).any():
        raise ValueError('Invalid exact U generating polynomial')
    return np.minimum(1,2*np.cumsum(coefficients)/float(comb(n+m,n)))


def mw_p(u,n,m,ties):
    n=np.broadcast_to(np.asarray(n),u.shape)
    total=n+m
    variance=n*m/12*(total+1-ties/(total*(total-1)))
    z=np.divide(np.maximum(u,n*m-u)-n*m/2-.5,np.sqrt(np.maximum(variance,0)),
                out=np.full_like(u,-np.inf,dtype=float),where=variance>0)
    p=np.minimum(1,2*norm.sf(z))
    exact=(np.minimum(n,m)<=8)&(ties==0)
    for size in np.unique(n[exact]):
        where=exact&(n==size)
        minimum=np.minimum(u[where],size*m-u[where])
        if not np.array_equal(minimum,np.rint(minimum)):
            raise ValueError('Untied U must be integral')
        p[where]=exact_lookup(int(size),int(m))[np.rint(minimum).astype(int)]
    return p


def stability(g,q,original_g):
    valid=np.isfinite(g)&np.isfinite(original_g)[None,:]&(original_g[None,:]!=0)
    same=valid&(np.sign(g)==np.sign(original_g)[None,:])
    frequency=(np.isfinite(q)&(q<.05)).mean(0)
    sign=same.mean(0)
    return frequency,sign,(frequency>=.8)&(sign>=.95)&valid.all(0)


def ladder(quotas,pools):
    stage=0
    configurations=[]
    while True:
        sizes={r:min(pools[r],quotas[r]*2**stage) if r!='European' else pools[r] for r in pools}
        if not configurations or sizes!=configurations[-1]:
            configurations.append(sizes)
        if all(sizes[r]==pools[r] for r in pools):
            return configurations
        stage+=1


def digest(weights):
    return hashlib.sha256(np.asarray(weights,dtype='<i4').tobytes()).hexdigest()


class Bank:
    def __init__(self,loaded,core):
        self.core=core
        self.disease,self.reference,self.features,self.definitions,self.audit,self.counts,self.primary,self.sampling,self.sources=loaded
        self.cells=self.counts.loc[self.counts.eligible].sort_values(['disease','canonical_race']).reset_index(drop=True)
        self.total=len(self.primary)
        self.races={}
        self.full_cache={}
        for race in core.RACES:
            ref=self.reference.loc[self.reference.canonical_race.eq(race)].sort_values('file').reset_index(drop=True)
            y=ref[self.features].to_numpy(float)
            cell_indices=self.cells.index[self.cells.canonical_race.eq(race)].to_numpy()
            cases=[self.disease.loc[self.disease.disease.eq(self.cells.disease.iloc[i])&self.disease.canonical_race.eq(race)]
                   .sort_values('patient_id')[self.features].to_numpy(float) for i in cell_indices]
            ns=np.array([len(x) for x in cases])
            contribution=np.empty((len(y),len(cases),len(self.features)),np.float32)
            group_ids=np.empty(y.shape,np.int32)
            offsets=[]
            count=0
            case_ties=np.zeros((len(cases),len(self.features)))
            overlap_rows=[]
            for j in range(len(self.features)):
                unique,inverse=np.unique(y[:,j],return_inverse=True)
                offsets.append(count)
                group_ids[:,j]=inverse+count
                for k,x in enumerate(cases):
                    ordered=np.sort(x[:,j])
                    left=np.searchsorted(ordered,y[:,j],side='left')
                    right=np.searchsorted(ordered,y[:,j],side='right')
                    contribution[:,k,j]=len(x)-(left+right)/2
                    values,frequencies=np.unique(x[:,j],return_counts=True)
                    case_ties[k,j]=np.sum(frequencies**3-frequencies)
                    positions=np.searchsorted(unique,values)
                    hit=positions<len(unique)
                    hit[hit]&=unique[positions[hit]]==values[hit]
                    overlap_rows.extend((k,j,count+pos,int(freq)) for pos,freq in zip(positions[hit],frequencies[hit]))
                count+=len(unique)
            center=y.mean(0)
            centered=y-center
            rows=np.array(overlap_rows,dtype=int).reshape(-1,4)
            self.races[race]={'ref':ref,'y':y,'cases':cases,'ns':ns,'cell_indices':cell_indices,
                'where':(cell_indices[:,None]*len(self.features)+np.arange(len(self.features))).ravel(),
                'contribution':contribution.reshape(len(y),-1),'groups':group_ids,'offsets':np.array(offsets),
                'group_count':count,'unique_reference_features':count==y.size,'case_ties':case_ties,'overlap':rows,
                'xmean':np.array([x.mean(0) for x in cases]),'xvar':np.array([x.var(0,ddof=1) for x in cases]),
                'center':center,'centered':centered,'centered_sq':centered**2}

    def evaluate(self,selections):
        p=np.ones(self.total)
        g=np.full(self.total,np.nan)
        means=[]
        sds=[]
        for race,data in self.races.items():
            indices=np.asarray(selections[race],int)
            if len(indices)<20 or (indices<0).any() or (indices>=len(data['y'])).any():
                raise ValueError('Invalid reference allocation')
            weights=np.bincount(indices,minlength=len(data['y']))
            full=np.all(weights==1)
            if full and race in self.full_cache:
                pp,gg,mean,var=self.full_cache[race]
            else:
                m=len(indices)
                u=(weights.astype(np.float32)@data['contribution']).reshape(len(data['cases']),len(self.features)).astype(float)
                sums=weights@data['centered']
                mean=data['center']+sums/m
                var=(weights@data['centered_sq']-sums**2/m)/(m-1)
                if data['unique_reference_features'] and not len(data['overlap']):
                    ties=data['case_ties']+np.sum(weights**3-weights)
                else:
                    take=np.flatnonzero(weights)
                    freq=np.bincount(data['groups'][take].ravel(),weights=np.repeat(weights[take],len(self.features)),minlength=data['group_count'])
                    ref_ties=np.add.reduceat(freq**3-freq,data['offsets'])
                    ties=data['case_ties']+ref_ties
                    if len(data['overlap']):
                        cell,feature,group,a=data['overlap'].T
                        b=freq[group]
                        np.add.at(ties,(cell,feature),3*a*b*(a+b))
                pp=mw_p(u,data['ns'][:,None],m,ties).ravel()
                gg=self.core.effect_moments(data['ns'][:,None],m,data['xmean'],data['xvar'],mean,var).ravel()
                if full:
                    self.full_cache[race]=(pp.copy(),gg.copy(),mean.copy(),var.copy())
            p[data['where']],g[data['where']]=pp,gg
            means.append(mean)
            sds.append(np.sqrt(np.maximum(var,0)))
        return p,g,np.asarray(means),np.asarray(sds)
