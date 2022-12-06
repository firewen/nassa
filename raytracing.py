#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 14:09:42 2021

@author: wenj
"""

import numpy as np

err = 1e-5

class VM:
    pass

def ray_out(vmodel,r0,sdep,wtype):
#compute the travel time and ray parameter for the first arrival phase
#vmodel : the velocity model
#r0 : the epicentral distance
#sdep : source depth
#wtype : phase type ('P' or 'S')
    [p,atime] = draw_ray(vmodel,r0,sdep,wtype)
    
    return(p,atime)
    
def model_info(modeldata):
    c = np.loadtxt(modeldata)
    nlayer = len(c)
    z = c[:,0]
    vp = c[:,2]
    vs = c[:,1]
    
    h = np.zeros(nlayer)
    h[0:nlayer-1] = np.diff(z)
    
    vmodel = VM()
    vmodel.z = z
    vmodel.h = h
    vmodel.vp = vp
    vmodel.vs = vs
    
    return(vmodel)
    
def draw_ray(vmodel,theta,zs,wtype):
    # for direct wave
    [p,atime] = direct(vmodel, theta, zs, wtype)
    # print(p,atime)
    # for refract wave
    [pn, antime] = refract(vmodel, theta, zs, wtype)
    # print(pn,antime)
    # for reflect wave
    # [pr, artime] = reflect(vmodel,theta,zs,wtype)
    # print(pr,artime)
    
    if pn.size == 0:
        fp = p[0]
        ftime = atime[0]
    else:
        ap = np.concatenate((p,pn))
        alltime = np.concatenate((atime,antime))
        
        ftime = np.min(alltime)
        idx = np.argmin(alltime)
        
        fp = ap[idx]
        
    return(fp,ftime)
#    return(p,atime)
    
def direct(vmodel, theta, zs, phase):
    s = sourcelayer(vmodel,zs)
    z = vmodel.z
    h = vmodel.h
    if phase == 'P' :
        v = vmodel.vp
    else:
        v = vmodel.vs
    
    # for source on the ground
    if s == 0 and zs == 0:
        p = np.nan
        atime = theta/v[0]
        return(p,atime)
        
    p = reflect_layered(vmodel,theta,zs,'d',phase)
    
    the = np.zeros((s+1,2))
    ttime = np.zeros(s+1)
    
    if s == 0:
        the[0,0] = theta
        the[0,1] = theta-p*v[0]*(zs-z[s])/np.sqrt(1-(p*v[0])**2)
        ttime[0] = np.sqrt((the[0,1]-the[0,0])**2+(zs-z[s])**2)/v[0]
    else:
        the[0,0] = theta
        the[0,1] = theta-p*v[0]*h[0]/np.sqrt(1-(p*v[0])**2)
        ttime[0] = np.sqrt((the[0,1]-the[0,0])**2+h[0]*h[0])/v[0]
        for i in range(1,s):
            the[i,0] = the[i-1,1]
            the[i,1] = the[i,0]-p*h[i]*v[i]/np.sqrt(1-(p*v[i])**2)
            ttime[i] = np.sqrt((the[i,1]-the[i,0])**2+h[i]*h[i])/v[i]
        the[s,0] = the[s-1,1]
#####################################################################
        #the[s,1] = the[s,0]-p*(zs-z[s])*v[s]/np.sqrt(1-(p*v[s])**2)
#####################################################################
        ttime[s] = np.sqrt((the[s,1]-the[s,0])**2+(zs-z[s])**2)/v[s]   
        
    atime = np.sum(ttime)
    
    p = np.array([p])
    atime = np.array([atime])
    return(p,atime)
    
    
def reflect(vmodel,theta,zs,phase):
    s = sourcelayer(vmodel,zs)
    z = vmodel.z
    h = vmodel.h
    n = len(h)
    
    if s == n:
        print('There is no reflected wave!')
        p = np.array([])
        atime = np.array([])
        return(p,atime)
        
    if phase == 'P' :
        # for P
        p = reflect_layered(vmodel,theta,zs,'r',phase)
        
        if p.size == 1 and p == 0.0:
            print('There are no reflected P waves for this source location')
            return(p,atime)
            
        v = vmodel.vp
        
    else:
        # for S
        p = reflect_layered(vmodel,theta,zs,'r',phase)
        
        if p.size == 1 and p == 0.0:
            print('There are no reflected S waves for this source location')
            return(p,atime)
        
        v = vmodel.vs
    
    atime = np.zeros(len(p))    
    for i in range(0,len(p)):
        m = s+1+i
        the = np.zeros((2*m-s,2))
        ttime = np.zeros(2*m-s)
        
        the[0,0] = theta
        the[0,1] = theta-p[i]*v[0]*h[0]/np.sqrt(1-(p[i]*v[0])**2)
        ttime[0] = np.sqrt((the[0,1]-the[0,0])**2+h[0]*h[0])/v[0]
        
        for k in range(1,m):
            the[k,0] = the[k-1,1]
            the[k,1] = the[k,0]-p[i]*h[k]*v[k]/np.sqrt(1-(p[i]*v[k])**2)
            ttime[k] = np.sqrt((the[k,1]-the[k,0])**2+h[k]*h[k])/v[k]
        
        j = m
        for k in range(m-1,s,-1):
            the[j,0] = the[j-1,1]
            the[j,1] = the[j,0]-p[i]*h[k]*v[k]/np.sqrt(1-(p[i]*v[k])**2)
            ttime[j] = np.sqrt((the[j,1]-the[j,0])**2+h[k]*h[k])/v[k]
            j = j+1
            
        the[2*m-s-1,0] = the[2*m-s-2,1]
        the[2*m-s-1,1] = the[2*m-s-2,1]-p[i]*(z[s+1]-zs)*v[s]/np.sqrt(1-(p[i]*v[s])**2)
        ttime[2*m-s-1] = np.sqrt((the[2*m-s-1,1]-the[2*m-s-1,0])**2+(z[s+1]-zs)**2)/v[s]
        
        atime[i] = np.sum(ttime)
        
    return(p,atime)
    
def refract(vmodel,theta,zs,phase):
    s = sourcelayer(vmodel,zs)
    z = vmodel.z
    h = vmodel.h
    n = len(h)
    ns = s+1
    
    if s+1 == n:
        print('There is no refracted wave!')
        pn = np.array([])
        traveltime = np.array([])
        return(pn,traveltime)
        
    if phase == 'P' :
        [p, lidx] = layered_2p_n(vmodel,theta,zs,phase)
        
        if p.size == 1 and p == 0.0:
            print('There are no refracted P waves for this source location')
            pn = np.array([])
            traveltime = np.array([])
            return(pn,traveltime)
            
        v = vmodel.vp
        
    else:
        [p, lidx] = layered_2p_n(vmodel,theta,zs,phase)
        
        if p.size == 1 and p == 0.0:
            print('There are no refracted S waves for this source location')
            pn = np.array([])
            traveltime = np.array([])
            return(pn,traveltime)
            
        v = vmodel.vs
        
    atime = np.zeros(len(p))
    flag = np.zeros(len(p))
    for i in range(0,len(p)):
        m = lidx[i]+1
        the = np.zeros((2*m-ns+2,2))
        ttime = np.zeros(2*m-ns+2)
        the[0,0] = theta
        the[0,1] = theta-p[i]*v[0]*h[0]/np.sqrt(1-(p[i]*v[0])**2)
        ttime[0] = np.sqrt((the[0,1]-the[0,0])**2+h[0]*h[0])/v[0]
        for k in range(1,m):
            the[k,0] = the[k-1,1]
            if p[i]*v[k] > 1:	# p *v should be less than 1
                flag[i] = 1
                break
            the[k,1] = the[k,0]-p[i]*v[k]*h[k]/np.sqrt(1-(p[i]*v[k])**2)
            ttime[k] = np.sqrt((the[k,1]-the[k,0])**2+h[k]*h[k])/v[k]
#            for given p, refract wave may cannot come back to the source
#            point. if the(k,2) is less than 0, it indicates that the
#            raytracing fails.
            if the[k,1] <= 0:
                flag[i] = 1
                break
            
        if flag[i] == 1:
            continue
        
        j = m
        for k in range(m-1,ns-1,-1):
            the[j,0] = the[j-1,1]
            the[j,1] = the[j,0]-p[i]*v[k]*h[k]/np.sqrt(1-(p[i]*v[k])**2)
            ttime[j] = np.sqrt((the[j,1]-the[j,0])**2+h[k]*h[k])/v[k]
            
            if the[j,1] <= 0:
                flag[i] = 1
                break
            
            j = j+1
            
        if flag[i] == 1:
            continue
        
        the[2*m-ns,0] = the[2*m-ns-1,1]
        the[2*m-ns,1] = the[2*m-ns-1,1]-p[i]*(z[s+1]-zs)*v[s]/np.sqrt(1-(p[i]*v[s])**2)
        ttime[2*m-ns] = np.sqrt((the[2*m-ns,1]-the[2*m-ns,0])**2+(z[s+1]-zs)**2)/v[s]
        
        if the[2*m-ns,1] <=0:
            flag[i] = 1
            continue
        
        l_refract = the[2*m-ns,1]
        for k in range(2*m-ns+1,m-1,-1):
            the[k,:] = the[k-1,:]-l_refract
            ttime[k] = ttime[k-1]
        the[m,0] = the[m-1,1]
        the[m,1] = the[m+1,0]
        ttime[m] = l_refract*p[i]
        
        atime[i] = np.sum(ttime)
        
    if len(p) == 1 and p == 0.0 :
        pn = np.array([])
        traveltime = np.array([])
    else:
        idx = np.where(flag == 0)
        pn = p[idx]
        traveltime = atime[idx]
        
    return(pn,traveltime)
    
def sourcelayer(vmodel,zs):
    z = vmodel.z
    nlayer = len(z)
    
    if zs == 0:
        idx = 0
        
    for i in range(1,nlayer):
        if zs > z[i-1] and zs <=z[i]:
            idx = i-1
            break
    
    if zs > z[nlayer-1]:
        idx = nlayer-1
    
    return(idx)

def reflect_layered(vmodel,theta,zs,rd,phase):
    s = sourcelayer(vmodel,zs)
    n = len(vmodel.h)
    
    if phase == 'P' :
        # for P wave
        if rd == 'd' :
            p = layered_2p_d(vmodel,theta,zs,phase)
        else:
            if n == 1 or s == n-1 :
                print('There are no reflected P waves for this source location')
                p = 0.0
                return(p)
                
            p = np.zeros(n-s-1)
            for i in range(s,n-1):
                p[i-s] = layered_2p_r(vmodel,i,theta,zs,phase)
                
    else:
        # for S wave
        if rd == 'd' :
            p = layered_2p_d(vmodel,theta,zs,phase)
        else:
            if n == 1 or s == n-1 :
                print('There are no reflected S waves for this source location')
                p = 0.0
                return(p)
                
            p = np.zeros(n-s)
            for i in range(s,n):
                p[i-s] = layered_2p_r(vmodel,i,theta,zs,phase)
                       
    return(p)
    
# for direct wave, n=s; for reflected wave, z(n) is the reflected plane
# string is either 'd' or 'r'
# for direct P/S wave    
def layered_2p_d(vmodel,theta,zs,phase):
    s = sourcelayer(vmodel,zs)
    if phase == 'P' :
        v = vmodel.vp
    else:
        v = vmodel.vs
    
    h = vmodel.h
    z = vmodel.z
    n = s+1
    
    hh = np.zeros(n)
    hh[0:s] = h[0:s]
    hh[s] = zs-z[s]
    
    M = 0
    for k in range(1,n):
        if v[k] > v[M]:
            M = k
            
    q = initial(n,M,v,hh,theta)
    err1 = np.abs(F(q,n,M,v,hh,theta))
    counter = 0
    while err1 > err :
        q = q-F(q,n,M,v,hh,theta)/F1(q,n,M,v,hh)
        err1 = np.abs(F(q,n,M,v,hh,theta))
        counter = counter+1
        if counter > 200:
            print('for direct wave, q is unstable')
            
    p = q/(v[M]*np.sqrt(1+q*q))
    
    return(p)
    
# for reflect P/S wave
def layered_2p_r(vmodel,i,theta,zs,phase):
    s = sourcelayer(vmodel,zs)
    if phase == 'P' :
        v = vmodel.vp
    else:
        v = vmodel.vs
    
    h = vmodel.h
    z = vmodel.z
    n = i+1
    
    hh = np.zeros(n)
    hh[0:s] = h[0:s]
    #################################
    # s -> s-1 in hh
    hh[s] = z[s+1]-zs+h[s]                        
    hh[s+1:n] = 2*h[s+1:n]
    #################################
    
    M = 0
    for k in range(1,n):
        if v[k] > v[M]:
            M = k
            
    q = initial(n,M,v,hh,theta)
    err1 = np.abs(F(q,n,M,v,hh,theta))
    counter = 0
    while err1 > err :
        q = q-F(q,n,M,v,hh,theta)/F1(q,n,M,v,hh)
        err1 = np.abs(F(q,n,M,v,hh,theta))
        counter = counter+1
        if counter > 200:
            print('for reflected wave, q is unstable')
            
    p = q/(v[M]*np.sqrt(1+q*q))
    
    return(p)

# for refract P/S wave    
def layered_2p_n(vmodel,theta,zs,phase):
    s = sourcelayer(vmodel,zs)
    
    h = vmodel.h
#    z = vmodel.z
    n = len(h)
    
    if n == 1 or s+1 == n:
        print('There are no refracted waves for this source location')
        p = 0.0
        lidx = 0
        return(p,lidx)
    
    if phase == 'P' :
        v = vmodel.vp
    else:
        v = vmodel.vs
        
    p = []
    lidx = []
    counter = 0
    for k in range(s,n-1):
        if v[k] < v[k+1]:
            counter = counter+1
            p.append(1/v[k+1])
            lidx.append(k)
            
    if counter == 0:
        print('There are no refracted waves for this source location')
        p = 0.0
        lidx = 0
    
    p = np.array(p)
    
    return(p,lidx)
    
    
def F(q,n,M,v,hh,theta):
    f = 0.0
    for k in range(0,n):
        f = f+hh[k]*v[k]/np.sqrt(v[M]*v[M]+q*q*(v[M]*v[M]-v[k]*v[k]))
    f = q*f-theta
    
    return(f)

def F1(q,n,M,v,hh):
    f = 0.0
    for k in range(0,n):
        f = f+hh[k]*v[k]/(v[M]*v[M]+q*q*(v[M]*v[M]-v[k]*v[k]))**1.5
    f = f*v[M]*v[M]
    
    return(f)
    
def initial(n,M,v,hh,theta):
    a = 0.0
    b = 0.0
    for k in range(0,n):
        etak = v[k]/v[M]
        a = a+etak*hh[k]/hh[M]
        
        if k != M:
            b = b+etak*hh[k]/np.sqrt(1-etak*etak)
    # print(a,b)            
    thetac = a*b/(a-1.0)
    if theta < thetac:
        q0 = theta/a
    else:
        q0 = theta-b
    
    return(q0)
    
