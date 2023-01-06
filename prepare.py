# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 09:24:35 2022

@author: wenj
"""

import numpy as np
import h5py
import random
import matplotlib.pyplot as plt
from obspy.geodetics import gps2dist_azimuth
from obspy import UTCDateTime, read
from scipy.signal import find_peaks
from scipy.ndimage import shift
from scipy.interpolate import griddata
from sklearn.cluster import DBSCAN
import os
import glob
import time
from pathos.multiprocessing import ProcessingPool as Pool
from raytracing import model_info, ray_out
import sys

def xydist(x,y):
    r0 = gps2dist_azimuth(x[1], x[0], y[1], y[0])[0]/1e3
    dist0 = np.sqrt(r0**2 + (x[2]-y[2])**2)
    return(dist0)

def getinfo(fnm):
    lb = []
    hb = []
    dh = []
    with open(fnm,'r') as reader:
        tmp = reader.readline()
        tmp = reader.readline()
        tmp = reader.readline()
        tt = tmp.strip().split()
        lb.append(float(tt[0]))
        hb.append(float(tt[1]))
        dh.append(float(tt[2]))
        tmp = reader.readline()
        tmp = reader.readline()
        tt = tmp.strip().split()
        lb.append(float(tt[0]))
        hb.append(float(tt[1]))
        dh.append(float(tt[2]))
        tmp = reader.readline()
        tmp = reader.readline()
        tt = tmp.strip().split()
        lb.append(float(tt[0]))
        hb.append(float(tt[1]))
        dh.append(float(tt[2]))
        tmp = reader.readline()
        stfnm = reader.readline().strip("\n")
        tmp = reader.readline()
        vnm = reader.readline().strip("\n")

    return(np.array(lb), np.array(hb), np.array(dh), stfnm, vnm)

def load_st(ddir,stfnm):
    # load station information
    with open(ddir + stfnm,'r') as reader:
        tmp = reader.readlines()
        
    nst = len(tmp)

    rx = np.zeros(nst)
    ry = np.zeros(nst)
    rz = np.zeros(nst)
    stnm = []
    vmodel = []
    for i,line in enumerate(tmp):
        tt = line.strip().split()
        stnm.append(tt[0])
        rx[i] = float(tt[1])
        ry[i] = float(tt[2])
        rz[i] = float(tt[3])
        vmodel.append(model_info(ddir + vnm))
    stcor = np.vstack((rx,ry,rz)).T
    return(nst,stnm,stcor,vmodel)

def datatrans(ddir, stnm, nt):
    nst = len(stnm)
    
    pall = np.zeros([nst, nt])
    # sall = np.zeros([nst, nt])
    
    record = np.zeros([nst, nt])

    for ist, st0 in enumerate(stnm):
        with h5py.File(ddir + 'probs/' + st0 + '.h5', 'r') as f:
            prob = f['prob']
            data = np.squeeze(prob[st0])
            pall[ist,:] = data[0:nt,1] 
            np.savetxt(ddir + 'probs/' + st0 + '.prob',data[0:nt,1])
            
        st = read(glob.glob(ddir + 'data/' + st0 + '.*.HHZ.sac')[0])
        st.detrend("demean")
        record[ist,:] = st[0].data[0:nt]

    return(pall,record)

def findNeighbor(j,X,eps):
    N = []
    for p in range(X.shape[0]):
        temp = np.sqrt(np.sum(np.square(X[j]-X[p])))
        if (temp<=eps):
            N.append(p)
    
    return(N)
            
def dbscan(X,eps,min_Pts):
    k = -1
    NeighborPts = []
    Ner_NeighborPts = []
    fil = []                # 初始时已访问对象列表为空
    gama = [x for x in range(len(X))]
    cluster = [-1 for y in range(len(X))]
    while len(gama) > 0:
        j = random.choice(gama)
        gama.remove(j)
        fil.append(j)
        
        NeighborPts = findNeighbor(j, X, eps)
        if len(NeighborPts) < min_Pts:
            cluster[j] = -1         #标记为噪点
        else:
            k = k+1
            cluster[j] = k
            for i in NeighborPts:
                if i not in fil:
                    gama.remove(i)
                    fil.append(i)
                    Ner_NeighborPts = findNeighbor(i, X, eps)
                    if len(Ner_NeighborPts) >= min_Pts:
                        for a in Ner_NeighborPts:
                            if a not in NeighborPts:
                                NeighborPts.append(a)
                    if (cluster[i]==-1):
                        cluster[i] = k
    
    return(cluster)

def result_plot(ddir,evid,nt,dt,twin1,twin2,pall,record,lb,hb,dh,stnm):

    llen = round(twin2/dt) - round(twin1/dt)                                                                  
    if twin1<0:                                                                                                   
        pc = np.zeros([nst,llen])                                                                                  
        # sc = np.zeros([nst,llen])                                                                                  
        record1 = np.zeros([nst,llen])                                                                            
        pc[:,-round(twin1/dt):] = pall[:,0:round(twin2/dt)]                                                        
    #    sc[:,-round(twin1/dt):] = sall[:,0:round(twin2/dt)]                                                        
        record1[:,-round(twin1/dt):] = record[:,0:round(twin2/dt)]                                                
    else:
        pc = pall[:,round(twin1/dt):round(twin2/dt)]                                                               
    #    sc = sall[:,round(twin1/dt):round(twin2/dt)]                                                               
        record1 = record[:,round(twin1/dt):round(twin2/dt)]

    starttime = time.time()
    data = np.load(ddir + 'outfig/na' + str(evid) + 'bright.npz')
    grid = data['grid']
    vm = data['vm']
    dist0 = data['dist']
    endtime = time.time()
    print('time for iteration: %f' % (endtime-starttime))
   
    # the directory for output 
    outfnm = ddir + 'outfig/' 
    
    idx = np.argmax(vm)
    sx = grid[idx,0]
    sy = grid[idx,1]
    sz = grid[idx,2]
    
    iidx=np.argwhere(vm>np.max(vm)*0.95)
    s_mean = np.mean(grid[iidx,:],axis=0)
    s_std = np.std(grid[iidx,:],axis=0)
    s_mean = np.reshape(s_mean,3)
    s_std = np.reshape(s_std,3)
    print('sx mean : %f, sy mean : %f, sz mean : %f\n' %(np.mean(grid[iidx,0]),np.mean(grid[iidx,1]),np.mean(grid[iidx,2])))
    print('sx std : %f, sy std : %f, sz std : %f\n' %(np.std(grid[iidx,0]),np.std(grid[iidx,1]),np.std(grid[iidx,2])))
    
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111,projection="3d")
    ax.scatter(grid[:,0],grid[:,1],grid[:,2],c=vm,cmap='jet',s=20)
    ax.invert_zaxis()
    plt.savefig(outfnm + 'na' + str(evid) +'_1.png')
    plt.close()
    
    x = np.arange(lb[0],hb[0],dh[0]*4)
    y = np.arange(lb[1],hb[1],dh[1]*4)
    z = np.arange(lb[2],hb[2],dh[2]*2)
    xv,yv,zv = np.meshgrid(x,y,z)
    nnp = np.size(xv)
    xi=np.concatenate([xv.reshape([nnp,1]),yv.reshape([nnp,1]),zv.reshape([nnp,1])],axis=1)
    gridv = griddata(grid,vm,xi,fill_value='nan')

    # approximately compute the confidence parameter
    dist = np.zeros(nnp)
    yy = [sx,sy,sz]
    F = lambda x : xydist(x,yy)
    print(yy,F(xi[0,:]))
    starttime = time.time()
    pool = Pool(8)
    dist = pool.map(F,xi)
    endtime = time.time()
    print('time: %f' %(endtime - starttime))
    
    lidx = np.isnan(gridv)
    iidx = []
    dist1 = []
    gridv1 = []
    for i,l in enumerate(lidx):
        if ~l :
            iidx.append(i)
            dist1.append(dist[i])
            gridv1.append(gridv[i])
    iidx = np.array(iidx)
    dist = np.array(dist1)
    gridv = np.array(gridv1)

    dist_min = np.min(dist)
    dist_max = np.max(dist)
    ndrange = 101
    drange = np.linspace(dist_min, dist_max, ndrange)
    v = np.zeros(ndrange)
    p = np.zeros(ndrange)
    pall = np.sum(np.exp(gridv))
    for i in range(1,ndrange):
        idx = np.argwhere(dist < drange[i])
        v[i] = np.size(idx)/len(dist)
        p[i] = np.sum(np.exp(gridv[idx]))/pall

    fig = plt.figure()
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    plt.subplot(221)
    plt.plot(drange,v,'b',label="P(dist)")
    plt.plot(drange,p,'g',label="V(dist)")
    plt.xlabel('Dist. (km)')
    plt.legend(loc='lower right')
    plt.subplot(223)
    plt.plot(drange,np.concatenate((np.zeros(1),np.diff(v))),'b',label="P(dist)")
    plt.plot(drange,np.concatenate((np.zeros(1),np.diff(p))),'g',label="V(dist)") 
    plt.xlabel('Dist. (km)')
    plt.legend(loc='upper right')
    plt.subplot(222)
    plt.plot(v,p,'b')
    plt.xlabel('V')
    plt.ylabel('P(V)')
    plt.text(0.5,0.5,'Pav = %5.3f'%np.trapz(p,v))
    plt.subplot(224)
    plt.plot(dist,gridv,'b.',label='interpolated')
    plt.plot(dist0,vm,'r.',label='original')
    plt.xlabel('Dist. (km)')
    plt.ylabel('Brightness')
    plt.legend(loc='upper right')
    #plt.show()   
    plt.savefig(outfnm + 'na' + str(evid) +'_4.png')
    plt.close()


    p0 = np.zeros(np.shape(pc))
    tt1 = np.zeros(nst)
    r0 = np.zeros(nst)
    for ist in range(nst):
        r0[ist] = gps2dist_azimuth(stcor[ist,1],stcor[ist,0],sy,sx)[0]/1e3
        pp,T1 = ray_out(vmodel[ist],r0[ist],sz-stcor[ist,2],'P')
        p0[ist,:] = shift(pc[ist,:], -round(T1/dt), cval=0)
        tt1[ist] = T1
    sump = np.sum(p0, axis=0)

    t0 = np.argmax(sump)*dt
    tt0 = t0 + twin1
    print('t0 : %f' %tt0)
    tdif = np.zeros(nst)
    
    dist_idx = np.argsort(r0)
    t = np.arange(0,llen)*dt + twin1

    plt.figure(figsize=(8,10))    
    ax1 = plt.subplot2grid((3,2), (0,0), colspan=2)
    ax1.plot(t, sump, 'b')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Brightness')
    
    ax2 = plt.subplot2grid((3,2), (1,0), rowspan=2)
    ax3 = plt.subplot2grid((3,2), (1,1), rowspan=2)

    for i,ist in enumerate(dist_idx):
        T1 = tt1[ist]
        ax2.plot(t,pc[ist,:]+i,'b')
        ax2.text(t[0],i+0.3,stnm[ist])
        iloc = round((t0+T1)/dt)

        peaks, _ = find_peaks(pc[ist,:],height=0.2,width=5)
        if np.size(peaks) > 0:
            piloc = np.argmin(abs(peaks-iloc))
            ploc = peaks[piloc]
        else:
            ploc = iloc - 100
    
        tdif[ist] = ploc*dt - (t0+T1)
        ax2.plot((tt0+T1)*np.array([1,1]),np.array([0,1])+i,'k--')
        ax3.plot(t,record1[ist,:]/np.max(abs(record1[ist,:]))+i,'b')
        ax3.plot((tt0+T1)*np.array([1,1]),np.array([0,1])+i,'k--')
        ax3.plot((twin1+ploc*dt)*np.array([1,1]),np.array([0,1])+i,'r--')
    
    ax2.set_xlabel('Time (s)')
    ax3.set_xlabel('Time (s)')
    #plt.show()
    plt.savefig(outfnm + 'na' + str(evid) +'_2.png')
    plt.close()
    
    plt.figure()
    plt.plot(r0,tdif,'b.')
    plt.xlabel('Distance (km)')
    plt.ylabel('Arrival Time Diff. (s)')
    #plt.show()
    plt.savefig(outfnm + 'na' + str(evid) +'_3.png')
    plt.close()

    return(np.array([sx,sy,sz]),tt0,s_mean,s_std,np.max(vm),np.trapz(p,x=v)) 


if __name__ == '__main__':
    
    if len(sys.argv)>1 :
        ddir = sys.argv[1]
    else:
        ddir = '2014.08.19.00.00.00/'

    # imput control file
    lb, hb, dh, stfnm, vnm = getinfo(ddir + "input.txt")

    if len(sys.argv) == 4:
        nt = int(sys.argv[2])
        dt = float(sys.argv[3])    
    else:
        nt = 8640000
        dt = 0.01
    
    tmp = ddir.split('/')[0].split('.')
    tb = tmp[0] + '-' + tmp[1] + '-' + tmp[2] + 'T' + tmp[3] + ':' + tmp[4] + ':' + tmp[5]
    tbeg = UTCDateTime(tb)
    print(tb)
    
    nst, stnm, stcor, vmodel = load_st(ddir , stfnm)
    
    p,record = datatrans(ddir, stnm, nt)
    
    tt = []
    st1 = []
    counter = 0
    for ist in range(nst):
        peaks, _ = find_peaks(p[ist,:],height=0.2,width=20)
        tt[counter:counter+len(peaks)] = peaks*dt
        st1[counter:counter+len(peaks)] = p[ist,peaks] + ist
        counter = counter + len(peaks)

    st = np.ones(len(tt))
    tt = np.array(tt)
    st1 = np.array(st1)
    X = np.concatenate((np.reshape(tt,(len(tt),1)),np.reshape(st,(len(tt),1))),axis=1)

    starttime = time.time()
    eps = 3                             # 控制两两台站之间的到时差，单位s
    min_Pts = 5                         # 每台周围的相邻点的个数
    #db = dbscan(X, eps, min_Pts)
    clustering = DBSCAN(eps=eps, min_samples=min_Pts).fit(X)
    db = clustering.labels_
    db = np.array(db)
    endtime = time.time()
    print("time for dbscan: %f"%(endtime-starttime))

    plt.figure()
    plt.scatter(X[:,0],st1,c=db)
    plt.show()

    k = np.max(db) + 1
    numk = np.zeros(k)
    # plt.figure()
    for i in range(k):
        idxt = np.argwhere(db==i)
        numk[i] = len(idxt)
        # plt.plot(X[idxt,0],st1[idxt],'o')
    # plt.show()

    idxk = np.argwhere(numk > 10)        # 每个事件至少有四个台拾取到了有效P初至
    nevent = len(idxk)

#    for i in range(0,nevent):
#        tmpt = X[db==idxk[i],0]
#        twin1 = np.min(tmpt) - 7       # when twin1 is less than 0, twin1 = 0
#        twin2 = np.max(tmpt) + 5
#        plt.plot([twin1,twin2,twin2,twin1,twin1],[-1,-1,nst+1,nst+1,-1],'k--')
#    plt.xlabel('Time (s)')
#    plt.ylabel('Station No.')
#    plt.show()

    print('The number of event condidates %d: ' %(nevent))
   
    tmin = np.zeros(nevent)
    for i in range(nevent):
        tmin[i] = np.min(X[db==idxk[i],0])
    itmp = np.argsort(tmin)
    idxk = idxk[itmp]
    
    #nevent = 1
    csvfile = open(ddir + 'outfig/eventna.csv','w')
    for i in range(0,nevent):
        tmpt = X[db==idxk[i],0]
        twin1 = np.min(tmpt) - 7       # when twin1 is less than 0, twin1 = 0
        twin2 = np.max(tmpt) + 5
        print('eventid is %d, twin: %f %f' %(i,twin1, twin2))
       
        tmp = './nassa/nassa ' + ddir + ' ' + str(twin1) + ' ' + str(twin2) + ' ' + str(nt) + ' ' + str(dt) + ' ' + str(i)
        print(tmp)
        r_v = os.system('./nassa/nassa ' + ddir + ' ' + str(twin1) + ' ' + str(twin2) + ' ' + str(nt) + ' ' + str(dt) + ' ' + str(i))   
        print(r_v) 
        source, tt0, s_mean, s_std, vm, confid =result_plot(ddir,i,nt,dt,twin1,twin2,p,record,lb,hb,dh,stnm)
        evttime = tbeg + tt0
        print('%9.4f %9.4f %9.4f ' % tuple(source))
        print('%9.4f ' % tt0)
        print('%9.4f %9.4f %9.4f' % tuple(s_mean))
        print('%9.4f %9.4f %9.4f' % tuple(s_std))
        csvfile.write('%s ' %(evttime.strftime('%Y %m %d %H %M %S %f')))
        csvfile.write('%9.4f %9.4f %9.4f ' % tuple(source))
        csvfile.write('%9.4f ' % tt0)
        csvfile.write('%9.4f ' % vm)
        csvfile.write('%9.4f ' % confid)
        csvfile.write('%9.4f %9.4f %9.4f' % tuple(s_mean))
        csvfile.write('%9.4f %9.4f %9.4f' % tuple(s_std))
        csvfile.writelines('\n')
    csvfile.close()
