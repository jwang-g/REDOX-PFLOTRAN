#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 20:56:35 2024

@author: jiaze
"""

import pandas as pd


def soils(tidenv='Little River',latgrd=1):
    
    pdens=2.65 #particle density g/cm3
    
    #read in soil data from maine.
    dis=latgrd*100/500 #lateral grid distance/resolutionn in cm with an arbitray scalinng factor1/500
    
    soil=pd.read_csv('./MEsoil/Vincent_and_Dionne_2023_depthseries.csv')
    
    #names=soil[soil.duplicated(['core_id']) == False]['core_id']

    siteID=soil[soil.duplicated(['site_id']) == False]['site_id']
    
    mshloc=['low','middle','high']
    
    tmp_bavg={}
    pfl_dblk={}
    pfl_por={}
    
    bdavg=[]
    
    ##estimation function between bulk density incretion and distance away from high marsh in cm  
    #Little River: dblk=BD[nm]+dis**2*5.27062778e-04-2.68739059e-03*dis+5.90699047e-01
    #Webhannet: dblk=BD[nm]+dis**2*3.57982576e-04+1.29780990e-03*dis+6.05921139e-01
    #Drakes Island: dblk=BD[nm]+dis**2*4.11859489e-05+1.00679304e-03*dis+6.15661153e-01
        

    for n in siteID:
        bdavg=[]
        for j in [0,5,10,15]:
            cnt=0
            bdtmp=0
            for i in range(0,100):
                inm=n+str(i+1)
                tmp=soil.loc[soil['core_id'] == inm]
                bdtmp+=tmp.loc[tmp['depth_min'] == j]['dry_bulk_density'].values[0]
                cnt=cnt+1
            bdavg+=[bdtmp/cnt,]
        if n not in tmp_bavg.keys():
            tmp_bavg[n]=bdavg
    if tidenv == 'Little River':
        pfl_dblk['high']=tmp_bavg[tidenv]
        incbd=[dis**2*5.27062778e-04-2.68739059e-03*dis+i for i in tmp_bavg[tidenv]]
        pfl_dblk['middle']=[m+n for m,n in zip(tmp_bavg[tidenv],incbd)] 
        pfl_dblk['low']=[m +n*2 for m,n in zip(tmp_bavg[tidenv],incbd)]
    elif tidenv == 'Webhannet':
        pfl_dblk['high']=tmp_bavg[tidenv]
        incbd=[dis**2*3.57982576e-04+1.29780990e-03*dis+i for i in tmp_bavg[tidenv]]
        pfl_dblk['middle']=[m+n for m,n in zip(tmp_bavg[tidenv],incbd)] 
        pfl_dblk['low']=[m+n*2 for m,n in zip(tmp_bavg[tidenv],incbd)]
    elif tidenv == 'Drakes Island':
        pfl_dblk['high']=tmp_bavg[tidenv]
        incbd=[dis**2*4.11859489e-05+1.00679304e-03*dis+i for i in tmp_bavg[tidenv]]
        pfl_dblk['middle']=[m+n for m,n in zip(tmp_bavg[tidenv],incbd)] 
        pfl_dblk['low']=[m +n*2 for m,n in zip(tmp_bavg[tidenv],incbd)]
    
    for i in mshloc:
        if i not in pfl_por.keys():
            pfl_por[i]=[1-m/pdens for m in pfl_dblk[i]]
    ##bulk density and porosity and accumulation rate is from bradley and morris et al., 1990
        #bd=[1.49, 0.79, 0.62, 0.45]
        #por=[43.4, 67.7, 77.0, 96.2]

        #fig, ax=plt.subplots(nrows=1,ncols=1,figsize=(35,45))

        #ax.plot(bd,por)

        #z = np.polyfit(bd, por, 2)

        #zval=[bd[i]**2*50.62-147.97*bd[i]+151.57 for i in range(len(bd))]   
        
    return pfl_dblk,pfl_por