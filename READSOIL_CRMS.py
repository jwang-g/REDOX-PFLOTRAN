
##print encoding for csv file to read in the right format
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def soils():
    pfl_bavg={}
    pfl_poravg={}
    pfl_satavg={}
    path="./MEsoil"     #"./soildata"
    #os.chdir(path)
    # iterate through all file
    for file in os.listdir(path):
        # Check whether file is in text format or not
        if file.endswith(".csv"):
            file_path = f"{path}/{file}"
        # call read text file function
            data=pd.read_csv(file_path,encoding = 'ISO-8859-1')
            stid=data['Station ID']
            day=data['Sample Date (mm/dd/yyyy)']

            watcon=data['Soil Moisture Content (%)']
            sal=data['Soil Salinity (ppt)']
            blkd=data['Bulk Density (g/cm3)']
            TN=data['Total Nitrogen (g/kg)']
            TC=data['Total Carbon (g/kg)']
            wetpH=data['Wet Soil pH (pH units)']
            wetvol=data['Wet Volume (cm3)']
            dryvol=data['Dry Volume (cm3)']
        #j=0
            pdens=2.65 #particle density g/cm3

            por=1-blkd/pdens
            std='2017-01-01'
            edd='2018-01-01'
            if stid[0][0:8]in['CRMS0147','CRMS0157','CRMS0118','CRMS0139']:
                corenm=['-S01','-S02','-S03']
            elif stid[0][0:8]in['CRMS0114','CRMS0131','CRMS0162']:
                corenm=['-S04','-S05','-S06']
            else:
                corenm=['-S01','-S02','-S03','-S04','-S05','-S06']
            stnm=[]
            for n in corenm:
                stnm.append(file[:8]+n)
            satura={}
            ssal={}
            sbulk={}
            sTN={}
            sTC={}
            swetpH={}
            spor={}
            for spec in stnm:
                tmpwt=[]
                tmpsal=[]
                tmpblk=[]
                tmpTN=[]
                tmpTC=[]
                tmppH=[]
                tmppor=[]
                for i in range(len(stid)):
                    if stid[i] == spec:
                        tmpwt+=[watcon[i],]
                        tmpsal+=[sal[i],]
                        tmpblk+=[blkd[i],]
                        tmpTN+=[TN[i],]
                        tmpTC+=[TC[i],]
                        tmppH+=[wetpH[i],]
                        #portm=1-blkd[i]/pdens
                        tmppor+=[por[i],]
                        satura[spec]=tmpwt
                        ssal[spec]=tmpsal
                        sbulk[spec]=tmpblk
                        sTN[spec]=tmpTN
                        sTC[spec]=tmpTC
                        swetpH[spec]=tmppH
                        spor[spec]=tmppor

#calculate soil property mean value for pflotran layers
            bavg=[]
            poravg=[]
            satavg=[]
            for j in range(len(stid)):
                if stid[j] == stnm[0]:
                    tmp1=0
                    tmp2=0
                    tmp3=0
                    for i in range(len(stnm)):#[0,1,2]:
                        tmp1+=blkd[j+i*6]
                        tmp2+=por[j+i*6]
                        tmp3+=watcon[j+i*6]/100
                    bavg+=[tmp1/len(stnm),]
                    poravg+=[tmp2/len(stnm),]
                    satavg+=[tmp3/len(stnm),]
            name=file[:8]
            if name not in pfl_bavg.keys():
                pfl_bavg[name]=bavg
                pfl_poravg[name]=poravg
                pfl_satavg[name]=satavg
    return pfl_bavg, pfl_poravg, pfl_satavg
#BD2pfl=pd.DataFrame(pfl_bavg)
#por2pfl=pd.DataFrame(pfl_bavg)
#sat2pfl=pd.DataFrame(pfl_bavg)
