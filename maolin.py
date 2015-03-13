#!/usr/bin/env python
from openeye.oechem import *
from openeye.oegraphsim import *
import sys
import os
import csv
import numpy as np
from sklearn import cluster
import math

# The functions based on OETKs
def OeReadfile(filename):
    """The function returns the ifs"""
    ifs = oemolistream()
    if not ifs.open(filename):
        OEThrow.Fatal("Unable to open %s" % filename)
    return ifs
def OeWritefile(filename):
    """The function returns the ofs"""
    ofs = oemolostream()
    if not ofs.open(filename):
        OEThrow.Fatal("Unable to create %s" % filename)
    return ofs
def OeManipulateFP(inputfile,outputfile,formats='sdf',fingerprint='false',addsame='false'):
    """Calculate the hex of binary MACCS and add an extra field
       inputfile can be any format
       outputfile should be sdf or csv
       formats is for outputfile
       fingerprint can be MACCSH or MACCSB
       addsame is the content of the extra field, this is for unsupervised learning. No specific meaning"""
    ifs=OeReadfile(inputfile)
    if formats=='sdf':
        ofs=OeWritefile(outputfile)
    if formats=='csv':
        csvfile=open(outputfile,'w')
        if fingerprint=='MACCSH':
            csvfile.write('MACCSH')
            if fingerprint=='MACCSB' or addsame!='false':
                csvfile.write(',')
        if fingerprint=='MACCSB':
            for i in range(0,166):
                csvfile.write(str(i+1))
                if i!=166:
                    csvfile.write(',')
                elif addsame!='false':
                    csvfile.write(',')
        if addsame!='false':
            csvfile.write(str(addsame))
        csvfile.write('\n')
    for mol in ifs.GetOEGraphMols():
        if fingerprint!='false':
            fp=OEFingerPrint()
            OEMakeFP(fp, mol, OEFPType_MACCS166)
            if fingerprint=='MACCSH':
                if formats=='sdf':
                    OESetSDData(mol, 'MACCS', str(fp.ToHexString()))
                if formats=='csv':
                    csvfile.write(str(fp.ToHexString()))
                    if fingerprint=='MACCSB' or addsame!='false':
                        csvfile.write(',')
            if fingerprint=='MACCSB':
                for bitnum in range(len(fp.ToHexString())-1):
                    four=str(mybin(int(fp.ToHexString()[bitnum],16)))
                    if bitnum!=(len(fp.ToHexString())-2):
                        for i in range(1,5):
                            if formats=='sdf':
                                OESetSDData(mol, str(bitnum*4+i),four[i-1])
                            if formats=='csv':
                                csvfile.write(str(four[i-1]))
                                csvfile.write(',')
                    else:
                        for i in range(3,5):
                            if formats=='sdf':
                                OESetSDData(mol, str(bitnum*4+i-2),four[i-1])
                            if formats=='csv':
                                csvfile.write(str(four[i-1]))
                                if (bitnum*4+i-2)!=166:
                                    csvfile.write(',')
                                elif addsame!='false':
                                    csvfile.write(',')
        if addsame!='false':
            if formats=='sdf':
                OESetSDData(mol, 'Added',str(addsame))
            if formats=='csv':
                csvfile.write(str(addsame))
        if formats=='sdf':
            OEWriteMolecule(ofs, mol)
        if formats=='csv':
            csvfile.write('\n')
    if formats=='csv':
        csvfile.close()
def OePickupsdf_from_list(alist,inputfile,outputfile):
    """Pick up the compounds from the given list"""
    ifs=OeReadfile(inputfile)
    ofs=OeWritefile(outputfile)
    index=0
    for mol in ifs.GetOEGraphMols():
        if index in alist:
            OEAddSDData(mol, 'index', str(index+1))
            OEWriteMolecule(ofs, mol)
        index+=1
def Probability_acid_amine(acidfile,aminefile,outputfile):
    """Reture the docked acid and amine whose distance is less than 2 for the C-N bond. The clash has been considered"""
    ifsacid=OeReadfile(acidfile)
    ifsamine=OeReadfile(aminefile)
    ofs=OeWritefile(outputfile)
    acidsmarts = "[CX3](=O)[OX1H0-,OX2H1]"
    qacid = OEQMol()
    OEParseSmarts(qacid, acidsmarts)
    acidss = OESubSearch()
    acidss.Init(qacid)
    aminesmarts = "[NH2;+0;!$(NC=O)]"
    #aminesmarts = "[NH20]"
    qamine = OEQMol()
    OEParseSmarts(qamine, aminesmarts)
    aminess = OESubSearch()
    aminess.Init(qamine)
    acidnum=0
    for acid in ifsacid.GetOEGraphMols():
        clash = OENearestNbrs(acid, float(1))
        acidnum+=1
        print('Progressing %d'% acidnum)
        for amine in ifsamine.GetOEGraphMols():
            nbrs = clash.GetNbrs(amine)
            if nbrs.IsValid():
                pass
            else:
                for acidmatch in acidss.Match(acid):
                    for aminematch in aminess.Match(amine):
                        for acidatom in acidmatch.GetTargetAtoms():
                            if acidatom.GetAtomicNum()==6:
                                for amineatom in aminematch.GetTargetAtoms():
                                    if amineatom.GetAtomicNum()==7:
                                        dist = OEGetDistance2(acid, acidatom, amine, amineatom)
                                        dist = math.sqrt(dist)
                                        if dist<2:                                 
                                            OEWriteMolecule(ofs, acid)
                                            OEWriteMolecule(ofs, amine)
        ifsamine=OeReadfile(aminefile)

# The functions based on scikit-learn
def LoadForScikit(inputfile):
    """The inputfile must be in csv format with header. """
    result=[]
    filetemp=open(inputfile,'r')
    i=0
    for line in filetemp:
        i+=1
    filetemp.close()
    data_file = csv.reader(open(inputfile))
    temp = next(data_file)
    n_samples = int(i-1)
    n_features = int(len(temp)-1)
    data = np.empty((n_samples, n_features))
    target = np.empty((n_samples,), dtype=np.int)
    for i, ir in enumerate(data_file):
        data[i] = np.asarray(ir[:-1], dtype=np.float)
        target[i] = np.asarray(ir[-1], dtype=np.int)
    result.append(data)
    result.append(target)
    return tuple(result)
def Cluster_Mkmeans(data,a):
    """data is from LoadForScikit(inputfile)
       a is the number of clusters
       init_size should be 3*number of compounds"""
    X = data[0]
    Y = data[1]
    M_kmeans = cluster.MiniBatchKMeans(n_clusters=a,init_size=40000)
    result=M_kmeans.fit_predict(X)
    return result

# The functions independent
def mybin(num):
    bstr = bin(num)
    l = (len(bstr) - 2) % 4
    if l > 0:
        bstr = ('0'*(4-l)) + bstr[2:]
    else:
        bstr = bstr[2:]
    return bstr
def pickup_from_clusters(clusters,inputfile,outputfile):
    """Pick up the compounds from the predicted clusters"""
    tempclu=[]
    tempindex=[]
    for i,ir in enumerate(clusters):
        if ir not in tempclu:
            tempindex.append(i)
            tempclu.append(ir)
    OePickupsdf_from_list(tempindex,inputfile,outputfile)
        
