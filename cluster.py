#!/usr/bin/env python
from maolin import *
if len(sys.argv) != 3:
    OEThrow.Usage("%s <input> <output>" % sys.argv[0])
fileinput=sys.argv[1]
file_maccs_class=fileinput[:-4]+'_maccs_class.csv'
fileoutput=sys.argv[2]
OeManipulateFP(fileinput,file_maccs_class,formats='csv',fingerprint='MACCSB',addsame='0')
data=LoadForScikit(file_maccs_class)
clusters=Cluster_Mkmeans(data,10)
pickup_from_clusters(clusters,fileinput,fileoutput)


