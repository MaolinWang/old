#!/usr/bin/env python
from maolin import *
if len(sys.argv) != 2:
    OEThrow.Usage("%s <input>" % sys.argv[0])
fileinput=sys.argv[1]
file_maccs_class=fileinput[:-4]+'_maccs_class.csv'

OeManipulateFP(fileinput,file_maccs_class,formats='csv',fingerprint='MACCSB',addsame='0')



