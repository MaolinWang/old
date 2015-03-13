#!/usr/bin/env python
from openeye.oechem import *
import sys
if len(sys.argv) != 4:
    OEThrow.Usage("%s <txt> <sdf> <output>" % sys.argv[0])
ifssdf = oemolistream()
ofs = oemolostream()
filetxt=open(sys.argv[1],'r')
txtnum=[line.strip() for line in filetxt]
if not ifssdf.open(sys.argv[2]):
    OEThrow.Fatal("Unable to open %s" % sys.argv[2])
if not ofs.open(sys.argv[3]):
    OEThrow.Fatal("Unable to create %s" % sys.argv[3])
for mol in ifssdf.GetOEGraphMols():
    num=OEGetSDData(mol,'num')
    if str(num) in txtnum:
        OEWriteMolecule(ofs, mol)
