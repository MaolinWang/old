#!/usr/bin/env python
from openeye.oechem import *
from openeye.oegraphsim import *
import sys
if len(sys.argv) != 3:
    OEThrow.Usage("%s <input> <output>. Only for active/inactive classification" % sys.argv[0])
ifs = oemolistream()
ofs = oemolostream()
if not ifs.open(sys.argv[1]):
    OEThrow.Fatal("Unable to open %s" % sys.argv[1])
if not ofs.open(sys.argv[2]):
    OEThrow.Fatal("Unable to open %s" % sys.argv[2])
if ofs.GetFormat() != OEFormat_CSV:
    OEThrow.Fatal("%s output file has to be an CSV file" % sys.argv[2])
file2=open(sys.argv[2],'w')
fp = OEFingerPrint()
fptype = OEGetFPType("Circular,ver=2.0.0,size=4096,radius=2-2,atype=AtmNum|Arom|Chiral|FCharge|HCount|EqHalo,btype=Order")
#fptype = OEGetFPType(OEFPType_Circular)
unique = True
active=[]
inactive=[]
allfrag=[]
activenum=0
inactivenum=0
for mol in ifs.GetOEGraphMols():
    value = int(OEGetSDData(mol, "classification"))
    if value:
        activenum+=1
        for abset in OEGetFPCoverage(mol, fptype, unique):
            newmol = OEGraphMol()
            adjustH = False
            RGroup = True
            OESubsetMol(newmol, mol, OEIsAtomMember(abset.GetAtoms()), adjustH, RGroup)
            active.append(OECreateIsoSmiString(newmol))
            allfrag.append(OECreateIsoSmiString(newmol))
    else:
        inactivenum+=1
        for abset in OEGetFPCoverage(mol, fptype, unique):
            newmol = OEGraphMol()
            adjustH = False
            RGroup = True
            OESubsetMol(newmol, mol, OEIsAtomMember(abset.GetAtoms()), adjustH, RGroup)
            inactive.append(OECreateIsoSmiString(newmol))
            allfrag.append(OECreateIsoSmiString(newmol))
file2.write("Fragments")
file2.write(",")
file2.write("Active")
file2.write(",")
file2.write("Inactive")
file2.write("\n")
freqactive=0.0
freqinactive=0.0
already=[]
for frag in allfrag:
    mol = OEGraphMol()
    freqactive=float(active.count(frag))/activenum
    freqinactive=float(inactive.count(frag))/inactivenum
    if frag not in already:
        file2.write(frag)
        file2.write(",")
        file2.write(str(freqactive))
        file2.write(",")
        file2.write(str(freqinactive))
        file2.write("\n")
        freqactive=0.0
        freqinactive=0.0
    already.append(frag)
    
    
