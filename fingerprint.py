#!/usr/bin/env python
from openeye.oechem import *
from openeye.oegraphsim import *
import sys
def OECosine(fpA, fpB):
    onlyA, onlyB, bothAB, neitherAB = OEGetBitCounts(fpA, fpB)
    yule = float(bothAB)
    yule /= float((onlyA + bothAB) * (onlyB + bothAB))
    return yule
def OEDice(fpA, fpB):
    onlyA, onlyB, bothAB, neitherAB = OEGetBitCounts(fpA, fpB)
    yule = float(2 * bothAB)
    yule /= float((onlyA + onlyB + 2 * bothAB))
    return yule
if len(sys.argv) != 4:
    OEThrow.Usage("%s <input> <reference> <output>" % sys.argv[0])
ifs = oemolistream()
ofs = oemolostream()
ifsr = oemolistream()
if not ifs.open(sys.argv[1]):
    OEThrow.Fatal("Unable to open %s" % sys.argv[1])
if not ifsr.open(sys.argv[2]):
    OEThrow.Fatal("Unable to open %s" % sys.argv[2])
if not ofs.open(sys.argv[3]):
    OEThrow.Fatal("Unable to create %s" % sys.argv[3])
if ofs.GetFormat() != OEFormat_SDF:
    OEThrow.Fatal("%s output file has to be an SDF file" % sys.argv[3])
fp = OEFingerPrint()
fpmaccsr = OEFingerPrint()
fplingor = OEFingerPrint()
fpcircularr = OEFingerPrint()
fppathr = OEFingerPrint()
fptreer = OEFingerPrint()
for mol in ifsr.GetOEGraphMols():
    OEMakeFP(fpmaccsr, mol, OEFPType_MACCS166)
    fptypestrmaccsr = fpmaccsr.GetFPTypeBase().GetFPTypeString()
    fphexdatamaccsr = fpmaccsr.ToHexString()
    
    OEMakeFP(fplingor, mol, OEFPType_Lingo)
    fptypestrlingor = fplingor.GetFPTypeBase().GetFPTypeString()
    fphexdatalingor = fplingor.ToHexString()

    OEMakeFP(fpcircularr, mol, OEFPType_Circular)
    fptypestrcircularr = fpcircularr.GetFPTypeBase().GetFPTypeString()
    fphexdatacircularr = fpcircularr.ToHexString()

    OEMakeFP(fppathr, mol, OEFPType_Path)
    fptypestrpathr = fppathr.GetFPTypeBase().GetFPTypeString()
    fphexdatapathr = fppathr.ToHexString()

    OEMakeFP(fptreer, mol, OEFPType_Tree)
    fptypestrtreer = fptreer.GetFPTypeBase().GetFPTypeString()
    fphexdatatreer = fptreer.ToHexString()
    
for mol in ifs.GetOEGraphMols():
    OEMakeFP(fp, mol, OEFPType_MACCS166)
    fptypestrmaccs = fp.GetFPTypeBase().GetFPTypeString()
    fphexdatamaccs = fp.ToHexString()
    tanimotomaccs=OETanimoto(fpmaccsr, fp)
    consinemaccs=OECosine(fpmaccsr, fp)
    dicemaccs=OEDice(fpmaccsr, fp)
    
    OEMakeFP(fp, mol, OEFPType_Lingo)
    fptypestrlingo = fp.GetFPTypeBase().GetFPTypeString()
    fphexdatalingo = fp.ToHexString()
    tanimotolingo=OETanimoto(fplingor, fp)
    consinelingo=OECosine(fplingor, fp)
    dicelingo=OEDice(fplingor, fp)

    OEMakeFP(fp, mol, OEFPType_Circular)
    fptypestrcircular = fp.GetFPTypeBase().GetFPTypeString()
    fphexdatacircular = fp.ToHexString()
    tanimotocircular=OETanimoto(fpcircularr, fp)
    consinecircular=OECosine(fpcircularr, fp)
    dicecircular=OEDice(fpcircularr, fp)

    OEMakeFP(fp, mol, OEFPType_Path)
    fptypestrpath = fp.GetFPTypeBase().GetFPTypeString()
    fphexdatapath = fp.ToHexString()
    tanimotopath=OETanimoto(fppathr, fp)
    consinepath=OECosine(fppathr, fp)
    dicepath=OEDice(fppathr, fp)

    OEMakeFP(fp, mol, OEFPType_Tree)
    fptypestrtree = fp.GetFPTypeBase().GetFPTypeString()
    fphexdatatree = fp.ToHexString()
    tanimototree=OETanimoto(fptreer, fp)
    consinetree=OECosine(fptreer, fp)
    dicetree=OEDice(fptreer, fp)
    
    OESetSDData(mol, fptypestrmaccs, fphexdatamaccs)
    OESetSDData(mol, fptypestrlingo, fphexdatalingo)
    OESetSDData(mol, fptypestrcircular, fphexdatacircular)
    OESetSDData(mol, fptypestrpath, fphexdatapath)
    OESetSDData(mol, fptypestrtree, fphexdatatree)
    OESetSDData(mol, "tanimotomaccs", str(tanimotomaccs))
    OESetSDData(mol, "consinemaccs", str(consinemaccs))
    OESetSDData(mol, "dicemaccs", str(dicemaccs))
    OESetSDData(mol, "tanimotolingo", str(tanimotolingo))
    OESetSDData(mol, "consinelingo", str(consinelingo))
    OESetSDData(mol, "dicelingo", str(dicelingo))
    OESetSDData(mol, "tanimotocircular", str(tanimotocircular))
    OESetSDData(mol, "consinecircular", str(consinecircular))
    OESetSDData(mol, "dicecircular", str(dicecircular))
    OESetSDData(mol, "tanimotopath", str(tanimotopath))
    OESetSDData(mol, "consinepath", str(consinepath))
    OESetSDData(mol, "dicepath", str(dicepath))
    OESetSDData(mol, "tanimototree", str(tanimototree))
    OESetSDData(mol, "consinetree", str(consinetree))
    OESetSDData(mol, "dicetree", str(dicetree))
    OEWriteMolecule(ofs, mol)
