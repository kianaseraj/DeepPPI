from utils import AAC, DPC, composition_descriptor, transition_descriptor, distribution_descriptor
from constants import *
import protpy
import torch
import numpy as np
import os



#filename = /PATH/to/file.npy
#dir_path = /PATH/to/directory to save features in tensor format


#Min-Max normalization in the range of 0,1
def norm(a):
    min_val = min(a)
    max_val = max(a)
    x_normalized = [(x - min_val) / (max_val - min_val) for x in a]
    return x_normalized


#Min-Max normalization in the range of -1,1
def norm_(x):
  min_val = min(x)
  max_val = max(x)
  normalized_composition = [2 * ((x - min_val) / (max_val - min_val)) - 1]
  return normalized_composition


#composition, distribution,transition descriptors
def CTD(sequence):
  ctd = []
  for i in physiochemical_properties:
    composition =composition_descriptor(i(sequence))
    ctd.extend(composition)
    distribution = distribution_descriptor(i(sequence))
    ctd.extend(distribution)
    transition = transition_descriptor(i(sequence))
    ctd.extend(transition)
  return ctd

#quasi sequence order
def QSOrder(protname):
  with open(f"/fs/pool/pool-schwille-user/Frohn_Bela/Seraj-Kiana/ProteinEvolution/ProteinEvolver/ppi/prot-feat/ProtFeat/src/feature_extraction_output/{protname}_QSOrder.txt") as f:
    qsorder = f.read()
    qsorder = [float(i) for i in qsorder.split("\t")[1:]]
  return qsorder 

#sequence coupling
def SOCNumber(protname):
  with open(f"/fs/pool/pool-schwille-user/Frohn_Bela/Seraj-Kiana/ProteinEvolution/ProteinEvolver/ppi/prot-feat/ProtFeat/src/feature_extraction_output/{protname}_SOCNumber.txt") as f:
    socnumber = f.read().strip()  # Remove leading/trailing whitespace
    soc = [float(x) for x in socnumber.split("\t")[1:]]  # Convert to float
    soc_schneider = soc[:30]  
    soc_grantham = soc[30:]  
    soc_grantham_normalized = norm(soc_grantham)
  return soc_schneider, soc_grantham_normalized


def APAAC(sequence):
  amphiphilic_composition = protpy.amphiphilic_pseudo_amino_acid_composition(sequence,lamda= 30, weight = 0.5)
  apaac = list(amphiphilic_composition.values[0])
  apaac = norm_(apaac)
  return apaac

tot = []
for filename in os.listdir("/fs/pool/pool-schwille-user/Frohn_Bela/Seraj-Kiana/ProteinEvolution/ProteinEvolver/ppi/Sequence"):
  sequence = (np.load(f"/fs/pool/pool-schwille-user/Frohn_Bela/Seraj-Kiana/ProteinEvolution/ProteinEvolver/ppi/Sequence/{filename}")).item()
  protname = filename.split(".")[0]
  aac = AAC(sequence)
  dpc = DPC(sequence)
  ctd = CTD(sequence)
  qsorder = QSOrder(protname)
  soc_schneider, soc_grantham_normalized = SOCNumber(protname)
  amphiphilic_composition = protpy.amphiphilic_pseudo_amino_acid_composition(sequence,lamda= 30, weight = 0.5)
  apaac = list(amphiphilic_composition.values[0])
  feature = aac + dpc + ctd + qsorder + soc_schneider + soc_grantham_normalized + apaac
  torch.save(feature,filename)
