
from utils import *
#physiochemical classes 
hydrophobicity = {"1": "RKEDQN", "2": "GASTPHY", "3": "CLVIMFW"}
Normalized_van_der_waals_vol = {"1" : "GASTPCD", "2" : "NVEQIL", "3":"MHKFRYW"}
polarity = {"1": "LIFWCMVY", "2" : "PATGS", "3" : "HQRKNED"}
polarizability = {"1" : "GASDT", "2": "CPNVEQIL", "3" : "KMHFRYW"}
charge = {"1" : "KR", "2" : "ANCQGHILMFPSTWYV", "3" : "DE"}
secondary_structure = {"1": "EALMQKRH", "2" : "VIYCWFT", "3" : "GNPSD"}
solvent_accessible = {"1": "ALFGCIVW", "2" : "PKQEND", "3" : "MPSTHY"}
surface_tension = {"1" : "GQDNAHR", "2" : "KTSEC", "3" : "ILMFPWYV"}
prot_prot_hotspot = {"1" : "DHIKNPRWY", "2" : "EQSTGAMF", "3" : "CLV"}
prot_prot_propensity = {"1" : "CDFMPQRWY", "2" : "AGHVLNST", "3" : "EIK"}
prot_dna_propensity_Schneider = {"1" : "GKNQRSTY", "2" : "ADEFHILVW", "3" : "CMP"}
prot_dna_propensity_Ahmad= {"1" : "GHKNQRSTY", "2" : "ADEFIPVW", "3" : "CLM"}
prot_RNA_propensity_Kim = {"1" : "HKMRY", "2" : "FGILNPQSVW", "3" : "CDEAT"}
prot_RNA_propensity_Ellis = {"1" : "HGKMRSYW", "2" : "AFINPQT", "3" : "CDELV"}
prot_RNA_propensity_Phipps = {"1" : "HKMQRS", "2" : "ADEFGLNPVY", "3" : "CITW"}
prot_ligand_propensity= {"1" : "CFHWY", "2" : "GILNMSTR", "3" : "AEDKPQV"}
prot_ligand_valid_propensity = {"1" : "CFHWYM", "2" : "DGILNSTV", "3" : "AEKPQR"}
prot_ligand_polar_propensity = {"1" : "DEHRY", "2" : "CFKMNQSTW", "3" : "AGILPV"}
molecular_weight = {"1" : "AGS", "2" : "CDEHIKLMNQPTV", "3" : "FRWY"}
cLogP = {"1" : "RKDNEQH", "2" : "PYSTGACV", "3" : "WMFLI"}
hydrogen_bond_donor = {"1" : "HKNQR", "2" : "DESTWY", "3" : "ACGFILMPV"}
hydrogen_bond_acceptor = {"1" : "DEHNQR", "2" : "KSTWY", "3" : "ACGFILMPV"}
water_Solubility = {"1" : "ACGKRT", "2" : "EFHILMNPQSVW", "3" : "DY"}
Amino_acid_flexibility = {"1" : "EGKNQS", "2" : "ADHIPRTV", "3" : "CFLMWY"}

#defined physiochemical functions in utils module
physiochemical_properties = [
hydrophobicity_descriptor,
Normalized_van_der_waals_vol_descriptor,
polarity_descriptor,
polarizability_descriptor,
charge_descriptor,
secondary_structure_descriptor,
solvent_accessible_descriptor,
surface_tension_descriptor,
prot_prot_hotspot_descriptor,
prot_prot_propensity_descriptor,
prot_dna_propensity_Schneider_descriptor,
prot_dna_propensity_Ahmad_descriptor,
prot_RNA_propensity_Kim_descriptor,
prot_RNA_propensity_Ellis_descriptor,
prot_RNA_propensity_Phipps_descriptor,
prot_ligand_propensity_descriptor,
prot_ligand_valid_propensity_descriptor,
prot_ligand_polar_propensity_descriptor,
molecular_weight_descriptor,
cLogP_descriptor,
hydrogen_bond_donor_descriptor,
hydrogen_bond_acceptor_descriptor,
water_Solubility_descriptor,
Amino_acid_flexibility_descriptor  
]