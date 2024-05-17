import math
import numpy as np
import protFeat
import protpy
import protFeat 
from protFeat.feature_extracter import extract_protein_feature

#feature1: Amino Acid Composition, 20 features
def AAC(seq_string):
    amino_acid_composition = protpy.amino_acid_composition(seq_string)
    return list((amino_acid_composition.values)[0])

#feeature2: Dipeptide Composition, 400 features
def DPC(seq_string):
    dipeptide_composition = protpy.dipeptide_composition(seq_string)
    return dipeptide_composition.values


#feature3: composition, transition, and distribution of physiochemical properties, 504 features
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


def hydrophobicity_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in hydrophobicity.items() if aa in val)
  return encoding

def Normalized_van_der_waals_vol_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in Normalized_van_der_waals_vol.items() if aa in val)
  return encoding

def polarity_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in polarity.items() if aa in val)
  return encoding

def polarizability_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in polarizability.items() if aa in val)
  return encoding

def charge_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in charge.items() if aa in val)
  return encoding

def secondary_structure_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in secondary_structure.items() if aa in val)
  return encoding

def solvent_accessible_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in solvent_accessible.items() if aa in val)
  return encoding

def surface_tension_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in surface_tension.items() if aa in val)
  return encoding

def prot_prot_hotspot_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in prot_prot_hotspot.items() if aa in val)
  return encoding
def prot_prot_propensity_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in prot_prot_propensity.items() if aa in val)
  return encoding
def prot_dna_propensity_Schneider_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in prot_dna_propensity_Schneider.items() if aa in val)
  return encoding
def prot_dna_propensity_Ahmad_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in prot_dna_propensity_Ahmad.items() if aa in val)
  return encoding
def prot_RNA_propensity_Kim_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in prot_RNA_propensity_Kim.items() if aa in val)
  return encoding
def prot_RNA_propensity_Ellis_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in prot_RNA_propensity_Ellis.items() if aa in val)
  return encoding
def prot_RNA_propensity_Phipps_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in prot_RNA_propensity_Phipps.items() if aa in val)
  return encoding

def prot_ligand_propensity_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in prot_ligand_propensity.items() if aa in val)
  return encoding


def prot_ligand_valid_propensity_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in prot_ligand_valid_propensity.items() if aa in val)
  return encoding

def prot_ligand_polar_propensity_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in prot_ligand_polar_propensity.items() if aa in val)
  return encoding

def molecular_weight_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in molecular_weight.items() if aa in val)
  return encoding

def cLogP_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in cLogP.items() if aa in val)
  return encoding

def hydrogen_bond_donor_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in hydrogen_bond_donor.items() if aa in val)
  return encoding

def hydrogen_bond_acceptor_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in hydrogen_bond_acceptor.items() if aa in val)
  return encoding


def water_Solubility_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in water_Solubility.items() if aa in val)
  return encoding

def Amino_acid_flexibility_descriptor(sequence):
  encoding = "".join(key for aa in sequence for key, val in Amino_acid_flexibility.items() if aa in val)
  return encoding

#composition(frequency of each defined class), 3 features
def composition_extractor(encoding):
  composition_embedd = [encoding.count("1")/len(encoding), encoding.count("2")/len(encoding), encoding.count("3")/len(encoding)]
  return composition_embedd

#Transition(shifts between two classes), 3 features
def transition_extractor(encoding):
    transitions_ = {}
    transitions = [
        "12",
        "13",
        "21",
        "23",
        "31",
        "32"
    ]

    for i in range(len(encoding) - 1):
        transition = encoding[i:i+2]
        if transition in transitions:
            transitions_[transition] = transitions_.get(transition, 0) + 1

    # Combine counts for transitions with the same group
    combined_transitions = {}
    for transition, count in transitions_.items():
        if transition in ["12", "21"]:
            group_key = "12"
        elif transition in ["13", "31"]:
            group_key = "13"
        else:
            group_key = "23"
        combined_transitions[group_key] = combined_transitions.get(group_key, 0) + count

    total_length = len(encoding)
    transition_embedd = [combined_transitions[key] / total_length for key in combined_transitions]

    return transition_embedd

#Distribution(the percentage of a sequence segment within which the first, 25%, 50%, 75%, 100% of a certain class appears), 15 features

def calculate_class_lengths(input_string, target_class):
    total_length = len(input_string)
    class_positions = [i for i, char in enumerate(input_string) if char == target_class]
    class_count = len(class_positions)
    lengths = {
        '0%': 0,
        '25%': 0,
        '50%': 0,
        '75%': 0,
        '100%': 0
    }

    if class_count > 0:
        lengths['0%'] = class_positions[0] + 1
        lengths['25%'] = class_positions[min(math.ceil(class_count * 0.25), class_count) - 1] + 1
        lengths['50%'] = class_positions[min(math.ceil(class_count * 0.50), class_count) - 1] + 1
        lengths['75%'] = class_positions[min(math.ceil(class_count * 0.75), class_count) - 1] + 1
        lengths['100%'] = class_positions[-1] + 1

    return lengths

#15 features
def distribution_extractor(input_string):
    class_lengths = {}
    classes = ["1", "2", "3"]
    for target_class in classes:
        class_lengths[target_class] = calculate_class_lengths(input_string, target_class)

    # Normalize lengths by the total length of the input string
    normalized_lengths = {key: {k: (v / len(input_string))*100 for k, v in lengths.items()} for key, lengths in class_lengths.items()}

    # Convert normalized lengths to a single list
    result = []
    for target_class in classes:
        result.extend(normalized_lengths[target_class].values())

    return result

#feature4:
#quasi sequence order, 100 dimensions
def QSOD(input_folder, fasta_file_name):
    descriptor = extract_protein_feature(QSOrder , input_folder, fasta_file_name, place_protein_id = 0)
    return descriptor

#sequence order coupling, 60 dimensions
def SOCD(input_folder, fasta_file_name):
    descriptor = extract_protein_feature(SOCNumber, 0, input_folder, fasta_file_name, place_protein_id = 0)
    return descriptor

#feeature5: Amphiphilic Pseudo Amino Acid Composition, 80 diminesions
def APAAC(seq_string):

    APAAC_embedding = list((PyPro.GetProDes(seq).GetAPAAC(lamda=30, weight=0.5)).values())
    return APAAC_embedding

