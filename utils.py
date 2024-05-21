
from constants import *
import math
import numpy as np
import protpy 
from protFeat.feature_extracter import extract_protein_feature

#feature1: Amino Acid Composition, 20 features
def AAC(sequence):
    amino_acid_composition = protpy.amino_acid_composition(sequence)
    aac = [i/100 for i in amino_acid_composition.values[0]]
    return aac

#feeature2: Dipeptide Composition, 400 features
def DPC(sequence):
    dipeptide_composition = protpy.dipeptide_composition(sequence)
    dpc = [i/100 for i in dipeptide_composition.values[0]]
    return dpc
#feature3: composition, transition, and distribution of physiochemical properties, 504 features
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
def composition_descriptor(encoding):
  composition_embedd = [encoding.count("1")/len(encoding), encoding.count("2")/len(encoding), encoding.count("3")/len(encoding)]
  return composition_embedd

#Transition(shifts between two classes), 3 featuresdef transition_descriptor(encoding):
def transition_descriptor(encoding):
    transitions_ = {}
    transitions = ["12", "13", "21", "23", "31", "32"]

    for i in range(len(encoding) - 1):
        transition = encoding[i:i+2]
        if transition in transitions:
            transitions_[transition] = transitions_.get(transition, 0) + 1

    # Combine counts for transitions with the same group
    combined_transitions = {"12": 0, "13": 0, "23": 0}
    for transition, count in transitions_.items():
        if transition in ["12", "21"]:
            combined_transitions["12"] += count
        elif transition in ["13", "31"]:
            combined_transitions["13"] += count
        else:  # transition in ["23", "32"]
            combined_transitions["23"] += count

    total_length = len(encoding)
    if total_length > 1:
        transition_embedd = [
            combined_transitions["12"] / (total_length - 1),
            combined_transitions["13"] / (total_length - 1),
            combined_transitions["23"] / (total_length - 1)
        ]
    else:
        transition_embedd = [0.0, 0.0, 0.0]

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
def distribution_descriptor(input_string):
    class_lengths = {}
    classes = ["1", "2", "3"]
    for target_class in classes:
        class_lengths[target_class] = calculate_class_lengths(input_string, target_class)

    # Normalize lengths by the total length of the input string
    normalized_lengths = {key: {k: (v / len(input_string)) for k, v in lengths.items()} for key, lengths in class_lengths.items()}

    # Convert normalized lengths to a single list
    result = []
    for target_class in classes:
        result.extend(normalized_lengths[target_class].values())

    return result

#feature4:
#quasi sequence order, 100 dimensions
def QSOD(input_folder, fasta_file_name):
    descriptor = extract_protein_feature("QSOrder" , input_folder, fasta_file_name, place_protein_id = 0)
    return descriptor

#sequence order coupling, 60 dimensions
def SOCD(input_folder, fasta_file_name):
    descriptor = extract_protein_feature("SOCNumber", 0, input_folder, fasta_file_name, place_protein_id = 0)
    return descriptor

#feature5: Amphiphilic Pseudo Amino Acid Composition, 80 diminesions
def APAAC(sequence):
  apaac = protpy.amphiphilic_pseudo_amino_acid_composition(sequence,lamda= 30, weight = 0.5)
  return apaac