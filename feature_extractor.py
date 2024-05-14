import protpy
from propy import PyPro
import math
import numpy as np


#feature1: Amino Acid Composition, 20 dimensional;
def AAC(seq_string):
    amino_acid_composition = protpy.amino_acid_composition(seq_string)
    return list((amino_acid_composition.values)[0])

#feeature2: Dipeptide Composition, 400 dimensional;
def DPC(seq_string):
    dipeptide_composition = protpy.dipeptide_composition(seq_string)
    return dipeptide_composition.values


#c1,c2,c3,t1,t2,t3,d1,d2,d3,....
#feature3: local sequence descriptors, 504 dimsional(84*6);
# composition, transition, and distribution of amino acids' hydrophobicity, normalized van der waals volume, polarity, and polarizability
hydrophobicity = {"1": "RKEDQN", "2": "GASTPHY", "3": "CLVIMFW"}
Normalized_van_der_waals_vol = {"1" : "GASTPCD", "2" : "NVEQIL", "3":"MHKFRYW"}
polarity = {"1": "LIFWCMVY", "2" : "PATGS", "3" : "HQRKNED"}
polarizability = {"1" : "GASDT", "2": "CPNVEQIL", "3" : "KMHFRYW"}

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



#totally 21 parameters
#composition:percentage of each class, 3 features;
def Composition_extractor(encoding):
  composition_embedd = [encoding.count("1")/len(encoding), encoding.count("2")/len(encoding), encoding.count("3")/len(encoding), *[0]*15]
  return composition_embedd



#T is the frequency of a certain class followed by another class(shifts bettween two classes), 3 features
def Transition_extractor(encoding):
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
    transition_embedd.extend([0] * (15))
    return transition_embedd

#D is the percentage of a sequence segment within which the first, 25%, 50%, 75%, 100% of a certain class appears

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

    # Append zeros to make the list length 80
    result.extend([0] * (75))

    return result


#partitioning sequences to 6 local segments:
def sequence_partition(sequence):
    length = len(sequence)
    part_length = round(length / 4)  # Integer division to get equal parts
    seq1_4 = sequence[:part_length]
    seq2_4 = sequence[part_length:2*part_length]
    seq3_4 = sequence[2*part_length:3*part_length]
    seq4_4 = sequence[3*part_length:]
    seq1_2 = sequence[0:round(length/2)]
    seq2_2 = sequence[round(length/2):length+1]
    return seq1_4, seq2_4, seq3_4, seq4_4, seq1_2, seq2_2


def local_descriptor(sequence):

  #encode each seq partition in 4 property groups:
  lc1_hydrophobicity = hydrophobicity_descriptor(sequence_partition(seq)[0])
  lc1_van_der_waals = Normalized_van_der_waals_vol_descriptor(sequence_partition(seq)[0])
  lc1_polarity = polarity_descriptor(sequence_partition(seq)[0])
  lc1_polarizability = polarizability_descriptor(sequence_partition(seq)[0])

  lc2_hydrophobicity = hydrophobicity_descriptor(sequence_partition(seq)[1])
  lc2_van_der_waals = Normalized_van_der_waals_vol_descriptor(sequence_partition(seq)[1])
  lc2_polarity = polarity_descriptor(sequence_partition(seq)[1])
  lc2_polarizability = polarizability_descriptor(sequence_partition(seq)[1])

  lc3_hydrophobicity = hydrophobicity_descriptor(sequence_partition(seq)[2])
  lc3_van_der_waals = Normalized_van_der_waals_vol_descriptor(sequence_partition(seq)[2])
  lc3_polarity = polarity_descriptor(sequence_partition(seq)[2])
  lc3_polarizability = polarizability_descriptor(sequence_partition(seq)[2])

  lc4_hydrophobicity = hydrophobicity_descriptor(sequence_partition(seq)[3])
  lc4_van_der_waals = Normalized_van_der_waals_vol_descriptor(sequence_partition(seq)[3])
  lc4_polarity = polarity_descriptor(sequence_partition(seq)[3])
  lc4_polarizability = polarizability_descriptor(sequence_partition(seq)[3])

  lc5_hydrophobicity = hydrophobicity_descriptor(sequence_partition(seq)[4])
  lc5_van_der_waals = Normalized_van_der_waals_vol_descriptor(sequence_partition(seq)[4])
  lc5_polarity = polarity_descriptor(sequence_partition(seq)[4])
  lc5_polarizability = polarizability_descriptor(sequence_partition(seq)[4])

  lc6_hydrophobicity = hydrophobicity_descriptor(sequence_partition(seq)[5])
  lc6_van_der_waals = Normalized_van_der_waals_vol_descriptor(sequence_partition(seq)[5])
  lc6_polarity = polarity_descriptor(sequence_partition(seq)[5])
  lc6_polarizability = polarizability_descriptor(sequence_partition(seq)[5])

  #
  

  lc1_composition = composition_extractor(lc1_hydrophobicity) + composition_extractor(lc1_van_der_waals) + composition_extractor(lc1_polarity) + composition_extractor(lc1_polarizability) 
  lc2_composition = composition_extractor(lc2_hydrophobicity) + composition_extractor(lc2_van_der_waals) + composition_extractor(lc2_polarity) + composition_extractor(lc2_polarizability)
  lc3_composition = composition_extractor(lc3_hydrophobicity) + composition_extractor(lc3_van_der_waals) + composition_extractor(lc3_polarity) + composition_extractor(lc3_polarizability)
  lc4_composition = composition_extractor(lc4_hydrophobicity) + composition_extractor(lc4_van_der_waals) + composition_extractor(lc4_polarity) + distribution_extractor(lc4_polarizability)

  lc1_transition = transition_extractor(lc1_hydrophobicity) + transition_extractor(lc1_van_der_waals) + transition_extractor(lc1_polarity) + transition_extractor(lc1_polarizability) 
  lc2_transition = transition_extractor(lc2_hydrophobicity) + transition_extractor(lc2_van_der_waals) + transition_extractor(lc2_polarity) + transition_extractor(lc2_polarizability)
  lc3_transition = transition_extractor(lc3_hydrophobicity) + transition_extractor(lc3_van_der_waals) + transition_extractor(lc3_polarity) + transition_extractor(lc3_polarizability)
  lc4_transition = transition_extractor(lc4_hydrophobicity) + transition_extractor(lc4_van_der_waals) + transition_extractor(lc4_polarity) + distribution_extractor(lc4_polarizability)

  lc1_distribution = distribution_extractor(lc1_hydrophobicity) + distribution_extractor(lc1_van_der_waals) + distribution_extractor(lc1_polarity) + distribution_extractor(lc1_polarizability) 
  lc2_distribution = distribution_extractor(lc2_hydrophobicity) + distribution_extractor(lc2_van_der_waals) + distribution_extractor(lc2_polarity) + distribution_extractor(lc2_polarizability)
  lc3_distribution = distribution_extractor(lc3_hydrophobicity) + distribution_extractor(lc3_van_der_waals) + distribution_extractor(lc3_polarity) + distribution_extractor(lc3_polarizability)
  lc4_distribution = distribution_extractor(lc4_hydrophobicity) + distribution_extractor(lc4_van_der_waals) + distribution_extractor(lc4_polarity) + distribution_extractor(lc4_polarizability)


import protFeat 
from protFeat.feature_extracter import extract_protein_feature
#feature4: Sequence Order
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

