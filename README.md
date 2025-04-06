
DeepPPI
=======
The paper : [DeepPPI: Boosting Prediction of Protein–Protein Interactions with Deep Neural Networks](https://pubs.acs.org/doi/full/10.1021/acs.jcim.7b00028/)

A deep learning-based framework for predicting protein-protein interactions (PPIs) using handcrafted biological features and fully connected neural networks.



### Project Overview

Protein-protein interactions are critical for understanding biological systems and disease mechanisms. DeepPPI builds a binary classifier for PPIs using a feature-rich representation of protein sequences. The model combines:

- Amino acid composition (AAC)
- Dipeptide composition (DPC)
- Physicochemical descriptors (composition, transition, distribution)
- Sequence-order features (QSO and SOC)
- Amphiphilic pseudo amino acid composition (APAAC)


### Project Structure



### Features Used for Sequence Encoding (Total = 1164 features)


| Feature Type                          | Dimensions |
|--------------------------------------|------------|
| Amino Acid Composition (AAC)         | 20         |
| Dipeptide Composition (DPC)          | 400        |
| CTD Descriptors (24 x 21 features)   | 504        |
| Quasi Sequence Order (QSO)           | 100        |
| Sequence Order Coupling (SOC)        | 60         |
| Amphiphilic Pseudo AAC (APAAC)       | 80         |
| **TOTAL**                            | **1164**   |


### Model Architecture


Each protein passes through its own 3-layer FC encoder:
- Linear → BatchNorm → LeakyReLU → Dropout (×3)
- Encoded vectors are concatenated
- Output passed through another FC + Sigmoid for binary classification

Layers:
- FC: Linear layers with decreasing size (1164 → 512 → 256 → 128)
- Dropout: 0.2 dropout rate
- Activation: LeakyReLU
- Output: Sigmoid on 2-class output


## How to Use


1. Clone the repo:
   git clone https://github.com/kianaseraj/DeepPPI.git
   cd DeepPPI

2. Install required packages:
   pip install numpy torch scikit-learn matplotlib tqdm protpy

3. Prepare the data:
   - Input `.npy` files for PPI pairs: [[prot1, prot2, label], ...]
   - Extract sequence features using `feature_generator.py`

4. Run training:
   python main.py

### Evaluation Metrics


All metrics are defined in `metrics.py`:

- Accuracy
- Precision
- Recall (Sensitivity)
- F1 Score

### Note on Dataset Privacy and Project Context

The DeepPPI model in this repository was implemented as part of a broader benchmarking project. The aim was to demonstrate how information leakage can occur in protein-protein interaction datasets, due to high sequence similarity between proteins across training and test sets. Our work highlights the importance of properly splitting datasets to avoid artificially inflated performance.

The code provided here benchmarks DeepPPI against other models using our custom dataset. However, due to privacy constraints, the dataset is not publicly available in this repository.



