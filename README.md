
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
- Specificity
- F1 Score
- Matthews Correlation Coefficient (MCC)
- AUROC (Area Under ROC Curve)
- AUPRC (Area Under Precision-Recall Curve)
- MSE (Mean Squared Error)


### Example Training Output


Epoch 5/30
Train loss: 0.0312 - Train accuracy: 91.34%
Val loss: 0.0279 - Val accuracy: 92.47%


### To-Do


- [x] Feature extraction pipeline
- [x] Model training & metrics
- [ ] Add `argparse` CLI interface
- [ ] Save model checkpoints
- [ ] ROC & PR Curve visualization
- [ ] Modular `src/` structure
- [ ] Write full `requirements.txt`


