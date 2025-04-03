# -*- coding: utf-8 -*-
"""
metrics.py

Custom evaluation metrics for binary classification.
Includes: accuracy, precision, recall, specificity, F1, MCC, AUC-ROC, AUPRC.

Author: Kiana Seraj
"""

import numpy as np
import math
from sklearn.metrics import roc_auc_score, average_precision_score


def pred_to_classes(predicted, threshold=0.5):
    """
    Converts probabilities to binary predictions.

    Args:
        predicted (list or np.array): Model predictions.
        threshold (float): Classification threshold.

    Returns:
        list: Binary predictions (0 or 1)
    """
    return [1 if p >= threshold else 0 for p in predicted]


def get_confusion_matrix(actual, predicted, threshold=0.5):
    """
    Computes TP, FP, TN, FN given a threshold.

    Returns:
        tuple: (TP, FP, TN, FN)
    """
    predicted_classes = pred_to_classes(predicted, threshold)
    tp = sum((a == 1 and p == 1) for a, p in zip(actual, predicted_classes))
    fp = sum((a == 0 and p == 1) for a, p in zip(actual, predicted_classes))
    tn = sum((a == 0 and p == 0) for a, p in zip(actual, predicted_classes))
    fn = sum((a == 1 and p == 0) for a, p in zip(actual, predicted_classes))
    return tp, fp, tn, fn


def get_accuracy(actual, predicted, threshold=0.5):
    """
    Accuracy = (TP + TN) / Total
    """
    predicted_classes = pred_to_classes(predicted, threshold)
    correct = sum([a == p for a, p in zip(actual, predicted_classes)])
    return correct / len(actual)


def precision(actual, predicted, threshold=0.5):
    """
    Precision = TP / (TP + FP)
    """
    tp, fp, _, _ = get_confusion_matrix(actual, predicted, threshold)
    return tp / (tp + fp) if (tp + fp) != 0 else 0.0


def recall(actual, predicted, threshold=0.5):
    """
    Recall (Sensitivity) = TP / (TP + FN)
    """
    tp, _, _, fn = get_confusion_matrix(actual, predicted, threshold)
    return tp / (tp + fn) if (tp + fn) != 0 else 0.0


def specificity(actual, predicted, threshold=0.5):
    """
    Specificity = TN / (TN + FP)
    """
    _, fp, tn, _ = get_confusion_matrix(actual, predicted, threshold)
    return tn / (tn + fp) if (tn + fp) != 0 else 0.0


def f1_score(actual, predicted, threshold=0.5):
    """
    F1 Score = 2 * (Precision * Recall) / (Precision + Recall)
    """
    p = precision(actual, predicted, threshold)
    r = recall(actual, predicted, threshold)
    return 2 * p * r / (p + r) if (p + r) != 0 else 0.0


def mcc(actual, predicted, threshold=0.5):
    """
    Matthews Correlation Coefficient (MCC)
    """
    tp, fp, tn, fn = get_confusion_matrix(actual, predicted, threshold)
    numerator = (tp * tn) - (fp * fn)
    denominator = math.sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
    return numerator / denominator if denominator != 0 else 0.0


def auroc(actual, predicted):
    """
    Area Under the ROC Curve (AUC-ROC)
    """
    return roc_auc_score(actual, predicted)


def auprc(actual, predicted):
    """
    Area Under the Precision-Recall Curve (AUPRC)
    """
    return average_precision_score(actual, predicted)


def get_mse(actual, predicted):
    """
    Mean Squared Error (MSE)
    """
    actual = np.array(actual)
    predicted = np.array(predicted)
    return np.mean((actual - predicted) ** 2)


def evaluate_metrics(actual, predicted, threshold=0.5):
    """
    Computes all relevant classification metrics.

    Returns:
        dict: A dictionary with all metrics.
    """
    return {
        "accuracy": get_accuracy(actual, predicted, threshold),
        "precision": precision(actual, predicted, threshold),
        "recall": recall(actual, predicted, threshold),
        "specificity": specificity(actual, predicted, threshold),
        "f1_score": f1_score(actual, predicted, threshold),
        "mcc": mcc(actual, predicted, threshold),
        "auroc": auroc(actual, predicted),
        "auprc": auprc(actual, predicted),
        "mse": get_mse(actual, predicted)
    }
