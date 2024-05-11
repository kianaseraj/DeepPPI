from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score
import math 
from math import sqrt


#mean squared error: 1/n∑(y-y^)^2
def get_mse(actual, predicted):
    loss = ((actual - predicted) ** 2).mean(axis=0)
    return loss
    

#computing acc = tp/total
def get_accuracy(actual, predicted, threshold):
    correct = 0
    predicted_classes = []
    for prediction in predicted :
      if prediction >= threshold :
        predicted_classes.append(1)
      else :
        predicted_classes.append(0)
    for i in range(len(actual)):
      if actual[i] == predicted_classes[i]:
        correct += 1
    return correct / float(len(actual)) * 100.0


#creating a list of predictions, having a value bigger than 0.5 will assing 1 otherwise 0!
def pred_to_classes(actual, predicted, threshold):
    predicted_classes = []
    for prediction in predicted :
      if prediction >= threshold :
        predicted_classes.append(1)
      else :
        predicted_classes.append(0)
    return predicted_classes
    
#True positive
def get_tp(actual, predicted, threshold):
    predicted_classes = pred_to_classes(actual, predicted, threshold)
    tp = 0
    for i in range(len(predicted_classes)):
      if predicted_classes[i] == 1 and actual[i] == 1:
       tp += 1
    return tp
    
     
#False positive    
def get_fp(actual, predicted, threshold):
    predicted_classes = pred_to_classes(actual, predicted, threshold)
    fp = 0
    for i in range(len(predicted_classes)):
      if predicted_classes[i] == 1 and actual[i] == 0:
       fp += 1
    return fp

#True negative

def get_tn(actual, predicted, threshold):
    predicted_classes = pred_to_classes(actual, predicted, threshold)
    tn = 0
    for i in range(len(predicted_classes)):
      if predicted_classes[i] == 0 and actual[i] == 0:
       tn += 1
    return tn

#False negative
def get_fn(actual, predicted, threshold):
    predicted_classes = pred_to_classes(actual, predicted, threshold)
    fn = 0
    for i in range(len(predicted_classes)):
      if predicted_classes[i] == 0 and actual[i] == 1:
       fn += 1
    return fn
    
    
#precision = TP/ (TP + FP)    
def precision(actual, predicted, threshold):
    prec = get_tp(actual, predicted, threshold) / (get_tp(actual, predicted, threshold) + get_fp(actual, predicted, threshold))
    return prec
    
    
    
   
# sensitivity, recall or true positive rate, recall = TP / (TP + FN)
def sensitivity(actual, predicted, threshold):
    sens = get_tp(actual, predicted, threshold)/ (get_tp(actual, predicted, threshold) + get_fn(actual, predicted, threshold))
    return sens
    

    
#Specificity = TN/(TN+FP)    
def specificity(actual, predicted, threshold):     
   spec =  get_tn(actual, predicted, threshold)/ (get_tn(actual, predicted, threshold) + get_fp(actual, predicted, threshold))
   return spec


#f1 score  = 2 / ((1/ precision) + (1/recall))   
def f_score(actual, predicted, threshold):
    f_sc = 2 / ( (1 / precision(actual, predicted, threshold)) + (1/ sensitivity(actual, predicted, threshold)))
    return f_sc

   
# Matthews correlation coefficient, mcc = (TP * TN - FP * FN) / sqrt((TN+FN) * (FP+TP) *(TN+FP) * (FN+TP)) 
def mcc(actual, predicted, threshold):
   tp = get_tp(actual, predicted, threshold) 
   tn = get_tn(actual, predicted, threshold)
   fp = get_fp(actual, predicted, threshold)
   fn = get_fn(actual, predicted, threshold)
   mcc_met = (tp*tn - fp*fn) / (sqrt((tn+fn)*(fp+tp)*(tn+fp)*(fn+tp)))
   return mcc_met
   
   
#auc-roc curve, a performance measurement for the classification problems at various threshold settings
def auroc(actual, predicted):
   return roc_auc_score(actual, predicted)
  

#The area under the precision-recall curve  
def auprc(actual, predicted):
   return average_precision_score(actual, predicted)
 
   