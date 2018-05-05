'''
compute_util.py
    - Contains utilities for preprocessing and special cost function
'''
import copy
# numpy
import numpy as np
# import Keras backend
from keras import backend as K

# sum up last dimension of ATAC-seq into one
def coarse_normalize_input(atac_seq):
    print "coarse normalizing..."
    o_shape = atac_seq.shape
    return atac_seq.sum(axis=2).reshape((o_shape[0], o_shape[1], 1))


# normalizes intervals to window-sized bp bins with summit at center
# non_inclusive
def normalize_interval(interval, window_size):
    normalized_interval = copy.deepcopy(interval)
    summit = int(interval.start) + int(interval[-1])
    normalized_interval.start = summit-(window_size-1)/2
    normalized_interval.end = summit+(window_size-1)/2 + 1
    if normalized_interval.start <= 0 or normalized_interval.end <= normalized_interval.start:
        return None
    return normalized_interval


# function that transforms output target 
def double_log_transform(d):
    return np.log(1.0+np.log(1.0+d))


# https://chat.stackoverflow.com/rooms/156491/discussion-between-julio-daniel-reyes-and-eleanora
# https://stackoverflow.com/questions/46619869/how-to-specify-the-correlation-coefficient-as-the-loss-function-in-keras
def pearson_loss(y_true, y_pred):
    x = y_true
    y = y_pred
    mx = K.mean(x)
    my = K.mean(y)
    xm, ym = x-mx, y-my
    r_num = K.sum(np.multiply(xm,ym))
    r_den = K.sqrt(np.multiply(K.sum(K.square(xm)), K.sum(K.square(ym))))
    r = r_num / r_den
    r = K.maximum(K.minimum(r, 1.0), -1.0)
    return 1 - K.square(r)
