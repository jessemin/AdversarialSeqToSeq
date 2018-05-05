'''
Baseline CNN trained on MSE

Author: Jesik Min
'''

import os, sys
sys.path.append("..")
# Custom file path package
from data import Data_Directories
# Custom util
from utils.compute_util import *
# Package for genomic data
from pybedtools import Interval, BedTool
from genomelake.extractors import ArrayExtractor, BigwigExtractor
# Package for plotting
import matplotlib.pyplot as plt
# Package for correlation
from scipy.stats.stats import pearsonr,spearmanr
# Tensorflow
import tensorflow as tf
# ArgParsing
import argparse
parser = argparse.ArgumentParser()

# Setup Arguments
parser.add_argument('-e',
                    '--num_epochs',
                    required=False,
                    type=int,
                    default=300,
                    dest="num_epochs",
                    help="Number of epochs for training")
parser.add_argument('-w',
                    '--window_size',
                    required=False,
                    type=int,
                    default=10001,
                    dest="window_size",
                    help="Window size for normalized intervals")
parser.add_argument('-d',
                    '--day',
                    required=False,
                    type=str,
                    default='day0',
                    dest="day",
                    help="Target day of data")
parser.add_argument('-f',
                    '--frag',
                    required=False,
                    type=str,
                    default='140',
                    dest="frag",
                    help="Target fragment length of data")
parser.add_argument('-o',
                    '--output',
                    required=False,
                    type=str,
                    default='H3K27ac',
                    dest="histone",
                    help="Target histone mark of data")
parser.add_argument('-m',
                    '--model_path',
                    required=False,
                    type=str,
                    default='',
                    dest="model_path",
                    help="Path to pre-trained model")
parser.add_argument('-save',
                    '--save_dir',
                    required=False,
                    type=str,
                    default='test',
                    dest="save_dir",
                    help="Where to save the best model, logs, and its predictions")
parser.add_argument('-sample_num',
                    '--sample_num',
                    required=False,
                    type=int,
                    default=10000,
                    dest="sample_num",
                    help="Total number of train sample; val sample is 0.2*train_num")
parser.add_argument('-cuda',
                    '--cuda',
                    required=True,
                    type=str,
                    dest="cuda",
                    help="Cuda visible devices")

args = parser.parse_args()

# Parse all arguments

window_size = args.window_size
day = args.day
frag = args.frag
histone = args.histone
model_path = args.model_path
save_dir = args.save_dir
sample_num = args.sample_num
cuda = args.cuda

os.environ["CUDA_VISIBLE_DEVICES"]=cuda

# Logging directories
model_dir = os.path.join("saved_models", save_dir)
log_dir = os.path.join("logs", save_dir)
srv_dir = os.path.join("/srv", "www", "kundaje", "jesikmin", "test_experiments", save_dir)
if not os.path.exists(model_dir):
    os.makedirs(model_dir)
if not os.path.exists(log_dir):
    os.makedirs(log_dir)
if not os.path.exists(srv_dir):
    os.makedirs(srv_dir)

# Train/val/test intervals
DATA_DIR = '/srv/scratch/jesikmin'
train_dir, val_dir, test_dir = os.path.join(DATA_DIR, 'train_interval'),\
                               os.path.join(DATA_DIR, 'val_interval'),\
                               os.path.join(DATA_DIR, 'test_interval')

print train_dir, val_dir, test_dir

# Get train/val/test intervals
train_intervals = list(BedTool(train_dir))
val_intervals = list(BedTool(val_dir))
test_intervals = list(BedTool(test_dir))
print '# of Train Intervals: {}'.format(len(train_intervals))
print '# of Val Intervals: {}'.format(len(val_intervals))
print '# of Test Intervals: {}'.format(len(test_intervals))

# Get input/output data directories
data = Data_Directories()
print data.intervals.keys()
print data.input_atac[day].keys()
print data.output_histone[day].keys()

# Extract input candidates
# Create an ArrayExtractor for ATAC-seq of a given day and specified fragment length
input_candidates = ArrayExtractor(data.input_atac[day][frag])
print 'Finished extracting bigwig for {}, {}bp'.format(day, frag)

# Extract output candiates
# Create a BigWigExtractor for histone mark of a given day
output_candidates = BigwigExtractor(data.output_histone[day][histone])
print 'Finished extracting bigwig for {}, {}'.format(day, histone)

# Normalize train intervals
normalized_train_intervals = [normalize_interval(interval, window_size) for interval in train_intervals if normalize_interval(interval, window_size)]
print 'Finished normalizing train intervals!'
# Normalize val intervals
normalized_val_intervals = [normalize_interval(interval, window_size) for interval in val_intervals if normalize_interval(interval, window_size)]
print 'Finished normalizing val intervals!'
# Normalize test intervals
normalized_test_intervals = [normalize_interval(interval, window_size) for interval in test_intervals if normalize_interval(interval, window_size)]
print 'Finished normalizing test intervals!'

# Fetch intervals of sample_num
normalized_train_intervals = normalized_train_intervals[:sample_num]
normalized_val_intervals = normalized_val_intervals[:int(sample_num*0.2)]
print 'Finished fethcing {} train set and {} val set'.format(sample_num, int(sample_num*0.2))

# Assertions of normalization step
assert (sample_num==len(normalized_train_intervals))
assert (int(sample_num*0.2)==len(normalized_val_intervals))
assert (len(test_intervals)==len(normalized_test_intervals))
# Examples of normalized intervals
print "Examples of original train intervals"
print [(int(_interval.start)+int(_interval[-1]), [int(_interval.start), int(_interval.end)])
       for _interval in train_intervals[:3]]
print "Examples of normalized train intervals with window size of {}".format(window_size)
print [([int(_interval.start), int(_interval.end)])
       for _interval in  normalized_train_intervals[:3]]

# Prune intervals that don's make sense
def prune_invalid_intervals(intervals, bigwig_file):
    for _interval in intervals[:]:
        try:
            bigwig_file([_interval])
        except:
            intervals.remove(_interval)
            pass

# Prune train intervals
print "Before pruning train: {}".format(len(normalized_train_intervals))
prune_invalid_intervals(normalized_train_intervals, input_candidates)
print "After pruning train: {}".format(len(normalized_train_intervals))
# Prune val intervals
print "Before pruning val: {}".format(len(normalized_val_intervals))
prune_invalid_intervals(normalized_val_intervals, input_candidates)
print "After pruning val: {}".format(len(normalized_val_intervals))
# Prune test intervals
print "Before pruning test: {}".format(len(normalized_test_intervals))
prune_invalid_intervals(normalized_test_intervals, input_candidates)
print "After pruning test: {}".format(len(normalized_test_intervals))


X_train = input_candidates(normalized_train_intervals)
print "Finished fetching X_train"
X_val = input_candidates(normalized_val_intervals)
print "Finished fetching X_val"
X_test = input_candidates(normalized_test_intervals)
print "Finished fetching X_test"
print X_train.shape, X_val.shape, X_test.shape

print "Dimension of ATAC-seq signal (input): {}".format(X_train[0].shape)

# Replace nan values with zeros
y_train = np.nan_to_num(output_candidates(normalized_train_intervals))
print "Finished fetching y_train"
y_val = np.nan_to_num(output_candidates(normalized_val_intervals))
print "Finished fetching y_val"
y_test = np.nan_to_num(output_candidates(normalized_test_intervals))
print "Finished fetching y_test"
print y_train.shape, y_val.shape, y_test.shape


y_train = np.expand_dims(y_train, axis=2)
y_val = np.expand_dims(y_val, axis=2)
y_test = np.expand_dims(y_test, axis=2)
print y_train.shape, y_val.shape, y_test.shape

print "Dimension of histone mark signal (output): {}".format(y_train[0].shape)

'''
Generative Adversarial Model for Genomics Seq-to-Seq
'''
# Import keras
from keras.layers import AveragePooling1D, Input, Dense, Conv1D, Dropout, BatchNormalization, Activation, ZeroPadding1D, Reshape, Flatten
from keras.layers.advanced_activations import LeakyReLU
from keras.models import Sequential, Model
from keras import optimizers
from keras import metrics
from keras import losses
from keras import backend as K
from keras.callbacks import Callback, TensorBoard, ReduceLROnPlateau, ModelCheckpoint
from keras.optimizers import Adam, SGD

'''
HYPERPARAMETERS
'''
# Dropout Rate
dropout_rate = 0.5
# First conv layer
hidden_filters_1 = 32
hidden_kernel_size_1 = window_size
# Second conv layer
output_filters = 1
output_kernel_size = 32
# Training
batch_size = 128
num_epochs = args.num_epochs


# Helper functions for writing the scores into bigwig file
from itertools import izip
from itertools import groupby
import subprocess

def interval_key(interval):
    return (interval.chrom, interval.start, interval.stop)

def merged_scores(scores, intervals, merge_type):
    # A generator that returns merged intervals/scores
    # Scores should have shape: #examples x #categories x #interval_size
    # Second dimension can be omitted for a 1D signal
    signal_dims = scores.ndim - 1
    assert signal_dims in {1, 2}

    # Only support max for now
    assert merge_type == 'max'
    score_first_dim = 1 if signal_dims == 1 else scores.shape[1]

    dtype = scores.dtype

    sort_idx, sorted_intervals =         zip(*sorted(enumerate(intervals),
                    key=lambda item: interval_key(item[1])))
    sorted_intervals = BedTool(sorted_intervals)

    # Require at least 1bp overlap
    # Explicitly convert to list otherwise it will keep opening a file when
    # retrieving an index resulting in an error (too many open files)
    interval_clust = list(sorted_intervals.cluster(d=-1))
    for _, group in groupby(izip(sort_idx, interval_clust),
                            key=lambda item: item[1].fields[-1]):
        idx_interval_pairs = list(group)
        group_idx, group_intervals = zip(*idx_interval_pairs)

        if len(idx_interval_pairs) == 1:
            yield group_intervals[0], scores[group_idx[0], ...]
        else:
            group_chrom = group_intervals[0].chrom
            group_start = min(interval.start for interval in group_intervals)
            group_stop = max(interval.stop for interval in group_intervals)

            # This part needs to change to support more merge_types (e.g. mean)
            group_score = np.full((score_first_dim, group_stop - group_start),
                                  -np.inf, dtype)
            for idx, interval in idx_interval_pairs:
                slice_start = interval.start - group_start
                slice_stop = slice_start + (interval.stop - interval.start)
                group_score[..., slice_start:slice_stop] = np.maximum(group_score[..., slice_start:slice_stop], scores[idx, ...])
            if signal_dims == 1:
                group_score = group_score.squeeze(axis=0)
            yield Interval(group_chrom, group_start, group_stop), group_score

def interval_score_pairs(intervals, scores, merge_type):
    return (izip(intervals, scores) if merge_type is None
            else merged_scores(scores, intervals, merge_type))

def _write_1D_deeplift_track(scores, intervals, file_prefix, merge_type='max',
                             CHROM_SIZES='/mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes'):
    assert scores.ndim == 2

    bedgraph = file_prefix + '.bedGraph'
    bigwig = file_prefix + '.bw'

    print 'Writing 1D track of shape: {}'.format(scores.shape)
    print 'Writing to file: {}'.format(bigwig)

    with open(bedgraph, 'w') as fp:
        for interval, score in interval_score_pairs(intervals, scores,
                                                    merge_type):
            chrom = interval.chrom
            start = interval.start
            for score_idx, val in enumerate(score):
                fp.write('%s\t%d\t%d\t%g\n' % (chrom,
                                               start + score_idx,
                                               start + score_idx + 1,
                                               val))
    print 'Wrote bedgraph.'

    try:
        output = subprocess.check_output(
            ['wigToBigWig', bedgraph, CHROM_SIZES, bigwig],
            stderr=subprocess.STDOUT)
        print 'wigToBigWig output: {}'.format(output)
    except subprocess.CalledProcessError as e:
        print 'wigToBigWig terminated with exit code {}'.format(
            e.returncode)
        print 'output was:\n' + e.output

    print 'Wrote bigwig.'

# Helper function for computing Pearson R in Keras
def pearson(y_true, y_pred):
    x = y_true
    y = y_pred
    mx = K.mean(x)
    my = K.mean(y)
    xm, ym = x-mx, y-my
    r_num = K.sum(np.multiply(xm,ym))
    r_den = K.sqrt(np.multiply(K.sum(K.square(xm)), K.sum(K.square(ym))))
    r = r_num / r_den
    r = K.maximum(K.minimum(r, 1.0), -1.0)
    return K.square(r)


print "Training the model..."

'''
Baseline CNN (based on KERAS functional API)
'''
inputs = Input(shape=(window_size, 5, ))
x = Conv1D(
        filters=hidden_filters_1,
        kernel_size=128,
        padding='same',
        activation='relu',
        strides=1,
        dilation_rate=10,
        name='gen_conv1d_1')(inputs)
x = Dropout(dropout_rate,
            name='gen_dropout_1')(x)

x = Conv1D(
        filters=hidden_filters_1,
        kernel_size=256,
        padding='same',
        activation='relu',
        strides=1,
        dilation_rate=5,
        name='gen_conv1d_2')(x)
x = Dropout(dropout_rate,
            name='gen_dropout_2')(x)

x = Conv1D(
        filters=hidden_filters_1,
        kernel_size=128,
        padding='same',
        activation='relu',
        strides=1,
        name='gen_conv1d_3')(x)
x = Dropout(dropout_rate,
            name='gen_dropout_3')(x)

outputs = Conv1D(
        filters=output_filters,
        kernel_size=output_kernel_size,
        padding='same',
        activation='linear',
        strides=1,
        name='gen_conv1d_output')(x)

model = Model(inputs=inputs, outputs=outputs)
if model_path:
    print "loading model from {}...".format(model_path)
    model.load_weights(model_path, by_name=True)


# setting an adam optimizer with 1.0 clip norm# setti 
adam = optimizers.Adam(lr=1e-3, clipnorm=1.)

print "Compiling a model with adam optimizer"
model.compile(loss=losses.mean_squared_error,
              optimizer=adam,
              metrics=[pearson, metrics.mse, metrics.mae])

# CallBack: reduce learning rate when validation loss meets plateau
reduce_lr = ReduceLROnPlateau(monitor='val_loss',
                              factor=0.96,
                              patience=5,
                              min_lr=1e-7)


# CallBack: store bigwig file for the best model
class SaveBigwig(Callback):
    def __init__(self, X_train, y_train, X_val, y_val, X_test, y_test):
        self.best_val_loss = float('Inf')
        self.best_epoch = -1
        self.X_train, self.y_train = X_train, y_train
        self.X_val, self.y_val = X_val, y_val
        self.X_test, self.y_test = X_test, y_test
        self.max_pearson = -1.0
        self.epochs = 0

    def on_epoch_end(self, batch, logs={}):
        self.epochs += 1

        predictions = model.predict(X_train).flatten()
        avg_pearson = pearsonr(predictions, y_train.flatten())[0]
        print "Pearson R on Train set: {}".format(avg_pearson)

        val_predictions = model.predict(X_val).flatten()
        avg_val_pearson = pearsonr(val_predictions, y_val.flatten())[0]
        print "Pearson R on Val set: {}".format(avg_val_pearson)

        cur_val_loss = logs['val_loss']
        if cur_val_loss < self.best_val_loss:
            print "Record low val loss..."
            self.best_val_loss = cur_val_loss
            self.best_epoch = self.epochs

            f = open(os.path.join(srv_dir, 'meta.txt'), 'wb')
            f.write(str(self.epochs) + " " + str(avg_pearson) + "  " + str(avg_val_pearson) + "\n")
            max_pearson = avg_val_pearson

            test_predictions = model.predict(X_test).flatten()
            avg_test_pearson = pearsonr(test_predictions, y_test.flatten())
            print "Pearson R on Test set: {}".format(avg_test_pearson)
            f.write("Test Pearson: " + str(avg_test_pearson))
            f.close()
            _write_1D_deeplift_track(test_predictions.reshape(X_test.shape[0], window_size),
                                                              normalized_test_intervals, os.path.join(srv_dir, 'test'))

        if self.epochs == num_epochs-1:
            _write_1D_deeplift_track(predictions.reshape(X_train.shape[0], window_size),
                                     normalized_train_intervals, os.path.join(srv_dir, 'last_train'))
            _write_1D_deeplift_track(val_predictions.reshape(X_val.shape[0], window_size),
                                     normalized_val_intervals, os.path.join(srv_dir, 'last_val'))



# CallBack: Save model checkpoint based on validation loss 
checkpointer = ModelCheckpoint(os.path.join(model_dir, "best_model.h5"),
                               monitor='val_loss',
                               verbose=1,
                               save_best_only=True,
                               mode='min')

model.fit(X_train, y_train,
          batch_size = batch_size,
          epochs = num_epochs,
          validation_data=(X_val, y_val),
          callbacks = [checkpointer,
                       SaveBigwig(X_train, y_train, X_val, y_val, X_test, y_test)])
