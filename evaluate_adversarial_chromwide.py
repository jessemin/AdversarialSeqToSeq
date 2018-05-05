'''
Chromosome wide evaluation of the trained adversarial model

This script is only hard-coded for Chromosome 22.

Author: Jesik Min
'''

import os, sys
sys.path.append("..")
# Custom file path package
from data import Data_Directories
# Custom utility package
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
from genomelake.backend import load_directory

# Setup Arguments
parser.add_argument('-e',
                    '--num_epochs',
                    required=False,
                    type=int,
                    default=300,
                    dest="num_epochs",
                    help="Number of epochs for training")
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
                    type=str,
                    default='saved_models/pretrained_gan/best_generator.h5',
                    dest="model_path",
                    help="Path to pre-trained model")
parser.add_argument('-save',
                    '--save_dir',
                    required=False,
                    type=str,
                    default='chromwide_evaluation',
                    dest="save_dir",
                    help="Where to save the best model, logs, and its predictions")
parser.add_argument('-cuda',
                    '--cuda',
                    required=True,
                    type=str,
                    dest="cuda",
                    help="Cuda visible devices")



args = parser.parse_args()

# Parse all arguments

day = args.day
frag = args.frag
histone = args.histone
model_path = args.model_path
save_dir = args.save_dir
cuda = args.cuda

os.environ["CUDA_VISIBLE_DEVICES"]=cuda

# Logging directories
srv_dir = os.path.join("/srv", "www", "kundaje", "jesikmin", "test_experiments", save_dir)
if not os.path.exists(srv_dir):
    os.makedirs(srv_dir)

data = Data_Directories()
X_test = load_directory(data.input_atac[day][frag], in_memory=True)['chr22']._arr
X_test = np.expand_dims(np.nan_to_num(X_test), axis=0)
print "Finished fetching X_test"
print X_test.shape
print "Dimension of ATAC-seq signal (input): {}".format(X_test[0].shape)

y_test = load_directory('/srv/scratch/jesikmin/output/bcolz/', in_memory=True)['chr22']
y_test = np.expand_dims(y_test, axis=0)
y_test = np.expand_dims(y_test, axis=2)
print "Finished fetching Y_test"
print y_test.shape
print "Dimension of ChIP-seq signal (output): {}".format(y_test[0].shape)

'''
Generator only
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
# GAN Discriminator
smooth_rate = 0.1
d_train_freq = 1
# Dropout Rate
dropout_rate = 0.5
# First conv layer
hidden_filters_1 = 32
# Second conv layer
output_filters = 1
output_kernel_size = 32
# Training
batch_size = 128


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

    sort_idx, sorted_intervals = zip(*sorted(enumerate(intervals),
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


class Generator():
    def __init__(self,
                 X_test,
                 y_test,
                 srv_dir):

        # Set train/val/test
        self.X_test, self.y_test = X_test, y_test

        self.srv_dir = srv_dir

        # Basic parameters
        # 1) Number of channels
        self.channels = 5
        # 2) Input and Output shape
        self.input_shape = (None, self.channels,)
        self.output_shape = (None, 1,)
        print "Input and Output Shape"
        print self.input_shape, self.output_shape

        # Build and compile the generator
        self.generator = self.build_generator()
        self.generator.compile(loss='binary_crossentropy', optimizer='rmsprop')

    def build_generator(self):
        # Generator
        # 1) 32 * window_size Conv1D layers with RELU and Dropout

        noise_shape = self.input_shape

        model = Sequential()

        model.add(Conv1D(hidden_filters_1,
                         128,
                         padding="same",
                         strides=1,
                         input_shape=(None, 5),
                         activation='relu',
                         dilation_rate=10,
                         name='gen_conv1d_1'))
        model.add(Dropout(dropout_rate,
                  name='gen_dropout_1'))

        model.add(Conv1D(hidden_filters_1,
                         256,
                         padding="same",
                         strides=1,
                         activation='relu',
                         dilation_rate=5,
                         name='gen_conv1d_2'))
        model.add(Dropout(dropout_rate,
                  name='gen_dropout_2'))

        model.add(Conv1D(hidden_filters_1,
                         128,
                         padding="same",
                         strides=1,
                         dilation_rate=1,
                         activation='relu',
                         name='gen_conv1d_3'))
        model.add(Dropout(dropout_rate,
                  name='gen_dropout_3'))

        # 2) 1 * 16 Conv1D layers with Linear
        # NOTE: All same padding
        model.add(Conv1D(output_filters,
                         output_kernel_size,
                         padding='same',
                         strides=1,
                         activation='linear',
                         name='gen_conv1d_output'))

        print "Generator"
        model.summary()

        noise = Input(shape=noise_shape)
        img = model(noise)

        # load weights for generator if specified
        # print "-"*50
        # print model.get_weights()
        model.load_weights(model_path, by_name=True)
        # print "-"*50
        # print model.get_weights()
        print "Model loading done!"

        return Model(noise, img)

    def evaluate(self):
        # ---------------------
        # Get generator's prediction and compute overall pearson on test set
        # ---------------------

        #51304566
        X_test_main = self.X_test[:,:51300000,:].reshape((10000, 5130, 5))

        X_test_left = self.X_test[:,51300000:,:]
        test_predictions_1 = self.generator.predict(X_test_main).reshape((51300000, 1))
        test_predictions_2 = self.generator.predict(X_test_left).reshape((4566, 1))
        print test_predictions_1.shape, test_predictions_2.shape

        full_x_test = np.concatenate([test_predictions_1, test_predictions_2], axis=0)
        full_x_test = full_x_test.flatten()

        #avg_test_pearson = pearsonr(full_x_test, self.y_test.flatten())
        avg_test_pearson = np.corrcoef(full_x_test, self.y_test.flatten())[0,1]
        print "Pearson R on Test set: {}".format(avg_test_pearson)

        f = open(os.path.join(self.srv_dir, 'meta.txt'), 'wb')
        f.write("Test Pearson: " + str(avg_test_pearson))
        f.close()

        _write_1D_deeplift_track(full_x_test.reshape(1, self.X_test.shape[1]),
                                 [Interval('chr22', 0, 51304566)],
                                 os.path.join(self.srv_dir, 'test'))


        print "Evaluation Complete!"


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
generator = Generator(X_test,
                y_test,
                srv_dir)
generator.evaluate()

