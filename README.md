#  A Sequence-to-sequence Regression of Genome-wide Chromatin Data Through Adversarial Training

This repository contains the implementation for my honors thesis paper: **A Sequence-to-sequence Regression of Genome-wide Chromatin Data Through Adversarial Training** We try to impute histone ChIP-seq signal from ATAC-seq signal using the adversarial training approach we propose in the paper.  For more details, please refer to the paper. Please also note that we cannot provide our ATAC-seq and ChIP-seq dataset.

## How to Run

### Implementation of adversarial model and CNN baseline model
 - [train_adversarial_model.py](https://github.com/jessemin/AdversarialSeqToSeq/blob/master/train_adversarial_model.py "train_adversarial_model.py")
 - [train_cnn.py](https://github.com/jessemin/AdversarialSeqToSeq/blob/master/train_cnn.py "train_cnn.py")
 
### Implementation of chromosome-wide evaluation
 - [evaluate_adversarial_chromwide.py](https://github.com/jessemin/AdversarialSeqToSeq/blob/master/evaluate_adversarial_chromwide.py "evaluate_adversarial_chromwide.py")

### Training scripts
 - [run_adversarial_training.py](https://github.com/jessemin/AdversarialSeqToSeq/blob/master/scripts/run_adversarial_training.py "run_adversarial_training.py")
 - [run_cnn_training.py](https://github.com/jessemin/AdversarialSeqToSeq/blob/master/scripts/run_cnn_training.py "run_cnn_training.py")
 - [run_chromwide_evaluation.py](https://github.com/jessemin/AdversarialSeqToSeq/blob/master/scripts/run_chromwide_evaluation.py "run_chromwide_evaluation.py")

## Dependencies

To install all dependencies, run:

    pip install -r requirements.txt

## Contact

If you have any questions, please contact:

 - Jesik (Jesse) Min: <jesikmin@stanford.edu>

## Collaborators

 - Johnny Israeli
 - Daniel Kim
 - Dr. Paul Khavari, MD
 - Dr. Mike Snyder
 - Dr. Anshul Kundaje

## References

[1] Partially adopted the implementation from:
<https://github.com/eriklindernoren/Keras-GAN/tree/master/gan>
