#  A Sequence-to-sequence Regression of Genome-wide Chromatin Data Through Adversarial Training

This repository contains the implementation for the paper: **A Sequence-to-sequence Regression of Genome-wide Chromatin Data Through Adversarial Training**.

We impute histone ChIP-seq signal from ATAC-seq signal using the adversarial training approach we propose in the paper. 

<p align="center">
    <img src="imgs/overview.png" width="30%" align="middle">
</p><br/>

Buenrostro et al.(2013)[3] showed that different chromatin states such as CTCF, TSS, enhancer sites, and repressed regions show clear fragment length distribution of ATAC-seq reads. In other words, the position and fragment length of ATAC-seq signal can be predictive
of different chromatin states and functionalities. This implies that ChIP-seq peaks may be correlatd with ATAC-seq peaks, because we can infer chromatin states using histone ChIP-seq signals. Given that ATAC-seq requires only 500 to 50,000 cells, while ChIP-Seq requires much larger biological samples, typically millions of cells, regressing ChIP-seq based on ATAC-seq data is a meaningful task.

While convolutional neural network has been widely used in the field of deep learning for genomics, we expect deep adversarial network to perform better in large-scale epigenomic signal imputation task than CNN. We introduce three major modifications to the vanilla generative adversarial network architecture.

1. Generator component of deep adversarial network so that the generator takes ATAC-seq signal instead of random noise and outputs ChIP-seq signal
2. Composite loss function that takes account of both mean squared error and adversarial loss
3. One-sided label smoothing as introduced in Salimans et al.(2016)[2]


<br/>
<p align="center">
    <img src="imgs/adv_architecture.png" align="middle">
</p>
<br/><br/>




We trained our model on GeForce GTX TITAN X for 300 epochs.

Here are two visualizations created with WashU Epigenome Browser.

The first visualization below shows how the generator trained through adversarial training outperforms the CNN-based model. We can observe that the adversarial model captures all three ChIP-seq peaks, while the CNN only capture one ChIP-seq peak in the middle. Although the leftmost ChIP-seq peak is not as clearly distinguishable as the other two peaks are, we still get a clear sense that there are three major peaks where ChIP-seq signal is signicant.

![Example 1](imgs/best_results.png)
The second visulization shows one of the ChIP-seq regression done by the generator trained through the adversarial training scheme. The position where the ATAC-seq peak occurs in the first channel represents nucleosome-free regions, and consequently, the ChIP-seq signal does not usually exist. The adversarial model correctly predicts that there would be no ChIP-seq signal at this particular position. The positions where the ATAC-seq peaks are present in the second and the third channels represent nucleosomal regions where the ChIP-seq peaks usually exist, and the adversarial model also very clearly predicts the two peaks in the nucleosomal regions. This result is promising and suggests that an adversarial model learns to capture biological implication between ATAC-seq and ChIP-seq signal.

![Example 2](imgs/perchannel_best2.png)
For more details, please refer to the paper. Note that we cannot provide our ATAC-seq and ChIP-seq dataset used in the paper.

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

We use Keras 2.0 with Tensorflow backend and [genomelake](https://pypi.org/project/genomelake/) for the implementation.

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

[1] We partially adopted the structure from:
<https://github.com/eriklindernoren/Keras-GAN/tree/master/gan>

[2] Salimans, T., Goodfellow, I., Zaremba, W., Cheung, V., Radford, A. and Chen, X. (2018). Improved Techniques for Training GANs. [online] Arxiv.org. Available at: https://arxiv.org/abs/1606.03498 [Accessed 6 May 2018].

[3] Jason D Buenrostro, Paul G Giresi, Lisa C Zaba, Howard Y Chang, and William J Greenleaf. Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, dna-binding proteins and nucleosome position. Nature Methods, 2013.
