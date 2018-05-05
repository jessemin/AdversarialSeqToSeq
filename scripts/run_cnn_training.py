import os

os.chdir('../')
os.system('python train_cnn.py\
          -w=10001\
          -save=cnn\
          -sample_num=100000\
          -cuda=5')
