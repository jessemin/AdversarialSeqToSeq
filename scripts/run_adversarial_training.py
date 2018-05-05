import os

os.chdir('../')
os.system('python train_adversarial_model.py\
          -w=10001\
          -save=gan\
          -sample_num=100000\
          -g_weight=0.4\
          -mse_weight=1.0\
          -g_lr=0.0005\
          -d_lr=0.0001\
          --smooth_rate=0.1\
          -cuda=0')
