import os

os.chdir('../')
os.system('python evaluate_adversarial_chromwide.py\
          -save=gan_genomwide\
          -cuda=0')
