from sys import argv
import os

config = f'''filename of reference sequence (FASTA) : test.{argv[1]}.{argv[2]}/_reference.fa
filname of the aligned reads (sam format) : test.{argv[1]}.{argv[2]}/reads_errFree.sam
paired-end (1 = true, 0 = false) : 0
SNV_thres : 0.001
reconstruction_start : 1
reconstruction_stop: 3000
min_mapping_qual : 60
min_read_length : 150
max_insert_length : 250
characteristic zone name : test
seq_err (assumed sequencing error rate(%)) : 0.01
MEC improvement threshold : 0.0395 '''

with open('abqr_config', 'w') as cfg:
    print(config, file=cfg)
os.system('../side-packages/aBayesQR-master/aBayesQR abqr_config')