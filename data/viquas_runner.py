from sys import argv
import os

d = f'test.{argv[1]}.{argv[2]}'
print(d)
os.system(f'samtools view -S -b {d}/reads_errFree.sam > {d}/reads_errFree.bam')
print('Built .bam')
os.chdir('../side-packages/ViQuaS1.3')
print(os.getcwd())
os.system(f'Rscript ViQuaS.R ../../data/{d}/_reference.fa ../../data/{d}/reads_errFree.bam 3 0.7 0 3000')
os.chdir('../../data')