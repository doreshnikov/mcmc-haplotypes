import os

def run(tname):
	f = list(map(str.strip, open('config').readlines()))
	ref = f[0].split()
	ref[-1] = f'../../data/{tname}/_reference.fa'
	reads = f[1].split()
	reads[-1] = f'../../data/{tname}/reads_bwa.sam'

	f[0] = ' '.join(ref)
	f[1] = ' '.join(reads)
	with open(f'config.{tname}', 'w') as cfg:
		for line in f:
			print(line, file=cfg)

	os.system(f'./aBayesQR config.{tname}')
	os.system(f'cp test_ViralSeq.txt ../../data/{tname}/abqr.txt')

for n in [500, 1000, 2500]:
	for mu in range(75, 40, -5):
		run(f'test.{n}.{mu}')
