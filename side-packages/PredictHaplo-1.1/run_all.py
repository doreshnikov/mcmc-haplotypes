import os

def run(tname):
	f = list(map(str.strip, open('config.txt').readlines()))
	f[4] = f'../../data/{tname}/_reference.fa'
	f[8] = f'../../data/{tname}/reads_bwa.sam'
	with open(f'config.{tname}.txt', 'w') as cfg:
		for line in f:
			print(line, file=cfg)

	os.system(f'./PredictHaplo config.{tname}.txt')
	os.system(f'cp phglobal_1_3000.fas ../../data/{tname}/ph.fas')

for n in [500, 1000, 2500]:
	for mu in range(75, 40, -5):
		run(f'test.{n}.{mu}')
