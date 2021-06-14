import os

def run(tname):
	f = list(map(str.strip, open('config.txt').readlines()))
	f[4] = f'../../data/{tname}/_reference.fa'
	f[8] = f'../../data/{tname}/reads_comb.sam'
	with open(f'config.{tname}.efree.txt', 'w') as cfg:
		for line in f:
			print(line, file=cfg)

	os.chdir('../../data')
	os.system(f'python3 ../src/main/kotlin/datagen/combine_reads.py {tname}')
	os.chdir('../side-packages/PredictHaplo-1.1')

	os.system(f'./PredictHaplo config.{tname}.efree.txt')
	os.system(f'cp phglobal_1_3000.fas ../../data/{tname}/ph_efree.fas')

for n in [500, 1000, 2500]:
	for mu in range(75, 40, -5):
		run(f'test.{n}.{mu}')
