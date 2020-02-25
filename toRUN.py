from MasterScript_v2 import RunFromQueue

exp_desc_file = '/home/peter/Desktop/MLFSOM/data_gaussian_N10_1.5A_10fr/exp-gaussian_N10_1.5A_10fr-desc.txt'
exp_queue_file = '/home/peter/Desktop/MLFSOM/data_gaussian_N10_1.5A_10fr/exp-gaussian_N10_1.5A_10fr-queue.txt'
RunFromQueue(exp_desc_file,exp_queue_file,threads=1)