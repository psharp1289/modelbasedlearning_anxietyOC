import os
import pandas as pd
import csv

path_to_self_report_data='/home/paul/MPC_Files/projects/Worry_model_2019/model_scripts'
path_to_data='/home/paul/MPC_Files/projects/Worry_model_2019/model_scripts/twostep_data_study2'
os.chdir(path_to_self_report_data)
sr_data=pd.read_csv('self_report_study2.csv')
os.chdir(path_to_data)
rel_files=[x+'.csv' for x in sr_data['subj'] if x.startswith('3')]
# rel_files.remove('3DH6GAKTYYPQQ6J3BR0D6CG2ANFZY7.csv')
cols=['trial','d1','d2','d3','d4','st1','a1','a1RT','i','st2','a2','s','a2RT','rew','irrel']
subs_to_old=[['old_file_name','new_sub_number']]
c=0
for file_name in rel_files:
	c+=1
	with open(file_name,'r') as f:
		r=csv.reader(f)
		ls=[l for l in r]
	found=0
	counter=1
	counter_blanks=1
	for line in range(len(ls)):
		if len(ls[line])==0:
			counter_blanks+=1
		h=ls[line]
		if 'instructionLoop' in h:
				start_index=counter
		counter+=1
	try:

		x=pd.read_csv(file_name,header=start_index-counter_blanks,names=cols)
		y=x[['a1','a2','s','rew']]
		y = y[y['s'] != -1]
		y = y[y['rew'] != -1]
		y.to_csv('sub_{}.csv'.format(c),header=False)
		subs_to_old.append([file_name,c])
	except:
		c-=1

with open('sub_conversion_file.csv','a') as g:
	w=csv.writer(g)
	w.writerows(subs_to_old)

