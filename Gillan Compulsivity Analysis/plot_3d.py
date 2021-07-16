import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from numpy.polynomial.polynomial import polyfit


from scipy.stats import shapiro
from scipy.stats import ttest_ind as tt
from scipy.stats import spearmanr as corrp
import numpy as np
from statsmodels.graphics.gofplots import qqplot

font = {'family' : 'sans-serif',
        'weight' : 'light',
        'size'   : 16}

matplotlib.rc('font', **font)
bad_indices=[]
sr_data=pd.read_csv('self_report_study2.csv') #load self-report data
mb_agnostic=pd.read_csv('mb_scores_rares_empirical_best.csv')
mb_scores=mb_agnostic['MB_behav']
state_t1=pd.read_csv('Gillan_TL_full_lrT.csv',header=None) #load state transition lrs
state_t=pd.read_csv('Gillan_Or_full_lrT_decay.csv',header=None) #load state transition lrs
print(len(state_t1))
r,p=corrp(state_t1[0],state_t[0])
print('CORREL ST TL both models : {}, p {}'.format(r,p))
it_mb=pd.read_csv('Gillan_Or_full_MB_decay.csv',header=None) #load MB beta
# it_mb=np.log(it_mb)
mf1=pd.read_csv('Gillan_Or_full_MF1_decay.csv',header=None)
mf2=pd.read_csv('Gillan_Or_full_mf2_decay.csv',header=None)
lr1=pd.read_csv('Gillan_Or_full_lr1_decay.csv',header=None)
# lr2=pd.read_csv('Gillan_O_full_lr2_decay.csv',header=None)
st=pd.read_csv('Gillan_Or_full_st_decay.csv',header=None)
mf1a=pd.read_csv('Gillan_Or_full_mf1a_decay.csv',header=None)
mf=mf1a+mf1

# mf1=pd.read_csv('mf1.csv',header=None)
# lr1=pd.read_csv('lr1.csv',header=None)
# lr2=pd.read_csv('lr2.csv',header=None)

# temp2=pd.read_csv('temp2nd.csv',header=None)
# stick=pd.read_csv('st.csv',header=None)

print('state transition LR mean: {}, sd: {}'.format(np.mean(state_t[0]),np.std(state_t[0])))
f2_low_lr=[]
f2_low_mb=[]
f2_high_lr=[]
f2_high_mb=[]
f2_low_mba=[]
f2_high_mba=[]
sr_data_r=sr_data
fac1=sr_data_r['Factor1']
fac2=sr_data_r['Factor2']
fac3=sr_data_r['Factor3']
print('mean MB agnostic score: {}, sd: {}'.format(np.mean(mb_scores),np.std(mb_scores)))
bad_rows=[]
low_mb=[]

# it_mb=np.log(it_mb)
# mf=np.log(mf)
# state_t=np.log(state_t)
# fac2=fac2+2.65
# fac2r=np.log(fac2+2.65)
# print('MINIMUM VALUE COMPULSIVITY : {}'.format(np.min(fac2r)))
high_lr_mbbeta=[]
high_lr_mbscores=[]
for i in range(len(fac2)):

#     if (it_mb[0][i]>15):
#         bad_rows.append(i)
        
    if (it_mb[0][i]<3):
        bad_rows.append(i)
        
    if (state_t[0][i])>0.35:
        high_lr_mbbeta.append(it_mb[0][i])
        bad_rows.append(i)
        high_lr_mbscores.append(mb_scores[i])

#     if mb_scores[i]<0.05:
#         bad_rows.append(i)
        
print('\ndistribution MB Betas for high TLs\n')
ax = sns.distplot(high_lr_mbbeta)
plt.show()

sns.scatterplot(x=high_lr_mbbeta,y=high_lr_mbscores)
plt.ylabel('mb scores (high TL)')
plt.xlabel('mb beta (high TL)')
plt.show()
# print(np.mode(it_mb[0]))
state_t = state_t.drop(labels=bad_rows)
state_t = state_t.reset_index(drop=True)
it_mb = it_mb.drop(labels=bad_rows)
it_mb = it_mb.reset_index(drop=True)
mf1 = mf1.drop(labels=bad_rows)
mf1 = mf1.reset_index(drop=True)
mf2 = mf2.drop(labels=bad_rows)
mf2 = mf2.reset_index(drop=True)
mf = mf.drop(labels=bad_rows)
mf = mf.reset_index(drop=True)
lr1 = lr1.drop(labels=bad_rows)
lr1 = lr1.reset_index(drop=True)
st = st.drop(labels=bad_rows)
st = st.reset_index(drop=True)
mf1a = mf1a.drop(labels=bad_rows)
mf1a = mf1a.reset_index(drop=True)
mb_scores = mb_scores.drop(labels=bad_rows)
mb_scores = mb_scores.reset_index(drop=True)
import scipy
# print(scipy.stats.beta.fit(state_t[0],floc=0,fscale=1))
# print(scipy.stats.gamma.fit(it_mb[0],floc=0,fscale=1))


fac2r=fac2.drop(labels=bad_rows)
fac2r = fac2r.reset_index(drop=True)
# fac2r=np.log(fac2r+2.65)
above_2_high=[]
above_2_low=[]
for i in range(len(fac2r)):

    if fac2r[i]>=1.0:
            f2_high_lr.append(state_t[0][i])
            f2_high_mb.append(it_mb[0][i])
            f2_high_mba.append(mb_scores[i])
            if state_t[0][i]>0.094:
                above_2_high.append(state_t[0][i])
    elif fac2r[i]<=-1:
            if state_t[0][i]>0.094:
                above_2_low.append(state_t[0][i])
            f2_low_lr.append(state_t[0][i])
            f2_low_mb.append(it_mb[0][i])
            f2_high_mba.append(mb_scores[i])
        

        
# it_mb=np.log(it_mb)

# state_t=np.log(state_t)
# lr1=np.log(lr1)

# state_t=state_t+1
print('low MB performers compulsivity scores: {}'.format(low_mb))
print('mean low comp: {}, mean high comp state LR: {}'.format(np.median(f2_low_lr),np.median(f2_high_lr)))
print('mean low comp: {}, mean high comp MB-beta: {}'.format(np.mean(f2_low_mb),np.mean(f2_high_mb)))
print('')
print('percentage high comp above mean TL : {}'.format(len(above_2_high)/len(f2_high_lr)))
print('percentage low above mean TL : {}'.format(len(above_2_low)/len(f2_low_lr)))
print('')
# fac2r=np.log(fac2r+2.65)
print('mean: {}, sd: {} of TL full sample'.format(np.mean(state_t[0]),np.std(state_t[0])))
print('')
print('mean: {}, sd: {} of compulsivity full sample'.format(np.mean(fac2),np.std(fac2)))
print('mean: {}, sd: {} of compulsivity reduced sample'.format(np.mean(fac2r),np.std(fac2r)))

t,p=tt(fac2,fac2r,equal_var=False)
print('difference in compulsivity before after t: {}, p:{}'.format(t,p))
print('here')
print(mf[0][1])
print('here')
print(len(it_mb[0]))
print('here')
ratio_mfmb=[(it_mb[0][i]-mf[0][i])/ (mf[0][i]+it_mb[0][i]) for i in range(len(it_mb[0]))]
print('B-MF median: {}'.format(np.median(mf1)))
print('B-MF mean: {}'.format(np.mean(mf1[0])))

fac1=fac1.drop(labels=bad_rows)
fac3=fac3.drop(labels=bad_rows)
sr_data= sr_data.drop(labels=bad_rows)
sr_data =sr_data.reset_index(drop=True)

# print('mean: {}, sd: {} of compulsivity small sample'.format(np.mean(fac2_r),np.std(fac2_r)))
r,p=corrp(mb_scores,state_t[0])
print('model agnostic scores and state_t: {}, pval: {}'.format(r,p))
r,p=corrp(mb_scores,it_mb[0])
print('model agnostic scores and MB beta: {}, pval: {}'.format(r,p))
r,p=corrp(mb_scores,lr1[0])
print('model agnostic scores and decay rate: {}, pval: {}'.format(r,p))
r,p=corrp(state_t[0],lr1[0])
print('state TL and decay rate: {}, pval: {}'.format(r,p))

r,p=corrp(it_mb[0],state_t[0])
print('MB beta and state_t: {}, pval: {}'.format(r,p))

TLMB=state_t*it_mb
r,p=corrp(state_t[0],TLMB[0])
print('TLMB and STATE TL: {}, pval: {}'.format(r,p))

lrt_r=state_t
corr_fac1=sr_data_r['Factor1'].corr(np.log(lrt_r[0]))
print(corr_fac1)
iq=sr_data_r['iq']
corr_fac2=sr_data_r['Factor2'].corr(np.log(lrt_r[0]))
print(corr_fac2)
corr_fac3=sr_data_r['Factor3'].corr(np.log(lrt_r[0]))
print(corr_fac3)

sns.scatterplot(x=fac2r,y=state_t[0])
plt.ylabel('State Transition Learning ')
plt.xlabel('Compulsivity')
plt.show()

sns.scatterplot(x=fac2r,y=mf[0])
plt.ylabel('MF  ')
plt.xlabel('Compulsivity')
plt.show()

sns.scatterplot(x=fac2r,y=lr1[0])
plt.ylabel('LR DECAY ')
plt.xlabel('Compulsivity')
plt.show()

sns.scatterplot(x=fac2r,y=ratio_mfmb)
plt.ylabel('ratio mfmb ')
plt.xlabel('Compulsivity')
plt.show()

sns.scatterplot(x=state_t[0],y=it_mb[0])
plt.ylabel('MB Beta ')
plt.xlabel('State T')
plt.show()

sns.scatterplot(x=fac2r,y=it_mb[0])
plt.ylabel('(log) MB Beta ')
plt.xlabel('Compulsivity')
plt.show()
sns.scatterplot(x=fac2r,y=mb_scores)
plt.ylabel('MB ma')
plt.xlabel('Compulsivity')

plt.show()

from mpl_toolkits.mplot3d import Axes3D

sns.set(style = "darkgrid")

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
df=pd.DataFrame()
df['compulsivity']=fac2r
df['model agnostic']=mb_scores
df['TL']=state_t
x =df['compulsivity']
y = df['model agnostic']
z = df['TL']

ax.set_xlabel("Compulsivity")
ax.set_ylabel("MA")
ax.set_zlabel("TL")

ax.scatter(x, y, z)

plt.show()


fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
df=pd.DataFrame()
df['compulsivity']=fac2r
df['mb beta']=it_mb
df['TL']=state_t
x =df['compulsivity']
y = df['mb beta']
z = df['TL']

ax.set_xlabel("Compulsivity")
ax.set_ylabel("MB_Beta")
ax.set_zlabel("TL")

ax.scatter(x, y, z)

plt.show()


sns.scatterplot(x=it_mb[0],y=mb_scores)
b, m = polyfit(it_mb[0], mb_scores, 1)
plt.plot(it_mb[0], mb_scores, '.')
plt.plot(it_mb[0], b + m * it_mb[0], '-')
plt.ylabel('MB MA')
plt.xlabel('MB beta')
plt.ylim(0.01,1)
plt.savefig('mbHigh_modelagnosticMB_mbBeta.png', dpi=300,bbox_inches='tight')
plt.show()
sns.scatterplot(x=state_t[0],y=mb_scores)
b, m = polyfit(state_t[0], mb_scores, 1)
plt.plot(state_t[0], mb_scores, '.')
plt.plot(state_t[0], b + m * state_t[0], '-')
plt.ylabel('Model Agnostic MB Scores')
plt.xlabel('LR')
plt.ylim(0.01,1)
plt.savefig('mbHigh_modelagnosticMB_learningRate.png', dpi=300,bbox_inches='tight')
plt.show()
stat,p=shapiro(lrt_r[0])
print("TL normality test: {} , p-val: {}".format(stat,p))
# q-q plot
print('qq plot lr transition')
qqplot(lrt_r[0], line='s')
plt.show()
print('qq plot mb beta')
qqplot(it_mb[0], line='s')
plt.show()
stat,p=shapiro(it_mb[0])
print("MB normality test: {} , p-val: {}".format(stat,p))
qqplot(mf1[0], line='s')
plt.show()
stat,p=shapiro(mf1[0])
print("mf1 normality test: {} , p-val: {}".format(stat,p))
qqplot(mf2[0], line='s')
plt.show()
stat,p=shapiro(mf2[0])
print("Mf2 normality test: {} , p-val: {}".format(stat,p))
qqplot(lr1[0], line='s')
plt.show()
stat,p=shapiro(lr1[0])
print("LR1 normality test: {} , p-val: {}".format(stat,p))
# qqplot(lr2[0], line='s')
# plt.show()
# stat,p=shapiro(lr2[0])
# print("LR2 normality test: {} , p-val: {}".format(stat,p))
qqplot(mf1a[0], line='s')
plt.show()
stat,p=shapiro(mf1a[0])
print("MF1a normality test: {} , p-val: {}".format(stat,p))
