%% Fit 2 step task
% Paul Sharp

% ------------ generate simulated data ------------------ %
rewardtrajectoriesgillan=readmatrix('reward_trajectories_gillan.csv');
nSubs=50;
% lr1=[0.065088, 0.70076, 0.78998, 0.46989, 0.62763, 0.5666899999999999, 0.5457, 0.052432000000000006, 0.28715999999999997, 0.23367, 0.61869, 0.17651, 0.46458, 0.6483800000000001, 0.28686999999999996, 0.65212, 0.37393000000000004, 0.50873, 0.89838, 0.9166799999999999, 0.56524, 0.4455, 0.8321700000000001, 0.37638, 0.83886, 0.5418, 0.57509, 0.066554, 0.083592, 0.7333, 0.70505, 0.94434, 0.6651699999999999, 0.53509, 0.53488, 0.39735, 0.29154, 0.35812, 0.21639, 0.92368, 0.35754, 0.78598, 0.32187, 0.38523, 0.45082, 0.5503399999999999, 0.56019, 0.8375299999999999, 0.20481, 0.63826, 0.62187, 0.15104, 0.24869000000000002, 0.6540600000000001, 0.6090800000000001, 0.8746, 0.3566, 0.25926, 0.72611, 0.21775999999999998, 0.50758, 0.44177, 0.79972, 0.28581999999999996, 0.76283, 0.024303, 0.70654, 0.5062, 0.3509, 0.74019, 0.7076600000000001, 0.17185, 0.39073, 0.8229299999999999, 0.1085, 0.6623100000000001, 0.88497, 0.42951999999999996, 0.14437, 0.28233, 0.86269, 0.72596, 0.5178, 0.11282, 0.48175, 0.078198, 0.9244899999999999, 0.66189, 0.77107, 0.35821, 0.28608, 0.7247399999999999, 0.81055, 0.4023, 0.32665, 0.5916100000000001, 0.18747, 0.08954, 0.47747, 0.41786, 0.77472, 0.69809, 0.70806, 0.51523, 0.26481, 0.62131, 0.83667, 0.040023, 0.9066299999999999, 0.8321299999999999, 0.39496, 0.34229, 0.57412, 0.066304, 0.75365, 0.44215, 0.36894, 0.7138899999999999, 0.81812, 0.13505999999999999, 0.19908, 0.36634, 0.54322, 0.86329, 0.41531999999999997, 0.18187, 0.83939, 0.96478, 0.25201, 0.2406, 0.55947, 0.13367, 0.43091, 0.6987800000000001, 0.6168, 0.26366, 0.3187, 0.58878, 0.60764, 0.26644, 0.10136, 0.12730999999999998, 0.6487, 0.42896000000000006, 0.44905, 0.8742700000000001, 0.46934, 0.8275299999999999, 0.37373, 0.60835, 0.5856899999999999, 0.55089, 0.10868, 0.28681999999999996, 0.22480999999999998, 0.57091, 0.48215, 0.031388, 0.59335, 0.892, 0.90727, 0.5806100000000001, 0.40894, 0.47206000000000004, 0.59737, 0.8524200000000001, 0.49686, 0.17342, 0.86721, 0.48733999999999994, 0.37722, 0.33765, 0.58794, 0.5834199999999999, 0.22154000000000001, 0.23653000000000002, 0.54727, 0.8645299999999999, 0.16718, 0.14627, 0.79023, 0.9107, 0.24070999999999998, 0.5576399999999999, 0.27275, 0.7048, 0.038018, 0.9150299999999999, 0.88734, 0.20035, 0.93174, 0.26314, 0.16187, 0.5767800000000001, 0.7552399999999999, 0.83734, 0.49401999999999996, 0.17406, 0.51466, 0.63999, 0.28102, 0.4753, 0.21203000000000002, 0.1311, 0.97575, 0.48373999999999995, 0.45361999999999997, 0.69184, 0.1177, 0.54373, 0.42513999999999996, 0.41731999999999997, 0.7398100000000001, 0.44791000000000003, 0.6774600000000001, 0.45003999999999994, 0.4183, 0.08302799999999999, 0.61835, 0.47093999999999997, 0.46458, 0.7481399999999999, 0.56297, 0.14589000000000002, 0.80913, 0.60263, 0.83714, 0.12425, 0.23824, 0.6434, 0.051629999999999995, 0.69329, 0.86816, 0.34807, 0.026586000000000002, 0.11160999999999999, 0.9445899999999999, 0.74073, 0.5363600000000001, 0.33726, 0.2603, 0.61126, 0.11824000000000001, 0.56067, 0.7286, 0.29268, 0.47198, 0.34609, 0.6096199999999999, 0.28597, 0.68733, 0.67038, 0.25376, 0.54769, 0.46528, 0.50474, 0.6913600000000001, 0.33692, 0.54845, 0.3925, 0.81793, 0.54056, 0.7279899999999999, 0.6624800000000001, 0.7235, 0.6554399999999999, 0.6378199999999999, 0.40529, 0.43141, 0.53092, 0.96225, 0.256, 0.26797, 0.043985, 0.52427, 0.28431, 0.47093, 0.67653, 0.5837399999999999, 0.23169, 0.26313000000000003, 0.49983999999999995, 0.43244, 0.76213, 0.71987, 0.73078, 0.8361700000000001, 0.73441, 0.62576, 0.11755, 0.32611999999999997, 0.23711, 0.38645, 0.35345, 0.4224, 0.52665, 0.49462, 0.32859, 0.7862100000000001, 0.44437, 0.27435, 0.40263000000000004, 0.9290799999999999, 0.30739, 0.13089, 0.07675900000000001, 0.75551, 0.28194, 0.88179, 0.5838899999999999, 0.23903000000000002, 0.7841899999999999, 0.18783, 0.35512, 0.47193, 0.27392, 0.2442, 0.6574300000000001, 0.4643, 0.34125, 0.75521, 0.83015, 0.28729, 0.10740999999999999, 0.292, 0.5193, 0.36361, 0.39864, 0.51012, 0.7957, 0.6147100000000001, 0.22179000000000001, 0.36493000000000003, 0.15975, 0.74344, 0.72336, 0.74545, 0.67356, 0.72569, 0.92655, 0.53846, 0.75035, 0.5150600000000001, 0.77278, 0.54371, 0.21706, 0.22824, 0.6545, 0.5901, 0.41181999999999996, 0.7482300000000001, 0.66347, 0.77688, 0.86884, 0.42751999999999996, 0.37041999999999997, 0.63002, 0.49935, 0.54357, 0.40819, 0.51184, 0.58395, 0.7291300000000001, 0.68137, 0.7008, 0.87641, 0.44369, 0.78275, 0.35602, 0.35129, 0.57512, 0.017658, 0.86534, 0.3482, 0.84483, 0.56707, 0.71233, 0.14386, 0.27965999999999996, 0.91795, 0.84893, 0.47342, 0.65789, 0.48683999999999994, 0.66351];
lr1=betarnd(1.76,0.57,[nSubs,1]);
decay=betarnd(1.50,3.50,[nSubs,1]);
lr_transition=betarnd(0.5,3.6,[nSubs,1]);
invtemp_mb=gamrnd(0.5,10,[nSubs,1])+2.5;
invtemp_mf=gamrnd(0.6,0.2,[nSubs,1]);
invtemp_mf2=gamrnd(1.02,1.32,[nSubs,1]);
invtemp_2nd=gamrnd(3.03,0.8,[nSubs,1]);
st=normrnd(0.84,0.73,[nSubs,1]);
nTrials=200;
% 
% %(lr1,lr2,invtemp_mb,invtemp_mf,invtemp_mf2,invtemp_2ndstage,lr_transition,st,ntrials,rew_probs)
i=1;
% % cd ../twostep_data_study2/
for counter=1:nSubs
[S,A,R,Tm] = twoStepTask_Simulation_learning_gillan_lr_decay(lr1(counter),invtemp_mb(counter),invtemp_mf(counter),invtemp_mf2(counter),invtemp_2nd(counter),lr_transition(counter),st(counter),decay(counter),nTrials,rewardtrajectoriesgillan);
        % d=zeros(200,1);
        % d(:)=invtemp_mb(counter);
        % tempd=[linspace(1,200,200)',A(:,1),A(:,2),S,R,d];
        % formatSpec = "%d_simulated.csv";
        % str = sprintf(formatSpec,counter);
        % writematrix(tempd,str);
        data(i).c1=A(:,1);
        data(i).c2=A(:,2);
        data(i).s=S;
        data(i).o=R;
        data(i).T=nTrials;
        i=i+1;
end
% cd ../Github_Repo_StateTransitionLearning_Paper/'Gillan Compulsivity Analysis'/
% writematrix(lr_transition,'LRT_sim_decay.csv')
% writematrix(invtemp_mb,'ITMB_sim_decay.csv')

% lev
% ------------ load data ------------------ %
% cd twostep_data_study2/
pwd
% nSubs=1413;

% for counter=1:nSubs
%         current_sub=sprintf('sub_%g.csv',counter);
%         mat_sub = readmatrix(current_sub);
%         data(counter).c1=mat_sub(:,2);
%         data(counter).c2=mat_sub(:,3);
%         data(counter).s=mat_sub(:,4);
%         data(counter).o=mat_sub(:,5);
%         data(counter).T=length(mat_sub(:,2));
%         
% 
%   
% end
cd 'Gillan 2-step task models'
pwd
N=length(data)
len_data=N
%% initialize models
% 




% 

model{1}.lik_func = @lik_MB_MF_gillan_orig_TL_seperate_decay;
model{1}.name = 'Gillan TL';
model{1}.spec.lrate.type = 'beta';
model{1}.spec.lrate.val = [1 1];
model{1}.spec.lr_transition.type = 'beta';
model{1}.spec.lr_transition.val = [1 1];
model{1}.spec.decay.type = 'beta';
model{1}.spec.decay.val = [1 1];
model{1}.spec.invtemp_mb.type = 'gamma';
model{1}.spec.invtemp_mb.val = [1 1];
model{1}.spec.invtemp_mf.type = 'gamma';
model{1}.spec.invtemp_mf.val = [1 1];
model{1}.spec.invtemp_mf2.type = 'gamma';
model{1}.spec.invtemp_mf2.val = [1 1];
model{1}.spec.invtemp_2ndstage.type = 'gamma';
model{1}.spec.invtemp_2ndstage.val = [1 1];
model{1}.spec.st.type = 'norm';
model{1}.spec.st.val = [1 1];
model{1}.bic = nan;

% 
% model{2}.lik_func = @lik_MB_MF_gillan_orig;
% model{2}.name = 'Gillan NO TL';
% model{2}.spec.lrate.type = 'beta';
% model{2}.spec.lrate.val = [1 1];
% model{2}.spec.invtemp_mb.type = 'gamma';
% model{2}.spec.invtemp_mb.val = [1 1];
% model{2}.spec.invtemp_mf.type = 'gamma';
% model{2}.spec.invtemp_mf.val = [1 1];
% model{2}.spec.invtemp_mf2.type = 'gamma';
% model{2}.spec.invtemp_mf2.val = [1 1];
% model{2}.spec.invtemp_2ndstage.type = 'gamma';
% model{2}.spec.invtemp_2ndstage.val = [1 1];
% model{2}.spec.st.type = 'norm';
% model{2}.spec.st.val = [1 1];
% model{2}.bic = nan;

%% fit models

S = 100000; % number of samples
for m = 1:length(model)

    improvement = nan;
    while ~(improvement < 0) % repeat until fit stops improving
        oldbic = model{m}.bic;
        
        for n = 1:N
            model{m} = mfUtil.randomP(model{m}, S); % sample random parameter values
            lik = model{m}.lik_func(model{m}.P, data(n)); % compute log-likelihood for each sample
            model{m} = mfUtil.computeEstimates(lik, model{m}, n); % resample parameter values with each sample weighted by its likelihoods
        end

        % fit prior to resampled parameters
        model{m} = mfUtil.fit_prior(model{m});

        % compute goodness-of-fit measures
        Nparams = 2*length(fieldnames(model{m}.spec)); % number of hyperparameters (assumes 2 hyperparameters per parameter)
        temp=[];
        for i=1:nSubs
            temp=[temp;data(i).c1];
        end
        Ndatapoints = numel(temp); % total number of samples
        model{m}.evidence = sum([model{m}.fit.evidence]); % total evidence
        model{m}.bic = -2*model{m}.evidence + Nparams*log(Ndatapoints); % Bayesian Information Criterion
        improvement = oldbic - model{m}.bic; % compute improvement of fit
        fprintf('%s - %s    old: %.2f       new: %.2f      \n', model{m}.name, 'bic', oldbic, model{m}.bic)
    end
end

% save('model_with_alldecayTLlapse','model');
%% display correspondence between model 1 parameters and original parameters
fits = [model{1}.fit.P];
invtempmb=[fits.invtemp_mb];
lrates=[fits.lrate];
invtempmf1=[fits.invtemp_mf];
invtempmf2=[fits.invtemp_mf2];
lr_transitions=[fits.lr_transition];
decays=[fits.decay];
invtemp=[fits.invtemp_2ndstage];
sts=[fits.st];% model{1}.spec.lr_transition.type = 'beta';
% model{1}.spec.lr_transition.val = [1 1];

correlation_mb_temp=corr([invtemp_mb [invtempmb.val]'],'Type','Pearson')
correlation_invtemps_2nd=corr([invtemp_2nd [invtemp.val]'],'Type','Pearson')
correlation_lrt=corr([lr_transition [lr_transitions.val]'],'Type','Pearson')
correlation_lrate1=corr([lr1 [lrates.val]'],'Type','Pearson')
correlation_decays=corr([decay [decays.val]'],'Type','Pearson')
correlation_stick=corr([st [sts.val]'],'Type','Pearson')
correlation_mf2=corr([invtemp_mf2 [invtempmf2.val]'],'Type','Pearson')
correlation_mf=corr([invtemp_mf [invtempmf1.val]'],'Type','Pearson')

figure; 
% subplot(1,2,1);
% scatter(lr_transition, [lr_transitions.val]);
% title('LR Transitions');
% xlabel('true value'); ylabel('fitted value');
% subplot(1,2,2);

yourvariables=[lr1 lr_transition invtemp_mf invtemp_mf2 invtemp_mb invtemp_2nd st decay [lrates.val]' [lr_transitions.val]' [invtempmf1.val]' [invtempmf2.val]' [invtempmb.val]' [invtemp.val]' [sts.val]' [decays.val]'];
yourlabelnames={'\alpha','\gamma','\beta_{MF0}','\beta_{MF1}','\beta_{MB}','\beta_{MF2}','\beta_{AP}','D','\alpha','\gamma','\beta_{MF0}','\beta_{MF1}','\beta_{MB}','\beta_{MF2}','\beta_{AP}','D'};
%% Code Snippet
x=corr(yourvariables,'Type','Pearson');
imagesc(x); % Display correlation matrix as an image
set(gca, 'XTick', 1:18); % center x-axis ticks on bins
set(gca, 'YTick', 1:18); % center y-axis ticks on bins
set(gca, 'XTickLabel', yourlabelnames); % set x-axis labels
set(gca, 'YTickLabel', yourlabelnames); % set y-axis labels
title('Corr Fitted & Real Params', 'FontSize', 16); % set title
colormap('jet'); % Choose jet or any other color scheme
exportgraphics(figure,'paramrecov_decay.png','Resolution',300)
