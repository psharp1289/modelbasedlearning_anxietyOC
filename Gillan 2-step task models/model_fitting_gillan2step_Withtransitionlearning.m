%% Fit 2 step task
% Paul Sharp
clear all

% ------------ generate simulated data ------------------ %
% nSubs=80;
% lr1=betarnd(2,6,[nSubs,1]);
% lr_transition=betarnd(2,2,[nSubs,1]);
% invtemp_mb=gamrnd(4,1,[nSubs,1]);
% invtemp_mf=gamrnd(2,1,[nSubs,1]);
% invtemp_mf2=gamrnd(2,1,[nSubs,1]);
% invtemp_2nd=gamrnd(4,1,[nSubs,1]);
% st=gamrnd(2,1,[nSubs,1]);
% %lr,invtemp,lr_transition,forgetting_rate,mbweight,ntrials
% nTrials=300;
% i=1;
% [S,A,R,Tm] = twoStepTask_Simulation_learning_gillan(lr1(counter),invtemp_mb(counter),invtemp_mf(counter),invtemp_mf2(counter),invtemp_2nd(counter),lr_transition(counter),st(counter),nTrials);
%         data(i).c1=A(:,2);
%         data(i).c2=A(:,1);
%         data(i).s=S;
%         data(i).o=R;
%         data(i).T=nTrials;
%         i=i+1;
% ------------ load data ------------------ %

nSubs=1410;

for counter=1:nSubs
        current_sub=sprintf('sub_%g.csv',counter);
        mat_sub = readmatrix(current_sub);
        data(counter).c1=mat_sub(:,2);
        data(counter).c2=mat_sub(:,3);
        data(counter).s=mat_sub(:,4);
        data(counter).o=mat_sub(:,5);
        data(counter).T=length(mat_sub(:,2));
end

N=length(data);
%% initialize models

% 
% % %Standard RL
model{1}.lik_func = @lik_MB_MF_gillan_noTL_nl_two;
model{1}.name = 'gillan no TL and 2 LR';
model{1}.spec.lrate.type = 'beta';
model{1}.spec.lrate.val = [1 1];
model{1}.spec.lrate2.type = 'beta';
model{1}.spec.lrate2.val = [1 1];
model{1}.spec.invtemp_mb.type = 'gamma';
model{1}.spec.invtemp_mb.val = [1 1];
model{1}.spec.invtemp_mf.type = 'gamma';
model{1}.spec.invtemp_mf.val = [1 1];
model{1}.spec.invtemp_mf2.type = 'gamma';
model{1}.spec.invtemp_mf2.val = [1 1];
model{1}.spec.invtemp_2ndstage.type = 'gamma';
model{1}.spec.invtemp_2ndstage.val = [1 1];
model{1}.spec.st.type = 'gamma';
model{1}.spec.st.val = [1 1];
model{1}.bic = nan;
% 
% % % %Standard RL
% model{1}.lik_func = @lik_MB_MF_gillan_noMF;
% model{1}.name = 'Gillan Model no MF ';
% model{1}.spec.lrate.type = 'beta';
% model{1}.spec.lrate.val = [1 1];
% model{1}.spec.lrate2.type = 'beta';
% model{1}.spec.lrate2.val = [1 1];
% model{1}.spec.lr_transition.type = 'beta';
% model{1}.spec.lr_transition.val = [1 1];
% model{1}.spec.invtemp_mb.type = 'gamma';
% model{1}.spec.invtemp_mb.val = [1 1];
% model{1}.spec.invtemp_2ndstage.type = 'gamma';
% model{1}.spec.invtemp_2ndstage.val = [1 1];
% model{1}.spec.st.type = 'gamma';
% model{1}.spec.st.val = [1 1];
% model{1}.bic = nan;
% % % % 

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

        % fit prior to resampled paramete                                                                                                                                                                                                                                                                                                                                                                                                                               rs
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

