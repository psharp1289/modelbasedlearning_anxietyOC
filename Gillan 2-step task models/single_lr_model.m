%% Fit 2 step task
% Paul Sharp
clear all
   
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
model{1}.lik_func = @lik_MB_MF_fixedTL;
model{1}.name = 'Bayesian TL Update';
model{1}.spec.lrate.type = 'beta';
model{1}.spec.lrate.val = [1 1];
model{1}.spec.lrate2.type = 'beta';
model{1}.spec.lrate2.val = [1 1]; 
model{1}.spec.invtemp_mf2.type = 'gamma';
model{1}.spec.invtemp_mf2.val = [1 1];
model{1}.spec.invtemp_mf.type = 'gamma';
model{1}.spec.invtemp_mf.val = [1 1];
model{1}.spec.invtemp_mb.type = 'gamma';
model{1}.spec.invtemp_mb.val = [1 1];
model{1}.spec.invtemp_2ndstage.type = 'gamma';
model{1}.spec.invtemp_2ndstage.val = [1 1];
model{1}.spec.st.type = 'gamma';
model{1}.spec.st.val = [1 1];
model{1}.bic = nan;
% 
%Standard RL
model{2}.lik_func = @lik_MB_MF_gillan_nl_two;
model{2}.name = 'Combined MB MF with Transition Learning';
model{2}.spec.lrate.type = 'beta';
model{2}.spec.lrate.val = [1 1];
model{2}.spec.lrate2.type = 'beta';
model{2}.spec.lrate2.val = [1 1]; 
model{2}.spec.lr_transition.type = 'beta';
model{2}.spec.lr_transition.val = [1 1];
model{2}.spec.invtemp_mf.type = 'gamma';
model{2}.spec.invtemp_mf.val = [1 1];
model{2}.spec.invtemp_mf2.type = 'gamma';
model{2}.spec.invtemp_mf2.val = [1 1];
model{2}.spec.invtemp_mb.type = 'gamma';
model{2}.spec.invtemp_mb.val = [1 1];
model{2}.spec.invtemp_2ndstage.type = 'gamma';
model{2}.spec.invtemp_2ndstage.val = [1 1];
model{2}.spec.st.type = 'gamma';
model{2}.spec.st.val = [1 1];
model{2}.bic = nan;

% 
%% fit models

S = 100000; % number of samples
improvement = nan;

while ~(improvement < 0) % repeat until fit stops improving
    oldbic = model{1}.bic;

    for n = 1:N
        model{2} = mfUtil.randomP(model{2}, S); % sample random parameter values
        lik = model{2}.lik_func(model{2}.P, data(n)); % compute log-likelihood for each sample
        model{2} = mfUtil.computeEstimates(lik, model{2}, n); % resample parameter values with each sample weighted by its likelihoods
    end
    mode(round(model{2}.P.lr_transition,2))
    lr_mode=repmat(mode(round(model{2}.P.lr_transition,2)),S,1);
    for n = 1:N
        model{1} = mfUtil.randomP(model{1}, S); % sample random parameter values
        lik = model{1}.lik_func(model{1}.P, data(n),lr_mode); % compute log-likelihood for each sample
        model{1} = mfUtil.computeEstimates(lik, model{1}, n); % resample parameter values with each sample weighted by its likelihoods
    end

    % fit prior to resampled paramete                                                                                                                                                                                                                                                                                                                                                                                                                               rs
    model{1} = mfUtil.fit_prior(model{1});

    % compute goodness-of-fit measures
    Nparams = 2*length(fieldnames(model{1}.spec)); % number of hyperparameters (assumes 2 hyperparameters per parameter)
    temp=[];
    for i=1:nSubs
        temp=[temp;data(i).c1];
    end
    Ndatapoints = numel(temp); % total number of samples
    model{1}.evidence = sum([model{1}.fit.evidence]); % total evidence
    model{1}.bic = -2*model{1}.evidence + Nparams*log(Ndatapoints); % Bayesian Information Criterion
    improvement = oldbic - model{1}.bic; % compute improvement of fit
    fprintf('%s - %s    old: %.2f       new: %.2f      \n', model{1}.name, 'bic', oldbic, model{1}.bic)
end


