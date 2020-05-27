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

%Standard RL
model{1}.lik_func = @lik_MB_MF_gillan_bayesian;
model{1}.name = 'Bayesian TL Update';
model{1}.spec.lrate.type = 'beta';
model{1}.spec.lrate.val = [1 1];
model{1}.spec.lrate2.type = 'beta';
model{1}.spec.lrate2.val = [1 1];
model{1}.spec.prior.type = 'beta';
model{1}.spec.prior.val = [1 1];
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
model{2}.lik_func = @lik_MB_MF_fixedTL;
model{2}.name = 'Bayesian TL Update';
model{2}.spec.lrate.type = 'beta';
model{2}.spec.lrate.val = [1 1];
model{2}.spec.lrate2.type = 'beta';
model{2}.spec.lrate2.val = [1 1]; 
model{2}.spec.invtemp_mf2.type = 'gamma';
model{2}.spec.invtemp_mf2.val = [1 1];
model{2}.spec.invtemp_mf.type = 'gamma';
model{2}.spec.invtemp_mf.val = [1 1];
model{2}.spec.invtemp_mb.type = 'gamma';
model{2}.spec.invtemp_mb.val = [1 1];
model{2}.spec.invtemp_2ndstage.type = 'gamma';
model{2}.spec.invtemp_2ndstage.val = [1 1];
model{2}.spec.st.type = 'gamma';
model{2}.spec.st.val = [1 1];
model{2}.bic = nan;
% 
% 
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


