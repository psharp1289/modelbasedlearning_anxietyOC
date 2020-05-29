function lik = lik_MB_MF_daw_lambda_tl(P,data) 
    
    %%  2-step task
    %
    %
    %
    %% Define Variables  
    S = size(P.lrate,1); % number of samples of parameters
    
    %Define transition matrices over all samples of distortion parameter
    for i=1:S
                             %s1a1 s1a2
    TransitionProbs(:,:,i) = [0.5  0.5; %s2
                              0.5  0.5];  %s3
      
    end
    lik = zeros(S,1);
    b1 = P.invtemp_1;
    b3 = P.invtemp_2;
    st = P.st; 
    lambda=P.lambda;
    lr1=P.lrate; %lrate value MF system 
    lr2=P.lrate2;
    lr_transition=reshape(P.lr_transition,1,1,S); %lr state transition
    
    w=P.mbweight;
    ntrials=data.T; 
    Qd1 = zeros(1,2);
    Qd1=repmat(Qd1,[S 1]);
    Qd2 = zeros(1,2);
    Qd2=repmat(Qd2,[S 1]);
    Qd3 = zeros(1,2);%3 states, 2 actions
    Qd3=repmat(Qd3, [S 1]);

    
    M = [0 0];
    M = repmat(M,[S 1]);
    %% loop through trials
    for t = 1:ntrials
        
        %load in data
        c1 = data.c1(t);
        if c1==1
            unchosen=2;
        else
            unchosen=1;
        end
        
        s = data.s(t);
        if s==2
            other_s=3;         
        else
            other_s=2;
        end
        c2 = data.c2(t);
        o = data.o(t);
  
        %first stage decision -- planning
        maxQ = [max(Qd2,[],2) max(Qd3,[],2)]; %max value of each action for each second level state state
        % model-based value Q(state1,action1) = P(state2|action1)*max(Q(state2,action1),Q(state2,action2)) + P(state3|action1)*max(...)
        Qm = squeeze(sum(reshape(maxQ',[2,1,S]).*TransitionProbs))';
        Qm_state1=w.*Qm+(1-w).*Qd1(1,:)+st.*M;
        
        
        lik=lik+b1.*Qm_state1(:,c1)- mfUtil.logsumexp(bsxfun(@times,b1,Qm_state1),2);
        if s==2
            lik=lik+b3.*Qd2(:,c2)- mfUtil.logsumexp(bsxfun(@times,b3,Qd2),2);
        elseif s==3
            lik=lik+b3.*Qd3(:,c2)- mfUtil.logsumexp(bsxfun(@times,b3,Qd3),2);
        end
        
        M=[0 0];
        M=repmat(M,[S 1]);  
        M(:,c1) = 1;
        %Update transition matrix
        TransitionProbs(s-1,c1,:) = TransitionProbs(s-1,c1,:) + lr_transition.*(1-TransitionProbs(s-1,c1,:));
        TransitionProbs(other_s-1,c1,:) = 1-TransitionProbs(s-1,c1,:);
%         if t==ntrials-1
%             tps=TransitionProbs
%         end
        %update q-values in model-free system
        if s==2
            Qd1(:,c1)=Qd1(:,c1)+lr1.*(Qd2(:,c2)-Qd1(:,c1));
            Qd1(:,c1)=Qd1(:,c1)+lambda.*lr1.*(o-Qd1(:,c1));
            Qd2(:,c2)=Qd2(:,c2)+lr2.*(o-Qd2(:,c2));
        elseif s==3
            Qd1(:,c1)=Qd1(:,c1)+lr1.*(Qd3(:,c2)-Qd1(:,c1));
            Qd1(:,c1)=Qd1(:,c1)+lambda.*lr1.*(o-Qd1(:,c1));
            Qd3(:,c2)=Qd3(:,c2)+lr2.*(o-Qd3(:,c2));
        end
   
    end
end
