function lik = lik_MB_MF_daw_orig(P,data) 
    
    %%  2-step task
    %
    %
    %
    %% Define Variables  
    S = size(P.lrate,1); % number of samples of parameters
    
    %Define transition matrices over all samples of distortion parameter
    for i=1:S
                             %s1a1 s1a2
    TransitionProbs1(:,:,i) = [0.7  0.3; %s2
                              0.3  0.7];  %s3
    
                             %s1a1 s1a2
    TransitionProbs2(:,:,i) = [0.3  0.7; %s2
                              0.7  0.3];  %s3
      
    end
    
    TransitionCounts=[0 0;
                      0 0];
    lik = zeros(S,1);
    b1 = P.invtemp_1;
    b3 = P.invtemp_2;
    st = P.st; %zeros(S,1)+0.1;
    w=P.mbweight;
    lr1=P.lrate; %lrate value MF system  
    lr2=P.lrate2;
    lambda=P.lambda;
%     lr_transition=zeros(S,1)+0.1;
%     lr_transition=reshape(P.lr_transition,1,1,S); %lr state transition
    
%     forgetting_rate=reshape(P.forgetting_rate,1,1,S);
%     w=P.mbweight;
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
        
        %determine which Tm to use
        one_sum=TransitionCounts(1,1)+TransitionCounts(2,2);
        two_sum=TransitionCounts(2,1)+TransitionCounts(1,2);
        if one_sum>=two_sum
            Tm=TransitionProbs1;
        else
            Tm=TransitionProbs2;
        end

        maxQ = [max(Qd2,[],2) max(Qd3,[],2)]; %max value of each action for each second level state state
        % model-based value Q(state1,action1) = P(state2|action1)*max(Q(state2,action1),Q(state2,action2)) + P(state3|action1)*max(...)
        Qm = squeeze(sum(reshape(maxQ',[2,1,S]).*Tm))';
        Qm_state1=w.*Qm+(1-w).*Qd1(1,:)+st.*M;
        
        
        lik=lik+b1.*Qm_state1(:,c1)- mfUtil.logsumexp(bsxfun(@times,b1,Qm_state1),2);
        
        if s==2
            lik=lik+b3.*Qd2(:,c2)- mfUtil.logsumexp(bsxfun(@times,b3,Qd2(:,:)),2);
        elseif s==3
            lik=lik+b3.*Qd3(:,c2)- mfUtil.logsumexp(bsxfun(@times,b3,Qd3(:,:)),2);
        end
        
        M=[0 0];
        M=repmat(M,[S 1]);  
        M(:,c1) = 1;
        
        %Update transition counts
        TransitionCounts(s-1,c1)=TransitionCounts(s-1,c1)+1;
        
        %update q-values in model-free system
        if s==2
            Qd1(:,c1)=Qd1(:,c1)+lr1.*(Qd2(:,c2)-Qd1(:,c1));
            Qd1(:,c1)=Qd1(:,c1)+lambda.*lr2.*(o-Qd1(:,c1));
            Qd2(:,c2)=Qd2(:,c2)+lr2.*(o-Qd2(:,c2));
        elseif s==3
            Qd1(:,c1)=Qd1(:,c1)+lr1.*(Qd3(:,c2)-Qd1(:,c1));
            Qd1(:,c1)=Qd1(:,c1)+lambda.*lr2.*(o-Qd1(:,c1));
            Qd3(:,c2)=Qd3(:,c2)+lr2.*(o-Qd3(:,c2));
        end
   
    end
end
