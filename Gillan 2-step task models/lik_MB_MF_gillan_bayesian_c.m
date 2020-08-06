function lik = lik_MB_MF_gillan_bayesian_c(P,data) 
    
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
    b1 = P.invtemp_mb;
    b2 = P.invtemp_mf;
    b2a = P.invtemp_mf2;
    b3 = P.invtemp_2ndstage;
    st = P.st; %zeros(S,1)+0.1;
    lr1=P.lrate; %lrate value MF system  
    lr2=P.lrate2;
    %Define prior parameters
    concentration=P.concentration;
    mode=zeros(S,1)+0.5;
    
    %Convert mode and concentration to alpha and beta parameters
    
    %Starting Alpha - Evidence in favor of T1
    alpha=(mode.*(concentration-2))+1;
    
    %Starting Beta - Evidence in favor of T2
    beta=((1-mode).*(concentration-2))+1;
    

    ntrials=data.T; 
    Qd1 = zeros(1,2);
    Qd1=repmat(Qd1,[S 1]);
    Qd1a = zeros(1,2);
    Qd1a=repmat(Qd1a,[S 1]);
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
        
        p1=alpha./(alpha+beta);
        p2=1-p1;

        

        maxQ = [max(Qd2,[],2) max(Qd3,[],2)]; %max value of each action for each second level state state
%         Qm_state1=maxQ;
        % model-based value Q(state1,action1) = P(state2|action1)*max(Q(state2,action1),Q(state2,action2)) + P(state3|action1)*max(...)
       
        Qm1 = squeeze(sum(reshape(maxQ',[2,1,S]).*TransitionProbs1))';
        Qm2 = squeeze(sum(reshape(maxQ',[2,1,S]).*TransitionProbs2))';
        Qm_state1=Qm1.*p1+Qm2.*p2;
        
        lik=lik+b1.*Qm_state1(:,c1)+b2.*Qd1(:,c1)+b2a.*Qd1a(:,c1)+st.*M(:,c1)- mfUtil.logsumexp([bsxfun(@times,b1,Qm_state1)+bsxfun(@times,b2,Qd1)+bsxfun(@times,b2a,Qd1a)+bsxfun(@times,st,M)],2);
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
        
        
        %Update evidence in favor of T1
        if s-1==1
            if c1==1
                alpha=alpha+1;
            else
                beta=beta+1;
            end
        else
            if c2==2
                alpha=alpha+1;
            else
                beta=beta+1;
            end
        end
        
        %update q-values in model-free system
        if s==2
            Qd1(:,c1)=Qd1(:,c1)+lr1.*(Qd2(:,c2)-Qd1(:,c1));
            Qd1a(:,c1)=Qd1a(:,c1)+lr1.*(o-Qd1a(:,c1));
            Qd2(:,c2)=Qd2(:,c2)+lr2.*(o-Qd2(:,c2));
        elseif s==3
            Qd1(:,c1)=Qd1(:,c1)+lr1.*(Qd3(:,c2)-Qd1(:,c1));
            Qd1a(:,c1)=Qd1a(:,c1)+lr1.*(o-Qd1a(:,c1));
            Qd3(:,c2)=Qd3(:,c2)+lr2.*(o-Qd3(:,c2));
        end
   
    end
end
