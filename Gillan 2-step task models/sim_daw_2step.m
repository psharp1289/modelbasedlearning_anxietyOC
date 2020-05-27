function [S,A,R,Tm,lr1,w,b1,b3,lr_transition,st] = sim_daw_2step(lr1,lr2,lambda,weight,invtemp_mb,invtemp_2ndstage,lr_transition,st,ntrials)


% x is params
lr1=lr1;
lr2=lr2;
lambda = lambda;     % eligibility trace decay
b1 = invtemp_mb;
b3 = invtemp_2ndstage;
st = st;
w=weight;

tr = 0.7; % transition prob common

%2-step classic set-up
nstates = 3;
nactions = 2;


% initialization
Qd = zeros(nstates,nactions); 	% Q(s,a): state-action value function for model-free Q-learning
S = zeros(ntrials,1); 			% second step states
A = zeros(ntrials,2); 			% chosen actions at first and second steps
R = zeros(ntrials,1); 			% second step rewards
% TransitionProbs = [0.5 0.5;0.5 0.5]; 			% transition counts
Tm = [...% s1a1 s1a2
            0.5 0.5; % s2
            0.5 0.5];% s3   % transition matrix

M = [0; 0];

rewards=[0.25 0.75; 
        0.75 0.25];
% probs transitioning to state 3
% from 	    s1a1   s1a2
pTrans3 = [(1-tr) tr];

% loop through trials
for t = 1:ntrials

	%  model-based and model-free share second-level action values (doesn't have to be that way, but is in this implementation)
	%    max value of Q for each second level state state
	%model-based value Q(state1,action1) = P(state2|action1)*max(Q(state2,action1),Q(state2,action2)) + P(state3|action1)*max(...)

	maxQ = max(Qd(2:3,:),[],2);
    Qm=Tm'*maxQ;
    Qm=Qm';
    Qweighted=w.*Qm+(1-w).*Qd(1,:)+st.*M;

% initialization
	if rand <= exp(b1.*Qweighted(1))./sum(exp(b1.*Qweighted))
 		a(1) = 1;
        unchosen = 2; %unchosen action
    else
        a(1) = 2;
        unchosen = 1; %unchosen action
	end

	% store first level choice (this is used for next trial, stickiness parameter)
	M = [0 0];
	M(a(1)) = 1;

	% transition to next state based on action chosen and T (to state 2 or 3)
	s = (rand<=pTrans3(a(1)))+2;
    if s==3
        other_s=2;
    else
        other_s=3;
    end

	% choose second level action
	if rand < exp(b3*Qd(s,1))/sum(exp(b3*Qd(s,:)))
		a(2) = 1; else a(2) = 2;
	end

	% update first choice Q value
	Qd(1,a(1)) = Qd(1,a(1)) + lr1*(Qd(s,a(2))-Qd(1,a(1)));

	% update transition matrix (RW rule, where sum of probabilities for either state for a given action = 1)
    Tm(s-1,a(1)) = Tm(s-1,a(1)) + lr_transition*(1-Tm(s-1,a(1)));
    Tm(other_s-1,a(1)) = 1-Tm(s-1,a(1)); 
    
    
    if mod(t,50)==0
        rewards=[rewards(:,2) rewards(:,1)];
    end
%     rewards(s-1,a(2));
    r = rand <= 0.5;

	% update second choice Q value
	Qd(s,a(2)) = Qd(s,a(2)) + lr2.*(r-Qd(s,a(2)));

	%  apply second level prediction error to first level choice based
 	Qd(1,a(1)) = Qd(1,a(1)) + lambda.*lr2.*(r-Qd(1,a(1)));

	% store trial results
	A(t,:) = a;
	R(t,1) = r;
	S(t,1) = s;
    
end
