function [neg2LL, VAF, neg2LL_gzero]=rbc_fxn(parms,modelnum, Ss_choices,nc,br,rewards,scale,lb,ub)

%v1.01; 14.08.19
%output paramters
%  neg2LL       -2LogLike
%  VAF          variance accounted for (percentage right choices using p=0.5
%               as threshold)
%  neg2LL_gzero -2LogLike for a model in which 1/n prob for every choice
%               (for n choices)
%
%input parameters (some of which ** are not used here)
%we need to jkeep the input params constant so we can call all model
%functions using the same code with just the function name changing
%  parms     model parameters in a row vector
%  modelnum  number code for model
%  Ss_choice subject's response choices 1, 2 etc
%  nc        number of response choices
%  br        biggest reward available **
%  reward    reward on a specific trial **
%  scale     have the paramters been (re-)scaled 1=yes; 0=no
%  lb        lower bounds of the parameters
%  ub        upper bounds of the parameters
%
%this fxn calls the following function(s) 
%
%loglike

if nc~=2
    %error message if you have more than 2 choices
    disp('Developed for 2 choices only at the moment.')
    error;
end;

if modelnum~=1
    %error message if you have more than 2 choices
    disp('Model number wrongly specified.')
    error;
end;


%convert scaled parms to alpha and beta within acceptable ranges
np=size(parms,2);
%        p_resp1=parms(1);
for k=1:np
    if scale==1
        parms(k)=(ub(k)-lb(k))/(1+exp(-parms(k)))+lb(k);
    end;
end;
%parms
%pause

%loop thru all the responses
nr=size(Ss_choices,1); %a is the vector of subject's choices on each trial; there are nr responses
smyhat=zeros(nr,1);  %initialise an array of predicted values under the softmax model

for t=1:nr %loops through the choices of each subject
   smyhat(t)=parms(1); %this is the probability of making choice 1
end

%compute -2*loglikelihood
%the next line uses the function loglike
%it matches choices of 1 (a==1, otherwise 0) against smyhat (the
%probability of making choice 1)
%smyhat
%pause;
neg2LL=-2.*loglike(Ss_choices==1,smyhat);  %need to minimize fit index, which is -2*log likelihood. This is same as maximizing log likelihood

%do this computation, of the variance accounted for, only if 2 or more output arguments are required
if nargout>1 
    %compute variance accounted for
    VAF=sum((smyhat >= 0.5) == (Ss_choices==1))./nr ;
end

%compute neg2LL for zero free param guessing model 
if nargout==3
    phalf=0.5.*ones(nr,1);
    neg2LL_gzero = -2.*loglike(Ss_choices==1,phalf);
end;
