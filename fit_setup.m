function [LBpick, UBpick, LowBpick, UppBpick, start_parms, numstarts] = fit_setup(modelnum, numpars, do_scaling)
%this completes the set up of the models for fitting
%@@@need to add description of fxn parameters here

n=numpars;

switch modelnum
    case 1
        %random biassed choice
        LB=0.01;
        UB=1;
        %set upper and lower bounds on fitted parameter
        if do_scaling==0
            LowB=0.01;
            UppB=1;
        elseif do_scaling==1
            LowB=-Inf;
            UppB=+Inf;
        end;
    case {2,3}
        %q learning with 2 parameters alpha(p) and beta plus alphan for
        %lower and upper bounds for alpha(p), beta and alphan respectively
        LB=[0.01 0.01 0.01];
        UB=[1     5    1];
        %set upper and lower bounds on fitted parameter
        if do_scaling==0
            LowB=[0.01 0.01 0.01];
            UppB=[1     5    1];
        elseif do_scaling==1
            LowB=[-Inf -Inf -Inf];
            UppB=[+Inf +Inf +Inf];
        end;     
    otherwise
        error('No model with this model number code exists yet');
end;

LBpick=LB(1:n);
UBpick=UB(1:n);
LowBpick=LowB(1:n);
UppBpick=UppB(1:n);

%start_parms=[-1.12 -3.00; -.02 0.2; 1.085 1.8]; %3 sets of random scaled starting values 
numstarts=3^numpars; %the number of starts values to do the fitting from
start_parms=zeros(numstarts,numpars);
param_bound=10;
for s=1:numstarts
    for k=1:numpars
        mysign=1;
        if rand<0.5
            mysign=-1;
        end;
        start_parms(s,k)=rand*param_bound*mysign;
    end
end;



