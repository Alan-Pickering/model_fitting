%getdata_and_callfitting_v3.m
%
%a data extractor and model fitter
%v3; 28.5.2020
%
%requires the following user created functions in the same folder as this
%programme
%
%   myfilesel.m      %%this is a simple file selector utility for reading and
%                      writing Excel files
%   dirchecker.m     %%checks the directory that a script is launched from to
%                      avoid path issues
%   fit_setup.m      %%this controls the setup of the model being fit
%   loglike.m        %%computes the loglikelihod of the data given the model
%   model_lister.m   %%lists the available models to fit and gives them
%                      numbers
%   BMC_v4.m         %%Bayesian Model comparison routine for use with the VBA
%                      toolbox
%
%   The various models in the model lister
%   QL_fxn.m         %%classic Qlearning with softmax (models 2 and 3)
%   rbc_fxn.m        %%a random biassed choice (model 1)
%
%Other important points
% -does fitting using MATLAB's functions fmincon or fminsearch
% -fits data stored in an Excel file with each trial per subject in a
%  separate row and certain key columns of data present (choice, reward etc)
%  all data other than column headers should be numeric
%  text column labels must be present to allow user to select the right columns
%  for fitting


%Plans for later versions
%************************
%a) Store best fitting parameters for each subject for each model
%b) Try to sort out the exit flag record (by allowing more iterations)
%c) Add further models to fit 
%d) Add choice of what to fit against 
%   eg add fit against Free energy 
%   rather than using simple MLE -- to do this we need to use 
%   the hessian and do the fit a number of times for each subject 
%   and or call the VBA toolbox / check against their results
%e) Add choices for options (eg fitting features)
%f) display pairwise model comparisons
%g) Add flagging of edge and other cases for exclusion
%h) Check the ability of fmincon to get same results
%i) Record Hessian from fits
%j) Add stimulus data in modelling, including stim sequences
%k) Add ability to model tasks without an overt response
%l) Add colourmap for pxp and routines to save EFmat, PXPMat
%
%@@@things to note/check/improve are marked like this below
%

clear variables;
rng('default');

clc;
disp('Your first task will be to select a performance dataset,');
disp('containing all participants with a single line per trial');
disp('for each participant.');
disp(' ');
disp('Hit any key to continue');
pause;
clc;

%these directories point to key folders
%@@@could replace these by uigetdir
homeDirectory='C:\Users\Alan Pickering\Google Drive\matlab code\work in EP docs Matlab\APfitting';
DataInDirectory='C:\Users\Alan Pickering\Google Drive\matlab code\work in EP docs Matlab\worthy models and data';
DataOutDirectory='C:\Users\Alan Pickering\Google Drive\matlab code\work in EP docs Matlab\APfitting\stored_fits';
currDir=pwd;
dirchecker(currDir,homeDirectory);

%here is where we select and then read in the Excel datafile to be fit
[myfilename, myinputpathname, filtindex] = myfilesel('*.xlsx', DataInDirectory, 'read');

if filtindex==1
    [rawdata, myheaders, ~]=xlsread(myfilename);
end

%read in some key data features that are required for fitting
[num_data_rows, num_inp_vars]=size(rawdata);

disp('Here are the column labels in your performance datafile:-');
for j=1:num_inp_vars
    disp([num2str(j) ' ... ' myheaders{j}]);
end

subcodes=0;
while (subcodes<1 || subcodes>num_inp_vars)
    subcodes=input('Type the column number used for unique subject code (+<Enter>) ');
end
choicecodes=0;
while (choicecodes<1 || choicecodes>num_inp_vars)
    choicecodes=input('Type the column number used for response choices (+<Enter>) ');
end
rewcodes=0;
while (rewcodes<1 || rewcodes>num_inp_vars)
    rewcodes=input('Type the column number used for rewards (+<Enter>) ');
end

%next we can compute how many subjects there are
%by counting the number of unique subject codes
nsubs=size(unique(rawdata(:,subcodes)),1);
scodelist=unique(rawdata(:,subcodes)); %list of the subject codes used

cd(homeDirectory); %need to be back in this directory for fitting in case path not set up correctly
warning('OFF', 'all'); %turn off some annoying non-serious warnings in Matlab that show up during fitting

clc;
%select how many models you want to fit
nm2fit=0;
maxnum_avail=3; %the maximum number of models currently available
    
while nm2fit<1 || nm2fit>maxnum_avail 
    disp(['You can fit between 1 and ' num2str(maxnum_avail) ' models.']);
    nm2fit=input('Enter number of models to fit (+<Enter>) ');
end
models2comp.BIC=zeros(nsubs, nm2fit+1); %add 1 for random guessing model
models2comp.names=cell(nm2fit+1,1);
models2comp.nparms=zeros(nm2fit+1,1);

for m=1:nm2fit
    clc;
    %select the model to be fit here from screen list
    [fitted_model_num, model_name, mod_fxn_name, numparms, parm_labs]  = model_lister(m);
    %fitting features that are fixed for all Ss
    %@@@select number of choices here
    nch=2; %number of resp choices available
    
    %@@@select how much data to store in the output file of fits
    n_fitted_vars=14;
    nvars=n_fitted_vars+numparms; %number of variables recorded for each subject
    wdata=zeros(nsubs,nvars,m); %this is where we will store the fitting data in one row per subject
    
    %@@@select various fitting options
    %@@@first the fitting method
    fitmeth=1; %1=max like estimation
    %@@@and then the initial reward values
    %initvaltype controls the way the initial reward prediction values are controled
    %  1 selects initial expected rewards to be half maximum possible
    %  2 selects zero as the initial expected rewards
    initvaltype=2;
    %@@@next choose which Matlab fitting routine to use
    con_or_srch=2; %controls which Matlab function minimisation routine is used 1=fmincon, 2=fminsearch
    %@@@and now the largest reward available (may often just be 1)
    bigrew=1; %this value is used for fittype=1
    
    %@@@next decide whether to (re)scale bounded values to unbounded
    if con_or_srch==2
        scale=1; %1= use scaling  0 = do not, scaling control has to be used for fminsearch
    elseif con_or_srch==1
        %free choice here
        scale=0; %1= use scaling  0 = do not
    end
    
    %now set up the fitting
    [LB, UB, LowB, UppB, parmmatrix, numstarts]=fit_setup(fitted_model_num, numparms, scale);
    warning('OFF', 'all'); %turn off some annoying non-serious warnings in Matlab that show up during fitting
    doshow=0; %print best fitting values subj by subj to screen; yes=1
    
    %go through the subject loop
    for s2=1:nsubs
        
        %read the subjects' data one subject at a time
        wdata(s2,1,m) = scodelist(s2); %records subject code
        
        %need to extract their choices (subchoices) and the rewards (rewardpts) they got for each trial
        %first extract the row numbers in the rawdata array for each subject code
        rownums=rawdata(:,subcodes)==scodelist(s2); %has value 1 for chosen subject
        rewardpts=rawdata(rownums==1,rewcodes);
        subchoices=rawdata(rownums==1,choicecodes);
        ntrials=sum(rownums); %number of trials for this subject
        
        disp(['Data extracted for subject code number ' num2str(scodelist(s2))]);
        disp('now fitting their data');
        %now here we try fitting models to the data
        bestfit=1000000; %arbitrarily high ie bad log(ML) value which we will improve upon by fitting
        % will fit the model once for each starting parameter combination
        for i=1:numstarts
            parm0=parmmatrix(i,:);
            if scale==1
                %leave parm0 alone
            else
                %parm0=[0.4 0.5]; %gives reasonable values, these values are alpha and beta respectively, but unscaled
                %or get the unscaled but bounded parms from the scaled values above
                for k=1:numparms
                    parm0(k)=(UB(k)-LB(k))/(1+exp(-parm0(k)))+LB(k);
                end
            end
            %parmbest=parm0;
            %now do fitting for each subject individually
            if con_or_srch==1
                %use fmincon for 2 different types of fitting (see above)
                if initvaltype==1
                    %Here we use MATLAB's constrained function minimisation function
                    %which returns parameter estimates into parmest and a fit
                    %value as fit. It fits the user supplied function QL2_fxn
                    %with initial parameters in parm0 and LowB and UppB as
                    %specified. [] are placeholders for unused options, defaults
                    %will apply. At the end of the list of parameters -- from fitted_model_num
                    %onwards -- are the parameters to pass to the function
                    %@@@replace function string with mod_fxn_name
                    [parmest, fit, exitflag]=fmincon(mod_fxn_name,parm0,[],[],[],[],LowB,UppB,[],[],fitted_model_num,subchoices,nch,bigrew,rewardpts,scale,LB,UB);
                elseif initvaltype==2
                    %note the diff value of the 3rd parameter passed
                    %@@@replace function string with mod_fxn_name
                    [parmest, fit, exitflag]=fmincon(mod_fxn_name,parm0,[],[],[],[],LowB,UppB,[],[],fitted_model_num,subchoices,nch,0,rewardpts,scale,LB,UB);
                end
            elseif con_or_srch==2
                %here we use fminsearch which has a slightly different syntax
                %and is unbounded, so we force bounding by scaling the
                %search range to the same bounds that we supply directly in fmincon
                %again 2 types of fit can be chosen, depending on initial
                %expectations assumptions
                if initvaltype==1
                    mystring=['[parmest, fit, exitflag]=fminsearch(@(x) ' mod_fxn_name '(x,fitted_model_num,subchoices,nch,bigrew,rewardpts,scale,LB,UB), parm0);'];
                    %[parmest, fit, exitflag]=fminsearch(@(x) QL_fxn(x,fitted_model_num,subchoices,nch,bigrew,rewardpts,scale,LB,UB), parm0);
                elseif initvaltype==2
                    mystring=['[parmest, fit, exitflag]=fminsearch(@(x) ' mod_fxn_name '(x,fitted_model_num,subchoices,nch,0,rewardpts,scale,LB,UB), parm0);'];
                    %[parmest, fit, exitflag]=fminsearch(@(x) QL_fxn(x,fitted_model_num,subchoices,nch,0,rewardpts,scale,LB,UB), parm0);
                end
                eval(mystring);
            end
            %recording the best fit value
            if fit<=bestfit
                bestfit=fit;
                parmbest=parmest;
                startcond=i;
            end
        end
        
        %computing fit indices etc based on best fitting parameters in order to
        %calculate VAF -- variance accounted for
        %and calculate -2LL for the model (fit)
        %and calculate -2LL for a zero-parameter guessing model (fit_g)
        lnlike=-0.5.*bestfit; %compute log-likelihood
        if initvaltype==1
            mystring=['[fit, vaf, fit_g]=' mod_fxn_name '(parmbest,fitted_model_num,subchoices,nch,bigrew,rewardpts,scale,LB,UB);'];
            %[fit, vaf, fit_g]=QL2_fxn(parmbest,subchoices,nch,bigrew,rewardpts,scale,LB,UB);
        elseif initvaltype==2
            mystring=['[fit, vaf, fit_g]=' mod_fxn_name '(parmbest,fitted_model_num,subchoices,nch,0,rewardpts,scale,LB,UB);'];
            %[fit, vaf, fit_g]=QL2_fxn(parmbest,subchoices,nch,0,rewardpts,scale,LB,UB);
        end
        eval(mystring);
        %check fit still same
        if fit ~=bestfit
            disp('uh oh; fitting result has done something very odd');
            error('You need to investigate');
        end
        %now compute BICs
        BIC_model=numparms.*log(ntrials) + fit;
        BIC_guess=fit_g;
        pmod_given_data=1./(1 +exp(-0.5.*(BIC_guess-BIC_model)));
        BFmg=exp(0.5.*(BIC_guess-BIC_model)); %Bayes factor for the model cf the guessing model
        
        if doshow==1
            disp('')
            disp('Best fitting parameters');
            disp(parmbest);
            disp('Log likelihood');
            disp(lnlike); %this is log-likelihood
            disp('%var accounted for by model across items');
            disp(vaf);
        end
        
        %convert scaled / raw parameters to alpha and beta
        if scale==1
            for k=1:numparms
                parmbest(k)=(UB(k)-LB(k))/(1+exp(-parmbest(k)))+LB(k);
            end
        end
        
        %record results of fitting into row of data for each subject
        wdata(s2,2,m)= fitted_model_num; %records model used
        wdata(s2,3,m)= lnlike;
        wdata(s2,4,m)= vaf;
        wdata(s2,5,m)= BIC_model;
        wdata(s2,6,m)= BIC_guess;
        wdata(s2,7,m)= BFmg;
        wdata(s2,8,m)= pmod_given_data;
        wdata(s2,9,m)= exitflag; %records the exit flag @@@ not working properly
        wdata(s2,10,m)= fitmeth; %records fitting method used
        wdata(s2,11,m)= startcond; %records which parameter starting condition worked best for this dataset
        wdata(s2,12,m)= initvaltype; %records the version of MLE fitting done
        wdata(s2,13,m)= con_or_srch; %records whether fmincon or fminsearch was used
        wdata(s2,14,m)= scale; %records whether scaling was used
        for k=1:numparms
            wdata(s2,n_fitted_vars+k,m)= parmbest(k);
        end
        
        %check
        checknvars=n_fitted_vars+numparms;
        if nvars~=checknvars  %nvars defined much earlier in the pgm
            disp('The number of variables doesn''t match');
            disp('Check the code carefully.');
        end
        
        %end of s1 loop
    end
    
    model_BIC=['BIC_' model_name];
    vlabs={'subcode','model_num','lnlike','vaf',model_BIC,'BIC_guess',...
        'Bayes_Factor_mg','p_mod_given_data', 'exitflag', 'fitting_method','start_cond',...
        'initval_type','con_or_search','scaled'};
    
    %and now add the labels for the parameters
    for k=1:numparms
        vlabs{n_fitted_vars+k} = parm_labs{k};
    end
    
    %write the data into a cell string array, alldata, with variable labels, vlabs, in the
    %first row.
    n_c=nvars;
    n_r=nsubs+1; %+1 for labels
    %set up empty cell string array, with the right dimensions
    alldata=cell(n_r,n_c);
    for r=1:n_r
        for c=1:n_c
            if r==1
                alldata{r,c}=vlabs{c};
            else
                alldata{r,c}=wdata(r-1,c,m);
            end
        end
    end
    
    %now select file and save fitting data to it
    [myoutfilename, myoutpathname, filtindex ] = myfilesel('savefits.xlsx', DataOutDirectory, 'save');
    
    if filtindex==1
        xlswrite(myoutfilename,alldata);
    end

    if m==1
        models2comp.BIC(:,1)=wdata(:,6,1);
        models2comp.names{1}='Random guess';
        models2comp.nparms(1)=0;
        models2comp.mn(1)=0; %model number 0 for random guess model
    end
    models2comp.BIC(:, m+1)=wdata(:,5,m); %BIC value
    models2comp.names{m+1}=model_name;
    models2comp.nparms(m+1)=numparms;
    models2comp.mn(m+1)=fitted_model_num;
    
%end of loop through various models being fitted    
end


disp(' ');
disp('Hit any key to compare pairs of models');
pause;
clc
%here we loop over set of models we want to compare pairwise
EFmat=zeros(nm2fit+1,nm2fit+1); %this is where we store the matrix of EFs
modnumgrid=zeros(nm2fit+1,nm2fit+1,2); %this is where we store the matrix of EFs
%@@@ could store a matrix for pxp too
for k=1:nm2fit
    for j=k+1:nm2fit+1
        %models2comp.mset=[models2comp.mn(k) models2comp.mn(j)]; %we can do a loop here to make all pairwise comparisons
        models2comp.mset=[k j]; %we can do a loop here to make all pairwise comparisons
        
        %and here we actually compare pairwise
        [my_py, my_oy] = BMC_v4(models2comp, 1, 'e', 0); %final zero suppresses graphic output
        
        modnumgrid(k,j,1:2)=[models2comp.mn(k) models2comp.mn(j)];
        EFmat(k,j)=my_oy.Ef(2); %stored the EF for the higher numbered model
         
    end
end

%duplicate the lower triangle of Efmat and display Efmat as colour matrix 
matsize=size(EFmat,1);
LowEFmat=tril(EFmat+(1-triu(EFmat)'))-(eye(matsize));
allEF=EFmat+LowEFmat;
qq=NaN(matsize); %to be used to set diagonal to white with pcolor
allmat=allEF+diag(diag(qq));
zz=[allmat zeros(matsize,1); zeros(1, matsize+1)]; %adjust so that right bit of matrix is graphed
figure;
shading('interp');
pcolor(zz);
axis ij;
axis square;
%the next quoted lines are one way to do this
%ax=gca;
%ax.XTick=[1.5:1:(nm2fit+1.5)];
%ax.XTickLabel=models2comp.names;
%ax.YTick=[1.5:1:nm2fit+1.5];
%ax.YTickLabel=models2comp.names;
%but here is another
set(gca,'XTick',1.5:1:(nm2fit+1.5));
set(gca,'YTick',1.5:1:(nm2fit+1.5));
set(gca,'XTicklabel',models2comp.names);
set(gca,'YTicklabel',models2comp.names);
colorbar;
title({'Expected frequencies of model vs model:', 'Column model vs. row model'});

%exit
warning('ON', 'all'); %turn warnings back on
cd(homeDirectory);