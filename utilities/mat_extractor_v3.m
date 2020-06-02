%mat_extractor_v3.m
%
%version 3; 1.6.20
%
%reads in multi-subject data files from .mat format
%and stores them in the right format 
%and layout for use with the AP fitting suite
%or the Piray fitting suite CBM
%
%this script uses the fxns
%
%dirchecker.m
%myfilesel_cs.m
%
%see lines of code which start %@@@ for more info

clc;
clear variables;

%% here should go a displayed description of what the pgm does currently
%-reads in multi-subject .mat files containing choices and actions
%-the data is stored in (a) cell array variable(s) in the
% .mat file
%each subject's data is in a different row or column element of the overall
% data cell array(s)
%-each subject must have the same number of trials
%-pgm saves the data to a file format appropriate for a specific model-fitting suite


%% you need to alter these fields for us on your particular set=up
homeDirectory='C:\Users\Alan Pickering\Google Drive\matlab code\work in EP docs Matlab\RL models\cbm_hbi\cbm-master\example_RL_task';
%homeDirectory='C:\Users\Alan Pickering\Google Drive\matlab code\work in EP docs Matlab\stolz_data';
%DataInDirectory='C:\Users\Alan Pickering\Google Drive\matlab code\work in EP docs Matlab\stolz_data\mat files';
%DataOutDirectory='C:\Users\Alan Pickering\Google Drive\matlab code\work in EP docs Matlab\stolz_data\excelfiles';
DataInDirectory='C:\Users\Alan Pickering\Google Drive\matlab code\work in EP docs Matlab\RL models\cbm_hbi\cbm-master\example_RL_task';
DataOutDirectory='C:\Users\Alan Pickering\Google Drive\matlab code\work in EP docs Matlab\RL models\cbm_hbi\cbm-master\example_RL_task';
currDir=pwd;
dirchecker(currDir,homeDirectory); 

%% @@@this is just to test the structure variable in selection
%kk.first=99;
%kk.second='hello';

%% here is where we select and then read in the .mat datafile to be
%  converted and saved in a different format
[myfilename, myinputpathname, filtindex] = myfilesel_cs('*.mat', DataInDirectory, 'read');
if filtindex==1
    load(myfilename);
end

%% @@@this is just to test other kinds of data
%first a cell array which has mutltiple rows and columns
%actplusout=cell(2,25);
%actplusout(1,:)=REW_choice(1,:);
%actplusout(2,:)=REW_feedback(1,:);
%and now to test a transposed cell array
%REW_choice=REW_choice'; 
%and now to set up cbm style data with multiple arrays per subject
%twodata=cell(40,2);
%twodata(:,1)=data(:,1);
%twodata(:,2)=data(:,1);

%% now we show the variables in the workspace
disp('Here is a list of variables in the Matlab workspace,');
disp('including those from the loaded file.');
xx=who;
nv=size(xx,1);
ctr=0;
vtype=zeros(nv,1);
for i=1:nv
    ctr=ctr+1;
    disp([num2str(ctr) '... ' xx{i}]);
    if isstruct(eval(xx{i}))
        vtype(ctr)=2; %denoting structure variable
        disp(['>>>Variable ' xx{i} ' is a structure with the following fields:-']);
        yy=fieldnames(eval(xx{i}));
        for k=1:size(yy,1)
            ctr=ctr+1;
            vtype(ctr)=3; %denoting structure field
            disp([num2str(ctr) '... ' yy{k}]);
        end
        disp('>>>End of structure variable');
    else
        vtype(ctr)=1; %denoting non-structure variable
    end
end

%% and now we select the variable which logs choices
disp('Please select the variable or structure field containing the choices made');
ch_num=0;
while ch_num<1 || ch_num>nv
    ch_num=input('Enter number of variable/field from the above list +<Enter> ');
end
disp(' ');
%now construct names for the choice (ch_name)
if vtype(ch_num)==1
    ch_name=xx{ch_num};
elseif vtype(ch_num)==2
    disp(char({'You selected a structure variable,', 'but the data are almost in a structure field!?', 'Please re-run pgm and make a different selection.'}));
    return;
elseif vtype(ch_num)==3
    %@@@need to do soemthing more complex for fields
end

%% checks on choices made for actions
[v1, v2]=size(eval(ch_name));
disp('The size of this choice array');
disp(['over all participants is: ' num2str(v1) ' (rows) & ' num2str(v2) ' (columns).']);
disp('Which dimension, row (#1) or column (#2), shows no. of participants?');
whichissub_ch=0;
while whichissub_ch~=1 && whichissub_ch~=2
    whichissub_ch=input('Make choice (1 or 2) followed by <Enter> ');
end
disp(' ');
nsubs=size(eval(ch_name),whichissub_ch); %extract number of subjects
if whichissub_ch==1
    rc='column';
else
    rc='row';
end

%% now check that the other cell array dimension is >1
if size(eval(ch_name),3-whichissub_ch)>1
    %other dimension is >1
    disp(['There are ' num2str(size(eval(ch_name),3-whichissub_ch)) ' ' rc 's in the data cell array you have selected.']);
    disp(['Select the ' rc ' that contains the choices made.']);
    whichrowcol_ch=0;
    while whichrowcol_ch<1 || whichrowcol_ch>size(eval(ch_name),3-whichissub_ch)
        whichrowcol_ch=input(['Enter no. of ' rc ' containing choices ']);
    end
    disp(' ');
else
    whichrowcol_ch=1;
end

%% next check what kind of variables in the cell array are numeric or structure
isnum_ch=isnumeric(eval([ch_name '{1,1}']));
if isnum_ch==1
    [v1, v2]=size(eval([ch_name '{1}']));
    disp('The size of the individual participant choice array');
    disp(['over trials is: ' num2str(v1) ' (rows) & ' num2str(v2) ' (columns).']);
    disp('Which dimension, row (#1) or column (#2), shows no. of trials?');
    whichistri_ch=0;
    while whichistri_ch~=1 && whichistri_ch~=2
        whichistri_ch=input('Make choice (1 or 2) followed by <Enter> ');
    end
    disp(' ');
    ntrials=size(eval([ch_name '{1}']),whichistri_ch); %extract number of trials
end
isstruck_ch=isstruct(eval([ch_name '{1}']));
ctr=0;
if isstruck_ch==1
    disp(['>>>Variable ' ch_name '{n} is a structure with the following fields:-']);
    yy=fieldnames(eval([ch_name '{1}']));
    for k=1:size(yy,1)
        ctr=ctr+1;
        disp([num2str(ctr) '... ' yy{k}]);
    end

    whichfield_ch=0;
    while whichfield_ch<1 || whichfield_ch>size(yy,1)
        whichfield_ch=input('Type number of field containing actions +<Enter> ');
    end
    [v1, v2]=size(eval([ch_name '{1}.' char(yy{whichfield_ch})]));

    disp(' ');
    disp('The size of the individual participant choice array');
    disp(['over trials is: ' num2str(v1) ' (rows) & ' num2str(v2) ' (columns).']);
    disp('Which dimension, row (#1) or column (#2), shows no. of trials?');
    whichistri_ch=0;
    while whichistri_ch~=1 && whichistri_ch~=2
        whichistri_ch=input('Make choice (1 or 2) followed by <Enter> ');
    end
    ntrials=size(eval([ch_name '{1}.' char(yy{whichfield_ch})]),whichistri_ch); %extract number of trials
    disp(' ');
end

%% now do choice for outcome/feedback (fb)
disp('Please select the variable or structure field containing the outcomes received ');
fb_num=0;
while fb_num<1 || fb_num>nv
    fb_num=input('Enter number of variable/field from the above list +<Enter> ');
end
disp(' ');
%now construct names for the outcome(fb_name)
if vtype(fb_num)==1
    fb_name=xx{fb_num};
elseif vtype(fb_num)==2
    disp(char({'You selected a structure variable,', 'but the data are almost in a structure field!?', 'Please re-run pgm and make a different selection.'}));
    return;
elseif vtype(fb_num)==3
    %@@@need to do soemthing more complex for fields
end

%% here is where we will duplicate the checks for the outcomes/feedback
[v1, v2]=size(eval(fb_name));
disp('The size of this outcome array');
disp(['over all participants is: ' num2str(v1) ' (rows) & ' num2str(v2) ' (columns).']);
disp('Which dimension, row (#1) or column (#2), shows no. of participants?');
whichissub_fb=0;
while whichissub_fb~=1 && whichissub_fb~=2
    whichissub_fb=input('Make choice (1 or 2) followed by <Enter> ');
end
disp(' ');
nsubs_chk=size(eval(fb_name),whichissub_fb);
if whichissub_fb==1
    rc='column';
else
    rc='row';
end

%% now check that the other cell array dimension is >1
if size(eval(fb_name),3-whichissub_fb)>1
    %other dimension is >1
    disp(['There are ' num2str(size(eval(fb_name),3-whichissub_fb)) ' ' rc 's in the data cell array you have selected.']);
    disp(['Select the ' rc ' that contains the outcomes experienced.']);
    whichrowcol_fb=0;
    while whichrowcol_fb<1 || whichrowcol_fb>size(eval(fb_name),3-whichissub_fb)
        whichrowcol_fb=input(['Enter no. of ' rc ' containing outcomes ']);
    end
    disp(' ');
else
    whichrowcol_fb=1;
end

%% check type of variable in the cell array as for actions
isnum_fb=isnumeric(eval([fb_name '{1,1}']));
if isnum_fb==1
    [v1, v2]=size(eval([fb_name '{1}']));
    disp('The size of the individual participant outcome array');
    disp(['over all trials is: ' num2str(v1) ' (rows) & ' num2str(v2) ' (columns).']);
    disp('Which dimension, row (#1) or column (#2), shows no. of trials?');
    whichistri_fb=0;
    while whichistri_fb~=1 && whichistri_fb~=2
        whichistri_fb=input('Make choice (1 or 2) followed by <Enter> ');
    end
    disp(' ');
    ntrials_chk=size(eval([fb_name '{1}']),whichistri_fb); %extract number of trials
end
isstruck_fb=isstruct(eval([fb_name '{1}']));
ctr=0;
if isstruck_fb==1
    disp(['>>>Variable ' fb_name '{n} is a structure with the following fields:-']);
    yy=fieldnames(eval([fb_name '{1}']));
    for k=1:size(yy,1)
        ctr=ctr+1;
        disp([num2str(ctr) '... ' yy{k}]);
    end
    whichfield_fb=0;
    while whichfield_fb<1 || whichfield_fb>size(yy,1)
        whichfield_fb=input('Type number of field containing outcomes +<Enter> ');
    end
    [v1, v2]=size(eval([fb_name '{1}.' char(yy{whichfield_fb})]));

    disp(' ');
    disp('The size of the individual participant outcomes array');
    disp(['over trials is: ' num2str(v1) ' (rows) & ' num2str(v2) ' (columns).']);
    disp('Which dimension, row (#1) or column (#2), shows no. of trials?');
    whichistri_fb=0;
    while whichistri_fb~=1 && whichistri_fb~=2
        whichistri_fb=input('Make choice (1 or 2) followed by <Enter> ');
    end
    ntrials_chk=size(eval([fb_name '{1}.' char(yy{whichfield_fb})]),whichistri_fb); %extract number of trials
    disp(' ');
end

%% now check that selected variables all match in size
if nsubs~=nsubs_chk
    disp('No. of subjects in stored data arrays for actions and outcomes don''t match');
    return;
end
if ntrials~=ntrials_chk
    disp('No. of trials in stored data arrays for actions and outcomes don''t match');
    return;
end

%% now choose the format appropriate to the fitting suite to be used
fsuites={'Alan P''s fitting suite','Piray''s CBM fitting suite'};
nforms=size(fsuites,2);
disp('Here is a list of fitting suite file formats we can use: -');
for i=1:nforms
    disp([num2str(i) ' ... ' char(fsuites(i))]);
end
whatformat=0;
while whatformat<1 || whatformat>nv
    whatformat=input('Enter number of fitting suite from the above list +<Enter> ');
end
disp(' ');

%% and now create and save the new datafile in the appropriate format
if whatformat==1
    %format for AP fitting suite
    %now set up arrays for all data combined
    numvars=4; %subj, trialnum, choice, outcome
    %and now fill the arrays
    nrows=nsubs.*ntrials +1; % +1 for header
    alldata=cell(nrows,numvars);
    %add headers
    alldata(1,:)={'subj','trial','choice','outcome'};
    for i=1:nsubs
        rowctr_s=1+(i-1).*ntrials;
        rowctr_e=rowctr_s+ntrials - 1;
        if isstruck_ch==1
            choices=eval([ch_name '{' num2str(i) ',' num2str(whichissub_ch) '}.' char(yy{whichfield_ch})]);
        elseif isnum_ch==1
            if whichissub_ch==1
                choices=eval([ch_name '{i,whichrowcol_ch}']); %note transpose, as trials arrays assumed (1,numtrials)
            else
                choices=eval([ch_name '{whichrowcol_ch,i}']); %note transpose etc
            end
        end
        if isstruck_fb==1
            outcomes=eval([fb_name '{' num2str(i) ',' num2str(whichissub_fb) '}.' char(yy{whichfield_fb})]);
        elseif isnum_fb==1
            if whichissub_fb==1
                outcomes=eval([fb_name '{i,whichrowcol_fb}']); %note transpose etc
            else
                outcomes=eval([fb_name '{whichrowcol_fb,i}']); %note transpose etc
            end
        end
        trials=1:ntrials;
        trials=trials';
        subj=i.*ones(ntrials,1);
        for j=1:ntrials
            alldata{rowctr_s+j,1}=subj(j);
            alldata{rowctr_s+j,2}=trials(j);
            alldata{rowctr_s+j,3}=choices(j);
            alldata{rowctr_s+j,4}=outcomes(j);
        end
    end
    
    %now save as .xlsx file
    [myoutfilename, myoutputpathname, filtindex] = myfilesel_cs('*.xlsx', DataOutDirectory, 'save');
    
    if filtindex==1
        xlswrite(myoutfilename,alldata);
    end
    
elseif whatformat==2
    %format for Piray CBM suite
end

cd(homeDirectory);  


