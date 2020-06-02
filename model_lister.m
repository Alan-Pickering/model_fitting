function [modnum, modname, fxn_name, npars, p_names]  = model_lister(mm)
%simply lists the models available
%and allows user to pick one
%here are the model details
mymodels={'mod_num', 'mod_name',      'num_params',  'fxn_name'; ...
            1,      'Biassed choice',     1,         'rbc_fxn'
            2,      'QLearn2softmax',     2,         'QL_fxn';...
            3,      'QLearn3softmax',     3,         'QL_fxn';...
};

totnmodels=size(mymodels,1)-1;

%here are the parameter labels for the models above
%null indicates no parameter
parm_names={
            'p_resp1', 'null', 'null'; ...
             'alpha', 'beta', 'null'; ...
             'alpha+','beta', 'alpha-'...
};

disp('Here are the models available:-');
for j=1:totnmodels
    disp([num2str(j) ' ... ' mymodels{j+1,2}]);
end;

modnum=0;

modlabs='ABCDEFGHIJKLMNOPQRSTUV';
while (modnum<1 || modnum>totnmodels)
    modnum=input(['Model ' modlabs(mm) ' to be fit: Select by entering number (+<Enter>) ']);
end;

modname=mymodels{modnum+1,2};
npars=mymodels{modnum+1,3};
fxn_name=mymodels{modnum+1,4};

p_names=cell(1,npars);
for k=1:npars
    p_names{k}=parm_names{modnum,k};
end;

end

