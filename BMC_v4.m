function [py, oy] = BMC_v4 (m2c, edgeflagged_inc, bicbase, wantfigs)
% //requires VBA toolbox //////////////////////////////////////////////////////////
%
% BMC_v4 ()
% reads in BICs from fitted models and compares them
% returns key stats in py and oy
%
% version 2; including edge fit removals 
% if edgeflagged_inc =1 then we included the edge fit cases, if edge fit ne 1 then
% they are removed
% if bicbase is a number (eg 2), then we have to do some conversion ('e' means bics used base = e)
% if wantfigs=1; then display figs; other values then no figs
% /////////////////////////////////////////////////////////////////////////

nmodels=2;
totalparms=sum(m2c.nparms);
disp('The pair of models being compared here are as follows:-');
for k=1:nmodels
    disp(['Model name ' m2c.names{m2c.mset(k)}]);
end;

% number of subjects (need to adjust code to take account of edge excluded
% subjects
nSubjects=size(m2c.BIC,1);
logEvidence_y=zeros(nmodels,nSubjects);

i=0;
for j=1:nSubjects 
    if edgeflagged_inc==1
        i=j;
        excludecase=0;
    else
        %@@@this needs rewriting as it doesn't work as currently set
        %now check all params in all models are not edge fits
        edgectr=0;
        for k=1:nmodels
            for p=1:nparms(k)
                edge=fitdata(j+1,model_notedge_col(k,p));
                if edge{1}==1
                    edgectr=edgectr+1;
                end;
            end;
        end;
        if edgectr==totalparms
            %all params not edges
            i=i+1;
            excludecase=0;
        else
            excludecase=1;
        end;
    end;
    %loop for models for included cases
    if excludecase==0
        for k=1:nmodels
            %read in Bic values in file to provide log model evidences
            if strcmp(bicbase,'e')
            else
                %now convert to a number from cell array and convert to base e
                m2c.BIC(i, m2c.mset(k))=log(bicbase)*m2c.BIC(i, m2c.mset(k));
            end;
            %now put them into the log model evidence arrays (one model per
            %row, can have more than 2)
            logEvidence_y(k, i) = -0.5*m2c.BIC(i, m2c.mset(k)); %convert to appropriate form for log model evidence
        end;
    end;
end

%display number of cases analysed and whether edge fits are excluded
if edgeflagged_inc==1
    disp('Edge fit cases were retained in this analysis');
else
    disp('Edge fit cases were excluded in this analysis');
end;
disp(['The number of cases analysed = ' num2str(i)]);
    
% display empirical histogram of log-Bayes factors
% -------------------------------------------------------------------------
if nmodels==2 && wantfigs==1
        plotBayesFactor (logEvidence_y(1:2,1:i),[], m2c.names); %can add additional sets to display as second parameter
end;

% perform model selection with the VBA
% =========================================================================
options.verbose = false;
if wantfigs ~=1
    options.DisplayWin = false; %put this line in if you want to suppress
end;
%figure
    
% perform group-BMS on the fits
[py, oy] = VBA_groupBMC (logEvidence_y(:,1:i), options);
%remove the next line if you want to suppress figure
if wantfigs ==1
    set (oy.options.handles.hf, 'name', 'group BMS: BIC fits')
end;

% display statistics
[~, idxWinner] = max(oy.Ef);
fprintf('The best model is the %s model: Ef = %4.3f (pxp = %4.3f)\n',m2c.names{m2c.mset(idxWinner)}, oy.Ef(idxWinner), oy.pxp(idxWinner));
disp(' ');
disp('Dirichlet stats for group model comparison')
disp('First alpha values for the 2 models, model 1 first')
disp(py.a);
%disp('and then exceedance probabilities, model 1 first')
%disp(py.r);
    
end


%% ########################################################################
% display subfunctions
% #########################################################################
function plotBayesFactor (mylogEvidence, ~, model_name) % second param might be mylogEvidence_2)
    %in this fxn I have commented out the plots for a second set of
    %evidences. Can do it better based on number of args supplied to function
    [n1, x1] = VBA_empiricalDensity ((mylogEvidence(1,:) - mylogEvidence(2, :))');
    %[n2, x2] = VBA_empiricalDensity ((mylogEvidence_2(1,:) - mylogEvidence_2(2 ,:))');
    hf = figure ('color' ,'w', 'name', 'Model Comparison: Distribution of log Bayes factors','NumberTitle','off');
    ha = axes ('parent', hf,'nextplot','add');
    plot (ha, x1, n1, 'color', 'r');
    %plot (ha, x2, n2, 'color', 'b');
    legend (ha, {['model 1= ' model_name{1} '; model 2= ' model_name{2}], 'label for plot 2'});
    xlabel (ha, 'log p(y|m1) - log p(y|m2)');
    ylabel (ha, 'proportion of participants');
end
