%%%%% PC-corr algorithm associate to PCA analysis%%%%%

%%% Released under MIT License
%%% Copyright (c) 16 Dec 2017 Sara Ciucci, Yan Ge, Claudio Durán and Carlo Vittorio Cannistraci

% Please cite:
% Enlightening discriminative network functional modules behind Principal Component Analysis separation in differential-omic science studies.
% Sara Ciucci, Yan Ge, Claudio Durán, Alessandra Palladini, Víctor Jiménez Jiménez, Luisa María Martínez Sánchez, 
% Yuting Wang, Susanne Sales, Andrej Shevchenko, Steven W. Poser, Maik Herbig, Oliver Otto, Andreas Androutsellis-Theotokis, 
% Jochen Guck, Mathias J. Gerl and Carlo Vittorio Cannistraci 
% Scientific Reports, 2017

% INPUT
%   x => (Numeric matrix MxN) Dataset with samples on the rows and features on the columns
%   sample_labels => (Cell array of character vectors Mx1) Labels of the samples
%   feat_names => (Cell array of character vectors Nx1) Names of the features
%   sample_names=> (Cell array of character vectors MX1) Names of the samples
%   dis=>'yes' or 'no', depending if you want to display the sample names
%   in the scatterplot. Default: 'no'

% OUTPUT
%   -For a single cut-off: 
%
%   Edges=>(Cell array Ox3) Edges in the PC-corr constructed network: 
%          -first/second columns: two nodes (features) that are connected by an edge 
%          -third column: PC-corr edge weight.
%   Nodes => (Cell array Px3) Nodes in the constructed network: 
%           -first column: nodes (features) in the network
%           -second column: corresponding node colors
%           -third column: correspondind loading value (V)
%
%   -For multiple cut-offs (say C cut-offs):
%
%   Edges=>(Cell array Cx2)
%          -first column: cut-offs
%          -second column: respective edge table in the PC-corr network for each cut-off 
%           (that is a cell array with the same structure as in the single cut-off case)
%   Nodes=>(Cell array Cx2)
%          -first column: cut-offs
%          -second column: respective node table in the PC-corr network  for each cut-off 
%           (that is a cell array with the same structure as in the single cut-off case)


function [Edges,Nodes]=PC_corr_v2(x,sample_labels,feat_names, sample_names,dis)

clc
%% initialisation and default options
if nargin<4, error('Not Enough Input Arguments'); end
if nargin<5, dis='no'; end

if sum(sum(isnan(x)))==1
    error('There is %d NaN value in your data matrix. \nPlease replace it.',sum(sum(isnan(x))));
elseif sum(sum(isnan(x)))>1
    error('There are %d NaN values in your data matrix. \nPlease replace them.',sum(sum(isnan(x))));
end

%%
url='https://www.nature.com/articles/srep43946';
sitename='Scientific Reports 7, Article number: 43946 (2017),doi:10.1038/srep43946';
fprintf('\nWelcome to PC-corr, a simple algorithm that associates to any PCA segregation a discriminative network of features.\n')
fprintf('For more details, see <a href = "%s">%s</a>\n',url,sitename);

labels = sample_labels; 

for i=1:size(labels,1)
    if ~ischar(class(labels(i)))
        labels(i)=cellfun(@num2str, labels(i,:),'UniformOutput',false);
    else
        labels(i)=labels(i);
    end
end

nameLabels = unique(labels); %name of the groups
numbLabels = size(nameLabels,1); %number of groups

for i=1:size(sample_names,1)
    if ~ischar(class(sample_names{i}))
        sample_names(i)=cellfun(@num2str, sample_names(i,:),'UniformOutput',false);
    else
        sample_names(i)=sample_names(i);
    end
end

for i=1:size(feat_names,1)
    if ~ischar(class(feat_names{i}))
        feat_names(i)=cellfun(@num2str, feat_names(i,:),'UniformOutput',false);
    else
        feat_names(i)=feat_names(i);
    end
end

%Remove features with same identical values across all the samples
x1=x;
x1=x1-repmat(mean(x1),size(x1,1),1);
ma=sum(x1==0,1)==size(x,1);
remov_feat=sum(ma); %number of removed fatures
x(:,ma)=[];
feat_names(ma')=[];


%% Type of labels
flag = 0;

while flag == 0
    
    fprintf('\nIs your data represented by\n[r] ranked labels(labels that are organized according to a progressive order. e.g. different stages of a disease, where Stage 1 < Stage 2 < Stage 3)\n[c] class labels (labels that are not necessary organized in a progressive order e.g. Condition A, Condition B, Condition C)? [r/c]:\n\n');
    u_lab = input('-> ','s');
    
    flag = 1;

    if ~strcmp(u_lab,'r') && ~strcmp(u_lab,'c')
        flag = 0;
        fprintf('\nPlease introduce either "r" for ranked labels or "c" for class labels\n')
    end
end

if strcmp(u_lab,'r')

    flag = 0;
    while flag == 0
        fprintf('\nAre the values of your ranked labels\n[d]   discrete (Stage 1 < Stage 2 < Stage 3)\n[con] continuous (different times of development of a cell line)? [d/con]:\n\n');
        u_lab = input('-> ','s');
        
        flag = 1;
        
        if ~strcmp(u_lab,'d') && ~strcmp(u_lab,'con')
            flag = 0;
            fprintf('\nPlease introduce either "d" for discrete labels or "con" for continuous labels\n')
        end
    end
    
end

if strcmp(u_lab,'d')
    % for correlation evaluators, labels are turned into numbers
    for i = 1:numbLabels
        labl_numb(ismember(labels,nameLabels{i}),1) = i;
    end
elseif strcmp(u_lab,'con')
    labl_numb = cell2mat(labels);
end

%% Normalizations of the dataset
flag = 0;

while flag == 0
    
    fprintf('\nThe analysis starts by normalizing or not the dataset.\n\nDo you want to apply:\n[1] no normalization \n[2] a prefered normalization \n[3] automatically all the set of available normalizations? [1/2/3]\n\n');
    u_norm_opt = input('-> ');
    
    flag = 1;

    if (u_norm_opt~=1)&&(u_norm_opt~=2) && (u_norm_opt~=3)
        flag = 0;
        fprintf('\nPlease introduce either 1 for no normalization, 2 for a prefered normalization or 3 for all the set of available normalizations. \n')
    else
        flag = 1;
    end
end

if u_norm_opt == 1
    norm{1}=x; norms{1} = '-'; %No normalization
elseif u_norm_opt == 3
    norm{1}=x./repmat(sum(x,1),size(x,1),1); norms{1} = 'DCS'; %dividing by the column sum
    norm{2}=x./repmat(sum(x,2),1,size(x,2)); norms{2} = 'DRS';%%dividing by the row sum
    norm{3}=log10(1+x); norms{3} = 'LOG';
    norm{4}=zscore(x); norms{4} = 'ZSCORE';
    norm{5}=quantilenorm(x')'; norms{5} = 'QUANTILE T';
    norm{6}=quantilenorm(x); norms{6} = 'QUANTILE';
    norm{7}=zscore(x')'; norms{7} = 'ZSCORE T';
    norm{8}=x+abs(min(min(x)));  norms{8} = 'PLUS(ABS(MIN))';
    
    xCenter=x-repmat(mean(x),size(x,1),1);
    norm{9}=xCenter./repmat(sqrt(std(xCenter)),size(xCenter,1),1); %pareto scaling
    norms{9} = 'PARETO SCALING';
    
    norm{10}=sqrt(x); norms{10} = 'SQRT';
    norm{11}=manorm(x);  norms{11} = 'MANORM';
    norm{12}=x; norms{12} = '-'; %No normalization
elseif u_norm_opt == 2
    
    flag = 0;
    
    while flag == 0
        norms_list_1= {'DCS','DRS','LOG','ZSCORE','QUANTILE T','QUANTILE','ZSCORE T','PLUS(ABS(MIN))','PARETO SCALING','SQRT','MANORM'};
        fprintf('\nInput a preferred normalization from the following list:\n');
        fprintf('DCS, DRS, LOG, ZSCORE, QUANTILE T, QUANTILE, ZSCORE T, PLUS(ABS(MIN)), PARETO SCALING, SQRT, MANORM\n');
        fprintf('(For detailed information on the type of normalization, see the User guide)\n\n');
        fprintf('Example: LOG\n\n');
        u_norm_choice = input('-> ','s');
        flag = 1;
        u_norm_choice1=find(strcmp(norms_list_1(:),u_norm_choice));
        if isempty(u_norm_choice1)==1
            flag = 0;
            fprintf('\nPlease introduce the exact name of the normalization.\n')
        end

    end
    
    norms{1}=u_norm_choice;
    switch u_norm_choice
        case 'DCS'
           norm{1}= x./repmat(sum(x,1),size(x,1),1); %dividing by the column sum
        case 'DRS'
            norm{1} = x./repmat(sum(x,2),1,size(x,2)); %dividing by the row sum
        case 'LOG'
            norm{1} = log10(1+x);
        case 'ZSCORE'
            norm{1} = zscore(x);
        case 'QUANTILE T'
            norm{1} = quantilenorm(x')';
        case 'QUANTILE'
            norm{1} = quantilenorm(x);
        case 'ZSCORE T'
            norm{1} = zscore(x')';
        case 'PLUS(ABS(MIN))'
            norm{1} = x+abs(min(min(x)));
        case 'PARETO SCALING'
            xCenter=x-repmat(mean(x),size(x,1),1);
            norm{1} = xCenter./repmat(sqrt(std(xCenter)),size(xCenter,1),1); %pareto scaling
        case 'SQRT'
            norm{1} = sqrt(x);
        case 'MANORM'
            norm{1} = manorm(x);
    end
    
    inf_norm=sum(sum(isinf(norm{1})))~=0;
    nan_norm=sum(sum(isnan(norm{1})))~=0;
    nonreal_norm=~isreal(norm{1});
    probl_norm=inf_norm+nan_norm+nonreal_norm;
        
end
    
    
for i=1:length(norm)
    inf_norm(i)=sum(sum(isinf(norm{i})))~=0;
    nan_norm(i)=sum(sum(isnan(norm{i})))~=0;
    nonreal_norm(i)=~isreal(norm{i});
end

probl_norm=inf_norm+nan_norm+nonreal_norm;
norm=norm(~probl_norm);
norms=norms(~probl_norm);
norms_list=norms;



% number of elements in each group 
for i=1:length(nameLabels)
    number_el_group(i)=sum(strcmp(labels,nameLabels{i}));
end

%Check if some groups have the same amount of elements
valCount = hist(number_el_group, unique(number_el_group));

flag = 0;
if (sum(valCount>1)==0) && (numbLabels>2)
    while flag == 0
        fprintf('\nFor the calculation of Area Under the ROC-Curve (AUC) and Area Under the Precision-Recall curve (AUPR), you need to provide a positive label for each pairwise group comparison.\n');
        fprintf('\nDo you want to calculate the AUC and AUPR values considering:\n[s] as positive label, the label of the smallest sample group in each pairwise group comparison\n[l] as positive label, the label of the largest sample group in each pairwise group comparison\n[r] a ranked list of possible positive labels \n\n');
            
        u_aupr = input('-> ','s');
        
        flag = 1;
        
        if ~strcmp(u_aupr,'s') && ~strcmp(u_aupr,'l') &&  ~strcmp(u_aupr,'r')
            flag = 0;

            fprintf('\nPlease introduce either "s" for the label of the smallest sample group as positive label, \n "l" for the labels of the largest sample group as positive label, \n "r" for a ranked list of positive labels.\n')
        end
        
    end
else
    flag = 1;
    u_aupr='r';
end


if strcmp(u_aupr,'r')

    flag = 0;
    while flag == 0
        if numbLabels==2 
            fprintf('\nInput the positive label for the calculation of AUC and AUPR values: \n');
        else
            fprintf('\nInput the ranked list of possible positive labels for the calculation of AUC and AUPR values: \n');
        end
        randperm_nameLabels=nameLabels(randperm(length(nameLabels)));
        if  length(nameLabels)>2
            fprintf('Example: {')
            for i=1:length(randperm_nameLabels)-2
                fprintf( '''%s'', ', randperm_nameLabels{i})
            end
            fprintf( '''%s''',randperm_nameLabels{length(randperm_nameLabels)-1})
            fprintf('}\n\n')
        elseif length(nameLabels)==2
            fprintf('Example: ''%s''\n\n',randperm_nameLabels{1})
        end
        u_aupr_r = input('-> ');
        
        flag = 1;
        
        if sum(ismember(u_aupr_r,nameLabels))==0
            flag = 0;
            if numbLabels > 2
                fprintf('\nPlease introduce a correct ranked list of positive labels.\n');
            elseif numbLabels == 2
                fprintf('\nPlease introduce a correct positive label.\n');
            end
        end
    end
    
end

if strcmp(u_aupr,'s') || strcmp(u_aupr,'l')
    u_aupr_r={};
end

for i = 1:length(norm)
    
    %%% non-centred PCA
    [~,snc,pc_nc{i}]=svd(norm{i},'econ');
    ncPCA{i}=norm{i}*pc_nc{i};
    
    %%% centred PCA
    normCenter=norm{i}-repmat(mean(norm{i}),size(norm{i},1),1);
    [~,sc,pc_c{i}]=svd(normCenter,'econ');
    cPCA{i}=normCenter*pc_c{i};
    
    latent_nc = (diag(snc).^2)/size(x,1); %component variance 
    explained_nc{i} = 100*latent_nc/sum(latent_nc);%explained variance (%)
    
    latent_c = (diag(sc).^2)/(size(x,1)-1); %component variance
    explained_c{i} = 100*latent_c/sum(latent_c);%explained variance (%)
    
    
    
    for k=1:size(ncPCA{i},2) %dimension
        
        if strcmp(u_lab,'c') || strcmp(u_lab,'d') 
            
            n = 1;
            m = 2;
            for j=1:nchoosek(numbLabels,2) %two-group comparison
                
                % Compute p-value of Mann-Whitney test
                mw_ncPCA{i,j,k} = ranksum(ncPCA{i}(ismember(labels,nameLabels{n}),k),ncPCA{i}(ismember(labels,nameLabels{m}),k));
                mw_cPCA{i,j,k} = ranksum(cPCA{i}(ismember(labels,nameLabels{n}),k),cPCA{i}(ismember(labels,nameLabels{m}),k));
                
                
                samp_lab = [labels(ismember(labels,nameLabels{n}));labels(ismember(labels,nameLabels{m}))];
                scores_nc = [ncPCA{i}(ismember(labels,nameLabels{n}),k);ncPCA{i}(ismember(labels,nameLabels{m}),k)];
                scores_c = [cPCA{i}(ismember(labels,nameLabels{n}),k);cPCA{i}(ismember(labels,nameLabels{m}),k)];
                

                possClass=positive_label_opt(u_aupr, nameLabels{n},nameLabels{m},labels,u_aupr_r);

                % Compute AUC
                [~,~,~,AUC_nc{i,j,k}] = perfcurve(samp_lab,scores_nc, possClass);
                [~,~,~,AUC_c{i,j,k}] = perfcurve(samp_lab,scores_c, possClass);
                
                if AUC_nc{i,j,k}<0.5
                    flag = 1;
                    AUC_nc{i,j,k} = 1-AUC_nc{i,j,k};
                    % Compute AUPR
                    flip_scores_nc=2*mean(scores_nc)-scores_nc;
                    AUPR_nc{i,j,k} = aupr_evaluation(samp_lab,flip_scores_nc,possClass);
                else
                    % Compute AUPR
                    AUPR_nc{i,j,k} = aupr_evaluation(samp_lab,scores_nc,possClass);
                end 
                
                if AUC_c{i,j,k}<0.5
                    flag = 1;
                    AUC_c{i,j,k} = 1-AUC_c{i,j,k};
                    % Compute AUPR
                    flip_scores_c=2*mean(scores_c)-scores_c;
                    AUPR_c{i,j,k} = aupr_evaluation(samp_lab,flip_scores_c,possClass);
                else
                    % Compute AUPR
                    AUPR_c{i,j,k} = aupr_evaluation(samp_lab,scores_c,possClass);
                end 
                

                
                m = m + 1;
                if(m > numbLabels)
                    n = n + 1;
                    m = n + 1;
                end
            end
            
            if strcmp(u_lab,'d')
                rank_pears_corr_ncPCA{i,k} = corr(ncPCA{i}(:,k),labl_numb,'Type','Pearson');
                rank_pears_corr_cPCA{i,k} = corr(cPCA{i}(:,k),labl_numb,'Type','Pearson');
                rank_spear_corr_ncPCA{i,k} = corr(ncPCA{i}(:,k),labl_numb,'Type','Spearman');
                rank_spear_corr_cPCA{i,k} = corr(cPCA{i}(:,k),labl_numb,'Type','Spearman');
            end
            
        else
            
            rank_pears_corr_ncPCA{i,k} = corr(ncPCA{i}(:,k),labl_numb,'Type','Pearson');
            rank_pears_corr_cPCA{i,k} = corr(cPCA{i}(:,k),labl_numb,'Type','Pearson');
            rank_spear_corr_ncPCA{i,k} = corr(ncPCA{i}(:,k),labl_numb,'Type','Spearman');
            rank_spear_corr_cPCA{i,k} = corr(cPCA{i}(:,k),labl_numb,'Type','Spearman');
            
        end
    end
end


%% Constructing the table of results
norms = repmat(norms',[size(ncPCA{1},2)*2 1]);
centred = repmat({'yes'},[size(ncPCA{1},2)*length(norm) 1]);
non_centred = repmat({'no'},[size(ncPCA{1},2)*length(norm) 1]);
centr = {non_centred{:} centred{:}};
dim = (1:size(ncPCA{1},2))';
dim = repmat(dim,[length(norm),1]);
dim = num2cell(sort(dim));
dim = {dim{:} dim{:}};
variance_c = [explained_c{:}];
variance_c = variance_c';
variance_c = num2cell(variance_c(:));
variance_nc = [explained_nc{:}];
variance_nc = variance_nc';
variance_nc = num2cell(variance_nc(:));
variance = vertcat(variance_nc,variance_c);



if numbLabels == 2 %two groups

    if strcmp(u_lab,'c')
    
        header = {'P-value','AUC','AUPR','Norm','Centering','Dim','expl Var'};
        header_xls = {'P-value','AUC','AUPR','Norm','Centering','Dim','expl Var',...
            'Trustworthiness(p-value)','Trustworthiness(AUC)','Trustworthiness(AUPR)'};
        pvals = {mw_ncPCA{:} mw_cPCA{:}};
        AUCs = {AUC_nc{:} AUC_c{:}};
        AUPRs = {AUPR_nc{:} AUPR_c{:}};
        results = horzcat(pvals',AUCs',AUPRs',norms,centr',dim',variance);
        results_xls=results;       
        
    elseif strcmp(u_lab,'d')
        
        header = {'P-value','AUC','AUPR','pears','spear','Norm','Centering','Dim','expl Var'};
        header_xls = {'P-value','AUC','AUPR','pears','spear','Norm','Centering','Dim','expl Var',...
             'Trustworthiness(p-value)','Trustworthiness(AUC)','Trustworthiness(AUPR)'};

        pvals = {mw_ncPCA{:} mw_cPCA{:}};
        AUCs = {AUC_nc{:} AUC_c{:}};
        AUPRs = {AUPR_nc{:} AUPR_c{:}};
        pears_corr = {rank_pears_corr_ncPCA{:} rank_pears_corr_cPCA{:}};
        spear_corr = {rank_spear_corr_ncPCA{:} rank_spear_corr_cPCA{:}};
        results = horzcat(pvals',AUCs',AUPRs',pears_corr',spear_corr',norms,centr',dim',variance);
        results_xls=results;
        
    elseif strcmp(u_lab,'con')
        
        header = {'pears','spear','Norm','Centering','Dim','expl Var'};
        pears_corr = {rank_pears_corr_ncPCA{:} rank_pears_corr_cPCA{:}};
        spear_corr = {rank_spear_corr_ncPCA{:} rank_spear_corr_cPCA{:}};
        results = horzcat(pears_corr',spear_corr',norms,centr',dim',variance);
        results_xls=results;
    end
    
    fprintf('\nThe best discrimination in a PCA result (combination of normalization and centering) and along one dimension (principal component, PCn) can be \nassessed by different evaluators.\n');
    
    flag = 0;
    
    while flag == 0
        
        if strcmp(u_lab,'c') || strcmp(u_lab,'d')
            
            if strcmp(u_lab,'c')
            
                fprintf('\nWould you like to rank the PCA results by \n[p]    P-value\n[auc]  AUC\n[aupr] AUPR? [p/auc/aupr]:\n\n');
                u_rank = input('-> ','s');           
            else
                fprintf('\nWould you like to rank the PCA results by \n[p]    P-value\n[auc]  AUC\n[aupr] AUPR\n[pc]   Pearson correlation\n[sc]   Spearman correlation? [p/auc/aupr/pc/sc]:\n\n');
                u_rank = input('-> ','s');
            end
            
            flag = 1;
            
            switch u_rank
                case 'p'
                    [~,idx] = sort([results{:,1}]);
                    results = results(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results([results{:,1}]<0.05,:))
                        fprintf('\nThe table below shows the best results (p-value<0.25), automatically ranked with respect to P-value from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');  
                        display(vertcat(header,results([results{:,1}]<0.25,:)))
                    else
                        fprintf('\nThe table below shows the best results (p-value<0.05), automatically ranked with respect to P-value from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                        display(vertcat(header,results([results{:,1}]<0.05,:)))
                    end
                    fprintf('\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results,\nautomatically ranked with respect to P-value.\n\n')
                    %Creation of the excel file with all the results
                    filename = 'results.xlsx';
                    results_xls = results_xls(idx,:);
                    if exist(filename,'file')~=0
                        delete(filename)
                        warning('off','MATLAB:xlswrite:AddSheet')
                        results_table = cell2table(results_xls, 'VariableNames', header);
                        writetable(results_table, filename, 'Sheet', 'PCA results');
                        % xlswrite(filename,vertcat(header,results_xls),'PCA results')
                        % RemoveSheet123(filename)
                    else
                        warning('off','MATLAB:xlswrite:AddSheet')                      
                        results_table = cell2table(results_xls, 'VariableNames', header);
                        writetable(results_table, filename, 'Sheet', 'PCA results');
                        % xlswrite(filename,vertcat(header,results_xls),'PCA results')
                        % RemoveSheet123(filename)
                    end
                    
                    
                case 'auc'
                    [~,idx] = sort([results{:,2}],'descend');
                    results = results(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results([results{:,2}]>=0.7,:))
                        fprintf('\nThe table below shows the best results (AUC>=0.6), automatically ranked with respect to AUC from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                        display(vertcat(header,results([results{:,2}]>=0.6,:)))
                    else
                        fprintf('\nThe table below shows the best results (AUC>=0.7), automatically ranked with respect to AUC from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                        display(vertcat(header,results([results{:,2}]>=0.7,:)))
                    end
                    fprintf('\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results,\nautomatically ranked with respect to AUC.\n\n')
                    %Creation of the excel file with all the results
                    filename = 'results.xlsx';
                    results_xls = results_xls(idx,:);
                    if exist(filename,'file')~=0
                        delete(filename)
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header,results_xls),'PCA results')
                        RemoveSheet123(filename)
                    else
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header,results_xls),'PCA results')
                        RemoveSheet123(filename)
                    end
                    
                case 'aupr'
                    [~,idx] = sort([results{:,3}],'descend');
                    results = results(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results([results{:,3}]>=0.7,:))
                        fprintf('\nThe table below shows the best results (AUPR>=0.6), automatically ranked with respect to AUPR from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                        display(vertcat(header,results([results{:,3}]>=0.6,:)))
                    else
                        fprintf('\nThe table below shows the best results (AUPR>=0.7), automatically ranked with respect to AUPR from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                        display(vertcat(header,results([results{:,3}]>=0.7,:)))
                    end
                    fprintf('\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results,\nautomatically ranked with respect to AUPR.\n\n')
                    %Creation of the excel file with all the results
                    filename = 'results.xlsx';
                    results_xls = results_xls(idx,:);
                    if exist(filename,'file')~=0
                        delete(filename)
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header,results_xls),'PCA results')
                        RemoveSheet123(filename)
                    else
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header,results_xls),'PCA results')
                        RemoveSheet123(filename)
                    end
                    
                case 'pc'
                    if strcmp(u_lab,'c')
                        flag = 0;
                        fprintf('\nPlease introduce either "p", "auc" or "aupr"\n')
                    else
                        [~,idx] = sort(abs([results{:,4}]),'descend');
                        results = results(idx,:);
                        %Only some results are shown on the screen
                        if isempty(results(abs([results{:,4}])>=0.6,:))
                            fprintf('\nThe table below shows the best results (|pears|>=0.5), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                            display(vertcat(header,results(abs([results{:,4}])>=0.5,:)))
                        else
                            fprintf('\nThe table below shows the best results (|pears|>=0.6), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                            display(vertcat(header,results(abs([results{:,4}])>=0.6,:)))
                        end
                        fprintf('\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results,\nautomatically ranked with respect to Pearson correlation.\n\n')
                        %Creation of the excel file with all the results
                        filename = 'results.xlsx';
                        results_xls = results_xls(idx,:);
                        if exist(filename,'file')~=0
                            delete(filename)
                            warning('off','MATLAB:xlswrite:AddSheet')
                            xlswrite(filename,vertcat(header,results_xls),'PCA results')
                            RemoveSheet123(filename)
                        else
                            warning('off','MATLAB:xlswrite:AddSheet')
                            xlswrite(filename,vertcat(header,results_xls),'PCA results')
                            RemoveSheet123(filename)
                        end
                        
                    end
                case 'sc'
                    if strcmp(u_lab,'c')
                        flag = 0;
                        fprintf('Please introduce either "p", "auc" or "aupr"\n')
                    else
                        [~,idx] = sort(abs([results{:,5}]),'descend');
                        results = results(idx,:);
                        %Only some results are shown on the screen
                        if isempty(results(abs([results{:,5}])>=0.6,:))
                            fprintf('\nThe table below shows the best results (|spear|>=0.5), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                            display(vertcat(header,results(abs([results{:,5}])>=0.5,:)))

                        else
                            fprintf('\nThe table below shows the best results (|spear|>=0.6), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                            display(vertcat(header,results(abs([results{:,5}])>=0.6,:)))
                        end
                        fprintf('\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results,\nautomatically ranked with respect to Spearman correlation.\n\n')
                        %Creation of the excel file with all the results
                        filename = 'results.xlsx';
                        results_xls = results_xls(idx,:);
                        if exist(filename,'file')~=0
                            delete(filename)
                            warning('off','MATLAB:xlswrite:AddSheet')
                            xlswrite(filename,vertcat(header,results_xls),'PCA results')
                            RemoveSheet123(filename)
                        else
                            warning('off','MATLAB:xlswrite:AddSheet')
                            xlswrite(filename,vertcat(header,results_xls),'PCA results')
                            RemoveSheet123(filename)
                        end
                        
                    end
                otherwise
                    if strcmp(u_lab,'c')
                        flag = 0;
                        fprintf('\nPlease introduce either "p", "auc" or "aupr"\n')
                    else
                        flag = 0;
                        fprintf('\nPlease introduce either "p", "auc", "aupr", "pc" or "sc"\n')
                    end
            end
            
        else
            
            fprintf('\nWould you like to rank the PCA results by \n[pc] Pearson correlation\n[sc] Spearman correlation? [pc/sc]:\n\n');
            u_rank = input('-> ','s');
            
            flag = 1;
            
            switch u_rank
                case 'pc'
                    [~,idx] = sort(abs([results{:,1}]),'descend');
                    results = results(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results(abs([results{:,1}])>=0.6,:))
                        fprintf('\nThe table below shows the best results (|pears|>=0.5), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                        display(vertcat(header,results(abs([results{:,1}])>=0.5,:)))
                    else
                        fprintf('\nThe table below shows the best results (|pears|>=0.6), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                        display(vertcat(header,results(abs([results{:,1}])>=0.6,:)))
                    end
                    fprintf('\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results,\nautomatically ranked with respect to Pearson correlation.\n\n')
                    %Creation of the excel file with all the results
                    filename = 'results.xlsx';
                    if exist(filename,'file')~=0
                        delete(filename)
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header,results),'PCA results')
                        RemoveSheet123(filename)
                    else
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header,results),'PCA results')
                        RemoveSheet123(filename)
                    end
                    
                case 'sc'
                    [~,idx] = sort(abs([results{:,2}]),'descend');
                    results = results(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results(abs([results{:,2}])>=0.6,:))
                        fprintf('\nThe table below shows the best results (|spear|>=0.5), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                        display(vertcat(header,results(abs([results{:,2}])>=0.5,:)))
                    else
                        fprintf('\nThe table below shows the best results (|spear|>=0.6), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                        display(vertcat(header,results(abs([results{:,2}])>=0.6,:)))
                    end
                    fprintf('\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results,\nautomatically ranked with respect to Spearman correlation.\n\n')
                    %Creation of the excel file with all the results
                    filename = 'results.xlsx';
                    if exist(filename,'file')~=0
                        delete(filename)
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header,results),'PCA results')
                        RemoveSheet123(filename)
                    else
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header,results),'PCA results')
                        RemoveSheet123(filename)
                    end
                    
                otherwise
                    flag = 0;
                    fprintf('\nPlease introduce either "pc" or "sc"\n');
            end
            
        end
    end
    
else %more than two groups
    
    if strcmp(u_lab,'c') || strcmp(u_lab,'d')
        
        pval1 = permute(mw_ncPCA,[1 3 2]);
        pval1 = reshape(pval1,length(norm)*size(mw_cPCA,3),nchoosek(numbLabels,2));
        pval2 = permute(mw_cPCA,[1 3 2]);
        pval2 = reshape(pval2,length(norm)*size(mw_cPCA,3),nchoosek(numbLabels,2));
        pvals = vertcat(pval1,pval2);
        avg = num2cell(mean(cell2mat(pvals),2));
        
        AUC1 = permute(AUC_nc, [1 3 2]);
        AUC1 = reshape(AUC1,length(norm)*size(mw_cPCA,3),nchoosek(numbLabels,2));
        AUC2 = permute(AUC_c, [1 3 2]);
        AUC2 = reshape(AUC2,length(norm)*size(mw_cPCA,3),nchoosek(numbLabels,2));
        AUCs = vertcat(AUC1,AUC2);
        avg_AUC = num2cell(mean(cell2mat(AUCs),2));
        
        AUPR1 = permute(AUPR_nc, [1 3 2]);
        AUPR1 = reshape(AUPR1,length(norm)*size(mw_cPCA,3),nchoosek(numbLabels,2));
        AUPR2 = permute(AUPR_c, [1 3 2]);
        AUPR2 = reshape(AUPR2,length(norm)*size(mw_cPCA,3),nchoosek(numbLabels,2));
        AUPRs = vertcat(AUPR1,AUPR2);
        avg_AUPR = num2cell(mean(cell2mat(AUPRs),2));
        
        o = 1;
        p = 2;
        for i = 1:nchoosek(numbLabels,2)
            group_head{i} = [nameLabels{o},' vs ',nameLabels{p}];
        
            p = p + 1;
            if(p > numbLabels)
                o = o + 1;
                p = o + 1;
            end
        end
        
        
        if strcmp(u_lab,'c')
            
            header = {'Avg P-val','Avg AUC','Avg AUPR','Norm','Centering','Dimension','expl Var'};
            results = horzcat(avg,avg_AUC,avg_AUPR,norms,centr',dim',variance);
            
        elseif strcmp(u_lab,'d')
            pears_corr = {rank_pears_corr_ncPCA{:} rank_pears_corr_cPCA{:}};
            spear_corr = {rank_spear_corr_ncPCA{:} rank_spear_corr_cPCA{:}};
            header = {'Avg Pval','Avg AUC','Avg AUPR','pears','spear','Norm','Centering','Dim','expl Var'};
            results = horzcat(avg,avg_AUC,avg_AUPR,pears_corr',spear_corr',norms,centr',dim',variance);
        end
        

        
        
        if strcmp(u_lab,'c')
            header_xls1 = {'Avg P-val',group_head{:},'Avg AUC',group_head{:},'Avg AUPR',group_head{:},'Norm','Centering','Dimension','explained Variance'};
            header_xls={'Avg P-val',group_head{:},'Avg AUC',group_head{:},'Avg AUPR',group_head{:},'Norm','Centering','Dimension','explained Variance',...
                'Trustworthiness(p-value)','Trustworthiness(AUC)','Trustworthiness(AUPR)'};
            results_xls = horzcat(avg,pvals,avg_AUC,AUCs,avg_AUPR,AUPRs,norms,centr',dim',variance);

        end
        if strcmp(u_lab,'d')
            header_xls1 = {'Avg P-val',group_head{:},'Avg AUC',group_head{:},'Avg AUPR',group_head{:},'pearson-correlation','spearman-correlation','Norm','Centering','Dimension','explained Variance'};
            header_xls = {'Avg P-val',group_head{:},'Avg AUC',group_head{:},'Avg AUPR',group_head{:},'pearson-correlation','spearman-correlation','Norm','Centering','Dimension','explained Variance',...
                          'Trustworthiness(p-value)','Trustworthiness(AUC)','Trustworthiness(AUPR)'};
            results_xls = horzcat(avg,pvals,avg_AUC,AUCs,avg_AUPR,AUPRs,pears_corr',spear_corr',norms,centr',dim',variance);
        end
        
        
    else
        header = {'pears','spear','Norm','Centering','Dimension','expl Var'};
        pears_corr = {rank_pears_corr_ncPCA{:} rank_pears_corr_cPCA{:}};
        spear_corr = {rank_spear_corr_ncPCA{:} rank_spear_corr_cPCA{:}};
        
        header_xls = {'pearson-correlation','spearman-correlation','Norm','Centering','Dimension','explained Variance'};
        results_xls = horzcat(pears_corr',spear_corr',norms,centr',dim',variance);
        results=results_xls;
    end
    
    
    flag = 0;
    while flag == 0
        
        if strcmp(u_lab,'c') || strcmp(u_lab,'d')
            
            if strcmp(u_lab,'c')
                
                fprintf('\nWould you like to rank the PCA results by \n[p]    P-value\n[auc]  AUC\n[aupr] AUPR? [p/auc/aupr]:\n\n');
                u_rank = input('-> ','s');
                
            else
                fprintf('\nWould you like to rank the PCA results by \n[p]    P-value\n[auc]  AUC\n[aupr] AUPR\n[pc]   Pearson correlation\n[sc]   Spearman correlation? [p/auc/aupr/pc/sc]:\n\n');
                u_rank = input('-> ','s');
            end
            
            flag = 1;
            
            switch u_rank
                case 'p'
                    [~,idx] = sort([results{:,1}]);
                    results = results(idx,:);
                    results_xls=results_xls(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results([results{:,1}]<0.05,:))
                        fprintf('\nThe table below show best results (p-value<0.25), automatically ranked with respect to P-value  from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.');
                        fprintf('\nSince there are more than two groups, only the average value (and not the values in each pairwise comparison) for each evaluator is present.\n\n')
                        display(vertcat(header,results([results{:,1}]<0.25,:)))
                    else
                        fprintf('\nThe table below show best results (p-value<0.05), automatically ranked with respect to P-value  from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.');
                        fprintf('\nSince there are more than two groups, only the average value (and not the values in each pairwise comparison) for each evaluator is present.\n\n')
                        display(vertcat(header,results([results{:,1}]<0.05,:)))
                    end
                    fprintf('\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results,\nautomatically ranked with respect to P-value.\n\n')
                    %Creation of the excel file with all the results
                    filename = 'results.xlsx';
                    if exist(filename,'file')~=0
                        delete(filename)
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header_xls1,results_xls),'PCA results')
                        RemoveSheet123(filename)
                    else
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header_xls1,results_xls),'PCA results')
                        RemoveSheet123(filename)
                    end
                    
                case 'auc'
                    [~,idx] = sort([results{:,2}],'descend');
                    results = results(idx,:);
                    results_xls=results_xls(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results([results{:,2}]>=0.7,:))
                        fprintf('\nThe table below show best results (AUC>=0.6), automatically ranked with respect to AUC from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.');
                        fprintf('\nSince there are more than two groups, only the average value (and not the values in each pairwise comparison) for each evaluator is present.\n\n')
                        display(vertcat(header,results([results{:,2}]>=0.6,:)))
                    else
                        fprintf('\nThe table below show best results (AUC>=0.7), automatically ranked with respect to AUC from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.');
                        fprintf('\nSince there are more than two groups, only the average value (and not the values in each pairwise comparison) for each evaluator is present.\n\n')
                        display(vertcat(header,results([results{:,2}]>=0.7,:)))
                    end
                    fprintf('\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results,\nautomatically ranked with respect to AUC.\n\n')
                    %Creation of the excel file with all the results
                    filename = 'results.xlsx';
                    if exist(filename,'file')~=0
                        delete(filename)
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header_xls1,results_xls),'PCA results')
                        RemoveSheet123(filename)
                    else
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header_xls1,results_xls),'PCA results')
                        RemoveSheet123(filename)
                    end
                    
                case 'aupr'
                    [~,idx] = sort([results{:,3}],'descend');
                    results = results(idx,:);
                    results_xls=results_xls(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results([results{:,3}]>=0.7,:))
                        fprintf('\nThe table below show best results (AUPR>=0.6), automatically ranked with respect to AUPR from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.');
                        fprintf('\nSince there are more than two groups, only the average value (and not the values in each pairwise comparison) for each evaluator is present.\n\n')
                        display(vertcat(header,results([results{:,3}]>=0.6,:)))
                    else
                        fprintf('\nThe table below show best results (AUPR>=0.7), automatically ranked with respect to AUPR from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.');
                        fprintf('\nSince there are more than two groups, only the average value (and not the values in each pairwise comparison) for each evaluator is present.\n\n')
                        display(vertcat(header,results([results{:,3}]>=0.7,:)))
                    end
                    fprintf('\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results,\nautomatically ranked with respect to AUPR.\n\n')
                    %Creation of the excel file with all the results
                    filename = 'results.xlsx';
                    if exist(filename,'file')~=0
                        delete(filename)
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header_xls1,results_xls),'PCA results')
                        RemoveSheet123(filename)
                    else
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header_xls1,results_xls),'PCA results')
                        RemoveSheet123(filename)
                    end
                    
                    
                case 'pc'
                    if strcmp(u_lab,'c')
                        flag = 0;
                        fprintf('\nPlease introduce either "p", "auc" or "aupr"\n')
                    else
                        [~,idx] = sort(abs([results{:,4}]),'descend');
                        results = results(idx,:);
                        results_xls=results_xls(idx,:);
                        %Only some results are shown on the screen
                        if isempty(results(abs([results{:,4}])>=0.6,:))
                            fprintf('\nThe table below show best results (|pears|>=0.5), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                            display(vertcat(header,results(abs([results{:,4}])>=0.5,:)))
                        else
                            fprintf('\nThe table below show best results (|pears|>=0.6), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                            display(vertcat(header,results(abs([results{:,4}])>=0.6,:)))
                        end
                        fprintf('\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results,\nautomatically ranked with respect to Pearson correlation.\n\n')
                        %Creation of the excel file with all the results
                        filename = 'results.xlsx';
                        if exist(filename,'file')~=0
                            delete(filename)
                            warning('off','MATLAB:xlswrite:AddSheet')
                            xlswrite(filename,vertcat(header_xls1,results_xls),'PCA results')
                            RemoveSheet123(filename)
                        else
                            warning('off','MATLAB:xlswrite:AddSheet')
                            xlswrite(filename,vertcat(header_xls1,results_xls),'PCA results')
                            RemoveSheet123(filename)
                        end
                        
                    end
                case 'sc'
                    if strcmp(u_lab,'c')
                        flag = 0;
                        fprintf('\nPlease introduce either "p", "auc" or "aupr"\n')
                    else
                        [~,idx] = sort(abs([results{:,5}]),'descend');
                        results = results(idx,:);
                        results_xls=results_xls(idx,:);
                        %Only some results are shown on the screen
                        if isempty(results(abs([results{:,5}])>=0.6,:))
                            fprintf('\nThe table below show best results (|spear|>=0.5), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                            display(vertcat(header,results(abs([results{:,5}])>=0.5,:)))
                        else
                            fprintf('\nThe table below show best results (|spear|>=0.6), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                            display(vertcat(header,results(abs([results{:,5}])>=0.6,:)))
                        end
                        fprintf('\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results,\nautomatically ranked with respect to Spearman correlation.\n\n')
                        %Creation of the excel file with all the results
                        filename = 'results.xlsx';
                        if exist(filename,'file')~=0
                            delete(filename)
                            warning('off','MATLAB:xlswrite:AddSheet')
                            xlswrite(filename,vertcat(header_xls1,results_xls),'PCA results')
                            RemoveSheet123(filename)
                        else
                            warning('off','MATLAB:xlswrite:AddSheet')
                            xlswrite(filename,vertcat(header_xls1,results_xls),'PCA results')
                            RemoveSheet123(filename)
                        end
                        
                    end
                otherwise
                    if strcmp(u_lab,'c')
                        flag = 0;
                        fprintf('\nPlease introduce either "p", "auc" or "aupr"\n')
                    else
                        flag = 0;
                        fprintf('\nPlease introduce either "p", "auc", "aupr", "pc" or "sc"\n')
                    end
            end
        else
            
            fprintf('\nWould you like to rank the PCA results by \n[pc] Pearson correlation \n[sc] Spearman correlation? [pc/sc]\n\n');
            u_rank = input('-> ','s');
            
            flag = 1;
            
            switch u_rank
                case 'pc'
                    [~,idx] = sort(abs([results{:,1}]),'descend');
                    results = results(idx,:);
                    results_xls=results_xls(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results(abs([results{:,1}])>=0.6,:))
                        fprintf('\nThe table below show best results (|pears|>=0.5), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                        display(vertcat(header,results(abs([results{:,1}])>=0.5,:)))
                    else
                        fprintf('\nThe table below show best results (|pears|>=0.6), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.\n\n');
                        display(vertcat(header,results(abs([results{:,1}])>=0.6,:)))
                    end
                    fprintf('\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results,\nautomatically ranked with respect to Pearson correlation.\n\n')
                    %Creation of the excel file with all the results
                    filename = 'results.xlsx';
                    if exist(filename,'file')~=0
                        delete(filename)
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header_xls,results_xls),'PCA results')
                        RemoveSheet123(filename)
                    else
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header_xls,results_xls),'PCA results')
                        RemoveSheet123(filename)
                    end
                    
                case 'sc'
                    [~,idx] = sort(abs([results{:,2}]),'descend');
                    results = results(idx,:);
                    results_xls=results_xls(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results(abs([results{:,2}])>=0.6,:))
                        fprintf('\nThe table below show best results (|spear|>=0.5), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.');
                        display(vertcat(header,results(abs([results{:,2}])>=0.5,:)))
                    else
                        fprintf('\nThe table below show best results (|spear|>=0.6), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, \nin order to assist the user in the creation of the network from the most discriminative PCn.');
                        display(vertcat(header,results(abs([results{:,2}])>=0.6,:)))
                    end
                    fprintf('\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results,\nautomatically ranked with respect to Spearman correlation.\n\n')
                    %Creation of the excel file with all the results
                    filename = 'results.xlsx';
                    if exist(filename,'file')~=0
                        delete(filename)
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header_xls,results_xls),'PCA results')
                        RemoveSheet123(filename)
                    else
                        warning('off','MATLAB:xlswrite:AddSheet')
                        xlswrite(filename,vertcat(header_xls,results_xls),'PCA results')
                        RemoveSheet123(filename)
                    end
                    
                otherwise
                    flag = 0;
                    fprintf('Please introduce either "pc" or "sc"\n')
            end
        end
    end
end


fprintf('\nAfter seeing the ranked results, you can choose any combination of normalization, centring,dimension (PCn) and cut-off \nfor the network to obtain the PC-corr network  according to your interest (or need).\n')
%% User interaction
if (u_norm_opt == 1) || (u_norm_opt ==2)
    u_norm=1;
    u_norm_n=norms_list{1};
else
    flag = 0;
    while flag == 0
        fprintf('\n\nSelect the normalization (for detailed information see the User guide):\n\n');
        poss_norm_case=norms_list(randperm(length(norms_list)));
        while strcmp(poss_norm_case{1},'-')
            poss_norm_case=norms_list(randperm(length(norms_list)));
        end
        fprintf('Examples: %s or - \n',poss_norm_case{1});
        fprintf('where - stands for no normalization.\n\n');
        u_norm_n = input('-> ','s');
        flag = 1;
        u_norm=find(strcmp(norms_list(:),u_norm_n));
        if isempty(u_norm)==1
            flag = 0;
            fprintf('\nPlease introduce the exact name of the normalization.\n')
        end
    end
end
flag = 0;
while flag == 0
    fprintf('\nCentering version?[y/n]\n[y] yes, centred PCA\n[n] no, non-centred PCA\n\n');
    u_cent = input('-> ','s');
    if strcmp(u_cent,'n') || strcmp(u_cent,'y')
        flag = 1;
    else
        fprintf('\nPlease introduce just "y" or "n"\n');
    end
end

flag = 0;
while flag == 0
    fprintf('\nSelect the dimension for generating the PC-corr network:\n\n');
    u_dim = input('-> ');
    if u_dim > 0 && u_dim <= size(ncPCA{1},2)
        flag = 1;
    else
        fprintf('\nPlease introduce an existing dimension.\n');
    end
end



if strcmp(u_cent,'y')
    PCA = cPCA;
    if strcmp(u_lab,'c') 
        switch u_rank
            case 'p'
                mw_PCA = mw_cPCA;
                mw_PCA1=AUC_c;
                mw_PCA2=AUPR_c;
                evaluat='P-value';
            case 'auc'
                mw_PCA = AUC_c;
                mw_PCA1=mw_cPCA;
                mw_PCA2=AUPR_c;
                evaluat='AUC';
            case 'aupr'
                mw_PCA = AUPR_c;
                mw_PCA1= mw_cPCA;
                mw_PCA2=AUC_c;
                evaluat='AUPR';
        end
    elseif strcmp(u_lab,'d') 
        switch u_rank
            case 'p'
                mw_PCA = mw_cPCA;
                mw_PCA1=AUC_c;
                mw_PCA2=AUPR_c;
                mw_PCA3=permute(rank_pears_corr_cPCA, [1 3 2]);
                mw_PCA4=permute(rank_spear_corr_cPCA, [1 3 2]);
                evaluat='P-value';
            case 'auc'
                mw_PCA = AUC_c;
                mw_PCA1=mw_cPCA;
                mw_PCA2=AUPR_c;
                mw_PCA3=permute(rank_pears_corr_cPCA, [1 3 2]);
                mw_PCA4=permute(rank_spear_corr_cPCA, [1 3 2]);
                evaluat='AUC';
            case 'aupr'
                mw_PCA = AUPR_c;
                mw_PCA1= mw_cPCA;
                mw_PCA2=AUC_c;
                mw_PCA3=permute(rank_pears_corr_cPCA, [1 3 2]);
                mw_PCA4=permute(rank_spear_corr_cPCA, [1 3 2]);
                evaluat='AUPR';
            case 'pc'
                mw_PCA = permute(rank_pears_corr_cPCA, [1 3 2]);
                mw_PCA1=mw_cPCA;
                mw_PCA2=AUC_c;
                mw_PCA3=AUPR_c;
                mw_PCA4=permute(rank_spear_corr_cPCA, [1 3 2]);
                evaluat='Pearson correlation';
            case 'sc'
                mw_PCA = permute(rank_spear_corr_cPCA, [1 3 2]);
                mw_PCA1=mw_cPCA;
                mw_PCA2=AUC_c;
                mw_PCA3=AUPR_c;
                mw_PCA4=permute(rank_pears_corr_cPCA, [1 3 2]);
                evaluat='Spearman correlation';
        end
    else
        if strcmp(u_rank,'pc')
            mw_PCA = permute(rank_pears_corr_cPCA, [1 3 2]);
            mw_PCA1 = permute(rank_spear_corr_cPCA, [1 3 2]);
            evaluat='Pearson correlation';
        else
            mw_PCA = permute(rank_spear_corr_cPCA, [1 3 2]);
            mw_PCA1= permute(rank_pears_corr_cPCA, [1 3 2]);
            evaluat='Spearman correlation';
        end
    end
    pc = pc_c;
    ttl = 'centred PCA';
    explained = explained_c;
elseif strcmp (u_cent,'n')
    PCA = ncPCA;
    if strcmp(u_lab,'c')
        switch u_rank
            case 'p'
                mw_PCA = mw_ncPCA;
                mw_PCA1=AUC_nc;
                mw_PCA2=AUPR_nc;
                evaluat='P-value';
            case 'auc'
                mw_PCA = AUC_nc;
                mw_PCA1=mw_ncPCA;
                mw_PCA2=AUPR_nc; 
                evaluat='AUC';
            case 'aupr'
                mw_PCA = AUPR_nc;
                mw_PCA1=mw_ncPCA;
                mw_PCA2=AUC_nc;
                evaluat='AUPR';
        end
    elseif strcmp(u_lab,'d')
        switch u_rank
            case 'p'
                mw_PCA = mw_ncPCA;
                mw_PCA1=AUC_nc;
                mw_PCA2=AUPR_nc;
                mw_PCA3=permute(rank_pears_corr_ncPCA, [1 3 2]);
                mw_PCA4=permute(rank_spear_corr_ncPCA, [1 3 2]);
                evaluat='P-value';
            case 'auc'
                mw_PCA = AUC_nc;
                mw_PCA1=mw_ncPCA;
                mw_PCA2=AUPR_nc;
                mw_PCA3=permute(rank_pears_corr_ncPCA, [1 3 2]);
                mw_PCA4=permute(rank_spear_corr_ncPCA, [1 3 2]);
                evaluat='AUC';
            case 'aupr'
                mw_PCA = AUPR_nc;
                mw_PCA1= mw_ncPCA;
                mw_PCA2=AUC_nc;
                mw_PCA3=permute(rank_pears_corr_ncPCA, [1 3 2]);
                mw_PCA4=permute(rank_spear_corr_ncPCA, [1 3 2]);
                evaluat='AUPR';
            case 'pc'
                mw_PCA = permute(rank_pears_corr_ncPCA, [1 3 2]);
                mw_PCA1=mw_ncPCA;
                mw_PCA2=AUC_nc;
                mw_PCA3=AUPR_nc;
                mw_PCA4=permute(rank_spear_corr_ncPCA, [1 3 2]);
                evaluat='Pearson correlation';
            case 'sc'
                mw_PCA = permute(rank_spear_corr_ncPCA, [1 3 2]);
                mw_PCA1=mw_ncPCA;
                mw_PCA2=AUC_nc;
                mw_PCA3=AUPR_nc;
                mw_PCA4=permute(rank_pears_corr_ncPCA, [1 3 2]);
                evaluat='Spearman correlation';
        end
        
    else
        if strcmp(u_rank,'pc')
            mw_PCA = permute(rank_pears_corr_ncPCA, [1 3 2]);
            mw_PCA1 = permute(rank_spear_corr_ncPCA, [1 3 2]);
            evaluat='Pearson correlation';
        else
            mw_PCA = permute(rank_spear_corr_ncPCA, [1 3 2]);
            mw_PCA1 = permute(rank_pears_corr_ncPCA, [1 3 2]);
            evaluat='Spearman correlation';
        end
    end
    pc = pc_nc;
    ttl = 'non-centred PCA';
    explained = explained_nc;
end




ind2 = [];
if strcmp(u_lab,'c') || strcmp(u_lab,'d')
    if strcmp(u_rank,'p')
        
        %P-value
        if numbLabels == 2
            [val,ind1] = min([mw_PCA{u_norm,1,:}]);
            if ind1 == u_dim
                [val,ind2] = min([mw_PCA{u_norm,1,[1:ind1-1,ind1+1:end]}]);
            end
        else
            [val,ind1] = min(mean(cell2mat(mw_PCA(u_norm,:,:))));
            if ind1 == u_dim
                [val,ind2] = min(mean(cell2mat(mw_PCA(u_norm,:,[1:ind1-1,ind1+1:end]))));
            end
        end
        
    elseif strcmp(u_rank,'pc')||strcmp(u_rank,'sc')
        
        %Pearson correlation and Spearman correlation
        [val,ind1] = max(abs([mw_PCA{u_norm,1,:}]));
        if ind1 == u_dim
            [val,ind2] = max(abs([mw_PCA{u_norm,1,[1:ind1-1,ind1+1:end]}]));
        end
        
    elseif strcmp(u_rank,'auc')||strcmp(u_rank,'aupr')
        if numbLabels == 2
            [val,ind1] = max([mw_PCA{u_norm,1,:}]);
            if ind1 == u_dim
                [val,ind2] = max([mw_PCA{u_norm,1,[1:ind1-1,ind1+1:end]}]);
            end
        else
            [val,ind1] = max(mean(cell2mat(mw_PCA(u_norm,:,:))));
            if ind1 == u_dim
                [val,ind2] = max(mean(cell2mat(mw_PCA(u_norm,:,[1:ind1-1,ind1+1:end]))));
            end
        end
    end
    
else
    [val,ind1] = max(abs([mw_PCA{u_norm,1,:}]));
    if ind1 == u_dim
        [val,ind2] = max(abs([mw_PCA{u_norm,1,[1:ind1-1,ind1+1:end]}]));
    end
end

if isempty(ind2)
    ind = ind1;
else
    if ind2 >= ind1
        ind = ind2+1;
    else
        ind = ind2;
    end
end

%% Scatter plot of the desired PCA

%Dimension 1: Best discriminating dimension according to the evaluator
%Dimension 2: Chosen dimension
if numbLabels==2
    if strcmp(u_rank,'p')
        if val<=mw_PCA{u_norm,1,u_dim}
            dim1=ind;
            dim2=u_dim;
        else dim1=u_dim;dim2= ind;
        end
    else
        if abs(val)>=abs(mw_PCA{u_norm,1,u_dim})
            dim1=ind;
            dim2=u_dim;
        else dim1=u_dim;dim2= ind;
        end
    end
else 
    if strcmp(u_rank,'p')
        if mean([mw_PCA{u_norm,:,ind}])<=mean([mw_PCA{u_norm,:,u_dim}])
            dim1=ind;
            dim2=u_dim;
        else dim1=u_dim;dim2= ind;
        end
    elseif  strcmp(u_rank,'auc')||strcmp(u_rank,'aupr')
        if mean([mw_PCA{u_norm,:,ind}])>=mean([mw_PCA{u_norm,:,u_dim}])
            dim1=ind;
            dim2=u_dim;
        else dim1=u_dim;dim2= ind;
        end
    else 
        if abs(val)>=abs(mw_PCA{u_norm,1,u_dim})
            dim1=ind;
            dim2=u_dim;
        else dim1=u_dim;dim2= ind;
        end
    end
end

%%Colours of the groups
%Group on the left of the most discriminative dimension in the PCA plot (median value) : black
%Group on the right of the most discriminative dimension in the PCA plot (median value) : red
%If there are more than two groups, the colours of the other groups are random
figure

for i=1:numbLabels
    a(i)=median(PCA{u_norm}(ismember(labels,nameLabels{i}),dim1));
end
[~,I]=sort(a);
col=cell(1,numbLabels);
col{I(1)}='k';
col{I(numbLabels)}='r';
if numbLabels>2
    for i=2:numbLabels-1
        col{I(i)}=0.15 + rand(1,3)*(0.85-0.15); %creating a vector of three random numbers in between 0.1 and 0.9
    end
end


%%Density plots
% x axis: density plot
subplot(3,3,[8,9])
hold on
for i = 1:numbLabels
    grp{i}=PCA{u_norm}(ismember(labels,nameLabels{i}),dim1);
    [f{i},xi{i}] = ksdensity(grp{i});
    plot(xi{i},f{i},'Color',col{i});
end
axs = [min([xi{:}]) max([xi{:}])];
xlim(axs)

% y axis: density plot

subplot(3,3,[1,4])
hold on
for i = 1:numbLabels
    grp1{i}=PCA{u_norm}(ismember(labels,nameLabels{i}),dim2);
    [f1{i},xi1{i}] = ksdensity(grp1{i});
    plot(xi1{i},f1{i},'Color',col{i});
end
axs1 = [min([xi1{:}]) max([xi1{:}])];
xlim(axs1)
view(90,-90)

% scatter plot

subplot(3,3,[2,3,5,6])

for i = 1:numbLabels
    scatter_plot(PCA{u_norm}(ismember(labels,nameLabels{i}),[dim1 dim2]),col{i},sample_names(ismember(labels,nameLabels{i})),dis)
end

%Labels in the axis
diground = @(x,d) round(x*10^d)/10^d; 
% 2 Groups of labels
if numbLabels == 2
    if strcmp(u_lab,'c')
        %class labels
        switch u_rank
            case 'p'
                %P-value
                if mw_PCA{u_norm,1,dim1}>=0.01
                    pdis1=num2str(diground(mw_PCA{u_norm,1,dim1},3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mw_PCA1{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mw_PCA1{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),' )']) 
                end
                if mw_PCA{u_norm,1,dim2}>=0.01
                    pdis2=num2str(diground(mw_PCA{u_norm,1,dim2},3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mw_PCA1{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mw_PCA1{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),' )'])
                end
            case 'auc'
                %AUC
                if mw_PCA1{u_norm,1,dim1}>=0.01
                    pdis1=num2str(diground(mw_PCA1{u_norm,1,dim1},3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),' )'])
                end
                if mw_PCA1{u_norm,1,dim2}>=0.01
                    pdis2=num2str(diground(mw_PCA1{u_norm,1,dim2},3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),' )'])
                end
                
            case 'aupr'
                %AUPR
                if mw_PCA1{u_norm,1,dim1}>=0.01
                    pdis1=num2str(diground(mw_PCA1{u_norm,1,dim1},3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),' )'])
                end
                if mw_PCA1{u_norm,1,dim2}>=0.01
                    pdis2=num2str(diground(mw_PCA1{u_norm,1,dim2},3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),' )'])
                end
        end
    elseif strcmp(u_lab,'d')
        %discrete labels
        switch u_rank
            case 'p'
                %P-value
                if mw_PCA{u_norm,1,dim1}>=0.01
                    pdis1=num2str(diground(mw_PCA{u_norm,1,dim1},3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mw_PCA1{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mw_PCA1{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                end
                if mw_PCA{u_norm,1,dim2}>=0.01
                    pdis2=num2str(diground(mw_PCA{u_norm,1,dim2},3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mw_PCA1{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mw_PCA1{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                end

            case 'auc'
                %AUC
                if mw_PCA1{u_norm,1,dim1}>=0.01
                    pdis1=num2str(diground(mw_PCA1{u_norm,1,dim1},3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                end
                if mw_PCA1{u_norm,1,dim2}>=0.01
                    pdis2=num2str(diground(mw_PCA1{u_norm,1,dim2},3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                end
                
            case 'aupr'
                %AUPR
                if mw_PCA1{u_norm,1,dim1}>=0.01
                    pdis1=num2str(diground(mw_PCA1{u_norm,1,dim1},3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                end
                if mw_PCA1{u_norm,1,dim2}>=0.01
                    pdis2=num2str(diground(mw_PCA1{u_norm,1,dim2},3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                end
            case 'pc'
                %Pearson correlation 
                if mw_PCA1{u_norm,1,dim1}>=0.01
                    pdis1=num2str(diground(mw_PCA1{u_norm,1,dim1},3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', pears = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', pears = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                end
                if mw_PCA1{u_norm,1,dim2}>=0.01
                    pdis2=num2str(diground(mw_PCA1{u_norm,1,dim2},3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', pears = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', pears = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                end
            case 'sc'
                %Spearman correlation
                if mw_PCA1{u_norm,1,dim1}>=0.01
                    pdis1=num2str(diground(mw_PCA1{u_norm,1,dim1},3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', pears = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim1},3)),', AUPR = ',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', pears = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),' )'])
                end
                if mw_PCA1{u_norm,1,dim2}>=0.01
                    pdis2=num2str(diground(mw_PCA1{u_norm,1,dim2},3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', pears = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mw_PCA2{u_norm,1,dim2},3)),', AUPR = ',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', pears = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),' )'])
                end
        end
    else
        switch u_rank
            case 'pc'
                xlabel(['PC',num2str(dim1),' ( pears = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA1{u_norm,1,dim1},3)),' )'])
                ylabel(['PC',num2str(dim2),' ( pears = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA1{u_norm,1,dim2},3)),' )'])
            case 'sc'
                xlabel(['PC',num2str(dim1),' ( pears = ',num2str(diground(mw_PCA1{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),' )'])
                ylabel(['PC',num2str(dim2),' ( pears = ',num2str(diground(mw_PCA1{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),' )'])
        end
    end
end

% More than 2 Groups of labels
if numbLabels > 2
    if strcmp(u_lab,'c')
        %class labels
        switch u_rank
            case 'p'
                %P-value
                if mean([mw_PCA{u_norm,:,dim1}])>=0.01
                    pdis1=num2str(diground(mean([mw_PCA{u_norm,:,dim1}]),3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA1{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA1{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),' )'])
                end
                if mean([mw_PCA{u_norm,:,dim2}])>=0.01
                    pdis2=num2str(diground(mean([mw_PCA{u_norm,:,dim2}]),3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA1{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA1{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),' )'])

                end
 
            case 'auc'
                %AUC
                if mean([mw_PCA1{u_norm,:,dim1}])>=0.01
                    pdis1=num2str(diground(mean([mw_PCA1{u_norm,:,dim1}]),3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),' )'])
                end
                if mean([mw_PCA1{u_norm,:,dim2}])>=0.01
                    pdis2=num2str(diground(mean([mw_PCA1{u_norm,:,dim2}]),3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),' )'])
                end
                
 
                
            case 'aupr'
                %AUPR
                if mean([mw_PCA1{u_norm,:,dim1}])>=0.01
                    pdis1=num2str(diground(mean([mw_PCA1{u_norm,:,dim1}]),3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA{u_norm,:,dim1}]),3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA{u_norm,:,dim1}]),3)),' )'])
                end
                if mean([mw_PCA1{u_norm,:,dim2}])>=0.01
                    pdis2=num2str(diground(mean([mw_PCA1{u_norm,:,dim2}]),3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA{u_norm,:,dim2}]),3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA{u_norm,:,dim2}]),3)),' )'])
                end
   
        end
        
    elseif strcmp(u_lab,'d')
        %discrete labels
        switch u_rank
            case 'p'
                %P-value
                if mean([mw_PCA{u_norm,:,dim1}])>=0.01
                    pdis1=num2str(diground(mean([mw_PCA{u_norm,:,dim1}]),3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA1{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA1{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                end
                if mean([mw_PCA{u_norm,:,dim2}])>=0.01
                    pdis2=num2str(diground(mean([mw_PCA{u_norm,:,dim2}]),3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA1{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA1{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                end
                

            case 'auc'
                %AUC
                if mean([mw_PCA1{u_norm,:,dim1}])>=0.01
                    pdis1=num2str(diground(mean([mw_PCA1{u_norm,:,dim1}]),3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),', pears =',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),', pears =',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                end
                if mean([mw_PCA1{u_norm,:,dim2}])>=0.01
                    pdis2=num2str(diground(mean([mw_PCA1{u_norm,:,dim2}]),3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),', pears =',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),', pears =',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                end
                
                
            case 'aupr'
                %AUPR
                if mean([mw_PCA1{u_norm,:,dim1}])>=0.01
                    pdis1=num2str(diground(mean([mw_PCA1{u_norm,:,dim1}]),3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA{u_norm,:,dim1}]),3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA{u_norm,:,dim1}]),3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                end
                if mean([mw_PCA1{u_norm,:,dim2}])>=0.01
                    pdis2=num2str(diground(mean([mw_PCA1{u_norm,:,dim2}]),3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA{u_norm,:,dim2}]),3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA{u_norm,:,dim2}]),3)),', pears = ',num2str(diground(mw_PCA3{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                end
                
                
            case 'pc'
                %Pearson correlation 
                if mean([mw_PCA1{u_norm,:,dim1}])>=0.01
                    pdis1=num2str(diground(mean([mw_PCA1{u_norm,:,dim1}]),3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA3{u_norm,:,dim1}]),3)),', pears = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA3{u_norm,:,dim1}]),3)),', pears = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),' )'])
                end
                if mean([mw_PCA1{u_norm,:,dim2}])>=0.01
                    pdis2=num2str(diground(mean([mw_PCA1{u_norm,:,dim2}]),3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA3{u_norm,:,dim2}]),3)),', pears = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA3{u_norm,:,dim2}]),3)),', pears = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),' )'])
                end
                
                
            case 'sc'
                %Spearman correlation
                if mean([mw_PCA1{u_norm,:,dim1}])>=0.01
                    pdis1=num2str(diground(mean([mw_PCA1{u_norm,:,dim1}]),3));
                    xlabel(['PC',num2str(dim1),' ( p = ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA3{u_norm,:,dim1}]),3)),', pears = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),' )'])
                else
                    pdis1='p < 0.01';
                    xlabel(['PC',num2str(dim1),' ( ',pdis1,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim1}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA3{u_norm,:,dim1}]),3)),', pears = ',num2str(diground(mw_PCA4{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),' )'])
                end
                if mean([mw_PCA1{u_norm,:,dim2}])>=0.01
                    pdis2=num2str(diground(mean([mw_PCA1{u_norm,:,dim2}]),3));
                    ylabel(['PC',num2str(dim2),' ( p = ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA3{u_norm,:,dim2}]),3)),', pears = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),' )'])
                else
                    pdis2='p < 0.01';
                    ylabel(['PC',num2str(dim2),' ( ',pdis2,', AUC = ',num2str(diground(mean([mw_PCA2{u_norm,:,dim2}]),3)),', AUPR = ',num2str(diground(mean([mw_PCA3{u_norm,:,dim2}]),3)),', pears = ',num2str(diground(mw_PCA4{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),' )'])
                end
                
        
        end
    else 
        switch u_rank
            case 'pc'
                xlabel(['PC',num2str(dim1),' ( pears = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA1{u_norm,1,dim1},3)),' )'])
                ylabel(['PC',num2str(dim2),' ( pears = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA1{u_norm,1,dim2},3)), ' )'])
            case 'sc'
                xlabel(['PC',num2str(dim1),' ( pears = ',num2str(diground(mw_PCA1{u_norm,1,dim1},3)),', spear = ',num2str(diground(mw_PCA{u_norm,1,dim1},3)),' )'])
                ylabel(['PC',num2str(dim2),' ( pears = ',num2str(diground(mw_PCA1{u_norm,1,dim2},3)),', spear = ',num2str(diground(mw_PCA{u_norm,1,dim2},3)),' )'])
        end
    end
end


title ([ttl,' for norm: ',u_norm_n])
legend(nameLabels{:})
xlim(axs)
ylim(axs1)



%% bar plots
figure
plot_1=subplot(2,2,[1 2]);
if numbLabels == 2
    if ~strcmp(u_rank,'pc') && ~strcmp(u_rank,'sc')
        if strcmp(u_rank,'p')
            plot_y= cell2mat(mw_PCA(u_norm,:));
            indexmin = find(min(plot_y) == plot_y);
            for i=1:size(ncPCA{1},2)
                if i~=indexmin
                    hold on
                    h(i)=bar(i,cell2mat(mw_PCA(u_norm,i)));
                else
                    h(i)=bar(i,cell2mat(mw_PCA(u_norm,i)),'facecolor','cyan');
                end
            end
            xlim1=get(gca,'xlim');
            hold on
            plot(xlim1,[0.05 0.05],'r--')
            text(size(x,1)+1,0.05,'0.05','Color','red','VerticalAlignment','middle')
            for i =1:length(indexmin)
                ymin= plot_y(indexmin(i));
                xmin=indexmin(i);  
                if ymin<0.01
                    text(xmin,ymin,'<0.01',...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','bottom','FontSize',9,'Color','red')
                else
                    text(xmin,ymin,num2str(ymin,'%0.2f'),...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','bottom','FontSize',9,'Color','red')
                end
            end
        else
            plot_y= cell2mat(mw_PCA(u_norm,:));
            indexmax = find(max(plot_y) == plot_y);
            for i=1:size(ncPCA{1},2)
                if i~=indexmax
                    hold on
                    h(i)=bar(i,cell2mat(mw_PCA(u_norm,i)));
                else
                    h(i)=bar(i,cell2mat(mw_PCA(u_norm,i)),'facecolor','cyan');
                end
            end
            xlim1=get(gca,'xlim');
            hold on
            plot(xlim1,[0.5 0.5],'r--')
            for i=1:length(indexmax)
                ymax = plot_y(indexmax(i));
                xmax=indexmax(i);
%                 if mw_PCA1{u_norm,1,indexmax(i)}<0.01
%                     str_pval='(<0.01*)';
%                 else
%                     str_pval=strcat('(',num2str(mw_PCA1{u_norm,1,indexmax(i)},'%0.2f'),'*)');
%                 end
                text(xmax,ymax,num2str(ymax,'%0.2f'),...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','bottom','FontSize',9,'Color','red')
            end
        end
    else
        a1=find(cell2mat(mw_PCA(u_norm,:))>=0);
        b1=find(cell2mat(mw_PCA(u_norm,:))<0);
        plot_y = abs(cell2mat(mw_PCA(u_norm,:)));
        indexmax = find(max(plot_y) == plot_y);
        for i=1:size(ncPCA{1},2)
            if ismember(i,a1)
               hold on
               h(i)=bar(i,abs(cell2mat(mw_PCA(u_norm,i))),'FaceColor','k');
            else
               h(i)=bar(i,abs(cell2mat(mw_PCA(u_norm,i))),'FaceColor',[0.5 0.5 0.5]);
            end
        end
        xlim1=get(gca,'xlim');
        hold on
        plot(xlim1,[0.5 0.5],'r--')
        for i=1:length(indexmax)
            ymax = plot_y(indexmax(i));
            xmax=indexmax(i);
            text(xmax,ymax,num2str(round(ymax,2),'%0.2f'),...
                'HorizontalAlignment','center',...
                'VerticalAlignment','bottom','FontSize',9,'Color','red')
        end
        if strcmp(u_rank,'pc')
            legend([h(a1(1)) h(b1(1))],'Pear. corr\geq0','Pear. corr< 0')
        else
            legend([h(a1(1)) h(b1(1))],'Spear. corr\geq0','Spear. corr< 0')
        end
    end
else
    if strcmp(u_lab,'c') || strcmp(u_lab,'d')
        if ~strcmp(u_rank,'pc') && ~strcmp(u_rank,'sc')
            if strcmp(u_rank,'p')
                plot_y= permute(mean(cell2mat(mw_PCA(u_norm,:,:))),[3 1 2]);
                indexmin = find(min(plot_y) == plot_y);
                for i=1:size(ncPCA{1},2)
                    if i~=indexmin
                        hold on
                        h(i)=bar(i,plot_y(i));
                    else
                        h(i)=bar(i,plot_y(i),'facecolor','cyan');
                    end
                end
                xlim1=get(gca,'xlim');
                hold on
                plot(xlim1,[0.05 0.05],'r--')
                text(size(x,1)+1,0.05,'0.05','Color','red','VerticalAlignment','middle')
                for i=1:length(indexmin)
                    ymin = plot_y(indexmin(i));
                    xmin=indexmin(i);
                    if ymin<0.01
                        text(xmin,ymin,'<0.01',...
                            'HorizontalAlignment','center',...
                            'VerticalAlignment','bottom','FontSize',9,'Color','red')
                    else
                        text(xmin,ymin,num2str(ymin,'%0.2f'),...
                            'HorizontalAlignment','center',...
                            'VerticalAlignment','bottom','FontSize',9,'Color','red')
                    end
                end
            else
                plot_y= permute(mean(cell2mat(mw_PCA(u_norm,:,:))),[3 1 2]);
                indexmax = find(max(plot_y) == plot_y);
                for i=1:size(ncPCA{1},2)
                    if i~=indexmax
                        hold on
                        h(i)=bar(i,plot_y(i));
                    else
                        h(i)=bar(i,plot_y(i),'facecolor','cyan');
                    end
                end
                xlim1=get(gca,'xlim');
                hold on
                plot(xlim1,[0.5 0.5],'r--')
                for i=1:length(indexmax)
                    ymax = plot_y(indexmax(i));
                    xmax=indexmax(i);
%                     if mean([mw_PCA1{u_norm,:,indexmax(i)}])<0.01
%                         str_pval='(<0.01*)';
%                     else
%                         str_pval=strcat('(',num2str(mean([mw_PCA1{u_norm,:,indexmax(i)}]),'%0.2f'),'*)');
%                     end
%                     text(xmax,ymax,strcat(num2str(ymax,'%0.2f'),str_pval),...
                    text(xmax,ymax,num2str(ymax,'%0.2f'),...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','bottom','FontSize',9,'Color','red')
                end
            end
        else
            a1=find(cell2mat(mw_PCA(u_norm,:))>=0);
            b1=find(cell2mat(mw_PCA(u_norm,:))<0);
            plot_y = abs(cell2mat(mw_PCA(u_norm,:)));
            indexmax = find(max(plot_y) == plot_y);
            for i=1:size(ncPCA{1},2)
                if ismember(i,a1)
                    hold on
                    h(i)=bar(i,abs(cell2mat(mw_PCA(u_norm,i))),'FaceColor','k');
                else
                    h(i)=bar(i,abs(cell2mat(mw_PCA(u_norm,i))),'FaceColor',[0.5 0.5 0.5]);
                end
            end
            xlim1=get(gca,'xlim');
            hold on
            plot(xlim1,[0.5 0.5],'r--')
            for i=1:length(indexmax)
                ymax = plot_y(indexmax(i));
                xmax=indexmax(i);
                text(xmax,ymax,num2str(round(ymax,2),'%0.2f'),...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','bottom','FontSize',9,'Color','red')
            end
            if strcmp(u_rank,'pc')
                legend([h(a1(1)) h(b1(1))],'Pear. corr\geq0','Pear. corr< 0')
            else
                legend([h(a1(1)) h(b1(1))],'Spear. corr\geq0','Spear. corr< 0')
            end
        end
    else
        a1=find(cell2mat(mw_PCA(u_norm,:))>=0);
        b1=find(cell2mat(mw_PCA(u_norm,:))<0);
        plot_y = abs(cell2mat(mw_PCA(u_norm,:)));
        indexmax = find(max(plot_y) == plot_y);
        for i=1:size(ncPCA{1},2)
            if ismember(i,a1)
                hold on
                h(i)=bar(i,abs(cell2mat(mw_PCA(u_norm,i))),'FaceColor','k');
            else
                h(i)=bar(i,abs(cell2mat(mw_PCA(u_norm,i))),'FaceColor',[0.5 0.5 0.5]);
            end
        end
        xlim1=get(gca,'xlim');
        hold on
        plot(xlim1,[0.5 0.5],'r--')
        for i=1:length(indexmax)
            ymax = plot_y(indexmax(i));
            xmax=indexmax(i);
            text(xmax,ymax,num2str(round(ymax,2),'%0.2f'),...
                'HorizontalAlignment','center',...
                'VerticalAlignment','bottom','FontSize',9,'Color','red')
        end
        if strcmp(u_rank,'pc')
            legend([h(a1(1)) h(b1(1))],'Pear. corr\geq0','Pear. corr< 0')
        else
            legend([h(a1(1)) h(b1(1))],'Spear. corr\geq0','Spear. corr< 0')
        end
    end
end


title([evaluat,'s over principal components in ',ttl,' for norm: ',u_norm_n])

ticknumb = 1:size(x,1);
set(gca,'XTick',ticknumb)
set(gca,'XTickLabelRotation',90)
xlabel('PC')
xlim([0 size(ncPCA{1},2)+1])
if strcmp(u_rank,'auc') || strcmp(u_rank,'aupr')
    ylim([0 1])
end
if ~strcmp(u_rank,'pc') && ~strcmp(u_rank,'sc')
    ylabel(evaluat)
else 
    ylabel(['|',evaluat,'|'])
end

box(plot_1,'on')

plot_2=subplot(2,2,[3 4]);
if strcmp(u_rank,'p') %p-value
    for i=1:size(ncPCA{1},2)
        if i~=indexmin
            hold on
            h1(i)=bar(i,explained{u_norm}(i));
        else
            h1(i)=bar(i,explained{u_norm}(i),'facecolor','c');
        end
    end
elseif strcmp(u_rank,'auc') || strcmp(u_rank,'aupr')
    for i=1:size(ncPCA{1},2)
        if i~=indexmax
            hold on
            h1(i)=bar(i,explained{u_norm}(i));
        else
            h1(i)=bar(i,explained{u_norm}(i),'facecolor','c');
        end
    end
else
    for i=1:size(ncPCA{1},2)
        if i~=indexmax
            hold on
            h1(i)=bar(i,explained{u_norm}(i));
        else
            h1(i)=bar(i,explained{u_norm}(i),'facecolor','c');
        end
    end
end

title('Explained variance for the respective principal components')
set(gca,'XTick',ticknumb)
set(gca,'XTickLabelRotation',90)
xlim([0 size(ncPCA{1},2)+1])
xlabel('PC')
ylabel('Explained Variance (%)')
box(plot_2,'on')



% User interaction for cut-off
flag = 0;
while flag == 0
    fprintf('\nSelect a cut-off or a set of cut-offs for generating the PC-corr network [number between 0 and 1]:\n\n');
    fprintf('Examples: 0.6 or [0.6 0.65 0.7]\n\n');
    cutoff = input('-> ');
    if max(cutoff)>= 0 && max(cutoff)<= 1
        flag = 1;
    else
        fprintf('\nPlease introduce a correct cut-off or set of cut-offs (in [0,1]).\n');
    end
end

%%%%% PC-corr %%%%%
if length(cutoff)==1
    [Edges,Nodes,pc_corr,x1,cutoff_f] = C_corr(norm{u_norm},pc{u_norm}(:,u_dim),feat_names,cutoff);
    
    
    [NodeColor,n1_f,n2_f]=match_V_samplCol(col,x1,labels, nameLabels, Nodes);
    Nodes(1:end,3)=Nodes(:,2);
    Nodes{1,2}='Colour';
    Nodes(2:end,2)=NodeColor;
    
    cut_off=cutoff_f;
    edges=Edges;
    nodes=Nodes;
    if (size(edges,1)-1)==1
        fprintf('\nAt cut-off %0.2f, the PC-corr network has %i nodes and %i edge.\n\n', cut_off,size(nodes,1)-1,size(edges,1)-1);
    else
        fprintf('\nAt cut-off %0.2f, the PC-corr network has %i nodes and %i edges.\n\n', cut_off,size(nodes,1)-1,size(edges,1)-1);
    end
else 
    for i=1:length(cutoff)

        [Edges{i,2},Nodes{i,2},pc_corr{i,2},x1,cutoff_f(i)] = C_corr(norm{u_norm},pc{u_norm}(:,u_dim),feat_names,cutoff(i));
        Edges{i,1}=cutoff_f(i);
        Nodes{i,1}=cutoff_f(i);
        pc_corr{i,1}=cutoff_f(i);
        x2{i,1}=cutoff_f(i);
        
        x2{i,2}=x1;
        
        [NodeColor,n1_f,n2_f]=match_V_samplCol(col,x2{i,2},labels, nameLabels, Nodes{i,2});
        Nodes{i,2}(1:end,3)=Nodes{i,2}(:,2);
        Nodes{i,2}{1,2}='Colour';
        Nodes{i,2}(2:end,2)=NodeColor;
        
        cut_off(i)=cutoff_f(i);
        edges=Edges{i,2};
        nodes=Nodes{i,2};
        if (size(edges,1)-1)==1
            fprintf('\nAt cut-off %0.2f, the PC-corr network has %i nodes and %i edge.\n\n', cut_off(i),size(nodes,1)-1,size(edges,1)-1);
        else
            fprintf('\nAt cut-off %0.2f, the PC-corr network has %i nodes and %i edges.\n\n', cut_off(i),size(nodes,1)-1,size(edges,1)-1);
        end
    end
end

filename = 'PC-corr_net_edges-nodes_results.xlsx';
if exist('PC-corr_net_edges-nodes_results.xlsx','file')~=0
    delete('PC-corr_net_edges-nodes_results.xlsx')
end

if length(cut_off) == 0
   disp(size(Edges{:,2}));
    Edges_data = Edges{:,2};  % Ensure Edges{:,2} is the correct 2D array
    if iscolumn(Edges_data)
        Edges_data = Edges_data';  % Convert to row if it's a column vector
    end
    
    Nodes_data = Nodes{:,2};  % Same for Nodes
    if iscolumn(Nodes_data)
        Nodes_data = Nodes_data';  % Convert to row if it's a column vector
    end
    
    % Create tables from Edges and Nodes
    Edges_table = cell2table(Edges_data, 'VariableNames', Edges{1,1});
    Nodes_table = cell2table(Nodes_data, 'VariableNames', Nodes{1,1});
    
    writetable(Edges_table, filename, 'Sheet', ['Edges-cutoff ', num2str(cut_off)]);
    writetable(Nodes_table, filename, 'Sheet', ['Nodes-cutoff ', num2str(cut_off)]);
    RemoveSheet123(filename);
else
    for i = 1:length(cut_off)
        % Create tables from Edges and Nodes for each cutoff
        Edges_table = cell2table(Edges(i,2), 'VariableNames', Edges(1,1));
        Nodes_table = cell2table(Nodes(i,2), 'VariableNames', Nodes(1,1));
        
        % Write each set of tables to separate sheets
        writetable(Edges_table, filename, 'Sheet', ['Edges-cutoff ', num2str(cut_off(i))]);
        writetable(Nodes_table, filename, 'Sheet', ['Nodes-cutoff ', num2str(cut_off(i))]);
    end
    % RemoveSheet123(filename);
end





fprintf('\nYou can find the tables of edge and node values of the PC-corr network in your current folder in an Excel file (in separate spreadsheet), named \n')
fprintf('PC-corr_net_edges-nodes_results.xlsx, so that you can easily visualize the graph with another network visualization program.\n\n')

flag = 0;
while flag == 0
    if length(cut_off)==1
        fprintf('\nDo you want to visualize the PC-corr network? [y/n]:\n\n');
    else
        fprintf('\nDo you want to visualize the PC-corr networks? [y/n]:\n\n');
    end
    u_vis_net = input('-> ','s');
    if strcmp(u_vis_net,'n') || strcmp(u_vis_net,'y')
        flag = 1;
    else
        fprintf('Please introduce just "y" or "n"\n');
    end
end

%% Graph plot for MATLAB 2015b and future releases

% g = graph(pc_corr);
% figure
% h = plot(g);
% labelnode(h,[1:size(nodes,1)-1],{nodes{2:end,1}})
% idx = [1:size(nodes,1)-1];
% highlight(h,[idx(sign([nodes{2:end,2}]) == -1)],'NodeColor','k')
% highlight(h,[idx(sign([nodes{2:end,2}]) == 1)],'NodeColor','r')
% highlight(h,[g.Edges.EndNodes(find(sign(g.Edges.Weight) == 1),1)],[g.Edges.EndNodes(find(sign(g.Edges.Weight) == 1),2)],'EdgeColor','r');
% highlight(h,[g.Edges.EndNodes(find(sign(g.Edges.Weight) == -1),1)],[g.Edges.EndNodes(find(sign(g.Edges.Weight) == -1),2)],'EdgeColor','k');

%% Biograph approach for generation of graph plot for MATLAB 2015a and previous releases
if strcmp(u_vis_net,'y')
    %Visualize the PC-corr networks
    for k=1:length(cut_off)
        if length(cut_off)==1
            pc_corr1=pc_corr;
            nodes1=Nodes;
        else
            pc_corr1=pc_corr{k,2};
            nodes1=Nodes{k,2};
        end
        
        if isempty(pc_corr1)
            continue
        else
            
            plot_graph(pc_corr1,nodes1,cut_off(k))
            
        end
    end
end




%% User interaction for the significance of measure of segregation 
if strcmp(u_lab,'c') || strcmp (u_lab, 'd')
    flag = 0;
    fprintf('\nAdditionally, you can compute the trustworthiness of p-value, AUC and AUPR results, that is a measure of significance of segregation.\n');
    fprintf('The trustworthiness (of p-value/AUC/AUPR results) evaluates if the results are discriminative because it is random (trustworthiness p-value>0.05)\nor it captures main variability (trustworthiness p-value<=0.05).\n');
    while flag == 0
        fprintf('\nWould you like to compute the trustworthiness of measure of segregation (p-value,AUC and AUPR)? [1/2/3]\n[1] no \n[2] yes, but just for the selected case (norm, dimension and centering)\n[3] yes, for all the cases\n\n');
        meas_seg=input('-> ');
        if (meas_seg~=1)&&(meas_seg~=2) && (meas_seg~=3)
            flag = 0;
            fprintf('\nPlease introduce either 1 for no calculation of trustworthiness, 2 for its calculation for the selected case (normalization, dimension and centering)\n or 3 for its calculation for all the cases. \n')
        else
            flag = 1;
        end
    end
end



if meas_seg == 3
    flag = 0;
    while flag == 0
        if (size(results,1)*nchoosek(numbLabels,2))> 400
            fprintf('\nAre you sure? Do you want to compute for all? It will take some time (hours). [y/n]');
            fprintf('\n(We suggest running it overnight)');
        else 
            fprintf('\nAre you sure? Do you want to compute for all? It will take some time (minutes). [y/n]');
        end
        fprintf('\n[y] yes, I want to compute for all the cases\n[n] no, I don’t want to compute for all the cases\n\n')
        meas_seg_all = input('-> ','s');
        if strcmp(meas_seg_all,'y') || strcmp(meas_seg_all,'n')
            flag = 1;
        else
            fprintf('\nPlease introduce just "y" or "n"\n');
        end
    end
end

if (meas_seg == 2) || (meas_seg == 3)
    labels_rp=random_permutation_labels(1000,labels);
end

%% Calculation of trustworthiness (pvalue,AUC,AUPR) for the selected options: norm, dimension and centred/non-centred version of PCA
if meas_seg == 2
    if numbLabels == 2 %two groups
        if strcmp(u_cent,'y')
            u_cent_n='yes';
        elseif strcmp(u_cent,'n')
            u_cent_n='no';
        end
        
        if strcmp(u_lab,'c')
            chosen_opt=strcmp(results_xls(:,4),u_norm_n)& strcmp(results_xls(:,5),u_cent_n)& (cell2mat(results_xls(:,6))==u_dim);
        elseif strcmp(u_lab,'d')
            chosen_opt=strcmp(results_xls(:,6),u_norm_n)& strcmp(results_xls(:,7),u_cent_n)& (cell2mat(results_xls(:,8))==u_dim);
        end
        
        idx_opt= find(chosen_opt);
        mw_opt=cell2mat(results_xls(idx_opt,1));
        auc_opt=cell2mat(results_xls(idx_opt,2));
        aupr_opt=cell2mat(results_xls(idx_opt,3));
        numb_rand = 1000;
        %Significance of segragation (pvalue,AUC,AUPR)
        pVal_opt=significance_segragation(numb_rand,PCA{u_norm},u_dim,labels_rp,nameLabels,numbLabels,[mw_opt auc_opt aupr_opt],'all',u_aupr,u_aupr_r);
        pVal_mw_opt=pVal_opt(1);
        pVal_auc_opt=pVal_opt(2);
        pVal_aupr_opt=pVal_opt(3);
        idx_excel= find(strcmp(header_xls,'Trustworthiness(p-value)'));
        results_xls(idx_opt,idx_excel:(idx_excel+2))=num2cell([pVal_mw_opt pVal_auc_opt pVal_aupr_opt]);
        xlswrite('results.xlsx',vertcat(header_xls,results_xls),'PCA results')
        
        
        
    else %more than two groups
        if strcmp(u_cent,'y')
            u_cent_n='yes';
        elseif strcmp(u_cent,'n')
            u_cent_n='no';
        end
        
        if strcmp(u_lab,'c')
            chosen_opt=strcmp(results(:,4),u_norm_n)& strcmp(results(:,5),u_cent_n)& (cell2mat(results(:,6))==u_dim);
        elseif strcmp(u_lab,'d')
            chosen_opt=strcmp(results(:,6),u_norm_n)& strcmp(results(:,7),u_cent_n)& (cell2mat(results(:,8))==u_dim);
        end
        idx_opt= find(chosen_opt);
        mw_opt=cell2mat(results(idx_opt,1));
        auc_opt=cell2mat(results(idx_opt,2));
        aupr_opt=cell2mat(results(idx_opt,3));
        numb_rand = 1000;
        %Significance of segragation (pvalue,AUC,AUPR)
        pVal_opt=significance_segragation(numb_rand,PCA{u_norm},u_dim,labels_rp,nameLabels,numbLabels,[mw_opt auc_opt aupr_opt],'all',u_aupr,u_aupr_r);
        pVal_mw_opt=pVal_opt(1);
        pVal_auc_opt=pVal_opt(2);
        pVal_aupr_opt=pVal_opt(3);
        idx_excel= find(strcmp(header_xls,'Trustworthiness(p-value)'));
        results_xls(idx_opt,idx_excel:(idx_excel+2))=num2cell([pVal_mw_opt pVal_auc_opt pVal_aupr_opt]);
        xlswrite('results.xlsx',vertcat(header_xls,results_xls),'PCA results')
         
        
    end
    
end
%% Calculation of trustworthiness (p-value,AUC and AUPR) for each case
if meas_seg == 3
    if strcmp(meas_seg_all,'y')
        
        if numbLabels == 2 %two groups
            
            if strcmp(u_lab,'c')
                tic
                idx_all=1:size(results_xls,1);
                [ia1,ia2]=unique(round(idx_all/length(idx_all),2)*100);
                %[ia1,ia2]=unique(round(idx_all/length(idx_all),1)*100);
                %ia2=ia2((mod(ia1,10)==0));
                fprintf('\n');
                textprogressbar('calculating trustworthiness: ');
                for i=1:length(idx_all)
                    numb_rand = 1000;
                    nm= find(strcmp(norms_list(:),results_xls{idx_all(i),4}));
                    if strcmp(results_xls{idx_all(i),5},'yes')
                        PCA_all=cPCA{nm};
                    elseif strcmp(results_xls{idx_all(i),5},'no')
                        PCA_all=ncPCA{nm};
                    end
                    
                    pVal_all(i,:)=significance_segragation(numb_rand,PCA_all,results_xls{idx_all(i),6},labels_rp,nameLabels,numbLabels,[results_xls{idx_all(i),1} results_xls{idx_all(i),2} results_xls{idx_all(i),3}],'all',u_aupr,u_aupr_r);
                    pVal_mw_all(i)=pVal_all(i,1);
                    pVal_auc_all(i)=pVal_all(i,2);
                    pVal_aupr_all(i)=pVal_all(i,3);
                    idx_excel_all= find(strcmp(header_xls,'Trustworthiness(p-value)'));
                    results_xls(idx_all(i),idx_excel_all:(idx_excel_all+2))=num2cell([pVal_mw_all(i) pVal_auc_all(i) pVal_aupr_all(i)]);
                  
                    if ismember(idx_all(i),ia2)
                        textprogressbar(round(idx_all(i)/length(idx_all),2)*100);
                    end
                end
                textprogressbar('done');
                xlswrite('results.xlsx',vertcat(header_xls,results_xls),'PCA results')
                toc
                
            elseif strcmp(u_lab,'d')
                tic
                idx_all=1:size(results_xls,1);
                [ia1,ia2]=unique(round(idx_all/length(idx_all),2)*100);
                %[ia1,ia2]=unique(round(idx_all/length(idx_all),1)*100);
                %ia2=ia2((mod(ia1,10)==0));
                fprintf('\n');
                textprogressbar('calculating trustworthiness: ');
                for i=1:length(idx_all)
                    numb_rand = 1000;
                    nm= find(strcmp(norms_list(:),results_xls{idx_all(i),6}));
                    if strcmp(results_xls{idx_all(i),7},'yes')
                        PCA_all=cPCA{nm};
                    elseif strcmp(results_xls{idx_all(i),7},'no')
                        PCA_all=ncPCA{nm};
                    end
                        pVal_all(i,:)=significance_segragation(numb_rand,PCA_all,results_xls{idx_all(i),8},labels_rp,nameLabels,numbLabels,[results_xls{idx_all(i),1} results_xls{idx_all(i),2} results_xls{idx_all(i),3}],'all',u_aupr,u_aupr_r);
                        pVal_mw_all(i)=pVal_all(i,1);
                        pVal_auc_all(i)=pVal_all(i,2);
                        pVal_aupr_all(i)=pVal_all(i,3);
                        idx_excel_all= find(strcmp(header_xls,'Trustworthiness(p-value)'));
                        results_xls(idx_all(i),idx_excel_all:(idx_excel_all+2))=num2cell([pVal_mw_all(i) pVal_auc_all(i) pVal_aupr_all(i)]);

                    if ismember(idx_all(i),ia2)
                        textprogressbar(round(idx_all(i)/length(idx_all),2)*100);
                    end
                end
                textprogressbar('done');
                xlswrite('results.xlsx',vertcat(header_xls,results_xls),'PCA results')
                toc
            end
            
        else %more than two groups
            
            if strcmp(u_lab,'c')
                tic
                idx_all=1:size(results,1);
                [ia1,ia2]=unique(round(idx_all/length(idx_all),2)*100);
                %[ia1,ia2]=unique(round(idx_all/length(idx_all),1)*100);
                %ia2=ia2((mod(ia1,10)==0));
                fprintf('\n');
                textprogressbar('calculating trustworthiness: ');
                for i=1:length(idx_all)
                    numb_rand = 1000;
                    nm= find(strcmp(norms_list(:),results{idx_all(i),4}));
                    if strcmp(results{idx_all(i),5},'yes')
                        PCA_all=cPCA{nm};
                    elseif strcmp(results{idx_all(i),5},'no')
                        PCA_all=ncPCA{nm};
                        
                    end
                    pVal_all(i,:)=significance_segragation(numb_rand,PCA_all,results{idx_all(i),6},labels_rp,nameLabels,numbLabels,[results{idx_all(i),1} results{idx_all(i),2} results{idx_all(i),3}],'all',u_aupr,u_aupr_r);
                    pVal_mw_all(i)=pVal_all(i,1);
                    pVal_auc_all(i)=pVal_all(i,2);
                    pVal_aupr_all(i)=pVal_all(i,3);
                    idx_excel_all= find(strcmp(header_xls,'Trustworthiness(p-value)'));
                    results_xls(idx_all(i),idx_excel_all:(idx_excel_all+2))=num2cell([pVal_mw_all(i) pVal_auc_all(i) pVal_aupr_all(i)]);

                    if ismember(idx_all(i),ia2)
                        textprogressbar(round(idx_all(i)/length(idx_all),2)*100);
                        %fprintf(strcat(' > Progress: ',{' '},num2str(round(idx_all(i)/length(idx_all),1)*100),'%%\n'))
                    end
                end
                textprogressbar('done');
                xlswrite('results.xlsx',vertcat(header_xls,results_xls),'PCA results')
                toc
            elseif strcmp(u_lab,'d')
                tic
                idx_all=1:size(results,1);
                [ia1,ia2]=unique(round(idx_all/length(idx_all),2)*100);
                %[ia1,ia2]=unique(round(idx_all/length(idx_all),1)*100);
                %ia2=ia2((mod(ia1,10)==0));
                fprintf('\n');
                textprogressbar('calculating trustworthiness: ');
                for i=1:length(idx_all)
                    numb_rand = 1000;
                    nm= find(strcmp(norms_list(:),results{idx_all(i),6}));
                    if strcmp(results{idx_all(i),7},'yes')
                        PCA_all=cPCA{nm};
                    elseif strcmp(results{idx_all(i),7},'no')
                        PCA_all=ncPCA{nm};
                        
                    end
                    
                    pVal_all(i,:)=significance_segragation(numb_rand,PCA_all,results{idx_all(i),8},labels_rp,nameLabels,numbLabels,[results{idx_all(i),1} results{idx_all(i),2} results{idx_all(i),3}],'all',u_aupr,u_aupr_r);
                    pVal_mw_all(i)=pVal_all(i,1);
                    pVal_auc_all(i)=pVal_all(i,2);
                    pVal_aupr_all(i)=pVal_all(i,3);
                    idx_excel_all= find(strcmp(header_xls,'Trustworthiness(p-value)'));
                    results_xls(idx_all(i),idx_excel_all:(idx_excel_all+2))=num2cell([pVal_mw_all(i) pVal_auc_all(i) pVal_aupr_all(i)]);
                    
                    if ismember(idx_all(i),ia2)
                        textprogressbar(round(idx_all(i)/length(idx_all),2)*100);
                    end
                end
                textprogressbar('done');
                xlswrite('results.xlsx',vertcat(header_xls,results_xls),'PCA results')
                toc
            end
        end
        
    end
end





function possClass=positive_label_opt(u_aupr, nameLabels1,nameLabels2, labels,u_aupr_r)
    
    if strcmp(u_aupr,'s')
        len_n=length(labels(ismember(labels,nameLabels1)));
        len_m=length(labels(ismember(labels,nameLabels2)));
        if len_n<len_m
            possClass=nameLabels1;
        elseif len_m<len_n
            possClass=nameLabels2;
        end
    elseif strcmp(u_aupr,'l')
        len_n=length(labels(ismember(labels,nameLabels1)));
        len_m=length(labels(ismember(labels,nameLabels2)));
        if len_n<len_m
            possClass=nameLabels2;
        elseif len_m<len_n
            possClass=nameLabels1;
        end
    elseif strcmp(u_aupr,'r')
        idx_n=find(strcmp(u_aupr_r,nameLabels1));
        idx_m=find(strcmp(u_aupr_r,nameLabels2));
        if isempty(idx_m)
            possClass=nameLabels1;
        elseif isempty(idx_n)
            possClass=nameLabels2;
        else
            if idx_n<idx_m
                possClass=nameLabels1;
            elseif idx_m<idx_n
                possClass=nameLabels2;
            end
        end
    end
    
function aupr = aupr_evaluation(samp_lab,scores,positiveLabel) 

[rec,prec,~,~] = perfcurve(samp_lab,scores,positiveLabel,'xCrit','reca','yCrit','prec');

% rec is the recall, prec is the precision.
% the first value of rec (at recall 0) is NaN (by definition, PRECISION = TP / (FP + TP))
% at recall 0 we have PRECISION = 0/(0 + 0) (we don't have neither TP nor FP)
% if at the first point of recall (prec(2)) the precision is 1, the NaN is changed
% for 1, in the contrary case (in the first point we have a FP), the NaN is changed for 0
if prec(2) == 1
    prec(1) = 1;
else
    prec(1) = 0;
end

aupr = trapz(rec,prec);
   

function scatter_plot(s,col,sample_names,dis)

n1=size(s,1); 

hold on

plot(s(:,1),s(:,2),'o',...
                'LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',col,...
                'MarkerSize',10); 
if strcmp(dis,'yes')
    for i=1:n1
        text('Interpreter','latex',...
            'String',[' ' sample_names(i,:)],... 
            'Position',[s(i,1) s(i,2)],...
            'FontSize',14,...
            'FontName','calibri',...
            'color','k')
    end
end


function [edges,nodes,pc_corr,x1,cutoff_f]=C_corr(x,V,feat_names,cutoff)
%x: normalized data matrix (NxM), that was given as an input of PCA. 
%The N rows contain the samples, while the M columns contain the features.
%V: loadings of one PCA dimension (i.e. one column of the third output matrix of SVD function)
%feat_names: the feature names  

if nargin<3, error('number of inputs must be at least three'),end
if nargin<4, cutoff=0.7; end


V=sign(V).*log10(1+abs(V)/mean(abs(V)));
V=sign(V).*((abs(V)-min(abs(V)))/(max(abs(V))-min(abs(V))));

%PC-corr on all the features
x_temp=x;
V_temp=V;
feat_names_temp=feat_names;

index=abs(V)>cutoff;
%index=index';
n=length(V)-sum(index);
if n>0
    if n==1
        display(strcat(num2str(n),' feature was deleted because |V_feature|<cutoff'));
    elseif n>1
        display(strcat(num2str(n),' features were deleted because |V_feature|<cutoff'));
    end
    x=x(:,index);
    V=V(index);
    feat_names=feat_names(index);
end

if length(V)>1
    % - Normal Pearson Corrleation calculation
    c=corr(x);
    % - Apply the C-Corr formula
    pc_corr=zeros(length(V));
    for i=1:length(V)-1
        for j=i+1:length(V)
            pc_corr(i,j)=sign(c(i,j))*min(abs([c(i,j),V(i),V(j)]));
        end
    end
else
    pc_corr=zeros(length(V));
end



if (max(abs(pc_corr(:)))<=cutoff) || (length(V)<=1)
    c_temp=corr(x_temp);
    pc_corr_temp=zeros(length(V_temp));
    for i=1:length(V_temp)-1
        for j=i+1:length(V_temp)
            pc_corr_temp(i,j)=sign(c_temp(i,j))*min(abs([c_temp(i,j),V_temp(i),V_temp(j)]));
        end
    end
    
    warning('With this cut-off , there are no edges that have |PC_corr(i,j)|>cutoff');

    
    max_cutoff = max(abs(pc_corr_temp(:)));
    diground = @(x,d) round(x*10^d)/10^d; 
    max_cutoff_d2= max_cutoff - diground(max_cutoff,2);
    if (sign(max_cutoff_d2) == 1)|| (sign(max_cutoff_d2) == 0)
        max_possibl_cutoff =diground(max_cutoff,2);
    else
        max_possibl_cutoff =diground(max_cutoff,2)-0.01;
    end
    
    warn = strcat('Try with another cut-off lower than ',{' '},num2str(max_possibl_cutoff,'%0.2f'));
    warning(warn{1});
    
    flag = 0;
    while flag == 0
        fprintf('\nThen select another cut-off for generating the PC-corr network [number between 0 and %0.2f]:\n\n',max_possibl_cutoff);
        cutoff = input('-> ');
        if cutoff>= 0 && cutoff<= max_possibl_cutoff
            flag = 1;
        else
            fprintf('Please introduce a correct cut-off (in [0,%0.2f]).\n',max_possibl_cutoff);
        end
    end
    
    if flag ==1
        index=abs(V_temp)>cutoff;
        %index=index';
        n=length(V_temp)-sum(index);
        if n>0
            if n==1
                display(strcat(num2str(n),' feature was deleted because |V_feature|<cutoff'));
            elseif n>1
                display(strcat(num2str(n),' features were deleted because |V_feature|<cutoff'));
            end
            x=x_temp(:,index);
            V=V_temp(index);
            feat_names=feat_names_temp(index);
        end
        
        
        % - Normal Pearson Corrleation calculation
        c=corr(x);
        
        % - Apply the C-Corr formula
        pc_corr=zeros(length(V));
        for i=1:length(V)-1
            for j=i+1:length(V)
                pc_corr(i,j)=sign(c(i,j))*min(abs([c(i,j),V(i),V(j)]));
            end
        end
    end

end


idx = pc_corr<-cutoff | pc_corr>cutoff;
pc_corr(idx==0)=0;

%Remove the singletons
pc_corr=pc_corr+pc_corr'; %symmetric
ma=sum(pc_corr==0,1)==size(pc_corr,1);
pc_corr(:,ma)=[];
pc_corr(ma,:)=[];
sig_PC_corr=pc_corr;
sig_PC_corr_Name=feat_names;
sig_PC_corr_Name(ma)=[];
%Datamatrix with only the features in the network.
x1=x;
x1(:,ma)=[];


cutoff_f=cutoff;
%Create the edges and nodes variables
%Edges
edges(1,:)={'Node i', 'Node j','PC-corr(i,j)'};
k=2;
for i=1:length(sig_PC_corr)-1
    for j=i+1:length(sig_PC_corr)
        if sig_PC_corr(i,j)~=0
            
            edges(k,1)=sig_PC_corr_Name(i);
            edges(k,2)=sig_PC_corr_Name(j);
            edges(k,3)=num2cell(sig_PC_corr(i,j));

        k=k+1;
        end    
    end
end

%Nodes
nodes(1,:)={'Node','Loading (V)'};
[~,Lcob]=ismember(sig_PC_corr_Name,feat_names);
for i=1:length(sig_PC_corr_Name)
    nodes(i+1,1)=sig_PC_corr_Name(i);
    nodes(i+1,2)=num2cell(V(Lcob(i)));    
end

function [NodeColor,n1_f,n2_f]=match_V_samplCol(col,x1,labels, nameLabels, Nodes)

n=find(strcmp(col,'r')==1);
m=find(strcmp(col,'k')==1);

nodes_feat=Nodes(2:end,:);
NodeColor=cell(size(nodes_feat,1),1);


for i=1:size(x1,2)
    men{i,1}=nodes_feat{i,1};
    men{i,2}=mean(x1(ismember(labels,nameLabels{n}),i)); %red group
    men{i,3}=mean(x1(ismember(labels,nameLabels{m}),i)); %black group
end
men_dif=zeros(size(x1,2),1);
men_dif=cell2mat(men(:,2))-cell2mat(men(:,3));
m1=men_dif>=0;
m2=men_dif<=0;


b=cell2mat(nodes_feat(:,2))<=0; % 1 where it is a negative loading
l1=sum(b==0);
l2=sum(b==1);
t1=sum(m1'.*(~b)')/sum(b==0);
t2=sum(m2'.*b')/sum(b==1); 
t3=sum(m1'.*b')/sum(b==1);
t4=sum(m2'.*(~b)')/sum(b==0);


if (l1==0)&&(t2>=t3)
    t=t2;
    NodeColor(b==1)={'Black'};
    n1_f=NaN;
    n2_f=1-t;
elseif (l1==0)&&(t3>=t2)
    t=t3;
    NodeColor(b==1)={'Red'};
    n1_f=NaN;
    n2_f=1-t;
elseif (l2==0)&&(t1>=t4)
    t=t1;
    NodeColor(b==0)={'Red'};
    n1_f=1-t;
    n2_f=NaN;
elseif  (l2==0)&&(t4>=t1)
    t=t4;
    NodeColor(b==0)={'Black'};
    n1_f=1-t;
    n2_f=NaN;
elseif (l2~=0)&&(l1~=0)&&(t1>t3)
    t=t1;
    NodeColor(b==0)={'Red'};
    NodeColor(b==1)={'Black'};
    n1_f=1-t;
    n2_f=1-t2;
elseif (l2~=0)&&(l1~=0)&&(t3>t1)
    t=t3;
    NodeColor(b==0)={'Black'};
    NodeColor(b==1)={'Red'};
    n1_f=1-t;
    n2_f=1-t4;
end

if (t1==0)&&(t1==t3)
    if sum(m2'.*b')>=sum(m2'.*(~b)')
        NodeColor(b==0)={'Red'};
        NodeColor(b==1)={'Black'};
        n1_f=1-t1;
        n2_f=1-t2;
    elseif sum(m2'.*(~b)')>sum(m2'.*b')
        NodeColor(b==0)={'Black'};
        NodeColor(b==1)={'Red'};
        n1_f=1-t3;
        n2_f=1-t4;
    end
end

if (t2==0)&&(t2==t4)
    if sum(m1'.*(~b)')>=sum(m1'.*b')
    NodeColor(b==0)={'Red'};
    NodeColor(b==1)={'Black'};
    n1_f=1-t1;
    n2_f=1-t2;
    elseif sum(m1'.*b')>sum(m1'.*(~b)')
        NodeColor(b==0)={'Black'};
        NodeColor(b==1)={'Red'};
        n1_f=1-t3;
        n2_f=1-t4;
    end
end




function plot_graph(pc_corr1,nodes1,cutoff)

G=sparse(triu(pc_corr1))
g = graph(G, nodes1(2:end, 1), 'upper');

num_nodes = length(g.Nodes.Name);
num_labels = length(nodes1(2:end, 1));
num_edges = length(g.Edges.Weight);
if num_nodes ~= num_labels
    error('Number of nodes in the graph does not match the number of labels in nodes1.');
end

% Set node labels and properties
for i = 1:(length(nodes1) - 1)
    len_label(i) = length(nodes1{i + 1, 1});
end
mean_len_label = mean(len_label);
if mean_len_label < 20
    size_lab_nodes = 20;
elseif 20 <= mean_len_label && mean_len_label < 30
    size_lab_nodes = 15;
elseif 30 <= mean_len_label && mean_len_label < 50
    size_lab_nodes = 10;
else
    size_lab_nodes = 8;
end

% Assign properties to nodes
node_colors = cell(1, length(g.Nodes.Name));
node_line_colors = cell(1, length(g.Nodes.Name));
for i = 1:length(g.Nodes.Name)
    node_colors{i} = [0 0 0];  % Default color (Black)
    node_line_colors{i} = [0 0 0];  % Default line color (Black)
    
    if strcmp(nodes1{i + 1, 2}, 'Red')
        node_colors{i} = [1 0 0];  % Red color
        node_line_colors{i} = [1 0 0];  % Red line color
    end
end

% Set node properties in graph
g.Nodes.Shape = repmat({'circle'}, length(g.Nodes.Name), 1);
g.Nodes.Size = repmat({[10, 10]}, length(g.Nodes.Name), 1);
g.Nodes.FontSize = size_lab_nodes * ones(length(g.Nodes.Name), 1);
g.Nodes.TextColor = repmat({[0.5 0.5 0.57]}, length(g.Nodes.Name), 1);

node_colors = cell(num_nodes, 1);

% Assign colors based on the 'Colour' column in nodes1
for i = 1:num_nodes
    if strcmp(nodes1{i + 1, 2}, 'Red')
        node_colors{i} = [1 0 0];  % Red color
    else
        node_colors{i} = [0 0 0];  % Black color
    end
end

g.Nodes.Color = node_colors;

% Edge colors based on weight
edge_weights = g.Edges.Weight;
max_edg_w = max(edge_weights);
min_edg_w = min(edge_weights);

% Assign edge colors
edge_line_colors = zeros(length(g.Edges.Weight), 3);
for i = 1:length(g.Edges.Weight)
    if sign(g.Edges.Weight(i)) == 1
        color_edg = 1 - g.Edges.Weight(i) / max_edg_w;
        edge_line_colors(i, :) = [1 color_edg color_edg];
    else
        color_edg = g.Edges.Weight(i) / abs(min_edg_w) + 1;
        edge_line_colors(i, :) = [color_edg color_edg color_edg];
    end
end
g.Edges.LineColor = edge_line_colors;

% Edge under frustration
edg_frust = 0;
for i = 1:length(g.Nodes.Name) - 1
    for j = i + 1:length(g.Nodes.Name)
        eh1 = findedge(g, i, j);
        if eh1 > 0
            n = strcmp(nodes1{i + 1, 2}, 'Red');
            m = strcmp(nodes1{j + 1, 2}, 'Red');
            if (n * m == 1) && (sign(g.Edges.Weight(eh1)) == -1)
                g.Edges.LineColor(eh1, :) = [0.9 0.9 0.9];
                edg_frust = edg_frust + 1;
            elseif (n + m == 0) && (sign(g.Edges.Weight(eh1)) == -1)
                g.Edges.LineColor(eh1, :) = [0.9 0.9 0.9];
                edg_frust = edg_frust + 1;
            elseif (n + m == 1) && (sign(g.Edges.Weight(eh1)) == 1)
                g.Edges.LineColor(eh1, :) = [0.9 0.9 0.9];
                edg_frust = edg_frust + 1;
            end
        end
    end
end

% For visualization
figure;
h = plot(g, 'NodeLabel', g.Nodes.Name, 'EdgeLabel', g.Edges.Weight);
h.NodeColor = [0 0 0];  % Default color (Black)
h.EdgeColor = g.Edges.LineColor;
h.NodeFontSize = size_lab_nodes;
h.NodeLabelColor = [0.5 0.5 0.57];  % Text color


%Replace the grey edges (edges under frustration) in the figure with dashed grey edges 
fr_edg_frust=edg_frust/num_edges;
axesHandle = gca;
plotHandle1 = findobj(axesHandle,'Type','line','Color',[1 0 0]);
if isempty(plotHandle1)
    plotHandle1=findobj(axesHandle,'Type','line','Color',[1 min_red min_red]);
end
plotHandle2 = findobj(axesHandle,'Type','line','Color','k','LineWidth',0.5);
if isempty(plotHandle2)
    plotHandle2=findobj(axesHandle,'Type','line','Color',[min_black min_black min_black]);
end
plotHandle3 = findobj(axesHandle,'Type','Patch','FaceColor','r','EdgeColor','r','LineStyle','-');
plotHandle4 = findobj(axesHandle,'Type','Patch','FaceColor','k','EdgeColor','k','LineStyle','-');
plotHandle = findobj(axesHandle,'Type','line','Color',[.9 .9 .9]);
for i=1:length(plotHandle)
    set(plotHandle(i),'Color',[.6 .6 .6],'Linestyle','--');
end
diground = @(x,d) round(x*10^d)/10^d; 
AxesH = axes('Parent', f, ...
  'Units', 'normalized', ...
  'Position', [0, 0, 1, 1], ...
  'Visible', 'off', ...
  'XLim', [0, 1], ...
  'YLim', [0, 1], ...
  'NextPlot', 'add');
frac_txt=['$\displaystyle\frac{',num2str(edg_frust),'}{',num2str(length(g.edges)),'}$'];
txt_figure={[],['   cut-off = ',num2str(cutoff)],[],['   frustration = ',frac_txt,' = ',num2str(diground(fr_edg_frust,3)*100),'\%']};

TextH = text(0,1, txt_figure, ...
  'HorizontalAlignment', 'left', ...
  'VerticalAlignment', 'top','Interpreter','latex','FontSize',14);

if isempty(plotHandle1) && isempty(plotHandle2) && ~isempty(plotHandle3) && ~isempty(plotHandle4) && ~isempty(plotHandle)
    h=legend([ plotHandle(1) plotHandle3(1) plotHandle4(1)],{'Edge under frustation', '$\uparrow$ Red sample group', '$\uparrow$ Black sample group'});
elseif isempty(plotHandle1) && isempty(plotHandle2) && isempty(plotHandle3) && ~isempty(plotHandle4) && ~isempty(plotHandle)
    h=legend([ plotHandle(1) plotHandle4(1)],{'Edge under frustation', '$\uparrow$ Black sample group'});
elseif isempty(plotHandle1) && isempty(plotHandle2) && ~isempty(plotHandle3) && isempty(plotHandle4) && ~isempty(plotHandle)
    h=legend([ plotHandle(1) plotHandle3(1)],{'Edge under frustation', '$\uparrow$ Red sample group'});    
elseif isempty(plotHandle4)&& isempty(plotHandle2) && ~isempty(plotHandle) && ~isempty(plotHandle1) && ~isempty(plotHandle3)
    h=legend([plotHandle1(1) plotHandle(1) plotHandle3(1)],{'PC-corr $> 0$', 'Edge under frustation', '$\uparrow$ Red sample group'});
elseif isempty(plotHandle1)&& isempty(plotHandle3) && ~isempty(plotHandle) && ~isempty(plotHandle2) && ~isempty(plotHandle4)
    h=legend([plotHandle2(1) plotHandle(1) plotHandle4(1)],{'PC-corr $< 0$', 'Edge under frustation', '$\uparrow$ Black sample group'});
elseif isempty(plotHandle2)&& ~isempty(plotHandle3) && ~isempty(plotHandle) && ~isempty(plotHandle1) && ~isempty(plotHandle4)
    h=legend([plotHandle1(1) plotHandle(1) plotHandle3(1) plotHandle4(1)],{'PC-corr $> 0$', 'Edge under frustation', '$\uparrow$ Red sample group', '$\uparrow$ Black sample group'});
elseif isempty(plotHandle1)&& ~isempty(plotHandle3) && ~isempty(plotHandle) && ~isempty(plotHandle2) && ~isempty(plotHandle4)
    h=legend([plotHandle2(1) plotHandle(1) plotHandle3(1) plotHandle4(1)],{'PC-corr $< 0$', 'Edge under frustation', '$\uparrow$ Red sample group', '$\uparrow$ Black sample group'});
elseif isempty(plotHandle3)&& isempty(plotHandle2) && ~isempty(plotHandle1) && ~isempty(plotHandle4)
    if isempty(plotHandle)
        h=legend([plotHandle1(1) plotHandle4(1)],{'PC-corr $> 0$','$\uparrow$ Black sample group'}); %only black nodes
    else
        h=legend([plotHandle1(1) plotHandle(1) plotHandle4(1)],{'PC-corr $> 0$', 'Edge under frustation','$\uparrow$ Black sample group'}); %only black nodes
    end
elseif isempty(plotHandle4)&& isempty(plotHandle2)&& ~isempty(plotHandle1) && ~isempty(plotHandle3)
    if isempty(plotHandle)
        h=legend([plotHandle1(1) plotHandle3(1)],{'PC-corr $> 0$','$\uparrow$ Red sample group'}); %only red nodes
    else
        h=legend([plotHandle1(1) plotHandle(1) plotHandle3(1)],{'PC-corr $> 0$','Edge under frustation','$\uparrow$ Red sample group'}); %only red nodes
    end
elseif isempty(plotHandle2) && ~isempty(plotHandle1) && ~isempty(plotHandle3) && ~isempty(plotHandle4)
    if isempty(plotHandle)
        h=legend([plotHandle1(1) plotHandle3(1) plotHandle4(1)],{'PC-corr $> 0$','$\uparrow$ Red sample group','$\uparrow$ Black sample group'}); %only positive PC-corr
    else
        h=legend([plotHandle1(1) plotHandle(1) plotHandle3(1) plotHandle4(1)],{'PC-corr $> 0$','Edge under frustation','$\uparrow$ Red sample group','$\uparrow$ Black sample group'}); %only positive PC-corr
    end
elseif isempty(plotHandle1) && ~isempty(plotHandle2) && ~isempty(plotHandle3) && ~isempty(plotHandle4)
    if isempty(plotHandle)
        h=legend([plotHandle2(1) plotHandle3(1) plotHandle4(1)],{'PC-corr $< 0$','$\uparrow$ Red sample group','$\uparrow$ Black sample group'}); %only negative PC-corr
    else
        h=legend([plotHandle2(1) plotHandle(1) plotHandle3(1) plotHandle4(1)],{'PC-corr $< 0$','Edge under frustation','$\uparrow$ Red sample group','$\uparrow$ Black sample group'}); %only negative PC-corr
    end
elseif ~isempty(plotHandle1)&& ~isempty(plotHandle2) && ~isempty(plotHandle3) && ~isempty(plotHandle4)
    if isempty(plotHandle)
        h= legend([plotHandle1(1) plotHandle2(1) plotHandle3(1) plotHandle4(1)],{'PC-corr $> 0$','PC-corr $< 0$','$\uparrow$ Red sample group','$\uparrow$ Black sample group'}); %no edges under frustation
    else
        h=legend([plotHandle1(1) plotHandle2(1) plotHandle(1) plotHandle3(1) plotHandle4(1)],{'PC-corr $> 0$','PC-corr $< 0$','Edge under frustation', '$\uparrow$ Red sample group','$\uparrow$ Black sample group'});
    end
end
set(h,'Interpreter','latex','FontSize',11)


% elseif ~isempty(plotHandle1) && ~isempty(plotHandle2) && ~isempty(plotHandle) && (sum(plotHandle2==plotHandle)==length(plotHandle)) && ~isempty(plotHandle3) && ~isempty(plotHandle4)
%     h=legend([plotHandle1(1) plotHandle(1) plotHandle3(1) plotHandle4(1)],{'PC-corr $> 0$','Edge under frustation', '$\uparrow$ Red sample group','$\uparrow$ Black sample group'});
% elseif ~isempty(plotHandle1) && ~isempty(plotHandle2) && ~isempty(plotHandle) && (sum(plotHandle1==plotHandle)==length(plotHandle)) && ~isempty(plotHandle3) && ~isempty(plotHandle4)
%     h=legend([plotHandle2(1) plotHandle(1) plotHandle3(1) plotHandle4(1)],{'PC-corr $< 0$','Edge under frustation', '$\uparrow$ Red sample group','$\uparrow$ Black sample group'});


function RemoveSheet123(excelFileName,sheetName)
% RemoveSheet123 - removes the sheets that are automatically added to excel
% file. 
% When Matlab writes data to a new Excel file, the Excel software
% automatically creates 3 sheets (the names are depended on the user
% languade). This appears even if the user has defined the sheet name to be
% added. 
%
% Usage:
% RemoveSheet123(excelFileName) - remove "sheet1", "sheet2","sheet3" from
% the excel file. excelFileName is a string of the Excel file name.
% RemoveSheet123(excelFileName,sheetName) - enables the user to enter the
% sheet name when the language is other than English.
% sheetName is the default sheet name, without the number.
%
%
%                       Written by Noam Greenboim
%                       www.perigee.co.il
%


%% check input arguments
if nargin < 1 || isempty(excelFileName)
    error('Filename must be specified.');
end

if ~ischar(excelFileName)
    error('Filename must be a string.');
end

try
    excelFileName = validpath(excelFileName);
catch 
    error('File not found.');
end


if nargin < 2
    sheetName = 'Sheet'; % EN: Sheet, DE: Tabelle, HE: âéìéåï , etc. (Lang. dependent)
else
    if ~ischar(sheetName)
        error('Default sheet name must be a string.');
    end
end

%%
% Open Excel file.
objExcel = actxserver('Excel.Application');
objExcel.Workbooks.Open(excelFileName); % Full path is necessary!
[~, sheets] = xlsfinfo(excelFileName);
n=sum(strcmp(sheets,'Sheet1')+strcmp(sheets,'Sheet2')+strcmp(sheets,'Sheet3'));
% Delete sheets.
try
    if n==1
      objExcel.ActiveWorkbook.Worksheets.Item([sheetName '1']).Delete;
    elseif n==2
        objExcel.ActiveWorkbook.Worksheets.Item([sheetName '1']).Delete;
        objExcel.ActiveWorkbook.Worksheets.Item([sheetName '2']).Delete;
    elseif n==3
        objExcel.ActiveWorkbook.Worksheets.Item([sheetName '1']).Delete;
        objExcel.ActiveWorkbook.Worksheets.Item([sheetName '2']).Delete;
        objExcel.ActiveWorkbook.Worksheets.Item([sheetName '3']).Delete;
    end

catch
    fprintf('\n')
    O=objExcel.ActiveWorkbook.Worksheets.get;
    if O.Count==1
        error('Can''t delete the last sheet. Excel file must containt at least one sheet.')
    else
      warning('Problem occured. Check excel file.'); 
    end
end
% Save, close and clean up.
objExcel.ActiveWorkbook.Save;
objExcel.ActiveWorkbook.Close;
objExcel.Quit;
objExcel.delete;

function filenameOut = validpath(filename)
    % VALIDPATH builds a full path from a partial path specification
    %   FILENAME = VALIDPATH(FILENAME) returns a string vector containing full
    %   path to a file. FILENAME is string vector containing a partial path
    %   ending in a file or directory name. May contain ..\  or ../ or \\. The
    %   current directory (pwd) is prepended to create a full path if
    %   necessary. On UNIX, when the path starts with a tilde, '~', then the
    %   current directory is not prepended.
    %
    %   See also XLSREAD, XLSWRITE, XLSFINFO.
    
    %   Copyright 1984-2012 The MathWorks, Inc.
    
    %First check for wild cards, since that is not supported.
    if strfind(filename, '*') > 0
        error(message('MATLAB:xlsread:Wildcard', filename));
    end
    
    % break partial path in to file path parts.
    [Directory, file, ext] = fileparts(filename);

    if ~isempty(ext)
        filenameOut = getFullName(filename);
    else
        extIn = matlab.io.internal.xlsreadSupportedExtensions;
        for i=1:length(extIn)
            try                                                                %#ok<TRYNC>
                filenameOut = getFullName(fullfile(Directory, [file, extIn{i}]));
                return;
            end
        end
        error(message('MATLAB:xlsread:FileDoesNotExist', filename));    
    end

function absolutepath=abspath(partialpath)
    
    % parse partial path into path parts
    [pathname, filename, ext] = fileparts(partialpath);
    % no path qualification is present in partial path; assume parent is pwd, except
    % when path string starts with '~' or is identical to '~'.
    if isempty(pathname) && strncmp('~', partialpath, 1)
        Directory = pwd;
    elseif isempty(regexp(partialpath,'(.:|\\\\)', 'once')) && ...
            ~strncmp('/', partialpath, 1) && ...
            ~strncmp('~', partialpath, 1);
        % path did not start with any of drive name, UNC path or '~'.
        Directory = [pwd,filesep,pathname];
    else
        % path content present in partial path; assume relative to current directory,
        % or absolute.
        Directory = pathname;
    end
    
    % construct absolute filename
    absolutepath = fullfile(Directory,[filename,ext]);
    
function filename = getFullName(filename)
    FileOnPath = which(filename);
    if isempty(FileOnPath)
        % construct full path to source file
        filename = abspath(filename);
        if isempty(dir(filename)) && ~isdir(filename)
            % file does not exist. Terminate importation of file.
            error(message('MATLAB:xlsread:FileDoesNotExist', filename));
        end
    else
        filename = FileOnPath;
    end
    
    function labels_rp=random_permutation_labels(numb_rand,labels)
        
        % This function returns numb_rand random permutation of the labels
        % INPUT
        %   numb_rand => (double number) number of permutations 
        %   labels => (cell array M x 1)labels of the samples 
        
        % OUTPUT
        %   labels_rp => (cell array M x numb_rand) numb_rand random permutation of the sample labels
        
        labels_rp=cell(size(labels,1),numb_rand);
        
        for pr=1:numb_rand
            
            labels_rp(:,pr)=labels(randperm(length(labels)));
        end
    
    
    function p_Value=significance_segragation(numb_rand,PCA,dim,labels_rp,nameLabels,numbLabels,segr_meas,type,u_aupr,u_aupr_r)
        
        % This function evaluates the significane of the segragation
        % measure (P-value, AUC or AUPR). Is the dimension discriminative
        % because it is random or capture the main variability?
        
        % INPUT
        %   numb_rand => number of permutation
        %   PCA=> PCA results either centred or not centred, with a certain
        %   normalization or not
        %   dim => dimension of sample segragation
        %   labels_rp => numb_rand random permutation of the sample labels
        %   nameLabels => unique name of labels
        %   numbLabels => number of unique sample labels
        %   segr_meas => segragation measure on the
        %   input dimension dim:
        %   - for 'p','auc','aupr', just one value: either segragation
        %   measure according to P-value, AUC or AUPR
        %   - for 'all': segragation measures according to P-value (in pos. 1), 
        %   AUC (in pos.2)and AUPR (in pos. 3)
        %   type => type of segragation measure, depending on the input
        %   segr_meas: 'p' (for MW p-value), 'auc' (for AUC), 'aupr' (for
        %   AUPR) or 'all' (if you want to compute for all: MW p-value, AUC
        %   and AUPR)
        
        
        % OUTPUT
        %   p_Value => significance (p-value) of the input segragation
        %   measure:
        %   - a number if computed for just one input segragation measure
        %   (type= 'p', 'auc' or 'aupr')
        %   - a numeric vector if computed for just one input segragation measure
        %   (type= 'all'): in idx 1 for P-value, in 2 for AUC and in 3 for
        %   AUPR
        
        rand_segr_meas=NaN(nchoosek(numbLabels,2),numb_rand);
        if strcmp(type,'all')
            rand_segr_meas_mw=NaN(nchoosek(numbLabels,2),numb_rand);
            rand_segr_meas_auc=NaN(nchoosek(numbLabels,2),numb_rand);
            rand_segr_meas_aupr=NaN(nchoosek(numbLabels,2),numb_rand);
        end
        for pr=1:numb_rand
            
            labels=labels_rp(:,pr);
            
            n = 1;
            m = 2;
            
            
            for j=1:nchoosek(numbLabels,2) %two-group comparison
                
                if strcmp(type,'p')
                    % Compute p-value of Mann-Whitney test
                    rand_segr_meas(j,pr) = ranksum(PCA(ismember(labels,nameLabels{n}),dim),PCA(ismember(labels,nameLabels{m}),dim));
                elseif strcmp(type,'auc')
                    samp_lab = [labels(ismember(labels,nameLabels{n}));labels(ismember(labels,nameLabels{m}))];
                    scores = [PCA(ismember(labels,nameLabels{n}),dim);PCA(ismember(labels,nameLabels{m}),dim)];
                   
                    possClass=positive_label_opt(u_aupr, nameLabels{n},nameLabels{m},labels,u_aupr_r);

                    % Compute AUC
                    [~,~,~,rand_segr_meas(j,pr)] = perfcurve(samp_lab,scores,possClass);
                    
                    if rand_segr_meas(j,pr)<0.5
                        rand_segr_meas(j,pr) = 1-rand_segr_meas(j,pr);
                    end
                    
                elseif strcmp(type,'aupr')
                    
                    samp_lab = [labels(ismember(labels,nameLabels{n}));labels(ismember(labels,nameLabels{m}))];
                    scores = [PCA(ismember(labels,nameLabels{n}),dim);PCA(ismember(labels,nameLabels{m}),dim)];
                    
                    possClass=positive_label_opt(u_aupr, nameLabels{n},nameLabels{m},labels,u_aupr_r);

                    % Compute AUC
                    [~,~,~,auc] = perfcurve(samp_lab,scores,possClass);
                    if auc<0.5
                        % Compute AUPR
                        flip_scores=2*mean(scores)-scores;
                        rand_segr_meas(j,pr)= aupr_evaluation(samp_lab,flip_scores,possClass);
                    else
                        % Compute AUPR
                        rand_segr_meas(j,pr) = aupr_evaluation(samp_lab,scores,possClass);
                    end
                    
                elseif strcmp(type,'all')
                    
                    % Compute p-value of Mann-Whitney test
                    rand_segr_meas_mw(j,pr) = ranksum(PCA(ismember(labels,nameLabels{n}),dim),PCA(ismember(labels,nameLabels{m}),dim));
                    % Compute AUC
                    samp_lab = [labels(ismember(labels,nameLabels{n}));labels(ismember(labels,nameLabels{m}))];
                    scores = [PCA(ismember(labels,nameLabels{n}),dim);PCA(ismember(labels,nameLabels{m}),dim)];

                    possClass=positive_label_opt(u_aupr, nameLabels{n},nameLabels{m},labels,u_aupr_r);

                    [~,~,~,rand_segr_meas_auc(j,pr)] = perfcurve(samp_lab,scores,possClass);
                    
                    if rand_segr_meas_auc(j,pr)<0.5
                        rand_segr_meas_auc(j,pr) = 1-rand_segr_meas_auc(j,pr);
                        % Compute AUPR
                        flip_scores=2*mean(scores)-scores;
                        rand_segr_meas_aupr(j,pr)= aupr_evaluation(samp_lab,flip_scores,possClass);
                    else
                        % Compute AUPR
                        rand_segr_meas_aupr(j,pr) = aupr_evaluation(samp_lab,scores,possClass);
                    end
                
                end
                
                m = m + 1;
                if(m > numbLabels)
                    n = n + 1;
                    m = n + 1;
                end
            end
        end
        
        if strcmp(type, 'p')|| strcmp(type,'auc') ||strcmp(type,'aupr')
            Rand_segr_meas=mean(rand_segr_meas,1);
        elseif strcmp(type,'all')
            Rand_segr_meas_mw=mean(rand_segr_meas_mw,1);
            Rand_segr_meas_auc=mean(rand_segr_meas_auc,1);
            Rand_segr_meas_aupr=mean(rand_segr_meas_aupr,1);
        end
        % Calculation of the significance (p-value) of the input segragation measure
        if strcmp(type,'p')
            p_Value =(sum(Rand_segr_meas<segr_meas)+1)/(numb_rand+1);
        elseif strcmp(type,'auc')||strcmp(type,'aupr')
            p_Value=(sum(Rand_segr_meas>segr_meas)+1)/(numb_rand+1);
        elseif strcmp(type,'all')
            p_Value(1) =(sum(Rand_segr_meas_mw<segr_meas(1))+1)/(numb_rand+1); % for p-value
            p_Value(2)=(sum(Rand_segr_meas_auc>segr_meas(2))+1)/(numb_rand+1); % for AUC
            p_Value(3)=(sum(Rand_segr_meas_aupr>segr_meas(3))+1)/(numb_rand+1); % for AUPR

        end

function textprogressbar(c)
% This function creates a text progress bar. It should be called with a 
% STRING argument to initialize and terminate. Otherwise the number correspoding 
% to progress in % should be supplied.
% INPUTS:   C   Either: Text string to initialize or terminate 
%                       Percentage number to show progress 
% OUTPUTS:  N/A
% Example:  Please refer to demo_textprogressbar.m

% Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
% Version: 1.0
% Changes tracker:  29.06.2010  - First version

% Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/

%% Initialization
persistent strCR;           %   Carriage return pesistent variable

% Vizualization parameters
strPercentageLength = 10;   %   Length of percentage string (must be >5)
strDotsMaximum      = 10;   %   The total number of dots in a progress bar

%% Main 

if isempty(strCR) && ~ischar(c),
    % Progress bar must be initialized with a string
    error('The text progress must be initialized with a string');
elseif isempty(strCR) && ischar(c),
    % Progress bar - initialization
    fprintf('%s',c);
    strCR = -1;
elseif ~isempty(strCR) && ischar(c),
    % Progress bar  - termination
    strCR = [];  
    fprintf([c '\n']);
elseif isnumeric(c)
    % Progress bar - normal progress
    c = floor(c);
    percentageOut = [num2str(c) '%%'];
    percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
    nDots = floor(c/100*strDotsMaximum);
    dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
    strOut = [percentageOut dotOut];
    
    % Print it on the screen
    if strCR == -1,
        % Don't do carriage return during first run
        fprintf(strOut);
    else
        % Do it during all the other runs
        fprintf([strCR strOut]);
    end
    
    % Update carriage return
    strCR = repmat('\b',1,length(strOut)-1);
    
else
    % Any other unexpected input
    error('Unsupported argument type');
end
        
   
        
        