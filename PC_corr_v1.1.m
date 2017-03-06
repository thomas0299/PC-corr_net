%%%%% PC-corr algorithm associate to PCA analysis%%%%%

%%% Released under MIT License
%%% Copyright (c) 2017 Sara Ciucci, Yan Ge, Claudio Durán and Carlo Vittorio Cannistraci

% Please cite:
% Enlightening discriminative network functional modules behind Principal Component Analysis separation in differential-omic science studies.
% Sara Ciucci, Yan Ge, Claudio Durán, Alessandra Palladini, Víctor Jiménez Jiménez, Luisa María Martínez Sánchez, 
% Yuting Wang, Susanne Sales, Andrej Shevchenko, Steven W. Poser, Maik Herbig, Oliver Otto, Andreas Androutsellis-Theotokis, 
% Jochen Guck, Mathias J. Gerl and Carlo Vittorio Cannistraci 
% Scientific Reports, 2017

% INPUT
%   x => (Numeric matrix MxN) Dataset with samples on the rows and features on the columns
%   sample_labels => (Cell array Mx1) Labels of the samples
%   feat_names => (Cell array Nx1) Names of the features
%   sample_names=> (Cell array NX1) Names of the samples
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


function [Edges,Nodes]=PC_corr(x,sample_labels,feat_names, sample_names,dis)

%% initialisation and default options
if nargin<4, error('Not Enough Input Arguments'); end
if nargin<5, dis='no'; end

if sum(sum(isnan(x)))==1
    error('There is %d NaN value in your data matrix. \nPlease replace it.',sum(sum(isnan(x))));
elseif sum(sum(isnan(x)))>1
    error('There are %d NaN values in your data matrix. \nPlease replace them.',sum(sum(isnan(x))));
end

%%

labels = sample_labels; 
nameLabels = unique(labels); %name of the groups
numbLabels = size(nameLabels,1); %number of groups


%Remove features with same identical values across all the samples
x1=x;
x1=x1-repmat(mean(x1),size(x1,1),1);
ma=sum(x1==0,1)==size(x,1);
remov_feat=sum(ma); %number of removed fatures
x(:,ma)=[];
feat_names(ma')=[];


%Normalizations of the dataset
norm{1}=x./repmat(sum(x,1),size(x,1),1); norms{1} = 'DSOR'; %dividing by the sum over the samples
norm{2}=x./repmat(sum(x,2),1,size(x,2)); norms{2} = 'DSOC';%%dividing by the sum over the features
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


flag = 0;

while flag == 0
    
    fprintf('Is your data represented by\n[r]ranked labels(labels that are organized according to a progressive order. e.g. different stages of a disease, where Stage 1 < Stage 2 < Stage 3)\n[c]class labels (labels that are not necessary organized in a progressive order e.g. Condition A, Condition B, Condition C)? [r/c]:\n\n');
    u_lab = input('-> ','s');
    
    flag = 1;

    if ~strcmp(u_lab,'r') && ~strcmp(u_lab,'c')
        flag = 0;
        fprintf('Please introduce either "r" for ranked labels or "c" for class labels\n')
    end
end

if strcmp(u_lab,'r')

    flag = 0;
    while flag == 0
        fprintf('Are the values of your ranked labels\n[d] discrete (Stage 1 < Stage 2 < Stage 3)\nor [con] continuous (different times of development of a cell line)? [d/con]:\n\n');
        u_lab = input('-> ','s');
        
        flag = 1;
        
        if ~strcmp(u_lab,'d') && ~strcmp(u_lab,'con')
            flag = 0;
            fprintf('Please introduce either "d" for discrete labels or "con" for continuous labels\n')
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
    
    
    
    for k=1:size(labels,1) %dimension
        
        if strcmp(u_lab,'c') || strcmp(u_lab,'d') 
            
            n = 1;
            m = 2;
            for j=1:nchoosek(numbLabels,2) %two-group comparison
                
                % Compute p-value of Mann-Whitney test
                mw_ncPCA{i,j,k} = ranksum(ncPCA{i}(ismember(labels,nameLabels{n}),k),ncPCA{i}(ismember(labels,nameLabels{m}),k));
                mw_cPCA{i,j,k} = ranksum(cPCA{i}(ismember(labels,nameLabels{n}),k),cPCA{i}(ismember(labels,nameLabels{m}),k));
                
                
                samp_lab = [sample_labels(ismember(labels,nameLabels{n}));sample_labels(ismember(labels,nameLabels{m}))];
                scores_nc = [ncPCA{i}(ismember(labels,nameLabels{n}),k);ncPCA{i}(ismember(labels,nameLabels{m}),k)];
                scores_c = [cPCA{i}(ismember(labels,nameLabels{n}),k);cPCA{i}(ismember(labels,nameLabels{m}),k)];
                % Compute AUC
                [~,~,~,AUC_nc{i,j,k}] = perfcurve(samp_lab,scores_nc,nameLabels{n});
                [~,~,~,AUC_c{i,j,k}] = perfcurve(samp_lab,scores_c,nameLabels{n});
                AUC_nc{i,j,k} = max([AUC_nc{i,j,k} 1-AUC_nc{i,j,k}]);
                AUC_c{i,j,k} = max([AUC_c{i,j,k} 1-AUC_c{i,j,k}]);
                
                % Compute AUPR
                [~,~,~,AUPR_nc{i,j,k}] = perfcurve(samp_lab,scores_nc,nameLabels{n},'xCrit','reca','yCrit','prec');
                [~,~,~,AUPR_c{i,j,k}] = perfcurve(samp_lab,scores_c,nameLabels{n},'xCrit','reca','yCrit','prec');
                AUPR_nc{i,j,k} = max([AUPR_nc{i,j,k} 1-AUPR_nc{i,j,k}]);
                AUPR_c{i,j,k} = max([AUPR_c{i,j,k} 1-AUPR_c{i,j,k}]);
                
                
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
norms = repmat(norms',[length(labels)*2 1]);
centred = repmat({'yes'},[length(labels)*12 1]);
non_centred = repmat({'no'},[length(labels)*12 1]);
centr = {non_centred{:} centred{:}};
dim = (1:length(labels))';
dim = repmat(dim,[12,1]);
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
    pvals = {mw_ncPCA{:} mw_cPCA{:}};
    AUCs = {AUC_nc{:} AUC_c{:}};
    AUPRs = {AUPR_nc{:} AUPR_c{:}};
    results = horzcat(pvals',AUCs',AUPRs',norms,centr',dim',variance);
   
    elseif strcmp(u_lab,'d')
        
        header = {'P-value','AUC','AUPR','pears','spear','Norm','Centering','Dim','expl Var'};
        pvals = {mw_ncPCA{:} mw_cPCA{:}};
        AUCs = {AUC_nc{:} AUC_c{:}};
        AUPRs = {AUPR_nc{:} AUPR_c{:}};
        pears_corr = {rank_pears_corr_ncPCA{:} rank_pears_corr_cPCA{:}};
        spear_corr = {rank_spear_corr_ncPCA{:} rank_spear_corr_cPCA{:}};
        results = horzcat(pvals',AUCs',AUPRs',pears_corr',spear_corr',norms,centr',dim',variance);
        
    elseif strcmp(u_lab,'con')
        
        header = {'pears','spear','Norm','Centering','Dim','expl Var'};
        pears_corr = {rank_pears_corr_ncPCA{:} rank_pears_corr_cPCA{:}};
        spear_corr = {rank_spear_corr_ncPCA{:} rank_spear_corr_cPCA{:}};
        results = horzcat(pears_corr',spear_corr',norms,centr',dim',variance);
        
    end
    
    
    
    flag = 0;
    
    while flag == 0
        
        if strcmp(u_lab,'c') || strcmp(u_lab,'d')
            
            if strcmp(u_lab,'c')
            
                fprintf('Would you like to rank the PCA results by P-value, AUC or AUPR [p/auc/aupr]:\n\n');
                u_rank = input('-> ','s');           
            else
                fprintf('Would you like to rank the PCA results by P-value, AUC, AUPR, Pearson correlation or Spearman correlation [p/auc/aupr/pc/sc]:\n\n');
                u_rank = input('-> ','s');
            end
            
            flag = 1;
            
            switch u_rank
                case 'p'
                    [~,idx] = sort([results{:,1}]);
                    results = results(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results([results{:,1}]<0.05,:))
                        display(vertcat(header,results([results{:,1}]<0.25,:)))
                    else
                        display(vertcat(header,results([results{:,1}]<0.05,:)))
                    end
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
                    
                    
                case 'auc'
                    [~,idx] = sort([results{:,2}],'descend');
                    results = results(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results([results{:,2}]>=0.7,:))
                        display(vertcat(header,results([results{:,2}]>=0.6,:)))
                    else
                        display(vertcat(header,results([results{:,2}]>=0.7,:)))
                    end
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
                    
                case 'aupr'
                    [~,idx] = sort([results{:,3}],'descend');
                    results = results(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results([results{:,3}]>=0.7,:))
                        display(vertcat(header,results([results{:,3}]>=0.6,:)))
                    else
                        display(vertcat(header,results([results{:,3}]>=0.7,:)))
                    end
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
                    
                case 'pc'
                    if strcmp(u_lab,'c')
                        flag = 0;
                        fprintf('Please introduce either "p", "auc" or "aupr"\n')
                    else
                        [~,idx] = sort(abs([results{:,4}]),'descend');
                        results = results(idx,:);
                        %Only some results are shown on the screen
                        if isempty(results(abs([results{:,4}])>=0.6,:))
                            display(vertcat(header,results(abs([results{:,4}])>=0.5,:)))
                        else
                            display(vertcat(header,results(abs([results{:,4}])>=0.6,:)))
                        end
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
                            display(vertcat(header,results(abs([results{:,5}])>=0.5,:)))
                        else
                            display(vertcat(header,results(abs([results{:,5}])>=0.6,:)))
                        end
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
                        
                    end
                otherwise
                    if strcmp(u_lab,'c')
                        flag = 0;
                        fprintf('Please introduce either "p", "auc" or "aupr"\n')
                    else
                        flag = 0;
                        fprintf('Please introduce either "p", "auc", "aupr", "pc" or "sc"\n')
                    end
            end
            
        else
            
            fprintf('Would you like to rank the PCA results by Pearson or Spearman correlation? [pc/sc]:\n\n');
            u_rank = input('-> ','s');
            
            flag = 1;
            
            switch u_rank
                case 'pc'
                    [~,idx] = sort(abs([results{:,1}]),'descend');
                    results = results(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results(abs([results{:,1}])>=0.6,:))
                        display(vertcat(header,results(abs([results{:,1}])>=0.5,:)))
                    else
                        display(vertcat(header,results(abs([results{:,1}])>=0.6,:)))
                    end
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
                        display(vertcat(header,results(abs([results{:,2}])>=0.5,:)))
                    else
                        display(vertcat(header,results(abs([results{:,2}])>=0.6,:)))
                    end
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
                    fprintf('Please introduce either "pc" or "sc"\n');
            end
            
        end
    end
    
else %more than two groups
    
    if strcmp(u_lab,'c') || strcmp(u_lab,'d')
        
        pval1 = permute(mw_ncPCA,[1 3 2]);
        pval1 = reshape(pval1,12*size(mw_cPCA,3),nchoosek(numbLabels,2));
        pval2 = permute(mw_cPCA,[1 3 2]);
        pval2 = reshape(pval2,12*size(mw_cPCA,3),nchoosek(numbLabels,2));
        pvals = vertcat(pval1,pval2);
        avg = num2cell(mean(cell2mat(pvals),2));
        
        AUC1 = permute(AUC_nc, [1 3 2]);
        AUC1 = reshape(AUC1,12*size(mw_cPCA,3),nchoosek(numbLabels,2));
        AUC2 = permute(AUC_c, [1 3 2]);
        AUC2 = reshape(AUC2,12*size(mw_cPCA,3),nchoosek(numbLabels,2));
        AUCs = vertcat(AUC1,AUC2);
        avg_AUC = num2cell(mean(cell2mat(AUCs),2));
        
        AUPR1 = permute(AUPR_nc, [1 3 2]);
        AUPR1 = reshape(AUPR1,12*size(mw_cPCA,3),nchoosek(numbLabels,2));
        AUPR2 = permute(AUPR_c, [1 3 2]);
        AUPR2 = reshape(AUPR2,12*size(mw_cPCA,3),nchoosek(numbLabels,2));
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
            header_xls = {'Avg P-val',group_head{:},'Avg AUC',group_head{:},'Avg AUPR',group_head{:},'Norm','Centering','Dimension','explained Variance'};
            results_xls = horzcat(avg,pvals,avg_AUC,AUCs,avg_AUPR,AUPRs,norms,centr',dim',variance);
        end
        if strcmp(u_lab,'d')
            pears_corr = {rank_pears_corr_ncPCA{:} rank_pears_corr_cPCA{:}};
            spear_corr = {rank_spear_corr_ncPCA{:} rank_spear_corr_cPCA{:}};
            
            header_xls = {'Avg P-val',group_head{:},'Avg AUC',group_head{:},'Avg AUPR',group_head{:},'pearson-correlation','spearman-correlation','Norm','Centering','Dimension','explained Variance'};
            results_xls = horzcat(avg,pvals,avg_AUC,AUCs,avg_AUPR,AUPRs,pears_corr',spear_corr',norms,centr',dim',variance);
        end
        
        
    else
        pears_corr = {rank_pears_corr_ncPCA{:} rank_pears_corr_cPCA{:}};
        spear_corr = {rank_spear_corr_ncPCA{:} rank_spear_corr_cPCA{:}};
        
        header_xls = {'pearson-correlation','spearman-correlation','Norm','Centering','Dimension','explained Variance'};
        results_xls = horzcat(pears_corr',spear_corr',norms,centr',dim',variance);
        
    end
    
    
    if strcmp(u_lab,'c')  
    
        header = {'Avg P-val','Avg AUC','Avg AUPR','Norm','Centering','Dimension','expl Var'};
        results = horzcat(avg,avg_AUC,avg_AUPR,norms,centr',dim',variance);
        
    elseif strcmp(u_lab,'d')
    
        header = {'Avg Pval','Avg AUC','Avg AUPR','pears','spear','Norm','Centering','Dim','expl Var'};
        results = horzcat(avg,avg_AUC,avg_AUPR,pears_corr',spear_corr',norms,centr',dim',variance);
        
    else
        header = {'pears','spear','Norm','Centering','Dimension','expl Var'};
        results = horzcat(pears_corr',spear_corr',norms,centr',dim',variance);
    end
    
    flag = 0;
    while flag == 0
        
        if strcmp(u_lab,'c') || strcmp(u_lab,'d')
            
            if strcmp(u_lab,'c')
                
                fprintf('Would you like to rank the PCA results by P-value, AUC or AUPR [p/auc/aupr]:\n\n');
                u_rank = input('-> ','s');
                
            else
                fprintf('Would you like to rank the PCA results by P-value, AUC, AUPR, Pearson correlation or Spearman correlation [p/auc/aupr/pc/sc]:\n\n');
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
                        display(vertcat(header,results([results{:,1}]<0.25,:)))
                    else
                        display(vertcat(header,results([results{:,1}]<0.05,:)))
                    end
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
                    
                case 'auc'
                    [~,idx] = sort([results{:,2}],'descend');
                    results = results(idx,:);
                    results_xls=results_xls(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results([results{:,2}]>=0.7,:))
                        display(vertcat(header,results([results{:,2}]>=0.6,:)))
                    else
                        display(vertcat(header,results([results{:,2}]>=0.7,:)))
                    end
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
                    
                case 'aupr'
                    [~,idx] = sort([results{:,3}],'descend');
                    results = results(idx,:);
                    results_xls=results_xls(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results([results{:,3}]>=0.7,:))
                        display(vertcat(header,results([results{:,3}]>=0.6,:)))
                    else
                        display(vertcat(header,results([results{:,3}]>=0.7,:)))
                    end
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
                    
                    
                case 'pc'
                    if strcmp(u_lab,'c')
                        flag = 0;
                        fprintf('Please introduce either "p", "auc" or "aupr"\n')
                    else
                        [~,idx] = sort(abs([results{:,4}]),'descend');
                        results = results(idx,:);
                        results_xls=results_xls(idx,:);
                        %Only some results are shown on the screen
                        if isempty(results(abs([results{:,4}])>=0.6,:))
                            display(vertcat(header,results(abs([results{:,4}])>=0.5,:)))
                        else
                            display(vertcat(header,results(abs([results{:,4}])>=0.6,:)))
                        end
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
                        
                    end
                case 'sc'
                    if strcmp(u_lab,'c')
                        flag = 0;
                        fprintf('Please introduce either "p", "auc" or "aupr"\n')
                    else
                        [~,idx] = sort(abs([results{:,5}]),'descend');
                        results = results(idx,:);
                        results_xls=results_xls(idx,:);
                        %Only some results are shown on the screen
                        if isempty(results(abs([results{:,5}])>=0.6,:))
                            display(vertcat(header,results(abs([results{:,5}])>=0.5,:)))
                        else
                            display(vertcat(header,results(abs([results{:,5}])>=0.6,:)))
                        end
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
                        
                    end
                otherwise
                    if strcmp(u_lab,'c')
                        flag = 0;
                        fprintf('Please introduce either "p", "auc" or "aupr"\n')
                    else
                        flag = 0;
                        fprintf('Please introduce either "p", "auc", "aupr", "pc" or "sc"\n')
                    end
            end
        else
            
            fprintf('Would you like to rank the PCA results by Pearson or Spearman correlation? [pc/sc]:\n\n');
            u_rank = input('-> ','s');
            
            flag = 1;
            
            switch u_rank
                case 'pc'
                    [~,idx] = sort(abs([results{:,1}]),'descend');
                    results = results(idx,:);
                    results_xls=results_xls(idx,:);
                    %Only some results are shown on the screen
                    if isempty(results(abs([results{:,1}])>=0.6,:))
                        display(vertcat(header,results(abs([results{:,1}])>=0.5,:)))
                    else
                        display(vertcat(header,results(abs([results{:,1}])>=0.6,:)))
                    end
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
                        display(vertcat(header,results(abs([results{:,2}])>=0.5,:)))
                    else
                        display(vertcat(header,results(abs([results{:,2}])>=0.6,:)))
                    end
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

%% User interaction
flag = 0;
while flag == 0
    fprintf('Select the normalization:\n\n');
    u_norm_n = input('-> ','s');
    flag = 1;
    switch u_norm_n
        case 'DSOR'
            u_norm = 1;
        case 'DSOC'
            u_norm = 2;
        case 'LOG'
            u_norm = 3;
        case 'ZSCORE'
            u_norm = 4;
        case 'QUANTILE T'
            u_norm = 5;
        case 'QUANTILE'
            u_norm = 6;
        case 'ZSCORE T'
            u_norm = 7;
        case 'PLUS(ABS(MIN))'
            u_norm = 8;
        case 'PARETO SCALING'
            u_norm = 9;
        case 'SQRT'
            u_norm = 10;
        case 'MANORM'
            u_norm = 11;
        case '-'
            u_norm = 12;
        otherwise
            flag = 0;
            fprintf('Please introduce the exact name of the normalization\n')
    end
end

flag = 0;
while flag == 0
    fprintf('Centering version? [y/n]:\n\n');
    u_cent = input('-> ','s');
    if strcmp(u_cent,'n') || strcmp(u_cent,'y')
        flag = 1;
    else
        fprintf('Please introduce just "y" or "n"\n');
    end
end

flag = 0;
while flag == 0
    fprintf('Select the dimension for generating the PC-corr network:\n\n');
    u_dim = input('-> ');
    if u_dim > 0 && u_dim <= length(labels);
        flag = 1;
    else
        fprintf('Please introduce an existing dimension.\n');
    end
end

flag = 0;
while flag == 0
    fprintf('Select a cut-off or a set of cut-offs for generating the PC-corr network [number between 0 and 1]:\n\n');
    fprintf('Example: [0.6 0.65 0.7]\n\n');
    cutoff = input('-> ');
    if max(cutoff)>= 0 && max(cutoff)<= 1
        flag = 1;
    else
        fprintf('Please introduce a correct cut-off or set of cut-offs (in [0,1]).\n');
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
        col{I(i)}=rand(1,3);
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
%% 2 Groups of labels
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

%% More than 2 Groups of labels
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



% bar plots
figure
subplot(2,2,[1 2])
if numbLabels == 2
    if ~strcmp(u_rank,'pc') && ~strcmp(u_rank,'sc')
        bar(cell2mat(mw_PCA(u_norm,:)))
    else
        a1=find(cell2mat(mw_PCA(u_norm,:))>=0);
        b1=find(cell2mat(mw_PCA(u_norm,:))<0);
        for i=1:size(labels,1)
            if ismember(i,a1)
               hold on
               h(i)=bar(i,abs(cell2mat(mw_PCA(u_norm,i))),'FaceColor','k');
            else
               h(i)=bar(i,abs(cell2mat(mw_PCA(u_norm,i))),'FaceColor',[0.5 0.5 0.5]);
            end
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
            bar(permute(mean(cell2mat(mw_PCA(u_norm,:,:))),[3 1 2]))
        else
            a1=find(cell2mat(mw_PCA(u_norm,:))>=0);
            b1=find(cell2mat(mw_PCA(u_norm,:))<0);
            for i=1:size(labels,1)
                if ismember(i,a1)
                    hold on
                    h(i)=bar(i,abs(cell2mat(mw_PCA(u_norm,i))),'FaceColor','k');
                else
                    h(i)=bar(i,abs(cell2mat(mw_PCA(u_norm,i))),'FaceColor',[0.5 0.5 0.5]);
                end
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
        for i=1:size(labels,1)
            if ismember(i,a1)
                hold on
                h(i)=bar(i,abs(cell2mat(mw_PCA(u_norm,i))),'FaceColor','k');
            else
                h(i)=bar(i,abs(cell2mat(mw_PCA(u_norm,i))),'FaceColor',[0.5 0.5 0.5]);
            end
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
xlim([0 size(x,1)+1])
xlabel('PC')
xlim([0 size(x,1)+1])
if ~strcmp(u_rank,'pc') && ~strcmp(u_rank,'sc')
    ylabel(evaluat)
else 
    ylabel(['|',evaluat,'|'])
end


subplot(2,2,[3 4])
bar(explained{u_norm})
title('Explained variance for the respective principal components')
set(gca,'XTick',ticknumb)
set(gca,'XTickLabelRotation',90)
xlim([0 size(x,1)+1])
xlabel('PC')
ylabel('Explained Variance (%)')


%%%%% PC-corr %%%%%
if length(cutoff)==1
    [Edges,Nodes,pc_corr,x1] = C_corr(norm{u_norm},pc{u_norm}(:,u_dim),feat_names,cutoff);
    
    
    [NodeColor,n1_f,n2_f]=match_V_samplCol(col,x1,labels, nameLabels, Nodes);
    Nodes(1:end,3)=Nodes(:,2);
    Nodes{1,2}='Colour';
    Nodes(2:end,2)=NodeColor;
    
    edges=Edges
    nodes=Nodes
else 
    for i=1:length(cutoff)
        Edges{i,1}=cutoff(i);
        Nodes{i,1}=cutoff(i);
        pc_corr{i,1}=cutoff(i);
        x2{i,1}=cutoff(i);
        [Edges{i,2},Nodes{i,2},pc_corr{i,2},x1] = C_corr(norm{u_norm},pc{u_norm}(:,u_dim),feat_names,cutoff(i));

        x2{i,2}=x1;
        
        [NodeColor,n1_f,n2_f]=match_V_samplCol(col,x2{i,2},labels, nameLabels, Nodes{i,2});
        Nodes{i,2}(1:end,3)=Nodes{i,2}(:,2);
        Nodes{i,2}{1,2}='Colour';
        Nodes{i,2}(2:end,2)=NodeColor;
        
        cut_off=cutoff(i)
        edges=Edges{i,2}
        nodes=Nodes{i,2}
    end
end

filename = 'PC-corr_net_edges-nodes_results.xlsx';
if exist('PC-corr_net_edges-nodes_results.xlsx','file')~=0
    delete('PC-corr_net_edges-nodes_results.xlsx')
end

% saving the PC-corr results in Excel spreadsheet file
if length(cutoff)==1
    xlswrite(filename,Edges,'Edges')
    xlswrite(filename,Nodes,'Nodes')
    RemoveSheet123(filename);
else
    for i=1:length(cutoff)
        warning('off','MATLAB:xlswrite:AddSheet')
        xlswrite(filename,Edges{i,2},['Edges-cutoff ',num2str(cutoff(i))])
        xlswrite(filename,Nodes{i,2},['Nodes-cutoff ',num2str(cutoff(i))])
    end
    RemoveSheet123(filename);
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
for k=1:length(cutoff)
    if length(cutoff)==1
        pc_corr1=pc_corr;
        nodes1=Nodes;
    else
        pc_corr1=pc_corr{k,2};
        nodes1=Nodes{k,2};
    end
    
    if isempty(pc_corr1)
       continue 
    else
        
        plot_graph(pc_corr1,nodes1,cutoff(k))

    end
end



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


function [edges,nodes,pc_corr,x1]=C_corr(x,V,feat_names,cutoff)
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

% - Normal Pearson Corrleation calculation
c=corr(x); 

% - Apply the C-Corr formula
pc_corr=zeros(length(V));
for i=1:length(V)-1
    for j=i+1:length(V)     
         pc_corr(i,j)=sign(c(i,j))*min(abs([c(i,j),V(i),V(j)]));       
    end
end

if max(abs(pc_corr(:)))<=cutoff
    c_temp=corr(x_temp);
    pc_corr_temp=zeros(length(V_temp));
    for i=1:length(V_temp)-1
        for j=i+1:length(V_temp)
            pc_corr_temp(i,j)=sign(c_temp(i,j))*min(abs([c_temp(i,j),V_temp(i),V_temp(j)]));
        end
    end
   warning('With this cut-off, there are no edges that have |PC_corr(i,j)|>cutoff');
   warn = strcat('Try with another cut-off lower than ',{' '},num2str(max(abs(pc_corr_temp(:)))));
   warning(warn{1});
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


if (l1==0)&&(t2>t3)
    t=t2;
    NodeColor(b==1)={'Black'};
    n1_f=NaN;
    n2_f=1-t;
elseif (l1==0)&&(t3>t2)
    t=t3;
    NodeColor(b==1)={'Red'};
    n1_f=NaN;
    n2_f=1-t;
elseif (l2==0)&&(t1>t4)
    t=t1;
    NodeColor(b==0)={'Red'};
    n1_f=1-t;
    n2_f=NaN;
elseif  (l2==0)&&(t4>t1)
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

G=sparse(triu(pc_corr1));
g=biograph(G,nodes1(2:end,1),'ShowArrows','off','EdgeType','straight', 'NodeAutoSize','off');


for i=1:length(g.nodes)
    g.nodes(i).Shape='circle';
    g.nodes(i).Size=[10,10];
    g.nodes(i).TextColor=[0.5 0.5 0.57];
    if strcmp(nodes1{i+1,2},'Red')
        g.nodes(i).Color=[1 0 0];
        g.nodes(i).LineColor=[1 0 0];
    else
        g.nodes(i).Color=[0 0 0];
        g.nodes(i).LineColor=[0 0 0];
    end
end

for i=1:length(g.edges)
    edg_w(i)=g.edges(i).Weight;
end

for i=1:length(g.edges)
    max_edg_w=max(edg_w);
    min_edg_w=min(edg_w);
    if sign(g.edges(i).Weight)==1
        color_edg=1-g.edges(i).Weight/max_edg_w;
        g.edges(i).LineColor=[1 color_edg color_edg];
    else
        color_edg=g.edges(i).Weight/abs(min_edg_w)+1;
        g.edges(i).LineColor=[color_edg color_edg color_edg];
    end
end

%Edge under frustation are in grey colour.
edg_frust=0;
for i=1:length(g.nodes)-1
    for j=i+1:length(g.nodes)
        eh1 = getedgesbynodeid(g,g.nodes(i).ID,g.nodes(j).ID);
        if isempty(eh1)
            continue
        else
            n=strcmp(nodes1{i+1,2},'Red');
            m=strcmp(nodes1{j+1,2},'Red');
            if (n*m==1)&&(sign(eh1.Weight)==-1)
                eh1.LineColor=[0.6 0.6 0.6];
                edg_frust=edg_frust+1;
            elseif (n+m==0)&&(sign(eh1.Weight)==-1)
                eh1.LineColor=[0.6 0.6 0.6];
                edg_frust=edg_frust+1;
            elseif (n+m==1)&& (sign(eh1.Weight)==1)
                eh1.LineColor=[0.6 0.6 0.6];
                edg_frust=edg_frust+1;
            end
        end
    end
end

val_red=[];
val_black=[];
for i=1:length(g.edges)
    if (g.edges(i).LineColor(1)==1) && (g.edges(i).LineColor(2)~=0)
        val_red=[val_red g.edges(i).LineColor(2)];
    else
         val_red=[val_red];
    end
    if g.edges(i).LineColor(1)==0
        val_black=[val_black g.edges(i).LineColor(1)];
    else
        val_black=[val_black];
    end
end
min_red=min(val_red);
min_black=min(val_black);

bg=biograph.bggui(g);
f1=get(bg.biograph.hgAxes,'Parent');
set(f1, 'HandleVisibility', 'on')
f = figure();
copyobj(bg.biograph.hgAxes,f);
set(f,'units','points');
set(f,'Position',f1.Position);
set(f,'Color',[1 1 1])
close(bg.hgFigure);


%Replace the grey edges (edges under frustration) in the figure with dashed grey edges 
fr_edg_frust=edg_frust/length(g.edges);
axesHandle = gca;
plotHandle = findobj(axesHandle,'Type','line','Color',[.6 .6 .6]);
plotHandle1 = findobj(axesHandle,'Type','line','Color',[1 0 0]);
if isempty(plotHandle1)
    plotHandle1=findobj(axesHandle,'Type','line','Color',[1 min_red min_red]);
end

plotHandle2 = findobj(axesHandle,'Type','line','Color','k','LineWidth',0.5);
if isempty(plotHandle2)
    plotHandle2=findobj(axesHandle,'Type','line','Color',[min_black min_black min_black]);
end
plotHandle3 = findobj(axesHandle,'Type','Patch','FaceColor','r');
plotHandle4 = findobj(axesHandle,'Type','Patch','FaceColor','k');
for i=1:length(plotHandle)
    set(plotHandle(i),'Color',[.6 .6 0.6],'Linestyle','--');
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
txt_figure={[],['   cut-off = ',num2str(cutoff)],[],['   frustation = ',frac_txt,' = ',num2str(diground(fr_edg_frust,3)*100),'\%']};

TextH = text(0,1, txt_figure, ...
  'HorizontalAlignment', 'left', ...
  'VerticalAlignment', 'top','Interpreter','latex','FontSize',11);

if isempty(plotHandle4)&& isempty(plotHandle2) && ~isempty(plotHandle) && ~isempty(plotHandle1) && ~isempty(plotHandle3)
    h=legend([plotHandle1(1) plotHandle(1) plotHandle3(1)],{'PC-corr $> 0$', 'Edge under frustation', '$\uparrow$ Red sample group'});
elseif isempty(plotHandle1)&& isempty(plotHandle3) && ~isempty(plotHandle) && ~isempty(plotHandle2) && ~isempty(plotHandle4)
    h=legend([plotHandle2(1) plotHandle(1) plotHandle4(1)],{'PC-corr $< 0$', 'Edge under frustation', '$\uparrow$ Black sample group'});
elseif isempty(plotHandle2)&& ~isempty(plotHandle3) && ~isempty(plotHandle) && ~isempty(plotHandle1) && ~isempty(plotHandle4)
    h=legend([plotHandle1(1) plotHandle(1) plotHandle3(1) plotHandle4(1)],{'PC-corr $> 0$', 'Edge under frustation', '$\uparrow$ Red sample group', '$\uparrow$ Black sample group'});
elseif isempty(plotHandle1)&& ~isempty(plotHandle3) && ~isempty(plotHandle) && ~isempty(plotHandle2) && ~isempty(plotHandle4)
    h=legend([plotHandle2(1) plotHandle(1) plotHandle3(1) plotHandle4(1)],{'PC-corr $< 0$', 'Edge under frustation', '$\uparrow$ Red sample group', '$\uparrow$ Black sample group'});
elseif isempty(plotHandle3)&& isempty(plotHandle2)
    if isempty(plotHandle)
        h=legend([plotHandle1(1) plotHandle4(1)],{'PC-corr $> 0$','$\uparrow$ Black sample group'}); %only black nodes
    else
        h=legend([plotHandle1(1) plotHandle(1) plotHandle4(1)],{'PC-corr $> 0$', 'Edge under frustation','$\uparrow$ Black sample group'}); %only black nodes
    end
elseif isempty(plotHandle4)&& isempty(plotHandle2)
    if isempty(plotHandle)
        h=legend([plotHandle1(1) plotHandle3(1)],{'PC-corr $> 0$','$\uparrow$ Red sample group'}); %only red nodes
    else
        h=legend([plotHandle1(1) plotHandle(1) plotHandle3(1)],{'PC-corr $> 0$','Edge under frustation','$\uparrow$ Red sample group'}); %only red nodes
    end
elseif isempty(plotHandle2)
    if isempty(plotHandle)
        h=legend([plotHandle1(1) plotHandle3(1) plotHandle4(1)],{'PC-corr $> 0$','$\uparrow$ Red sample group','$\uparrow$ Black sample group'}); %only positive PC-corr
    else
        h=legend([plotHandle1(1) plotHandle(1) plotHandle3(1) plotHandle4(1)],{'PC-corr $> 0$','Edge under frustation','$\uparrow$ Red sample group','$\uparrow$ Black sample group'}); %only positive PC-corr
    end
elseif isempty(plotHandle1)
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