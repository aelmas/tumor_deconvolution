clear 
%% 1. load data
load data_MC3.mat
load data_mRNA.mat
load data_PRO.mat
% common samples between the transcriptome and proteome data sets
[common_samples_mRNA_PRO, is_mRNA, is_PRO] = intersect(data_mRNA.Properties.VariableNames, data_PRO.Properties.VariableNames);
% common genetic elements between the transcriptome and proteome data sets
[common_genes_mRNA_PRO, ig_mRNA, ig_PRO] = intersect(data_mRNA.Properties.RowNames, data_PRO.Properties.RowNames);

%% 2. clean data
raw_values_mRNA = table2array(data_mRNA);
raw_values_PRO = table2array(data_PRO);
flag_normalize_mRNA = false;
flag_normalize_PRO = false;
% 2.a. impute data
if ~isempty(find(isnan(raw_values_mRNA), 1))
    tic; values_mRNA = knnimpute(raw_values_mRNA); toc;
    flag_normalize_mRNA = true;
else
    values_mRNA = raw_values_mRNA;
end
if ~isempty(find(isnan(raw_values_PRO), 1))
    tic; values_PRO = knnimpute(raw_values_PRO); toc;
    flag_normalize_PRO = true;
else
    values_PRO = raw_values_PRO;
end

% 2.b. remove outlier samples
c_dispersion = 3; %dispersion constant
non_outliers_mRNA = iqr(values_mRNA) <= mean(values_mRNA)+c_dispersion*std(values_mRNA);
values_mRNA = values_mRNA(:,non_outliers_mRNA);
non_outliers_PRO = iqr(values_PRO) <= mean(values_PRO)+c_dispersion*std(values_PRO);
values_PRO = values_PRO(:,non_outliers_PRO);

% 2.c. normalize samples
% mRNA | normalized across samples at 75%% quantile | log2 normalized
if flag_normalize_mRNA
    values_mRNA = values_mRNA-repmat(mean(values_mRNA),size(values_mRNA,1),1); % 0-mean
    values_mRNA = values_mRNA./repmat(std(values_mRNA),size(values_mRNA,1),1); % 1-SD
    log_mean_norm_error_mRNA = log(mean(abs(mean(values_mRNA)-0)))
    log_std_norm_error_mRNA = log(mean(abs(std(values_mRNA)-1)))
end
% PRO | samples normalized | 0-mean, 1-SD | 
if flag_normalize_PRO
    values_PRO = values_PRO-repmat(mean(values_PRO),size(values_PRO,1),1); % 0-mean
    values_PRO = values_PRO./repmat(std(values_PRO),size(values_PRO,1),1); % 1-SD
    log_mean_norm_error_PRO = log(mean(abs(mean(values_PRO)-0)))
    log_std_norm_error_PRO = log(mean(abs(std(values_PRO)-1)))
end
% check distributions
figure(1);
subplot(211); boxplot(values_mRNA); ylabel('mRNA');
subplot(212); boxplot(values_PRO); ylabel('Proteome');


%% 3. deconvolution methods
%% 3.a. tumor transcriptome
% xCell, performs cell type enrichment analysis from gene expression data for 64 
%   immune and stroma cell types by reducing associations between closely related cell types.
% xCell produces enrichment scores, not percentages.
% xCell uses the expression levels ranking and not the actual values, 
%   thus normalization does not have an effect.
% installation\\$ devtools::install_github('dviraran/xCell')
%
% pipe data to xCell
data_mRNA_cleaned = array2table(values_mRNA);
data_mRNA_cleaned.Properties.VariableNames = data_mRNA.Properties.VariableNames;
data_mRNA_cleaned.Properties.RowNames = data_mRNA.Properties.RowNames;
writetable(data_mRNA_cleaned,'BRCA_mRNA_formatted_normalized_cleaned.txt','Delimiter','\t','WriteRowNames',true);
% script to call xCell 
fid = fopen('call_xCell_mRNA.R','w');
fprintf(fid,[
    'library(xCell)\n' ...
    'expression_matrix = read.table("BRCA_mRNA_formatted_normalized_cleaned.txt",header=TRUE,row.names=1, as.is=TRUE)\n' ... 
    'xCell_result = xCellAnalysis(expression_matrix'...
    ', parallel.type = "FORK"' ... % faster in unix
    ', cell.types.use = NULL)\n' ... % use all cell types
    'write.table(xCell_result,file = "xCell_result_BRCA_mRNA_formatted_normalized_cleaned.txt",sep="\t",quote=FALSE)\n' ... 
    ]);
fclose(fid); 
% grand user permission of file executions
system('chmod u+x call_xCell_mRNA.R xCell-master/R/xCell.R') 
% call xCell
system('/usr/local/bin/Rscript call_xCell_mRNA.R');
% read xCell result
data_mRNA_cleaned_xCell_result = readtable('xCell_result_BRCA_mRNA_formatted_normalized_cleaned.txt','ReadRowNames',true);
% xCell has a bug that assigns the first variable to the row names | fix this
data_mRNA_cleaned_xCell_result.Properties.VariableNames = data_mRNA_cleaned.Properties.VariableNames; 
% display cell enrichments
figure(2);
heatmap(data_mRNA_cleaned_xCell_result.Properties.VariableNames,data_mRNA_cleaned_xCell_result.Properties.RowNames,table2array(data_mRNA_cleaned_xCell_result))

%% 3.b. tumor proteome
% a quick observation shows that there are ~10.5K common genetic elements between the analyzed transcriptome and proteome data sets:
length(common_genes_mRNA_PRO)
% a deconvolution method designed for transcriptome data may worth studying
%   the proteome data that contains many common genetic signatures with transcriptome
%
% pipe data to xCell
data_PRO_cleaned = array2table(values_PRO);
data_PRO_cleaned.Properties.VariableNames = data_PRO.Properties.VariableNames;
data_PRO_cleaned.Properties.RowNames = data_PRO.Properties.RowNames;
writetable(data_PRO_cleaned,'BRCA_PRO_formatted_normalized_cleaned.txt','Delimiter','\t','WriteRowNames',true);
% script to call xCell 
fid = fopen('call_xCell_PRO.R','w');
fprintf(fid,[
    'library(xCell)\n' ... 
    'expression_matrix = read.table("BRCA_PRO_formatted_normalized_cleaned.txt",header=TRUE,row.names=1, as.is=TRUE)\n' ... 
    'xCell_result = xCellAnalysis(expression_matrix'...
    ', parallel.type = "FORK"' ... % faster in unix
    ', cell.types.use = NULL)\n' ... % use all cell types
    'write.table(xCell_result,file = "xCell_result_BRCA_PRO_formatted_normalized_cleaned.txt",sep="\t",quote=FALSE)\n' ... 
    ]);
fclose(fid); 
% grand user permission of file executions
system('chmod u+x call_xCell_PRO.R xCell-master/R/xCell.R') 
% call xCell
system('/usr/local/bin/Rscript call_xCell_PRO.R');
% read xCell result
data_PRO_cleaned_xCell_result = readtable('xCell_result_BRCA_PRO_formatted_normalized_cleaned.txt','ReadRowNames',true);
% xCell has a bug that assigns the first variable to the row names | fix this
data_PRO_cleaned_xCell_result.Properties.VariableNames = data_PRO_cleaned.Properties.VariableNames; 
% display cell enrichments
figure(3);
heatmap(data_PRO_cleaned_xCell_result.Properties.VariableNames,data_PRO_cleaned_xCell_result.Properties.RowNames,table2array(data_PRO_cleaned_xCell_result))

[common_genes_mRNA_PRO_cleaned, ig_mRNA_cleaned, ig_PRO_cleaned] = intersect(data_mRNA_cleaned.Properties.RowNames, data_PRO_cleaned.Properties.RowNames);
