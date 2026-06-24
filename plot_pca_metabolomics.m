function [SCR, LDS] = plot_pca_metabolomics(X,y, colors, legend_labels, text_label,C1, C2, trainset, var_select)
%%

%% meancenter data
if nargin <= 7
    trainset = true(size(X,1),1);
end
X_mean = mean(X(trainset,:));
X_std = std(X(trainset,:));
X_asc = (X - X_mean);%./X_std;
X_asc(:, any(isnan(X_asc))) = [];

% X_asc = X_log;
%% calculate PCA

if nargin <= 8
    [~, LDS, VAR] = pcafunction(X_asc(trainset,:));
else
    [~, LDS, VAR] = pcafunction(X_asc(trainset,var_select));
end
SCR = X_asc*LDS;
%% create figure
figure
g = gscatter(SCR(:,C1), SCR(:,C2), y,colors, 'ox+^ds', 16);
def = max(SCR(:,1)) - min(SCR(:,1));
set(g, 'linewidth', 2)
if ~isempty(text_label)
    text(SCR(:,C1)+def/100, SCR(:,C2), num2cell(text_label), 'fontsize', 16)
end
legend(legend_labels)
set(gca,'fontsize', 20)
xlabel(['PC ' num2str(C1) '(' num2str(VAR(C1)) '% variance explained)'])
ylabel(['PC ' num2str(C2) '(' num2str(VAR(C2)) '% variance explained)'])
axis('square')