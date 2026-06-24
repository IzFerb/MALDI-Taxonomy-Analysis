dirinfo = dir();
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
dirinfo(ismember({dirinfo.name}, 'resources')) = [];
dirinfo(ismember({dirinfo.name}, 'icoshift3p0')) = [];
dirinfo(ismember({dirinfo.name}, '.git')) = [];
filenames = [];
subdirinfo = cell(length(dirinfo));
count = 0;
y = [];
for K = 1 : length(dirinfo)
    thisdir = dirinfo(K).name;
    tmp = dir(fullfile(thisdir, '*.txt'));
    filenames = [filenames {tmp.name}];
    count = count + 1;
    y(length(y)+1:length(y)+length(tmp)) = count;
end
filenames(5) = [];
y(5) = [];
reference_samples = (y == 3) | (y == 4) | (y == 6);
reference_samples = find(reference_samples);
%%
xq = 800:0.01:1300;
vq2 = zeros(length(filenames), length(xq));
for l1 = 1:length(filenames)
    x_tmp = readtable(filenames{l1});
    x_tmp = table2array(x_tmp);
    vq2(l1,:) = interp1(x_tmp(:,1),x_tmp(:,2),xq,'spline');
end
%%
[X,ints,ind,target] = icoshift('median',vq2,'whole',34);
figure
plot(xq,X)
%% SVN
% Q = X./median(X(reference_samples,:));
% X = X./repmat(median(Q,2), 1,size(X,2));
X = (X - mean(X,2))./std(X,[],2); 
figure
plot(xq,X)
%% peak picking
for l1 = 1:size(X,1)
    [pks{l1},locs{l1},w{l1},p{l1}] = findpeaks(X(l1,:),xq,'MinPeakProminence',std(X(l1,:))*10);
    [~,locs_ref{l1},~,~] = findpeaks(X(l1,:),xq,'MinPeakProminence',std(X(l1,:))*10);
end
%% peak screening
maximum_width = mean([w{reference_samples}]) + std([w{reference_samples}])*3;
all_peaks = unique([locs_ref{reference_samples}]);
l1 = 0;
X_new = [];
while ~isempty(all_peaks)
    l1 = l1 + 1;
    selected_peaks = all_peaks > all_peaks(1) - maximum_width & all_peaks < all_peaks(1) + maximum_width;
    VariableNames(l1) = median(all_peaks(selected_peaks));
    all_peaks(selected_peaks) = []; 
    for l2 = 1:size(X,1)
        a = locs{l2};
        a = find(a(a > VariableNames(l1) - maximum_width & a < VariableNames(l1) + maximum_width));
        if isempty(a)
            X_new(l2,l1) = 0;
        elseif length(a) == 1
            X_new(l2,l1) = p{l2}(a);
        else
            X_new(l2,l1) = max(p{l2}(a));
        end
    end
end
%% PQN
VariableNames = VariableNames(~any(X_new(1:4,:)));
Z = X_new(:,~any(X_new(1:4,:)));
%%
Z2 = Z(any(Z,2) & y' ~= 7,:);
filenames_kept = filenames(any(Z,2) & y' ~= 7);
filenames_removed = filenames(~any(Z,2) | y' == 7);
y2 = y(any(Z,2) & y' ~= 7); 
%% Acacia
Z_acacia = Z2(y2 == 3,:);
filenames_acacia  = filenames_kept(y2 == 3); 
k = strfind(filenames_acacia , 'Acacia');
class = categorical(1); 
for l1 = 1:length(filenames_acacia )
    k1 = strfind(filenames_acacia {l1}, '_');
    k2 = k1(k1 > k{l1});
    class(l1) = filenames_acacia {l1}(k2(1)+1:k2(2)-1);
end
class(39) = 'seyal';
class(50) = 'stenocarpa'; 
%%
class(20) = [];
Z_acacia(20,:) = [];
filenames_acacia(20) = []; 
D = [];
for l1 = 1:length(class)
    for l2 = 1:length(class)
        D(l1,l2) = simbin(Z_acacia(l1,:), Z_acacia(l2,:), 'Pei1');
    end
end

%%
colors = distinguishable_colors(22,[0 0 0; 1 1 1]);
[SCR, LDS] = plot_pca_metabolomics(D,cellstr(class)', colors, unique(class, 'stable'), 1:length(class),1,2);
[SCR, LDS] = plot_pca_metabolomics(D,cellstr(class)', colors, unique(class, 'stable'), 1:length(class),1,3);
[SCR, LDS] = plot_pca_metabolomics(D,cellstr(class)', colors, unique(class, 'stable'), 1:length(class),1,4);
[SCR, LDS] = plot_pca_metabolomics(D,cellstr(class)', colors, unique(class, 'stable'), 1:length(class),2,3);
[SCR, LDS] = plot_pca_metabolomics(D,cellstr(class)', colors, unique(class, 'stable'), 1:length(class),2,4);
[SCR, LDS] = plot_pca_metabolomics(D,cellstr(class)', colors, unique(class, 'stable'), 1:length(class),3,4);
%%
Z_historical = Z(y == 7,:);
filenames_historical  = filenames(y == 7); 
filenames_historical([6 9 10]) = [];
Z_historical([6 9 10],:) = [];
class2 = class;
for l1 = 1:length(filenames_historical)
    for l2 = 1:length(filenames_acacia)
        D_historical(l1,l2) = simbin(Z_historical(l1,:) , Z_acacia(l2,:), 'Pei1');
    end
    class2(end+1) = 'historical';
end
%%
trainset = [true(56,1); false(11,1)];

[SCR, LDS] = plot_pca_metabolomics([D; D_historical],cellstr(class2)', colors, unique(class2, 'stable'), [],1,2,trainset);
%% project pigmented references
Z_pigmented = Z(y == 5,:);
filenames_pigmented  = filenames(y == 5); 
class3 = class;
for l1 = 1:length(filenames_pigmented)
    for l2 = 1:length(filenames_acacia)
        D_pigmented(l1,l2) = simbin(Z_pigmented(l1,:) , Z_acacia(l2,:), 'Pei1');
    end
    class3(end+1) = 'pigmented'; %cannot run historical
end
%%
trainset = [true(56,1); false(11,1)];

[SCR, LDS] = plot_pca_metabolomics([D; D_pigmented],cellstr(class3)', colors, unique(class3, 'stable'), [],1,2,trainset);

