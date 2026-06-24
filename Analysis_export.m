% Analysis_export.m — export reference results for R validation
%
% Prerequisites:
%   - icoshift toolbox on MATLAB path
%   - sample_data/ extracted from "Gums files.zip"
%
% Usage (from project root):
%   cd('/path/to/art_institute_code')
%   Analysis_export
%
% Outputs saved to validation/matlab/

function Analysis_export()
    export_dir = fullfile('validation', 'matlab');
    if ~exist(export_dir, 'dir')
        mkdir(export_dir);
    end

    data_dir = 'sample_data';
    if ~exist(data_dir, 'dir')
        error('sample_data/ not found. Extract Gums files.zip first.');
    end

    dirinfo = dir(data_dir);
    dirinfo(~[dirinfo.isdir]) = [];
    exclude = {'.', '..', 'resources', 'icoshift3p0', '.git', 'validation'};
    dirinfo(ismember({dirinfo.name}, exclude)) = [];
    [~, ord] = sort({dirinfo.name});
    dirinfo = dirinfo(ord);

    filenames = {};
    y = [];
    count = 0;
    group_names = {};
    for K = 1:length(dirinfo)
        thisdir = dirinfo(K).name;
        tmp = dir(fullfile(data_dir, thisdir, '*.txt'));
        [~, tord] = sort({tmp.name});
        tmp = tmp(tord);
        for t = 1:length(tmp)
            filenames{end+1} = fullfile(data_dir, thisdir, tmp(t).name); %#ok<AGROW>
        end
        count = count + 1;
        group_names{count} = thisdir; %#ok<AGROW>
        y(end+1:end+length(tmp)) = count; %#ok<AGROW>
    end

    filenames(5) = [];
    y(5) = [];

    reference_groups = [4];
    acacia_group = 4;
    exclude_group = 2;

    reference_samples = find(ismember(y, reference_groups));

    xq = 800:0.01:1300;
    vq2 = zeros(length(filenames), length(xq));
    for l1 = 1:length(filenames)
        x_tmp = readtable(filenames{l1});
        x_tmp = table2array(x_tmp);
        vq2(l1,:) = interp1(x_tmp(:,1), x_tmp(:,2), xq, 'spline');
    end

    [X, ~, ~, ~] = icoshift('median', vq2, 'whole', 34);
    X = (X - mean(X,2)) ./ std(X, [], 2);

    for l1 = 1:size(X,1)
        [pks{l1}, locs{l1}, w{l1}, p{l1}] = findpeaks(X(l1,:), xq, 'MinPeakProminence', std(X(l1,:))*10); %#ok<ASGLU>
        [~, locs_ref{l1}, ~, ~] = findpeaks(X(l1,:), xq, 'MinPeakProminence', std(X(l1,:))*10);
    end

    ref_widths = [w{reference_samples}];
    maximum_width = mean(ref_widths) + std(ref_widths) * 3;
    all_peaks = unique([locs_ref{reference_samples}]);
    l1 = 0;
    X_new = [];
    VariableNames = [];
    while ~isempty(all_peaks)
        l1 = l1 + 1;
        selected_peaks = all_peaks > all_peaks(1) - maximum_width & all_peaks < all_peaks(1) + maximum_width;
        VariableNames(l1) = median(all_peaks(selected_peaks)); %#ok<AGROW>
        all_peaks(selected_peaks) = [];
        for l2 = 1:size(X,1)
            a = locs{l2};
            a = find(a > VariableNames(l1) - maximum_width & a < VariableNames(l1) + maximum_width);
            if isempty(a)
                X_new(l2,l1) = 0;
            elseif length(a) == 1
                X_new(l2,l1) = p{l2}(a);
            else
                X_new(l2,l1) = max(p{l2}(a));
            end
        end
    end

    keep = ~any(X_new(1:4,:));
    Z = X_new(:, keep);

    Z2 = Z(any(Z,2) & y' ~= exclude_group, :);
    y2 = y(any(Z,2) & y' ~= exclude_group);

    Z_acacia = Z2(y2 == acacia_group, :);
    [~, fn_acacia] = fileparts(filenames(y2 == acacia_group));

    class = categorical(1);
    for l1 = 1:length(fn_acacia)
        k = strfind(fn_acacia{l1}, 'Acacia');
        k1 = strfind(fn_acacia{l1}, '_');
        k2 = k1(k1 > k{l1});
        class(l1) = fn_acacia{l1}(k2(1)+1:k2(2)-1);
    end
    class(39) = 'seyal';
    class(50) = 'stenocarpa';

    class(20) = [];
    Z_acacia(20,:) = [];
    fn_acacia(20) = [];

    D = zeros(length(class));
    for l1 = 1:length(class)
        for l2 = 1:length(class)
            D(l1,l2) = simbin(Z_acacia(l1,:), Z_acacia(l2,:), 'Pei1');
        end
    end

    writematrix(D, fullfile(export_dir, 'similarity_matrix.csv'));
    writematrix(Z_acacia, fullfile(export_dir, 'acacia_feature_matrix.csv'));
    writematrix(Z, fullfile(export_dir, 'feature_matrix.csv'));
    writetable(table(string(class)), fullfile(export_dir, 'acacia_classes.csv'));

    short_names = cellfun(@(f) erase(f, [data_dir filesep]), filenames, 'UniformOutput', false);
    manifest = table((1:length(filenames))', string(short_names)', y', ...
        string(group_names(y))', 'VariableNames', ...
        {'index', 'filename', 'group_id', 'group_name'});
    writetable(manifest, fullfile(export_dir, 'sample_manifest.csv'));

    fid = fopen(fullfile(export_dir, 'summary.txt'), 'w');
    fprintf(fid, 'n_samples_total=%d\n', length(filenames));
    fprintf(fid, 'n_features=%d\n', size(Z, 2));
    fprintf(fid, 'n_acacia_samples=%d\n', length(class));
    fprintf(fid, 'n_acacia_features=%d\n', size(Z_acacia, 2));
    fprintf(fid, 'similarity_dim=%dx%d\n', size(D, 1), size(D, 2));
    fprintf(fid, 'similarity_mean=%.8f\n', mean(D(:)));
    fprintf(fid, 'similarity_sd=%.8f\n', std(D(:)));
    fclose(fid);

    fprintf('MATLAB exports saved to %s\n', export_dir);
end
