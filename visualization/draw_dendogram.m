function draw_dendogram(X, profile_names,region_name, corr_type)

[X,profile_names] = flatten_cell(X,profile_names, region_name);
figure;
if ~exist('corr_type');
    corr_type = 'pearson';
end

corr_X = corr(X, 'type',corr_type);

distances = 1 - corr_X;
Z = linkage(distances);
[H,T,outperm] = dendrogram(Z);

names_ordered = profile_names(outperm);
ax = gca;
ax.XTickLabel = names_ordered;
ax.XTickLabelRotation	=45;

figure;
imagesc(corr_X(outperm,outperm));colormap(jet);colorbar;
ax = gca;
ax.XTick = 1:length(names_ordered);
ax.XTickLabel = names_ordered;
ax.XTickLabelRotation	=45;

end

function [flat_X,flat_names] = flatten_cell(X,profile_names,regions)



    flat_X = [];
    flat_names = {};
    for i=1:length(X)   
        flat_X = [flat_X , X{i}];
        profile_str = strcat(profile_names{i},'/',regions{i});
        flat_names = [flat_names ; profile_str(:)];
    end
    
    
end