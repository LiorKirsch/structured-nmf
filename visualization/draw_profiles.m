function draw_profiles(profiles,proportions,parms)

region_name = parms.relation_regions;
[ k, num_genes] = size(profiles{1});

all_profiles = cat(1,profiles{:});
region_name = repmat(region_name,1,k)';
region_name = region_name(:);

[wcoeff,score,latent,tsquared,explained] = pca(all_profiles);

regions = unique(region_name,'stable');
[~,ids] = ismember(region_name, regions);


% figure(); 
hold on;
gscatter(score(:,1),score(:,2),ids,[],[],23);
xlabel('1st Principal Component');
ylabel('2nd Principal Component');

% legend(strrep(regions, '_', ' '));

ind = 0;
samples2D = [];
sample_group = [];
for i = 1:length(regions)
    centers = score(ind+1: ind+3, 1:2);
    
    samples_2d = proportions{i}' * centers;
    group = i*ones(size(proportions{i} ,2),1);
    
    samples2D = [samples2D; samples_2d];
    sample_group = [sample_group; group];
    ind = ind +3;
end

% figure();
gscatter(samples2D(:,1),samples2D(:,2),sample_group);
xlabel('1st Principal Component');
ylabel('2nd Principal Component');

legend(strrep(regions, '_', ' '), 'location', 'best');

end

function draw_the_proportions(centers,proportions)

      samples_2d = proportion_matrix * centers;
     scatter(samples_2d(:,1), samples_2d(:,2),23,sample_group,'filled');
end