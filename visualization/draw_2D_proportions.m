function draw_2D_proportions(proportion_matrix, sample_group)
    figure; hold on;
    [num_samples, k] = size(proportion_matrix);
    
    if ~exist('sample_group','var')
        sample_group = ones(num_samples,1);
    end

    
    
    centers = get_circle_n_points(k);
    
    samples_2d = proportion_matrix * centers;

    
    scatter(samples_2d(:,1), samples_2d(:,2),23,sample_group,'filled');
    scatter(centers(:,1), centers(:,2),30,0*ones(size(centers,1),1),'filled');
end

function centers = get_circle_n_points(n)
    degree = 2*pi / n;
    
    x = cos(degree *(1:n)');
    y = sin(degree *(1:n)');
    centers = [x,y];
end