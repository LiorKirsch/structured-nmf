function draw_structure_matrix(strucutre_matrix, regions)

    figure;imagesc(strucutre_matrix);colorbar; colormap(jet);
    ax = gca;
    ax.XTick = 1:length(regions);
    ax.XTickLabel = regions;
    ax.YTick = 1:length(regions);
    ax.YTickLabel = regions;
    ax.XTickLabelRotation	=45;
    
end