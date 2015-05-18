function draw_imagesc_with_seperation(matrix_to_draw, seperation_indices, cell_type_order)
     num_celltype = size(matrix_to_draw,1);
    figure; hold on;
    imagesc(matrix_to_draw);  
    colormap(jet); 
    h=colorbar;
    set(h,'fontsize',20);
    
     for i =1 :(length(seperation_indices)-1)
        x = [0.5 ,size(matrix_to_draw,2)+0.5];
        y = [seperation_indices(i) ,seperation_indices(i)];
        plot(x,y,'Color','w','LineStyle','-');
        plot(x,y,'Color','k','LineStyle',':');
     end
    hold off;
    
    axis ij;
    ylim([0.5 size(matrix_to_draw,1)+0.5]);
    
    ax = gca;
    z = [1,seperation_indices(1:end-1) ; seperation_indices];
    ax.YTick = round(mean(z,1));
%     ax.YTick = [ round((1 + sep_index1)/2), round((sep_index1 + sep_index2)/2), round((sep_index2 + num_okaty_celltype)/2)];
    ax.YTickLabel = cell_type_order;
    
%      matrix_to_draw2 = [matrix_to_draw(1:sep_index1,:); -10*ones(3,3); matrix_to_draw( (sep_index1+1):(sep_index2 ) ,:); -10*ones(3,3); matrix_to_draw((sep_index2+1): end,:) ];
%      imagesc(matrix_to_draw2);
%      
end