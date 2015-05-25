function parms =load_markers(dataset, gene_info, num_genes,num_regions,parms)

    switch dataset
        case {'kang_regions','Human6_selected_regions', 'brainspan_rnaseq_DFC_OFC'}
            neuro_markers = {'STMN2','BRUNOL4','BRUNOL5','BRUNOL6','SYT1'};
            astro_markers = { 'GFAP' , 'ALDH1L1'};
            oligo_markers = {'OLIG2' ,'MBP'};
        case 'Zapala_isocortex_medulla_striatum_cerebellum'
            neuro_markers = {'Stmn2','Celf4','Syt1','Gria3','Dlg3','Dlg4','Tubb3','Map2'};
            astro_markers = { 'Gfap' , 'Aqp4' ,'Fgfr3' ,'Slc1a2' ,'Gjb6' };
            oligo_markers = {'Mbp' ,'Sox10' ,'Mag' ,'Mog'};
        otherwise
            error('unkown dataset %s - cannot load markers', dataset)
    end


    H_markers = false(parms.num_types,num_genes);
    H_markers(1, ismember(gene_info.gene_symbols,neuro_markers) ) = true;
    H_markers(2, ismember(gene_info.gene_symbols,astro_markers) ) = true;
    H_markers(3, ismember(gene_info.gene_symbols,oligo_markers) ) = true;

    % TODO - change so each region has its own markers
    parms.H_markers = repmat({H_markers}, num_regions,1);
    fprintf('Using cell type specific markers for dataset %s\n', dataset);
end