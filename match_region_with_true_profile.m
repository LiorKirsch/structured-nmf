function [celltype_expression, celltypes_used] = match_region_with_true_profile(mouse_cell_types, regions)
    celltype_expression = cell(length(regions),1);
    celltypes_used = cell(length(regions),1);
    for i = 1:length(regions)
        curr_region = regions{i};
        
%         {'Cerebral_cortex';'Motor Cortex';'Entorhinal Cortex';'Perirhinal Cortex';...
%     'Hippocampus';'Dentate Gyrus';'Hippocampus CA1';'Hippocampus CA3';...
%     'Striatum';'Cerebellum';'Medulla'};

%         {'H_lva_cortex';'H_lvb_cortex';'H_l6_cortex';'H_cck+_cortex';'H_mn_brainstem';'H_chol_basal_fb';...
%             'H_chol_spinal';'H_chol_corp_striatum';'H_cort+_cortex';'H_drd1+_msn_striatum';'H_drd2+_msn_striatum';...
%             'H_golgi_cerebellum';'H_uni_brush_cerebellum';'H_stell_basket_cerebellum';'H_granule_cerebellum';...
%             'H_molig_cerebellum';'H_molig_cortex';'H_mixed_olig_cerebellum';'H_mixed_olig_cortex';...
%             'H_purkinje_cerebellum';'H_pnoc+_cortex';'H_bergmann_cerebellum';'H_astro_cerebellum';'H_astro_cortex'}

   all_strct = {'Nervous_system','Telencephalon','Diencephalon','Mesencephalon',...
        'Metencephalon','cerebellar cortex','primary auditory (A1) cortex',...
    'hippocampus','posterior inferior parietal cortex',...
    'primary motor (M1) cortex',...
    'primary somatosensory (S1) cortex', 'primary visual (V1) cortex',...
    'temporal cortex', 'superior temporal cortex','inferior temporal cortex',...
    'prefrontal cortex', 'medial prefrontal cortex','orbital prefrontal cortex',...
    'ventrolateral prefrontal cortex','dorsolateral prefrontal cortex',...
    'Cerebral Nuclei','striatum','amygdala',...
    'mediodorsal nucleus of the thalamus'};

        switch curr_region
            case {'Cerebral_cortex','Motor Cortex','Entorhinal Cortex','Perirhinal Cortex',...
                  'Hippocampus','Dentate Gyrus','Hippocampus CA1','Hippocampus CA3',...
                  'primary auditory (A1) cortex',...
                'hippocampus','posterior inferior parietal cortex',...
                'primary motor (M1) cortex',...
                'primary somatosensory (S1) cortex', 'primary visual (V1) cortex',...
                'temporal cortex', 'superior temporal cortex','inferior temporal cortex',...
                'prefrontal cortex', 'medial prefrontal cortex','orbital prefrontal cortex',...
                'ventrolateral prefrontal cortex','dorsolateral prefrontal cortex',...
                'Frontal Lobe','Occipital Lobe','Parietal Lobe','Temporal Lobe'};
                neuro_id = {'H_lva_cortex';'H_lvb_cortex';'H_l6_cortex';'H_cck+_cortex';...
                            'H_cort+_cortex'; 'H_pnoc+_cortex'};
                oligo_id = {'H_molig_cortex';'H_mixed_olig_cortex'};
                astro_id = {'H_astro_cortex'};
            case {'Cingulate gyrus'}
                neuro_id = {'S_yfp_cg';'S_bcg';'S_g30_cg';'S_gin_cg';'S_g42_cg'} ;
                oligo_id = {};
                astro_id = {};
            case {'Mesencephalon'}
                neuro_id = {'I_A9';'I_A10'} ;
                oligo_id = {};
                astro_id = {};
            case {'hippocampal formation'}
                neuro_id = {'S_yfp_h';'S_gin_h'} ;
                oligo_id = {};
                astro_id = {};
            case {'Basal Forebrain'}
                neuro_id = {'H_chol_basal_fb'} ;
                oligo_id = {};
                astro_id = {};
            case {'Thalamus'}
                neuro_id = {'S_g42_lg'} ;
                oligo_id = {};
                astro_id = {};
            case {'Striatum', 'striatum'}
                neuro_id = { 'H_chol_corp_striatum';'H_drd1+_msn_striatum';'H_drd2+_msn_striatum'} ;
                oligo_id = {};
                astro_id = {};
            case {'amygdala','Amygdala'}
                neuro_id = {'S_yfp_a';'S_g30_a'};
                oligo_id = {};
                astro_id = {};
            case {'Cerebellum','Cerebellar Cortex', 'cerebellar cortex'}
                neuro_id = { 'H_golgi_cerebellum';'H_stell_basket_cerebellum';'H_granule_cerebellum';...
                    'H_purkinje_cerebellum'};
                oligo_id = {'H_molig_cerebellum';'H_mixed_olig_cerebellum'};
                astro_id = {'H_uni_brush_cerebellum';'H_bergmann_cerebellum';'H_astro_cerebellum'};
            case {'Medulla','Myelencephalon'}
                neuro_id = { 'H_mn_brainstem'};
                oligo_id = {};
                astro_id = {};
            otherwise
                error('unknown region %s' ,curr_region)
        end
        
        expression = [];
        filterneuron = ismember(mouse_cell_types.cell_type_id, neuro_id);
        sample2type = mouse_cell_types.sample2type(:,filterneuron);
        if ~isempty(sample2type)
            all_samples = any(sample2type,2);
            mean_expression = mouse_cell_types.expression * double(all_samples) / sum(all_samples);
            expression = cat(2,expression, mean_expression);
        end
  
        filterneuron = ismember(mouse_cell_types.cell_type_id, astro_id);
        sample2type = mouse_cell_types.sample2type(:,filterneuron);
        if ~isempty(sample2type)
            all_samples = any(sample2type,2);
            mean_expression = mouse_cell_types.expression * double(all_samples) / sum(all_samples);
            expression = cat(2,expression, mean_expression);
        end
        
        filterneuron = ismember(mouse_cell_types.cell_type_id, oligo_id);
        sample2type = mouse_cell_types.sample2type(:,filterneuron);
        if ~isempty(sample2type)
            all_samples = any(sample2type,2);
            mean_expression = mouse_cell_types.expression * double(all_samples) / sum(all_samples);
            expression = cat(2,expression, mean_expression);
        end
        
        celltype_expression{i} = expression;
        types.neuron = neuro_id;
        types.astro = astro_id;
        types.oligo = oligo_id;
        celltypes_used{i} = types;
    end


end