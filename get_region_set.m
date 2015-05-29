function [region_set, species] = get_region_set(dataset)
%
% regionset_name is ignoreed at this point.
     switch dataset
        case 'kang2011',
          species = 'human';
          region_set =  {'A1C','AMY','CBC','DFC','HIP', 'IPC','ITC', ...
                         'M1C','MFC','OFC','S1C', 'STC','V1C','VFC'};

        case 'brainspan2014',
          species = 'human';
          region_set =  {'A1C','AMY','CBC','DFC','HIP', 'IPC','ITC', ...
                         'M1C','MFC','OFC','S1C', 'STC','V1C','VFC'};
          region_set =  {'DFC','OFC'};

        case 'zapala2005',
          species = 'mouse';
          region_set = {'Isocortex';'Motor Cortex';'Entorhinal Cortex';'Perirhinal Cortex';...
                        'Hippocampus';'Dentate Gyrus';'Hippocampus CA1';'Hippocampus CA3';...
                        'Striatum';'Cerebellum';'Medulla'};

        case 'human6',
          species = 'human';
          region_set = {'Frontal Lobe';'Cingulate gyrus';'hippocampal formation';...
                         'Occipital Lobe';'Parietal Lobe';'Temporal Lobe';'Amygdala';...
                         'Basal Forebrain';...
                         'Striatum';'Thalamus';'Mesencephalon';'Cerebellar Cortex';...
                         'Myelencephalon'};

        otherwise 
          error('invalid dataset = [%s]\n', dataset);
    end
end