if exist('inited','var')
    fprintf('Already inited, skip.\n');
    return
end
fprintf('Start initialization: begin init.m\n');

addpath('/cortex/code/cellmix/nmf/');
addpath('/home/lab/lior/Projects/load datasets/');
addpath('visualization/');
addpath('evaluation/');
addpath('baselines/');
addpath('dataset_creation/');
addpath('infrastructure/');

parms.dummy=1;

[do_save, parms] = take_from_struct(parms, 'do_save', false);




delete(gcp('nocreate'));parpool local;

inited = 1;
