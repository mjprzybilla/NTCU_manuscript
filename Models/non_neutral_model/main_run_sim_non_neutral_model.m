%% This script runs simulations of the non neutral model, used to produce
% to reproduce the theory curve shown in Figure 2G
%
% This script runs simulatiosn for the neutral competition model.
%
%% ------------------------------------------------------------------------
%
% 1) Directory with images of initial conditions for the numerical simulations
ic_path = 'data_sim_initial_condition\';

% 2) Directory where simulation results will be saved
work_path = 'output';

imfile = {};
imfile{1} = [ic_path,'Control_pre_treatment_Dorsal_1.jpg'];
imfile{2} = [ic_path,'Control_pre_treatment_Ventral_1.jpg'];
imfile{3} = [ic_path,'Control_pre_treatment_Dorsal_2.jpg'];
imfile{4} = [ic_path,'Control_pre_treatment_Ventral_2.jpg'];

% Simulation parameters
run_sim = 1;
save_files = 1;

nrepeats = 1;
fracs_prolif = [0.1]; % fraction of mutant cells
t_max = 15; % sim time
Nrows = 100;
Ncols = 300; %2mm
% here rate_divisions and rate_divisions_fast allows a slower division rate
% to one population
rate_divisions_fast = 0.1;%
rate_divisions = rate_divisions_fast; %
probs_duplicate = [1]; % duplication probability

factor = 1; % record interval

param_combos = [];
save_paths = {};
ns = [];
npars = 0;

% we run though all parameters and create folders to save files.
if run_sim == 1
    for nrdf = 1:length(rate_divisions_fast)
        rate_division_fast = rate_divisions_fast(nrdf);
        rate_division = rate_division_fast;
        for fprof = 1:length(fracs_prolif)
            frac_prolif = fracs_prolif(fprof);
            for pdup = 1:length(probs_duplicate)
                prob_duplicate = probs_duplicate(pdup);

                npars = npars + 1;
                param_combos.rate_division(npars) = rate_division;
                param_combos.rate_division_fast(npars) = rate_division_fast;
                param_combos.prob_duplicate(npars) = prob_duplicate;
                param_combos.frac_prolif(npars) = frac_prolif;

                sim_params = containers.Map;
                sim_params('nrepeats') = nrepeats;
                sim_params('t_max') = t_max;
                sim_params('factor') = factor;
                sim_params('Nrows') = Nrows;
                sim_params('Ncols') = Ncols;
                sim_params('rate_division') = rate_division;
                sim_params('rate_division_fast') = rate_division_fast;
                sim_params('prob_duplicate') = prob_duplicate;
                sim_params('frac_prolif') = frac_prolif;



                if save_files == 1
                    [exi,save_path] = create_save_folder(work_path,sim_params);
                    if exi == 0
                        save_paths{npars} = save_path;
                        save([save_path,'sim_params.mat'],'sim_params')
                    end
                    ns(end+1) = npars;
                end
            end
        end
    end
end
%%
% run simulations
if run_sim == 1
    % we run though each set of parameters
    for np=1:length(ns)
        n = ns(np);

        rate_division = param_combos.rate_division(n);
        rate_division_fast = param_combos.rate_division_fast(n);
        prob_duplicate = param_combos.prob_duplicate(n);
        frac_prolif = param_combos.frac_prolif(n);

        % codegen -args {300,0.5,0,100,100,0.1,0,0.25,1} sim_lattice_c
        sim_params = containers.Map;
        sim_params('nrepeats') = nrepeats;
        sim_params('t_max') = t_max;
        sim_params('factor') = factor;
        sim_params('Nrows') = Nrows;
        sim_params('Ncols') = Ncols;
        sim_params('rate_division') = rate_division;
        sim_params('rate_division_fast') = rate_division_fast;
        sim_params('prob_duplicate') = prob_duplicate;
        sim_params('frac_prolif') = frac_prolif;

        save_path = save_paths{n};
        if save_files == 1
            parsave_sim_params([save_path,'sim_params.mat'],sim_params)
        end

        % we run each realization in paralel
        for n_rep = 1:nrepeats
            im = [];

            if n_rep < nrepeats/4
                im = imread(imfile{1});
                im = double(im2bw(im))+1;
            elseif n_rep < nrepeats/2
                im = imread(imfile{2});
                im = double(im2bw(im))+1;
            elseif n_rep < 3*nrepeats/4
                im = imread(imfile{3});
                im = double(im2bw(im))+1;
            else
                im = imread(imfile{4});
                im = double(im2bw(im))+1;
            end

            % compile sim_non_neutral_model with codegen for improved speed
            [time,sim_grid_t] = sim_non_neutral_model(t_max,Nrows,Ncols,frac_prolif,rate_division,rate_division_fast,prob_duplicate,factor,im);

            if save_files == 1
                parsave([save_path,'time_',num2str(n_rep.','%04d'),'.mat'],time)
                parsave([save_path,'sim_grid_t_',num2str(n_rep.','%04d'),'.mat'],sim_grid_t)
            end
        end
    end
end

%% ANALYSIS
% Here we load each output and compute void sizes.
fold_paths = find_folders(work_path);
if ~ispc; slsh = '/'; else; slsh = '\'; end

for nfold = 1:length(fold_paths)
    sim_params = load([fold_paths{nfold},'sim_params.mat']); sim_params = sim_params.sim_params;
    nrepeats = sim_params('nrepeats');
    %     nrepeats = 10;

    factor = sim_params('factor');
    t_max = sim_params('t_max');
    N = ceil(t_max/factor);
    data = [];

    for nrep = 1:nrepeats

        if exist([fold_paths{nfold},'sim_grid_t_',num2str(nrep.','%04d'),'.mat']) ~= 0

            sim_grid_t = load([fold_paths{nfold},'sim_grid_t_',num2str(nrep.','%04d'),'.mat']);
            sim_grid_t = sim_grid_t.x;
            time = load([fold_paths{nfold},'time_',num2str(nrep.','%04d'),'.mat']);
            time = time.x;

            for npars = 1:N
                BW =sim_grid_t(:,:,npars)<2;
                stats = regionprops(BW,'area');
                data(nrep).time(npars,:) = time(npars);
                avg_area = mean([stats.Area]);
                data(nrep).nholes(npars,:) = length(stats);
                data(nrep).avgarea(npars,:) = avg_area;
                data(nrep).area{npars} = [stats.Area];
                nclones = length(unique(sim_grid_t(:,:,npars)))-1;
                data(nrep).nclones(npars,:) = nclones;
            end
        end
    end
    parsave([fold_paths{nfold},'data_holes.mat'],data)
end
%% UTILS
function fold_paths = find_folders(work_path)
% function save_path = create_save_folder(work_path,parmas)
dinfo = dir(work_path);
dinfo(ismember( {dinfo.name}, {'.', '..','.DS_Store'})) = [];
fold_paths{length(dinfo)} = [];

if ~ispc; slsh = '/'; else; slsh = '\'; end

if ~isempty(dinfo)
    n_strPadded = sprintf( '%05d', 1);
    save_path = [work_path,slsh,n_strPadded,slsh];

    for i = 1:length(dinfo)
        fold = dinfo(i).name;
        fold_paths{i} = [work_path,slsh,fold,slsh];
    end
end
end

function [exi,save_path] = create_save_folder(work_path,newparams)
exi = 0;
save_path =[];
dinfo = dir(work_path);
dinfo(ismember( {dinfo.name}, {'.', '..','.DS_Store'})) = [];

if ~ispc; slsh = '/'; else; slsh = '\'; end

if length(dinfo) == 0
    n_strPadded = sprintf( '%05d', 1);
    save_path = [work_path,slsh,n_strPadded,slsh];
    mkdir(save_path)
else
    for i = 1:length(dinfo)
        fold = dinfo(i).name;
        load([work_path,slsh,fold,slsh,'sim_params.mat']);
        ks = keys(sim_params);
        check_if_exists = [];
        for j = 1:length(ks)
            if isstr(sim_params(ks{j}))
                check_if_exists(j) = strcmp(sim_params(ks{j}),newparams(ks{j}));
            else
                check_if_exists(j) = all(sim_params(ks{j}) == newparams(ks{j}));
            end
        end

        if all(check_if_exists)
            disp(['Simulation data exists for these parameters: ',dinfo(i).name]);
            exi = 1;
            break
        end
    end

    n = 0;
    if exi == 0
        for i = 1:length(dinfo)
            n(i) = str2num(dinfo(i).name);
        end
        n_strPadded = sprintf( '%05d', max(n) + 1 );
        save_path = [work_path,slsh,n_strPadded,slsh];
        mkdir(save_path)
    end
end
end

function parsave(fname, x)
save(fname, 'x')
end

function parsave_sim_params(fname, sim_params)
save(fname, 'sim_params')
end