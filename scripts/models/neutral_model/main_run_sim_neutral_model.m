%% This script runs simulations of the neutral model, used to produce
% to reproduce the theory curve shown in Figure 2F
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
t_max = 17;
factor = 1; % record time
Nrows = 100;
Ncols = 300;
rate_divisions = 0.05; 
rate_divisions_fast = rate_divisions;
probs_duplicate = 1;
ps = 0.25; % prob of moving to the right: 0.25 is unbiased rw

param_combos = [];
save_paths = {};
ns = [];
npars = 0;
if run_sim == 1
    for pdup = 1:length(probs_duplicate)
        prob_duplicate = probs_duplicate(pdup);
        for rdiv = 1:length(rate_divisions)
            rate_division = rate_divisions(rdiv);
            for nrdf = 1:length(rate_divisions_fast)
                rate_division_fast = rate_divisions_fast(nrdf);
                for np = 1:length(ps)
                    p = ps(np);
                    npars = npars + 1;
                    param_combos.rate_division(npars) = rate_division;
                    param_combos.rate_division_fast(npars) = rate_division_fast;
                    param_combos.p(npars) = p;
                    param_combos.prob_duplicate(npars) = prob_duplicate;

                    sim_params = containers.Map;
                    sim_params('nrepeats') = nrepeats;
                    sim_params('t_max') = t_max;
                    sim_params('factor') = factor;
                    sim_params('Nrows') = Nrows;
                    sim_params('Ncols') = Ncols;
                    sim_params('rate_division') = rate_division;
                    sim_params('rate_division_fast') = rate_division_fast;
                    sim_params('p') = p;
                    sim_params('prob_duplicate') = prob_duplicate;

                    [exi,save_path] = create_save_folder(work_path,sim_params);
                    save_paths{npars} = save_path;
                    if exi == 0
                        if save_files == 1
                            save([save_path,'sim_params.mat'],'sim_params')
                        end
                        ns(end+1) = npars;
                    end
                    %     end
                    % end
                end
            end
        end
    end
end

if run_sim == 1
    for np=1:length(ns)
        n = ns(np);
        rate_division = param_combos.rate_division(n);
        rate_division_fast = param_combos.rate_division_fast(n);
        p = param_combos.p(n);
        prob_duplicate = param_combos.prob_duplicate(n);
        sim_params = containers.Map;
        sim_params('nrepeats') = nrepeats;
        sim_params('t_max') = t_max;
        sim_params('factor') = factor;
        sim_params('Nrows') = Nrows;
        sim_params('Ncols') = Ncols;
        sim_params('rate_division') = rate_division;
        sim_params('rate_division_fast') = rate_division_fast;
        % sim_params('rate_decay') = rate_decay;
        sim_params('p') = p;
        sim_params('prob_duplicate') = prob_duplicate;
        save_path = save_paths{n};
        if save_files == 1
            parsave_sim_params([save_path,'sim_params.mat'],sim_params)
        end
 
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

            [time,sim_grid_t] = sim_neutral_model(t_max,Nrows,Ncols,rate_division,rate_division_fast,p,factor,im);
                             
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
   factor = sim_params('factor');
    t_max = sim_params('t_max');
    N = ceil(t_max/factor);
    data = [];
    nrep = 0;
    for nn = 1:nrepeats
        nrep = nrep+1
        if exist([fold_paths{nfold},'sim_grid_t_',num2str(nn.','%04d'),'.mat'],'file') ~= 0

            sim_grid_t = load([fold_paths{nfold},'sim_grid_t_',num2str(nn.','%04d'),'.mat']);
            sim_grid_t = sim_grid_t.x;
            time = load([fold_paths{nfold},'time_',num2str(nn.','%04d'),'.mat']);
            time = time.x;

            for npars = 1:N
                BW =sim_grid_t(:,:,npars)<2;
                stats = regionprops(BW,'area');
                data(nrep).time(npars,:) = time(npars);
                avg_area = mean([stats.Area]);
                data(nrep).nholes(npars,:) = length(stats);
                data(nrep).avgarea(npars,:) = avg_area;
                data(nrep).area{npars} = [stats.Area];
            end
        end
    end
    parsave([fold_paths{nfold},'data_holes.mat'],data)
end

%% UTILS
function fold_paths = find_folders(work_path)
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
% function save_path = create_save_folder(work_path,parmas)
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

function Im = apply_blur(Im,sigma)
Im = imgaussfilt(Im,sigma);
end

function im=im_2_bw(im)
T = graythresh(im);
im = im2bw(im,T);
end