pc = 1;

% set directory to save simulation results
work_path = 'D:\...';
% Save output files? toggle to 0 if not.
save_files = 0;

%% Set simulation parameters
% Number of mutant clones
nclones = 80; 
% fraction of mutant clones that have driver mutations
frac_mutant_w_driver = [0.2]; 

% if include_driver == 1 then one mutant clone is considered to have driver
% mutations
include_driver = 1;
% p_imbalance in [0,2]: 
% p_imbalance = 0: mutant nad mutant+drivers are equivalent.
% p_imbalance = 2: driver mutant replace neighbouring mutant with prob. 1
p_imbalance = [0.1];
% system size
Nrows = 350;
Ncols = 100; %2mm
rate_division = 0.1;
t_max = 100;
factor = 1; % data record intervals

param_combos = [];
save_paths = {};
ns = [];
npars = 0;

%% RUN simulations

im = ones(Nrows,Ncols);

[time,sim_grid_t] = sim_non_neutral_w_driver(t_max,Nrows,Ncols,nclones,frac_mutant_w_driver,include_driver,p_imbalance,rate_division,factor,im);

if save_files == 1
    % STORE params
    sim_params = containers.Map;
    sim_params('t_max') = t_max;
    sim_params('factor') = factor;
    sim_params('p_imbalance') = p_imbalance;
    sim_params('Nrows') = Nrows;
    sim_params('Ncols') = Ncols;
    sim_params('rate_division') = rate_division;
    sim_params('include_driver') = include_driver;
    sim_params('frac_prolif') = frac_mutant_w_driver;
    sim_params('nclones') = nclones;
    parsave_sim_params([save_path,'sim_params.mat'],sim_params)
    % save output
    parsave([work_path,'time_',num2str(n_rep.','%04d'),'.mat'],time)
    parsave([work_path,'sim_grid_t_',num2str(n_rep.','%04d'),'.mat'],sim_grid_t)
end

%% UTILS

function parsave(fname, x)
save(fname, 'x')
end

function parsave_sim_params(fname, sim_params)
save(fname, 'sim_params')
end