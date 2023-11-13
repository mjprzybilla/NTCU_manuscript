% Script to run simulations of the neutral voter model dynamics.
% For efficiency, generate C-compiled mex file:
%% 1) RUN:
codegen -args {300,100,100,0.3,0.1,0.1,0.5,1} script_voter_model

%% 2) SET INPUTS
%
t_max = 170; % maximum simulation time (days)
label_fraction = 0.45; % fraction of initially labelled cells
% Each voxel represents a cell:
Nrows = 250; % system length
Ncols = 750; % system width
prob_duplicate = 1; % probability of dubplication in [0,1]
rate_division = 3/110;%1/11; % 1/days
rate_division_fast = rate_division; % 1/days
%    
record_interval = 1; % save state time interval (days)

%% 3) RUN SIMULATIONS
[time,sim_grid_t,inds_label,inds_unlabel] = script_voter_model_mex(t_max,Nrows,Ncols,label_fraction,rate_division,rate_division_fast,prob_duplicate,record_interval);

%% 4) PLOT TIMELAPSE
ns = [];
figure;
subplot(3,1,1)
imagesc(sim_grid_t(:,:,1)); title('Time: 1 day'); axis off; axis equal
for i = 1:size(sim_grid_t,3)-1
   subplot(3,1,2)
   imagesc(sim_grid_t(:,:,i)); title({'Clones',['Time: ',num2str(i),' days']}); axis off
   axis equal
   subplot(3,1,3)
   imshow(sim_grid_t(:,:,i)>1); title({'Labelled/unlabelled',['Time: ',num2str(i),' days']}); axis off
   axis equal
   pause(0.1)
end

%% 5) PLOTE VOID SIZE DISTRIBUTION FOR A SINGLE REALIZATION
ns = [];
figure;
for i = 1:size(sim_grid_t,3)-1
    BW =sim_grid_t(:,:,i)<2;
    stats = regionprops(BW,'area');
    areas = [stats.Area];
    [cdf,bincents] = plot_cdf(areas);
    set(gca,'yscale','log'); box on;
    ylabel('Cumulative probability')
    xlabel('Void size (cells)')
    title(['Time: ',num2str(i),' days'])
    xlim([0,200])
    pause(0.1)
end

function [cdf,bincents] = plot_cdf(areas)
binedges = 0.5:1:max(areas)+0.5;
bincents = (binedges(2:end) + binedges(1:end - 1))/2;
% compute pdf
Nnorm = histcounts(areas,binedges,'Normalization','pdf');
% compute cdf
cdf = 1 - cumsum(Nnorm)./sum(Nnorm);
plot(bincents,cdf,'-k')
end