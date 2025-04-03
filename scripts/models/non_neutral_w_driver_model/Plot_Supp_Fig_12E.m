%% PLOT SUPPLEMENTARY FIGURE 12E
figure;
load('sim_output_350x100_n16_mutant_n1_w_driver_mutation.mat')
sim = x;

% prepare data to plot
x = zeros(size(sim));
x(sim>1)=2; % sim values >1 correspond to mutant clones
x(sim==-4) = 5; % sim values equal to -4 correspond to mutant with driver mutations
   
subplot(1,3,1)
imagesc(x(:,:,1)); axis equal; axis tight
% adjust colormap
cmap = gray(256); cmap(9,:) = [1,0,0];
colormap(cmap);

subplot(1,3,2)
imagesc(x(:,:,400)); axis equal; axis tight
% adjust colormap
cmap = gray(256); cmap(9,:) = [1,0,0];
colormap(cmap); % Apply the custom colormap.

subplot(1,3,3)
imagesc(x(:,:,800)); axis equal; axis tight
% adjust colormap
cmap = gray(256); cmap(9,:) = [1,0,0];
colormap(cmap); 