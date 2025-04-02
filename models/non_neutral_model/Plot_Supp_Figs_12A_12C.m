%% PLOT SUPPLEMENTARY FIGURES 12A and 12C
figure;
load('data/sim_output_350x100_n80_mutant_clones.mat')
  
subplot(1,3,1)
imagesc(x(:,:,1)); axis equal; axis tight
% adjust colormap
cmap = gray(10); 
colormap(cmap);

subplot(1,3,2)
imagesc(x(:,:,400)); axis equal; axis tight
% adjust colormap
cmap = gray(256); 
colormap(cmap); % Apply the custom colormap.

subplot(1,3,3)
imagesc(x(:,:,800)); axis equal; axis tight
% adjust colormap
cmap = gray(256); 
colormap(cmap); 