%% PLOT SUPPLEMENTARY FIGURES 12D

load('data/clone_count.mat')

figure;
hold on

lbls = {'16 initial clones','40 initial clones','80 initial clones'};
for n = 1:3
    cs = [];
    plot(1:100,mean(clone_count{n},2))
end

xlabel('Distance from proximal side (sites - cell sizes)')
ylabel('Number of clones')
legend(lbls)

box on

