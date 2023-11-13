% Plot Figure 2
% Loads experimental void sizes and void sizes obtained from numerical
% simulations of the Voter model and plots the cumulative distributions for
% each conditions at t= 168 days.

dset = {};
dset{end+1} = '24 weeks - control'; 
dset{end+1} = '24 weeks - short term'; 
dset{end+1} = '24 weeks - treated'; 

% Load experimental void sizes
load('data_expvoid_size.mat') % areas
% Load simulated void sizes
load('data_votervoid_size.mat') % areas

regs = {'dorsal','ventral'};

figure;
min_hole_size = 3;

for nset = 1:length(dset) % run for three conditions
    lgds = [];
    for nsamp = 1:size(areas{nset},1) % for all sample in condition

    subplot(length(dset),1,nset)
    csa = areas{nset}{nsamp,1};
    csa(csa<min_hole_size) = [];

    plot_cdf(csa,'-',[0 0 0]+0.5*(nsamp-1));
    hold on
  
    csb = areas{nset}{nsamp,2};
    csb(csb<min_hole_size) = [];

    plot_cdf(csb,'-',[0 0 0]+0.5*(nsamp-1));
    xlim([1,10000])
    end
    title(dset{nset})
    
    xlabel('Void size (cells)')
    ylabel('Cumulative probability')
    ylim([10e-4,1])
    set(gca,'yScale','log')
    set(gca,'xScale','log')

    x = [10,100]; y = 7./x;
    plot(x,y,'k-','linewidth',2)

    legend([regs,regs,'powerlaw'])
    legend boxoff
end
set(gcf,'Color','w')

% plot voter model on top

time = 168;
for nset = 1:length(dset)
    subplot(length(dset),1,nset)

    hold on
        nrepeats = length(voter_voids);

        cdfs = [];
        for nrep = 1:nrepeats
            avg_area = voter_voids(nrep).avgarea(time,:);
            hole_areas = voter_voids(nrep).area{time};

            hole_areas(hole_areas<min_hole_size) = [];
            [cdf,bincents] = calc_cdf(hole_areas);
            cdfs(nrep,:) = cdf;
        end
        m = mean(cdfs);
        stdev = std(cdfs,1);
        plot(bincents,m,'-.k')
        plot(bincents,m+stdev,':k','handlevisibility','off')
        plot(bincents,m-stdev,':k','handlevisibility','off')
        set(gca,'xscale','log')
        set(gca,'yscale','log')
end

% functions
% ---------------------------------------------------------------------
function [cdf,bincents] = calc_cdf(areas,min_size)
if nargin == 2
binedges = (min_size-0.5):1:max(areas)+0.5;
else
binedges = 0.5:1:10^6+0.5;
end
bincents = (binedges(2:end) + binedges(1:end - 1))/2;
% compute pdf
Nnorm = histcounts(areas,binedges,'Normalization','pdf');
% compute cdf
cdf = 1 - cumsum(Nnorm)./sum(Nnorm);
end

function plot_cdf(areas,line,col)
[cdf,bincents] = calc_cdf(areas);
plot(bincents,cdf,line,'color',col,'LineWidth',1)
end