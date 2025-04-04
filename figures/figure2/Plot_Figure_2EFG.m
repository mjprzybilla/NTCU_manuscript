%% Plot Figure 2 E,F,G
% directory with void data

root_path = 'data';

% datasets
dset = {}; lstyle = {};
dset{end+1} = 'Tam Ctrl (4 days)'; lstyle{end+1} = '--';
dset{end+1} = 'Tam Ctrl (24 weeks)'; lstyle{end+1} = '-';
dset{end+1} = 'Tam +NTCU (24 weeks)'; lstyle{end+1} = '.-';
% order
regs = {'dorsal','ventral'};
% load data
load([root_path,'/area_hole_cells.mat'])
hole_sizes = areas;

figure;
min_hole_size = 2;
all_hole_areas_exp = {}; 
cdfs_exps = {};
for nset = 1:3
    lgds = [];
    ccdfsa = [];
    ccdfsb = [];
    for nsamp = 1:size(areas{nset},1)

    csa = hole_sizes{nset}{nsamp,1};

    csa(csa<min_hole_size) = [];

%     histogram(cs,0:100:10^5,'normalization','pdf');
% csa = csa./mean(csa);
    binedges = min_hole_size-0.5:1:10^5+0.5; bincents = (binedges(2:end) + binedges(1:end - 1))/2; 
    Nnorm = histcounts(csa,binedges,'Normalization','pdf');
    cdf = cumsum(Nnorm)./sum(Nnorm);
    ccdf = 1 - cdf;

    ccdfsa(nsamp,:) = ccdf;

    csb = hole_sizes{nset}{nsamp,2};
    csb(csb<min_hole_size) = [];
    binedges = min_hole_size-0.5:1:10^5+0.5; bincents = (binedges(2:end) + binedges(1:end - 1))/2; 
    Nnorm = histcounts(csb,binedges,'Normalization','pdf');
    cdf = cumsum(Nnorm)./sum(Nnorm);
    ccdf = 1 - cdf;

    ccdfsb(nsamp,:) = ccdf;
    
    end
    cdfs_exps{nset,1} = ccdfsa;
    cdfs_exps{nset,2} = ccdfsb;
    subplot(1,3,nset)
    hold on
    m = mean(ccdfsa);
    stdev = std(ccdfsa,1);
    n = size(ccdfsa,1);
    minvals = m-stdev;
    plot(bincents,m,'-k')
    plot(bincents,m+stdev,'--k','handlevisibility','off')
    plot(bincents,minvals,'--k','handlevisibility','off')
    box on

    m = mean(ccdfsb);
    stdev = std(ccdfsb,1);
    plot(bincents,m,'Color',[0.5,0.5,0.5])
    plot(bincents,m+stdev,':','Color',[0.5,0.5,0.5],'handlevisibility','off')
    plot(bincents,m-stdev,':','Color',[0.5,0.5,0.5],'handlevisibility','off')
    set(gca,'xscale','log')
    set(gca,'yscale','log')

    xlim([1,100000])
    title(dset{nset})
    legend(regs)
    legend boxoff
    
    xlabel('Void size (cells)')
    ylabel('Cumulative probability')
    ylim([10e-4,1])
    set(gca,'yScale','log')
    set(gca,'xScale','log')

    if nset<3
        x = [70,150]; y = 9./x;
        plot(x,y,'k-','linewidth',2,'handlevisibility','off')
    end

    if nset == 2
        plot_theory([root_path,'/neutral_sim_params.mat'],[root_path,'/neutral_data_holes.mat'],93)
    elseif nset ==3
        plot_theory([root_path,'/non_neutral_sim_params.mat'],[root_path,'/non_neutral_data_holes.mat'],130)

    end
end
set(gcf,'Color','w')

function plot_theory(sim_params_file,data_holes_file,tplot)
        min_hole_size = 2;
        sim_params = load(sim_params_file); sim_params = sim_params.sim_params;
        data_holes = load(data_holes_file); data_holes = data_holes.x;

        nrepeats = sim_params('nrepeats');
        factor = sim_params('factor');
        t_max = sim_params('t_max');
        N = ceil(t_max/factor);

        cdfs = [];
        all_hole_areas = [];
        for nrep = 1:nrepeats
            hole_areas = data_holes(nrep).area{tplot};
            hole_areas(hole_areas<min_hole_size) = [];
            [cdf,bincents] = compute_ccdf(hole_areas,min_hole_size);
            cdfs(nrep,:) = cdf;
        end
        m = mean(cdfs);
        stdev = std(cdfs,1);
        tplot = size(cdfs,1);
        minvals = m-stdev;
        hold on
        plot(bincents,m)
        plot(bincents,m+stdev,'-k','handlevisibility','off')
        plot(bincents,m-stdev,'-k','handlevisibility','off')
        hold off
end

function [cdf,bincents] = compute_ccdf(clone_sizes,min_hole_size)
binedges = min_hole_size-0.5:1:5*10^5;
bincents = (binedges(2:end) + binedges(1:end - 1))/2;
% compute pdf
Nnorm = histcounts(clone_sizes,binedges,'Normalization','pdf');
% compute cdf
cdf = 1 - cumsum(Nnorm)./sum(Nnorm);
end