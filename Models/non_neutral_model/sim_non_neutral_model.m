% simulation of non neutral model
% codegen -args {300,100,100,0.3,0.3,0.1,0.1,0.5,0.5,0.1,0.25,1,zeros(100)} sim_lattice_trachea_aw_eden_rev1
% sim_lattice_hierarchy_aw %#codegen pragma
function [time,sim_grid_t] = sim_non_neutral_model(t_max,Nrows,Ncols,frac_prolif,rate_division,rate_division_label,prob_duplicate,factor,sim_grid) %#codegen pragma

N = ceil(t_max/factor)+1;
time = zeros(N,1);

inds_label = cell(N,1);
inds_unlabel = cell(N,1);
inds_long_lived = cell(N,1);
inds_label_prolif = cell(N,1);
inds_unlabel_prolif = cell(N,1);

for i = 1:N % ensures arrays are full (for codegen)
    %     area{i,1} = 0;
    inds_label{i,1} = 0;
    inds_unlabel{i,1} = 0;
    inds_long_lived{i,1} = 0;
end
sim_grid_t = zeros(Nrows,Ncols,N+1);

% sim_grid = ones(Nrows,Ncols); % unlabelled are 1; long lived 0; labelled > 1
Ntotal = Nrows*Ncols;
all_inds = (1:Ntotal)';

Nlong = 0;
inds_long = zeros(0,1);
active_inds = all_inds;

Nactive = Ntotal-Nlong;

label_inds = find(sim_grid>1); % find all labeled cells
label_inds = intersect(active_inds,label_inds);
Nlabel = numel(label_inds); % number of labeled active cells

in = randperm(Nlabel);
in(floor(Nlabel*frac_prolif)+1:end) = [];

label_inds_prolif = label_inds(in); % select labeled proliferative cells
Nlabel_prolif = length(in);
% of those labelled, only a fraction are proliferative (i.e. mutant)

unlabeled_inds = find(sim_grid==1); % find all labeled cells
unlabeled_inds = intersect(active_inds,unlabeled_inds);
Nunlabeled = numel(unlabeled_inds); % number of labeled active cells

in = randperm(Nunlabeled);

in(floor(Nunlabeled*frac_prolif)+1:end) = [];
unlabeled_inds_prolif = unlabeled_inds(in); % select unlabeled proliferative cells
Nunlabeled_prolif = length(in);
% of those unlabelled, only a fraction are proliferative (i.e. mutant)
%%
t = 0;
n = 0;
avg_area = 0;
while t<=t_max && Nlabel<Ntotal && Nunlabeled<Ntotal
    all_inds_prolif = [label_inds_prolif;unlabeled_inds_prolif];
    all_inds_prolif = sort(all_inds_prolif);
    all_inds_nonprolif = setdiff(active_inds,all_inds_prolif);

    w0 = Nlabel_prolif*rate_division_label;
    w1 = Nunlabeled_prolif*rate_division;
    w = w0+w1;
    r=rand();
    if w > 0
        if r<=w0/w
            rd = rand();
            if rd < prob_duplicate
                % 1
                [sim_grid,Nlabel_prolif,label_inds_prolif,Nunlabeled_prolif,unlabeled_inds_prolif] = move_rand(sim_grid,Nrows,Ncols,Nlabel_prolif,label_inds_prolif,Nunlabeled_prolif,unlabeled_inds_prolif,all_inds_prolif,all_inds_nonprolif);
            else
                if Nlabel_prolif>1
                    il = randsample(length(label_inds_prolif),1);
                else
                    il = 1;
                end
                inds_long = [inds_long;label_inds_prolif(il)];
                Nlong = Nlong + 1;
                label_inds_prolif(il) = [];
                Nlabel_prolif = Nlabel_prolif - 1;
            end
        else
            rd = rand();
            if rd < prob_duplicate
                [sim_grid,Nunlabeled_prolif,unlabeled_inds_prolif,Nlabel_prolif,label_inds_prolif] = move_rand(sim_grid,Nrows,Ncols,Nunlabeled_prolif,unlabeled_inds_prolif,Nlabel_prolif,label_inds_prolif,all_inds_prolif,all_inds_nonprolif);
            else
                if Nunlabeled_prolif>1
                    il = randsample(length(unlabeled_inds_prolif),1);
                else
                    il = 1;
                end
                inds_long = [inds_long;Nunlabeled_prolif(il)];
                Nlong = Nlong + 1;
                unlabeled_inds_prolif(il) = [];
                Nunlabeled_prolif = Nunlabeled - 1;
                %             end
            end
        end
        t = t - (1/w)*log(1-rand());

        if t > n*factor && n<=N
            n = n + 1;
            time(n) = n*factor;
            sim_grid_t(:,:,n) = sim_grid;
        end
    end

end
end

%% FUNCTIONS

function [ind_closest,dmin] = find_closest(ind_ref,ind_others,Nrows,Ncols)
[xi,yi] = ind2sub([Nrows,Ncols],ind_ref);
[xu,yu] = ind2sub([Nrows,Ncols],ind_others);

Du = pdist2([xi,yi],[xu,yu]);
Du(Du == 0) = nan;
% min(Du)
logico =  (abs(Du - min(Du,[],'all'))<1e-4)';
Dumind = find(logico);

if length(Dumind)>1
    ind_closest = randsample(Dumind,1);
elseif length(Dumind)==1
    ind_closest = Dumind;
else
    disp('ERROR!!')
    ind_closest = 0;
end
dmin = min(Du,[],2);
ind_closest = ind_closest(1);
end

function [sim_grid,nactive,active_inds,nother,other_inds] = move_rand(sim_grid,Nrows,Ncols,nactive,active_inds,nother,other_inds,all_inds_prolif,all_inds_nonprolif)
% select an active cell
i = randsample(length(active_inds),1);
ind = active_inds(i);
% check if any of the closest cells is nonproliferative
% [iclose,dclose] = find_closest(ind,ind_non_lbl,Nrows,Ncols); % find closest 
if isempty(all_inds_nonprolif)
    dclose1 = inf;
    iclose1 = 0;
else
    [iclose1,dclose1] = find_closest(ind,all_inds_nonprolif,Nrows,Ncols); % find closest non prolif
end
if isempty(all_inds_prolif)
    dclose2 = inf;
    iclose2 = 0;
else
    [iclose2,dclose2] = find_closest(ind,all_inds_prolif,Nrows,Ncols); % find closest prolif
end

if dclose1<=dclose2 % closest cell in nonpoliferative
    indactive = all_inds_nonprolif(iclose1);
    active_inds = [active_inds;indactive];
    nactive = nactive+1;
else
    indactive = all_inds_prolif(iclose2);
    a = sim_grid(indactive)>1;
    b = sim_grid(ind)>1;
    if a ~= b
        nactive = nactive+1;
        active_inds = [active_inds;indactive];
        other_inds(other_inds == indactive) = [];
        nother = nother - 1;
    end
end
sim_grid(indactive) = sim_grid(ind);
end % move_rand

