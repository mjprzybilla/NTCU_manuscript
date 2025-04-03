% Runtime can be reduced by compilation using codegen:
% codegen -args {300,100,100,0.3,0.3,0.1,0.1,0.5,0.5,0.1,0.25,1,zeros(100)} sim_non_neutral_w_driver
% and by adding #codegen pragma at the end of the function definition:
function [time,sim_grid_t] = sim_non_neutral_w_driver(t_max,Nrows,Ncols,nclones,frac_prolif,include_driver,p_imbalance,rate_division,factor,sim_grid) %#codegen pragma

% define savetime
N = ceil(t_max/factor)+1;
time = zeros(N,1);

% initialize arrays to store data 
inds_label_prolif = cell(N,1);
inds_unlabel_prolif = cell(N,1);
% drivers
inds_label_prolif_driver = cell(N,1);
inds_unlabel_prolif_driver = cell(N,1);

for i = 1:N % ensures arrays are full (for codegen)
    inds_label_prolif{i,1} = 0;
    inds_unlabel_prolif{i,1} = 0;
    inds_label_prolif_driver{i,1} = 0;
    inds_unlabel_prolif_driver{i,1} = 0;
end
sim_grid_t = zeros(Nrows,Ncols,N+1);
% cell types are assigned according to:
% sim_grid = ones(Nrows,Ncols); % unlabelled are 1;
% labelled > 1; unlab prolif = -1; lab prolif = -2; unlab driver = -3; lab driver = -4 


Ntotal = Nrows*Ncols;
all_inds = (1:Ntotal)';
active_inds = all_inds;
Nactive = Ntotal;

% generate initial condition with clones at the proximal side
[sim_grid,~,~,maxcloneind] = label_clusters_left(sim_grid,Nrows,nclones,frac_prolif);

% fraction of labeled and unlabeled cells that are atively proliferating
label_inds = find(sim_grid>1); % find all labeled cells
label_inds = intersect(active_inds,label_inds);
Nlabel = numel(label_inds); % number of labeled active cells

% if include_driver == 1 then one mutant clone is considered to have driver
% mutations
if include_driver>0
    label_inds_driver = maxcloneind;
    label_inds_prolif = setdiff(label_inds,label_inds_driver);

    Nlabel_driver = length(label_inds_driver);
    Nlabel_prolif = length(label_inds_prolif);
    sim_grid(label_inds_driver) = -4; % labeled driver
else  
    label_inds_driver = [];
    Nlabel_driver = 0;
    label_inds_prolif = label_inds;
    Nlabel_prolif = length(label_inds_prolif);
end

unlabeled_inds = find(sim_grid==1); % find all unlabeled cells
unlabeled_inds = intersect(active_inds,unlabeled_inds);
Nunlabeled = numel(unlabeled_inds); % number of unlabeled active cells

unlabel_inds_driver = [];
Nunlabel_driver = 0;

unlabeled_inds_prolif = []; % select unlabeled proliferative cells
Nunlabeled_prolif = 0;

%% RUN MODEL
t = 0;
n = 0;
avg_area = 0;
while t<=t_max && Nlabel<Ntotal && Nunlabeled<Ntotal
    all_inds_prolif = [label_inds_prolif;label_inds_driver;unlabeled_inds_prolif;unlabel_inds_driver];
    all_inds_prolif = sort(all_inds_prolif);
    all_inds_nonprolif = setdiff(active_inds,all_inds_prolif);

    if (Nlabel_prolif+Nunlabeled_prolif) ~= (numel(label_inds_prolif) + numel(unlabeled_inds_prolif)); disp('breaking2'); break; end;
    % write propensity functions for labelled, unlabelled, w/ driver mutations and long lived
    w0 = (Nlabel_prolif+Nunlabeled_prolif)*rate_division;
    w1 = (Nlabel_driver+Nunlabel_driver)*rate_division;
    w = w0+w1;
    r=rand();
    if w > 0
        if r<=w0/w % a (nondriver) mutant duplicates
            [sim_grid, label_inds_prolif, Nlabel_prolif, unlabeled_inds_prolif, Nunlabeled_prolif, label_inds_driver, Nlabel_driver, unlabel_inds_driver, Nunlabel_driver] = move_rand_drivers(sim_grid,Nrows,Ncols, Nlabel_prolif, Nunlabeled_prolif, label_inds_prolif, unlabeled_inds_prolif, Nlabel_driver, label_inds_driver, Nunlabel_driver, unlabel_inds_driver,all_inds_nonprolif,active_inds,p_imbalance);
        else % a driver mutant duplicates
            [sim_grid, label_inds_prolif, Nlabel_prolif, unlabeled_inds_prolif, Nunlabeled_prolif, label_inds_driver, Nlabel_driver, unlabel_inds_driver, Nunlabel_driver] = move_rand_drivers(sim_grid,Nrows,Ncols, Nlabel_driver, Nunlabel_driver, label_inds_driver, unlabel_inds_driver, Nlabel_prolif, label_inds_prolif, Nunlabeled_prolif, unlabeled_inds_prolif,all_inds_nonprolif,active_inds,p_imbalance);
        end
    end
    t = t - (1/w)*log(1-rand());

    % lets save intermediate states
    if t > n*factor && n<=N
        n = n + 1;
        time(n) = n*factor; % times
        sim_grid_t(:,:,n) = sim_grid; % full simulation grid
    end
end
end

%% FUNCTIONS

function [ind_closest,dmin] = find_closest(ind_ref,ind_others,Nrows,Ncols)
% find closest active sites to given site (including BCs)
[xi,yi] = ind2sub([Nrows,Ncols],ind_ref);
[xu,yu] = ind2sub([Nrows,Ncols],ind_others);

xu((xi-xu)>Nrows/2) = xu((xi-xu)>Nrows/2) + Nrows;
xu((xi-xu)<-Nrows/2) = xu((xi-xu)<-Nrows/2) - Nrows;

Du = pdist2([xi,yi],[xu,yu]);
Du(Du == 0) = nan;
logico =  (abs(Du - min(Du,[],'all'))<1e-4)';
Dumind = find(logico);

if nargin == 5
    dmin = min(Du,[],2);
    ind_closest = Dumind;
else
    if length(Dumind)>1
        ind_closest = randsample(Dumind,1);
    elseif length(Dumind)==1
        ind_closest = Dumind;
    else
        disp('ERROR: no sites found by find closest')
        ind_closest = 0;
    end
    dmin = min(Du,[],2);
    ind_closest = ind_closest(1);
end
end

function [sim_grid,Nlabel,label_inds,maxcloneind] = label_clusters_left(sim_grid,Ny,nclones,frac_lbl)
% label clones in the proximal side.
if nclones>0
    ind = Ny;
    n = 2:(nclones+1); % 1 is reserved for unlabelled; 
    [x,y] = ind2sub(size(sim_grid),1:numel(sim_grid(:)));

    a = exprnd(8,nclones,1);
    a = floor(Ny*a./sum(a)); a(a ==0) = 1;
    a = cumsum(a);
    a = [a-a(1)+1;Ny];
    clone_size = zeros(nclones,1);
    for i = 1:nclones
        if rand()<frac_lbl
        cent = mean(a(i):a(i+1));
        r = abs(a(i)-a(i+1))/2;

        dx = abs(x'-cent);
        dy = abs(y'-1);
        d = sqrt(dx.^2+dy.^2);
        sim_grid(d<=r) = n(i);
        clone_size(i) = sum(sim_grid(d<=r),'all');
        end
    end
    Nlabel = sum(sim_grid(:)>1);
    label_inds = find(sim_grid(:)>1);

    [maxval,maxi] = max(clone_size);
    maxcloneind = find(sim_grid == n(maxi));
else
    Nlabel = 0;
    label_inds = zeros(0,1);
    clone_size = 0;
    maxcloneind = [];
end
end

function [sim_grid,ind_lbl_prolif,n_lbl_prolif,ind_unlbl_prolif,n_unlbl_prolif,ind_lbl_driv,n_lbl_driv,ind_unlbl_driv,n_unlbl_driv] = move_rand_drivers(sim_grid,Nrows,Ncols,nactive1,nactive2,active_inds1,active_inds2,nother1,other_inds1,nother2,other_inds2,all_inds_nonprolif,active_inds,p_imbalance)
% duplicate a mutant considering advantage (p_imbalance) of cells with driver mutations
i = randsample(nactive1+nactive2,1);
if i <= nactive1
    ind = active_inds1(i);
    set = 1;
else
    ind = active_inds2(i-nactive1);
    set = 2;
end
indtype = sim_grid(ind);

% lets organise input data
if any(sim_grid(active_inds1) > 1) % active are proliferative: label -2, unlabel -1
    ind_lbl_prolif = active_inds1;
    n_lbl_prolif = nactive1; % = [other_inds1;other_inds2];
    ind_unlbl_prolif = active_inds2;
    n_unlbl_prolif = nactive2;

    ind_lbl_driv = other_inds1;
    n_lbl_driv = nother1;
    ind_unlbl_driv = other_inds2;
    n_unlbl_driv = nother2; % nactive1,nactive2,active_inds1,active_inds2,nother1,other_inds1,nother2,other_inds2,all_inds_prolif,all_inds_n
else % active are drivers
    ind_lbl_prolif = other_inds1;
    n_lbl_prolif = nother1; % = [other_inds1;other_inds2];
    ind_unlbl_prolif = other_inds2;
    n_unlbl_prolif = nother2;

    ind_lbl_driv = active_inds1;
    n_lbl_driv = nactive1;
    ind_unlbl_driv = active_inds2;
    n_unlbl_driv = nactive2; % nactive1,nactive2,active_inds1,active_inds2,nother1,other_inds1,nother2,other_inds2,all_inds_prolif,all_inds_n
end

all_prolif = [ind_lbl_prolif;ind_unlbl_prolif];
all_drivers = [ind_lbl_driv;ind_unlbl_driv];

% check if any of the closest cells is nonproliferative or driver
[iclose,dclose] = find_closest(ind,active_inds,Nrows,Ncols); % check 4 neighbours
inds = active_inds(iclose);
indNNnonprolif = inds(ismember(inds,all_inds_nonprolif));
nNNnonprolif = numel(indNNnonprolif);
indNNprolif = inds(ismember(inds,all_prolif));
nNNprolif = numel(indNNprolif);
indNNdriver = inds(ismember(inds,all_drivers));
nNNdriver = numel(indNNdriver);

if nNNnonprolif>0% closest cell in nonpoliferative
    if nNNnonprolif>1
        ind_closest = randsample(indNNnonprolif,1);
    else
        ind_closest = indNNnonprolif;
    end
    [ind_lbl_prolif,n_lbl_prolif,ind_unlbl_prolif,n_unlbl_prolif,ind_lbl_driv,n_lbl_driv,ind_unlbl_driv,n_unlbl_driv] = add_to_set(ind,ind_closest,ind_lbl_prolif,n_lbl_prolif,ind_unlbl_prolif,n_unlbl_prolif,ind_lbl_driv,n_lbl_driv,ind_unlbl_driv,n_unlbl_driv);
    sim_grid(ind_closest) = sim_grid(ind);
else % if neighbour is proliferative or driver
    ntot = nNNprolif+nNNdriver;
    wprolif = (1+p_imbalance/2)*nNNprolif;
    wdriver = (1-p_imbalance/2)*nNNdriver;

    w = wprolif + wdriver;
    if w>0
    r = rand();
    if r < wprolif/w % replace proliferative
        if nNNprolif > 1
            ind_closest = randsample(indNNprolif,1);
        else
            ind_closest = indNNprolif;
        end
        % disp('yup')
    else % replace driver
        if nNNdriver > 1
            ind_closest = randsample(indNNdriver,1);
        else
            ind_closest = indNNdriver;
        end
        % indtype
    end
    [ind_lbl_prolif,n_lbl_prolif,ind_unlbl_prolif,n_unlbl_prolif,ind_lbl_driv,n_lbl_driv,ind_unlbl_driv,n_unlbl_driv] = rm_from_set(ind_closest,ind_lbl_prolif,n_lbl_prolif,ind_unlbl_prolif,n_unlbl_prolif,ind_lbl_driv,n_lbl_driv,ind_unlbl_driv,n_unlbl_driv);     
    [ind_lbl_prolif,n_lbl_prolif,ind_unlbl_prolif,n_unlbl_prolif,ind_lbl_driv,n_lbl_driv,ind_unlbl_driv,n_unlbl_driv] = add_to_set(ind,ind_closest,ind_lbl_prolif,n_lbl_prolif,ind_unlbl_prolif,n_unlbl_prolif,ind_lbl_driv,n_lbl_driv,ind_unlbl_driv,n_unlbl_driv);
    sim_grid(ind_closest) = sim_grid(ind);
    end
end
end % move_rand_drivers

function [ind_lbl_prolif,n_lbl_prolif,ind_unlbl_prolif,n_unlbl_prolif,ind_lbl_driv,n_lbl_driv,ind_unlbl_driv,n_unlbl_driv] = add_to_set(ind_ref,ind_to_add,ind_lbl_prolif,n_lbl_prolif,ind_unlbl_prolif,n_unlbl_prolif,ind_lbl_driv,n_lbl_driv,ind_unlbl_driv,n_unlbl_driv)
% find to which cell compartment the ind_ref belong to and add ind_to_add to the
% corresponding dataset
if ismember(ind_ref,ind_lbl_prolif)
    ind_lbl_prolif = [ind_lbl_prolif;ind_to_add];
    n_lbl_prolif = n_lbl_prolif + 1;
elseif ismember(ind_ref,ind_unlbl_prolif)
    ind_unlbl_prolif = [ind_unlbl_prolif;ind_to_add];
    n_unlbl_prolif = n_unlbl_prolif + 1;
elseif ismember(ind_ref,ind_lbl_driv)
    ind_lbl_driv = [ind_lbl_driv;ind_to_add];
    n_lbl_driv = n_lbl_driv + 1;
elseif ismember(ind_ref,ind_unlbl_driv)
    ind_unlbl_driv = [ind_unlbl_driv;ind_to_add];
    n_unlbl_driv = n_unlbl_driv + 1;
end
end

function [ind_lbl_prolif,n_lbl_prolif,ind_unlbl_prolif,n_unlbl_prolif,ind_lbl_driv,n_lbl_driv,ind_unlbl_driv,n_unlbl_driv] = rm_from_set(ind_ref,ind_lbl_prolif,n_lbl_prolif,ind_unlbl_prolif,n_unlbl_prolif,ind_lbl_driv,n_lbl_driv,ind_unlbl_driv,n_unlbl_driv)
% find to which dataset the ind_ref belong to and remove ind_to_add from it
if ismember(ind_ref,ind_lbl_prolif)
    ind_lbl_prolif(ind_lbl_prolif == ind_ref) = [];
    n_lbl_prolif = n_lbl_prolif - 1;
elseif ismember(ind_ref,ind_lbl_driv)
    ind_lbl_driv(ind_lbl_driv == ind_ref) = [];
    n_lbl_driv = n_lbl_driv - 1;
end
end