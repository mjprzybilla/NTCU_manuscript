%% STOCHASTIC SIMULATIONS OF THE VOTER MODEL ON A SQUARE LATTICE
% see Supplementary Note
%
% INPUTS
% t_max : maximum simulation time (days)
%
% label_fraction : fraction of initially labelled cells
% Nrows : system length
% Ncols : system width
% rate_division : rate of division of unlabelled cells
% rate_division_label : rate of division of labelled cells
% prob_duplicate : probability of dubplication in [0,1]
%
% record_interval : save state time interval
%
%--------------------------------------------------------------------------
function [time,sim_grid_t,inds_label,inds_unlabel] = script_voter_model(t_max,Nrows,Ncols,label_fraction,rate_division,rate_division_label,prob_duplicate,factor) %#codegen pragma
N = ceil(t_max/factor)+1;
time = zeros(N,1);

inds_label = cell(N,1);
inds_unlabel = cell(N,1);
% inds_long_lived = cell(N,1);

for i = 1:N % ensures arrays are full (for codegen compilation)
    inds_label{i,1} = 0;
    inds_unlabel{i,1} = 0;
end
sim_grid_t = zeros(Nrows,Ncols,N+1);

sim_grid = ones(Nrows,Ncols); % unlabelled are 1; labelled > 1
Ntotal = Nrows*Ncols;
active_inds = (1:Ntotal)';

[sim_grid,Nlabel,label_inds] = label_random(sim_grid,Ntotal,active_inds,label_fraction);

Nunlabeled = Ntotal - Nlabel;
unlabeled_inds = setdiff(active_inds,label_inds);

t = 0;
n = 0;
while t<=t_max
    % write propensity functions for labelled and unlabelled cells
    w0 = Nlabel*rate_division_label;
    w1 = Nunlabeled*rate_division;
    w = w0+w1;
    r=rand();
    if w > 0
        if r<=w0/w % A -> A+A
            rd = rand();
            if rd < prob_duplicate
                [sim_grid,Nlabel,label_inds,Nunlabeled,unlabeled_inds] = move_rand(sim_grid,Nrows,Ncols,Nlabel,label_inds,Nunlabeled,unlabeled_inds);
            else % A - > 0
                [sim_grid,Nlabel,label_inds,Nunlabeled,unlabeled_inds] = kill_rand(sim_grid,Nrows,Ncols,Nlabel,label_inds,Nunlabeled,unlabeled_inds);
            end
        else %if r<=(w0+w1)/w 
            rd = rand();
            if rd < prob_duplicate % B -> B+B
                [sim_grid,Nunlabeled,unlabeled_inds,Nlabel,label_inds] = move_rand(sim_grid,Nrows,Ncols,Nunlabeled,unlabeled_inds,Nlabel,label_inds);
            else % B -> 0
                [sim_grid,Nunlabeled,unlabeled_inds,Nlabel,label_inds] = kill_rand(sim_grid,Nrows,Ncols,Nunlabeled,unlabeled_inds,Nlabel,label_inds);
            end
     
        end
    end
    t = t - (1/w)*log(1-rand());
    if t > n*factor && n<=N
        n = n + 1;
        time(n) = n*factor;
        inds_label{n} = label_inds;
        inds_unlabel{n} = unlabeled_inds;
        sim_grid_t(:,:,n) = sim_grid;
    end
end
end

function [sim_grid,nactive,active_inds,nother,other_inds,inds] = move_rand(sim_grid,Nrows,Ncols,nactive,active_inds,nother,other_inds)
% select an active cell
i = randsample(length(active_inds),1);
ind = active_inds(i);
% find indices of the active cell
[row,col] = ind2sub([Nrows,Ncols],ind);
% initialize "move" array with all movement possibilities for the active cell
moves_tmp = [-1,0;1,0;0,-1;0,1];

moves = [];
% find neighbours of the active particle
neighs = [row,col] + moves_tmp;

moves = moves_tmp;

% check boundaries

pos_new_tmp = [row,col] + moves;
% if new pos is outside of system, then do nothing
[i,~] = find(pos_new_tmp(:,1)<1 | pos_new_tmp(:,2)<1 | pos_new_tmp(:,1) > Nrows | pos_new_tmp(:,2) > Ncols);
moves(i,:) = [];

% number of allowed moves
nm = size(moves,1);
if nm>1
    imove = randsample(nm,1);
    move = moves(imove,:);
elseif nm == 1
    move = moves;
else
    move = [0,0];
end
pos_new = [row,col] + move;

indactive = sub2ind([Nrows,Ncols],pos_new(1),pos_new(2));
sim_grid(indactive) = sim_grid(ind);
if  all(active_inds ~= indactive) && any(active_inds == ind) % if 
    nactive = nactive+1;
    active_inds = [active_inds;indactive];
    nother = nother - 1;
    other_inds(other_inds == indactive) = [];
end
inds = [indactive,ind];
end % end move_rand

function [sim_grid,Nlabel,label_inds] = label_random(sim_grid,N,active_inds,label_fraction)
if label_fraction>0
    ind = randsample(length(active_inds),ceil(N*label_fraction));
    n = randperm(length(ind)+1); n(n==1) = []; % 1 is reserved for unlabelled;
    for i = 1:length(ind)
        sim_grid(active_inds(ind(i))) = n(i);
    end
    Nlabel = length(ind);
    label_inds = sort(active_inds(ind));
else
    Nlabel = 0;
    label_inds = zeros(0,1);
end
end % end label_random

function [sim_grid,nactive,active_inds,nother,other_inds,inds] = kill_rand(sim_grid,Nrows,Ncols,nactive,active_inds,nother,other_inds)
% select an active cell
i = randsample(length(active_inds),1);
ind = active_inds(i);
% find indices of the active cell
[row,col] = ind2sub([Nrows,Ncols],ind);
% initialize "move" array with all movement possibilities for the active cell
moves_tmp = [-1,0;1,0;0,1;0,-1];

moves = [];
% find neighbours of the active particle
neighs = [row,col] + moves_tmp;
moves = moves_tmp;

% check boundaries
pos_new_tmp = [row,col] + moves;
% if new pos is outside of system, then do nothing
[i,~] = find(pos_new_tmp(:,1)<1 | pos_new_tmp(:,2)<1 | pos_new_tmp(:,1) > Nrows | pos_new_tmp(:,2) > Ncols);
moves(i,:) = [];

% number of allowed moves
nm = size(moves,1);
if nm>1
    imove = randsample(nm,1);
    move = moves(imove,:);
elseif nm == 1
    move = moves;
else
    move = [0,0];
end
pos_new = [row,col] + move;

indactive = sub2ind([Nrows,Ncols],pos_new(1),pos_new(2));
sim_grid(ind) = sim_grid(indactive); % replace sellected cell by neighbour
if  ~any(active_inds == ind) && any(active_inds == indactive)
    nactive = nactive+1;
    active_inds = [active_inds;indactive];
    nother = nother - 1;
    other_inds(other_inds == indactive) = [];
else

end
inds = [indactive,ind];
end % end kill_rand