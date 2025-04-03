% simulation of neutral model
% compile with codegen for better performance
function [time,sim_grid_t,inds_label,inds_unlabel] = sim_neutral_model(t_max,Nrows,Ncols,rate_division,rate_division_label,p,factor,sim_grid) %#codegen pragma
N = ceil(t_max/factor)+1;
time = zeros(N,1);

inds_label = cell(N,1);
inds_unlabel = cell(N,1);
for i = 1:N % ensures arrays are full (for codegen)
    inds_label{i,1} = 0;
    inds_unlabel{i,1} = 0;
end
sim_grid_t = zeros(Nrows,Ncols,N+1);
% cell types are encoded as
% sim_grid = ones(Nrows,Ncols); % unlabelled are 1; labelled > 1
Ntotal = Nrows*Ncols;
active_inds = (1:Ntotal)';
Nactive = Ntotal;

label_inds = find(sim_grid);
Nlabel = numel(label_inds);

Nunlabeled = Nactive - Nlabel;

unlabeled_inds = setdiff(active_inds,label_inds);
t = 0;
n = 0;
avg_area = 0;
while t<=t_max
    % write propensity functions for labelled and unlabelled cells
    w0 = Nlabel*rate_division_label;
    w1 = Nunlabeled*rate_division;
    w = w0+w1;
    r=rand();
    if w > 0
        if r<=w0/w
            [sim_grid,Nlabel,label_inds,Nunlabeled,unlabeled_inds] = move_rand(sim_grid,Nrows,Ncols,Nlabel,label_inds,Nunlabeled,unlabeled_inds,p);
        else
            [sim_grid,Nunlabeled,unlabeled_inds,Nlabel,label_inds] = move_rand(sim_grid,Nrows,Ncols,Nunlabeled,unlabeled_inds,Nlabel,label_inds,p);
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

function [sim_grid,nactive,active_inds,nother,other_inds,inds] = move_rand(sim_grid,Nrows,Ncols,nactive,active_inds,nother,other_inds,p)
% the code allows for introducing a bias in the direction to which a cell
% divides, this is included as:
% p: prob of moving to the right
% q: prob of moving anywhere else
% such that p+3q = 1 [q = (1-p)/3]
x = 3*p/(1-p); % relation between moving to the right or elsewhere
% then we perform an action
i = randsample(numel(active_inds),1);
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
    % check boundary conditions
    if any(moves(:,2) == 1) % if any move is right (+1 column), then the move is biased according to p and q
        r = rand();
        % prob of moving to the right (pp: prob of p), normalized to the
        % number of allowed moves.
        pp = 1/(1+(nm-1)/x);
        if r<pp % move to the right
            move = [0,1];
        else % move elsewhere
            imove = randsample(nm-1,1);
            move = moves(imove,:);
        end
    else % if move right is not an option, then choose randomly
        imove = randsample(nm,1);
        move = moves(imove,:);
    end
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
end % move_rand