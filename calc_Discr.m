function [d] = calc_Discr(x, subj_id)
% Performs Discriminability calculation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   REQUIRED INPUTS
%        x                    efc vectorized x scans
%        subj_id              Index of subjects per scan
%        distance metric      'euc', ... (default = Euclidean distance)
%
%   OUTPUTS
%        d                    Discriminability (Bridgeford et al., 2020)
%
%   Example:
%        [d]=calc_Discr(x, subj_id, 'euc');
%
%   References:
%        If you use this script, please cite:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% calculate discriminability based on Euclidean distance
N = length(subj_id);
dist_mtrx = zeros(N, N);% NxN of distance between scans
for dist = 1:size(x, 2)
    dist_mtrx(:, dist) = sum(abs(x - x(:, dist)));
end

nsubj = max(subj_id);
nrep = sum(subj_id == 1);

d0 = repmat(1:nsubj,[nrep,1]);
d0 = sort(d0(:));
d0 = agreement(d0);

g = sum(d0(:));% number of self comparisons

% f = number of between subject distances smaller than within (calculated for each subject and then summed)
f_subj = zeros(nsubj, 1);

for s = 1:nsubj
    nullmask = zeros(N,N);
    nullmask(nrep*(s-1)+1:nrep*s, nrep*(s-1)+1:nrep*s) = 1;
    nullmask = nullmask .* d0;
    
    %within_subj_dist = mean(nonzeros(dist_mtrx .* nullmask));
    within_subj_dist = max(nonzeros(dist_mtrx .* nullmask));
    btw_subj_dist = triu(dist_mtrx, 1) .* ~nullmask;
    f_subj(s) = sum(nonzeros(btw_subj_dist(2*(s-1)+1:2*s, :)) < within_subj_dist);
end

f = sum(f_subj)*2;

d = 1 - (f/(N^2 - g));% fraction of times across subjects are smaller than within subject distances
end