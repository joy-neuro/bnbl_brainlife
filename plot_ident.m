function [d] = plot_ident(x, subj_id, niter, idiff, d);
% Plot identifiability measures versus subject randomization
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

Idiff_rand = zeros(niter,1);
discr_rand = zeros(niter,1);
N = length(subj_id);
nsubj = max(subj_id);
nrep = sum(subj_id == 1);
d0 = repmat(1:nsubj, [nrep, 1]);
d0 = sort(d0(:));
d0 = agreement(d0);
R = corr(x);
d0corr = d0(all(~isnan(R)),all(~isnan(R)));
g = sum(d0(:));

cols = size(x, 2);
for i = 1:niter
    % create randomized Iself
    rp = randperm(cols);
    fc_all_tmp = x(:, rp);
    R_tmp = corr(fc_all_tmp);
    Idiff_rand(i) = mean(R_tmp(triu(d0corr == 1, 1))) - mean(R_tmp(triu(d0corr == 0, 1)));

    % randomized discr
    dist_mtrx = zeros(N, N);% NxN of distance between scans
    for dist = 1:size(fc_all_tmp, 2)
        dist_mtrx(:, dist) = sum(abs(fc_all_tmp - fc_all_tmp(:, dist)));
    end
    
    f_subj_rand = zeros(nsubj, 1);

    for s = 1:nsubj
        nullmask = zeros(N,N);
        nullmask(nrep*(s-1)+1:nrep*s, nrep*(s-1)+1:nrep*s) = 1;
        nullmask = nullmask .* d0;

        %within_subj_dist = mean(nonzeros(dist_mtrx .* nullmask));
        within_subj_dist = max(nonzeros(dist_mtrx .* nullmask));
        btw_subj_dist = triu(dist_mtrx, 1) .* ~nullmask;
        f_subj_rand(s) = sum(nonzeros(btw_subj_dist(nrep*(s-1)+1:nrep*s, :)) < within_subj_dist);
    end

    f_rand = sum(f_subj_rand)*2;

    discr_rand(i) = 1 - (f_rand/(N^2 - g));% fraction of times across subjects are smaller than within subject distances
end

z_Idiff = (idiff - mean(Idiff_rand)) / std(Idiff_rand);
z_discr = (d - mean(discr_rand)) / std(discr_rand);

h1 = histogram(zscore(Idiff_rand));
hold on;
h2 = histogram(zscore(discr_rand));
xline(z_Idiff, '--', 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
xline(z_discr, '--', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
legend({'Rand Iself','Rand Discr', 'Iself', 'Discr'});
%breakxaxis([15 45]);
%set(gca, 'xscale', 'log');
axis square;
hold off;
end
