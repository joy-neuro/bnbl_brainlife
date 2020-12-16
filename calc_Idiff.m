function [idiff, iself, iothers] = calc_Idiff(x, subj_id)
% Performs Idiff calculation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   REQUIRED INPUTS
%        x                    efc vectorized x scans
%        subj_id              Index of subjects per scan
%
%   OUTPUTS
%        idiff                Differential identifiability (Iself - Iothers)
%        iself (opt)          Scan self-similarity
%        iothers (opt)        Scan between-subject similarity
%
%   Example:
%        [idiff]=calc_Idiff(x, subj_id);
%
%   References:
%        If you use this script, please cite:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsubj = max(subj_id);
nrep = sum(subj_id == 1);
d0 = repmat(1:nsubj, [nrep, 1]);
d0 = sort(d0(:));
d0 = agreement(d0);
R = corr(x);
d0corr = d0(all(~isnan(R)), all(~isnan(R)));

iself = mean(R(triu(d0corr == 1, 1)));
iothers = mean(R(triu(d0corr == 0, 1)));
idiff = iself - iothers;
end