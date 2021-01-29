function [efc, subj_id] = create_efc(files, nnode, nrep)
% Performs EFC generation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   REQUIRED INPUTS
%        files                .mat files for conversion to x (nodes x time series x scans)
%        nnode                Number of nodes in parcellation
%        nrep (Optional)      Scans per subject
%
%   OUTPUTS
%        efc                  Edge functional connectivity
%        subj_id              Subject IDs
%
%   Example:
%        [efc,ets]=create_efc(x);
%
%   References:
%        If you use this script, please cite:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsubj = length(files);
count = 0;
nedge = nnode*(nnode-1)/2;
utmask = triu(ones(nedge), 1) > 0;

efc_all = zeros(sum(utmask(:)), nsubj*nrep);
sub_all = zeros(nsubj*nrep,1);
for isubj = 1:nsubj
    data = load(files(isubj).name);
    nscan = length(data.parcel_time);
    scans = data.scans;
    
    %% loop scans
    for iscan = 1:2:nscan
        ts1 = data.parcel_time{iscan};
        ts2 = data.parcel_time{iscan+1};
        ts = vertcat(ts1, ts2);
        zi = zscore(ts);
        
        count = count + 1;
        ets = fcn_edgets(zi);
        efc = fcn_edgets2edgecorr(ets);
        efc_all(:, count) = efc(utmask);
        sub_all(count) = isubj;
    end
end

efc = efc_all;
subj_id = sub_all;
end