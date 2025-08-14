function [EMD] = envm_emd_metrics(imf,w)

for i=1:length(imf)
    EMD.pow_imf(i) = sum(abs(imf{i}));
    EMD.mu_w(i) = nanmean(w{i}); %#ok<*NANMEAN>
    EMD.var_w(i) = nanvar(w{i}); %#ok<*NANVAR>
    EMD.sd_w(i) = nanstd(w{i}); %#ok<NANSTD>
end

if numel(imf)>1
    EMD.imf_ratio21 = EMD.pow_imf(2)./EMD.pow_imf(1);
else
    EMD.imf_ratio21 = nan;
end

end