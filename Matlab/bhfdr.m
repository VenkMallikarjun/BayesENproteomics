function fdr = bhfdr(p)
% Compute Benjamini-Hochberg FDRs
nanidx = isnan(p);

m = numel(p(~nanidx));
[p_ord, idx] = sort(p);
fdr_ord =  min(1, p_ord(:) .* (m./[(1:m),NaN(1,sum(nanidx,1))])');
fdr(idx) = fdr_ord;
fdr = fdr';
end
