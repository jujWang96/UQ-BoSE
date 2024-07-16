function [] = make_hist_areaflux(flux,area,edge_flux,edge_area)
nrow = length(flux);
for i = 1:nrow
    subplot(nrow,2,2*i-1)
    histogram(flux{i},edge_flux)
    subplot(nrow,2,2*i)
    histogram(area{i},edge_area)
end
