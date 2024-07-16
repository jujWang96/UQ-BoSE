function flux = get_adjust_flux(flux,area,n, total_area)
if nargin==3
    total_area = 1;
end
if length(n)<length(flux)
    n = n*ones(length(flux),1);
end
if length(total_area)<length(flux)
    total_area = total_area*ones(length(flux),1);
end

for ii = 1:length(flux)
    flux(ii) = round(flux(ii)- area(ii)*(n(ii)-flux(ii))/(total_area(ii)-area(ii)));
end
%get the flux adjusted for the background intensity 
