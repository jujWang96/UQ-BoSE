addpath(genpath('../../G-SRG'))
addpath(genpath('../../Astro_sim')) 
load X_circle.mat
load check_X_pois_020.mat
rej=0;
l=0;
for b = 1:B
    if isempty(Contour{b})
        continue;
    end
    l=l+1;
    if test_similar(X,check_X_set{b})
        rej = rej+1;
    end
end
power = rej/l