estimate_intensity =  @(x,y) 3*(0 <= x) & (x <= 1) & (0 <= y) & (y <= 1);
true_intensity =  @(x,y) 2*(0 <= x) & (x <= 1) & (0 <= y) & (y <= 1);
MISE = calc_MISE(estimate_intensity,true_intensity);