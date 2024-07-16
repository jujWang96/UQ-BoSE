function MISE = calc_MISE(estimate_intensity,true_intensity)
fun = @(x,y) (estimate_intensity-true_intensity).^2;
MISE =  integral2(fun,0,1,0,1);
