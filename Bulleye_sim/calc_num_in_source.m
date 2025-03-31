function num = calc_num_in_source(sz,unit_area,contrast,background_intensity,SNR,contrast_def)

%calculate the expected number in each source 
% input: size: the ratio between area of object
%unit area: the unit area of objects
% contrast: contrast between extended sources from inner to outmost 
% background_intensity: intensity of background
% SNR: signal to noise ratio between the source and background. 
% contrast_def: diffierent definition of contrast 
%       1.contrast = (B1-B2)/B1
%       2.contrast = ((S1/A1)/(S1+S2)/(A1+A2))
syms x
eqn = x^2 == (x+background_intensity*sum(sz)*unit_area)*SNR^2;
SC = solve(eqn,x);
SC = round(SC(SC>0))
if contrast_def==1
syms x
% for i =1:length(size)
%     if i==1
%         intensity(end) = x;
%     else
%         intensity(end-i+1) = (intensity(end-i+2)+contrast*background_intensity)/(1-contrast);
%     end
%     
% end
for i =1:length(sz)
    if i==1
        intensity(1) = x;
    else
        intensity(i) = (intensity(i-1)+contrast*background_intensity)/(1-contrast);
    end
    
end

eqn = dot(flip(intensity),sz)*unit_area == SC;

inten_out= solve(eqn,x);
for i =1:length(sz)
    if i==1
        intensity(1) = inten_out;
    else
        intensity(i) = (intensity(i-1)+contrast*background_intensity)/(1-contrast);
    end
end
intensity = flip(intensity);
for i =1:length(sz)
    num(i) = double(round(intensity(i)*sz(i)*unit_area));
end
elseif contrast_def==2
    syms x
    for i =1:length(sz)
        if i==1
            count(1) = x;
        else
            count(i) = (count(i-1)/contrast/sz(i-1)*(sz(i-1)+sz(i)))-count(i-1);
        end
        
    end
    eqn = sum(count) == SC;

    S1= round(solve(eqn,x));
    for i =1:length(sz)
        if i==1
            count(1) = S1;
        else
            count(i) = (count(i-1)/contrast/sz(i-1)*(sz(i-1)+sz(i)))-count(i-1);
        end
        
    end
    num = double(round(count));
    num(end) = SC-sum(num(1:end-1)); 
else
    num = zeros(1,length(sz));
    disp("the contrast definition is not spcified correctly")
end
