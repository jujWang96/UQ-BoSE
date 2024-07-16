function y = calc_sec(x,aboxwdt,boxmin,boxcen,darcsecperpix,offset,Axis)
y = x-offset;
if Axis==1
    y = -y/darcsecperpix;
else
    y = y/darcsecperpix;
end
y = y+boxcen;
y = y-boxmin;
y = y/aboxwdt;