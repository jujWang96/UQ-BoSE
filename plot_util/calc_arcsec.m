function y = calc_arcsec(x,boxwdt,boxmin,boxcen,arcsecperpix,offset,Axis)
%convert X-Y coordinate to arcsec coordinate centered on Chandra field
%Axis: 1  transform x to RA (reverted)
%      0  transform y to DEC (nonreverted)
%
y = x*boxwdt;
y = y+boxmin;
y = y-boxcen;
if Axis==1
    y = -y*arcsecperpix;
else
    y = y*arcsecperpix;
end
y = y+offset;