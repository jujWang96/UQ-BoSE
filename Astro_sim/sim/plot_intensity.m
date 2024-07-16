function [] = plot_intensity( cx,cy, selected, contour, region_intensity,pt_size)
%plot the 2D histogram based k  intensity of each region 
%
%
%
%
% r=1000;
% hist = zeros(r,r);
% 
% for i = contour.regionId
%     if any(contour.background==i)
%         
%     end
%     disp(i)
%     [x,y] = get_curve(contour.contourV{i},true);
%     hist = hist+region_intensity(i)*double(poly2mask(x*r,y*r,r,r));
%     
% end
% intense_image = flip(hist/max(region_intensity));
% intense_image = rot90(intense_image,3);
% figure
% histogram2('XBinEdges',0:r,'YBinEdges',0:r,'BinCounts',intense_image,'DisplayStyle','tile')
% axis equal
% colorbar
colormap jet;
p = jet(1000);
L = length(p);
a = (L-1)/(log(max(region_intensity(contour.regionId)))-log(min(region_intensity(contour.regionId))));
b = ((1+L)-a*(log(max(region_intensity(contour.regionId)))+log(min(region_intensity(contour.regionId)))))/2;
for seg = reshape(contour.regionId,1,[])
      scatter(cx(selected{seg}), cy(selected{seg}), pt_size,  p(min(max(round(a*log(region_intensity(seg))+b),1),length(p)),:), '.')
      hold on 
end
%legend(string(contour.regionId))
 
% for i = contour.regionId
%      contourVcurr = contour.contourV{i};
%     x1 = contourVcurr(:,1);
%     x2 = contourVcurr(:,2);
%     y1 = contourVcurr(:,3);
%     y2 = contourVcurr(:,4);
%     plot([x1';x2'], [y1';y2'],'Color','k');         
% end
cb=colorbar;
set(cb,'position',[0.65 .1 .02 .15])
caxis([min(region_intensity(contour.regionId)), max(region_intensity(contour.regionId))]);
axis equal
axis([0 1 0 1])
