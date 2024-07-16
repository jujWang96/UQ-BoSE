function plot_arc(loc_ring, radius_in, radius_out,lw,colr)

if nargin == 4
    colr = [0.6 0.6 0.6];
end 
ang(loc_ring, radius_out, [-pi/2 pi/2], colr, lw);
ang(loc_ring, radius_in, [-pi/2 pi/2], colr, lw);
plot([loc_ring(1) loc_ring(1)], [loc_ring(2)-radius_out loc_ring(2)-radius_in], 'Color', colr, 'LineWidth', lw)
plot([loc_ring(1) loc_ring(1)], [loc_ring(2)+radius_in loc_ring(2)+radius_out], 'Color', colr, 'LineWidth', lw)

end