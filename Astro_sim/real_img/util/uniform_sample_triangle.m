function [sx, sy] = uniform_sample_triangle(vx, vy, n)
% Sample n points uniformly in a triangle.
% Reference: https://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle
%
% Args:
%   vx: A 1-by-3 vector giving the x-coordinates of the three vertices.
%   vy: A 1-by-3 vector giving the y-coordinates of the three vertices.
%   n: Number of points to sample.
%
% Returns:
%   sx: A 1-by-n vector giving the x-coordinates of the sampled points.
%   sy: A 1-by-n vector giving the y-coordinates of the sampled points.

r1 = rand(1, n);
r2 = rand(1, n);

coef1 = 1 - sqrt(r1);
coef2 = sqrt(r1) .* (1 - r2);
coef3 = r2 .* sqrt(r1);

sx = coef1 * vx(1) + coef2 * vx(2) + coef3 * vx(3);
sy = coef1 * vy(1) + coef2 * vy(2) + coef3 * vy(3);

end