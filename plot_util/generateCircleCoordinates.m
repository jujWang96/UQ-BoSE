function [x, y] = generateCircleCoordinates(radius, centerX, centerY, numPoints)
    % Generates x and y coordinates of a circle.
    % Input:
    %   radius: Radius of the circle.
    %   centerX: X-coordinate of the circle's center.
    %   centerY: Y-coordinate of the circle's center.
    %   numPoints: Number of points to generate on the circle.
    % Output:
    %   x: X-coordinates of the points on the circle.
    %   y: Y-coordinates of the points on the circle.
    
    if nargin < 4
        numPoints = 100; % Default number of points
    end
    
    % Generate equally spaced angles between 0 and 2*pi
    angles = linspace(0, 2 * pi, numPoints);
    
    % Calculate x and y coordinates using the parametric equation of a circle
    x = radius * cos(angles) + centerX;
    y = radius * sin(angles) + centerY;
    
end