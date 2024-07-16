function [distinctRegionsX, distinctRegionsY] = partitionPolygons(polygonsX, polygonsY)
    % Partitions a region into distinct subregions using polygons.
    % Input:
    %   polygonsX: Cell array of x-coordinates of vertices of polygons.
    %   polygonsY: Cell array of y-coordinates of vertices of polygons.
    % Output:
    %   distinctRegionsX: Cell array of x-coordinates of distinct regions.
    %   distinctRegionsY: Cell array of y-coordinates of distinct regions.

    numPolygons = numel(polygonsX);
    
    % Initialize the result arrays
    distinctRegionsX = cell(numPolygons, 1);
    distinctRegionsY = cell(numPolygons, 1);
    
    % Iterate through each polygon
    for i = 1:numPolygons
        % Initialize the current polygon with the i-th polygon
        currentX = polygonsX{i};
        currentY = polygonsY{i};
        
        % Initialize the result with the i-th polygon
        resultX = currentX;
        resultY = currentY;
        
        % Perform boolean intersection with other polygons
        for j = (i+1):numPolygons
            [resultX, resultY] = polybool('intersection', currentX, currentY, polygonsX{j}, polygonsY{j});
            % Store the distinct region in the result arrays
            distinctRegionsX{i} = resultX;
            distinctRegionsY{i} = resultY;
        end
    end

    distinctRegionsDiffX = cell(numPolygons, 1);
    distinctRegionsDiffY = cell(numPolygons, 1);
    
    % Iterate through each polygon
    for i = 1:numPolygons
        % Initialize the current polygon with the i-th polygon
        currentX = polygonsX{i};
        currentY = polygonsY{i};
        
        % Initialize the result with the i-th polygon
        resultX = currentX;
        resultY = currentY;
        
        % Perform boolean intersection with other polygons
        for j = 1:numPolygons
            if j ~= i
                [resultX, resultY] = polybool('subtraction', currentX, currentY, polygonsX{j}, polygonsY{j});
                 % Store the distinct region in the result arrays
                distinctRegionsDiffX{i} = resultX;
                distinctRegionsDiffY{i} = resultY;
            end
        end
       
    end
    distinctRegionsX = [distinctRegionsX; distinctRegionsDiffX];
    distinctRegionsY = [distinctRegionsY; distinctRegionsDiffY];
    nonEmptyCells = cellfun(@(x) ~isempty(x), distinctRegionsX);
    distinctRegionsX = distinctRegionsX(nonEmptyCells);
    distinctRegionsY = distinctRegionsY(nonEmptyCells);
end