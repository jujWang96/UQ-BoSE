function gcf = set_image_margin(gcf)
    leftMargin = 10;   % Left margin
    rightMargin = 1.0;  % Right margin
    topMargin = 0.5;    % Top margin
    bottomMargin = 10; % Bottom margin
    figPos = get(gcf, 'Position');
        
    % Calculate the new figure position and size to set the margins
    newWidth = figPos(3) + leftMargin + rightMargin;
    newHeight = figPos(4) + topMargin + bottomMargin;
    newLeft = figPos(1) - leftMargin;
    newBottom = figPos(2) - bottomMargin;
    
    % Set the new figure position and size
    set(gcf, 'Position', [newLeft, newBottom, newWidth, newHeight]);