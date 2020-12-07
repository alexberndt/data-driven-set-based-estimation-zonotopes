function maxDist = getMaxDist(zonoObj)
    
    plotdata = plot(zonoObj,[1 2],'r-');
    x = plotdata.XData;
    y = plotdata.YData;
    
    numberOfCoords = length(x);

    for k = 1 : numberOfCoords
      % Label the kth point in the plot.
      text(x(k)+0.01, y(k), num2str(k));
    end
    maxDistance = zeros(1, numberOfCoords);
    indexOfMax = zeros(1, numberOfCoords, 'int32');
    %-----------------------------------------------
    % Main engine:
    % Find the furthest away points.
    for k = 1 : numberOfCoords
      distances = sqrt((x-x(k)).^2 + (y-y(k)).^2);
      [maxDistance(k), indexOfMax(k)] = max(distances);
    end
    %-----------------------------------------------
    % Done!  Now show out results in the command window.
    % Display in command window.
    for k = 1 : numberOfCoords
      thisDistance = maxDistance(k);
      thisIndex = indexOfMax(k);
    %   fprintf('Point #%d at (%.1f, %.1f) is farthest away from point #%d at (%.1f, %.1f) and is at a distance of %f\n',...
    %     thisIndex, x(thisIndex), y(thisIndex),...
    %     k, x(k), y(k), thisDistance);
    end
    
    maxDist = max(maxDistance);

end