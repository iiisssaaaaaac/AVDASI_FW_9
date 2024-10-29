function [Ixx, Iyy, centroid, area, lengths, midpoints, area_disc, parent_element] = ...
    x_section_anal(nodes, connectivity, t, disc_density, plot_yn)
%%%% Function for finding the structural attributes of a shape
% Inputs:
%   nodes: coordinates of the nodes defining the cross-section
%   connectivity: list of node pairs representing elements
%   t: thickness of elements
%   disc_density: discretisation density (number of divisions per element length)
%   plot_yn: flag to indicate whether to plot the shape (1 for yes, 0 for no)
%
% Outputs:
%   Ixx: moment of inertia about the x-axis
%   Iyy: moment of inertia about the y-axis
%   centroid: coordinates of the centroid of the cross-section
%   area: total area of the cross-section
%   lengths: lengths of individual elements
%   midpoints: midpoints of discretised elements
%   area_disc: area of each discrete section
%   parent_element: corresponding parent element for each discrete section

    % Set default values for optional inputs
    if isempty(plot_yn)
        plot_yn = 1; % Default is to plot the shape
    end
    if isempty(disc_density)
        disc_density = 200; % Default discretization density
    end

    % Get the number of elements and nodes
    n_element = length(connectivity);
    n_nodes = length(nodes);
    
    % Ensure thickness array matches the number of elements, assume uniform if not
    if length(t) < n_element
        t = ones(n_element) * t(1);
        warning('Fewer thickness values provided than required. Thickness assumed uniform')
    end

    % Initialize arrays for lengths and vectors of elements
    lengths = zeros(n_element,1);
    vector = zeros(n_element,2);

    % Calculate lengths and vectors for each element
    for n = 1:n_element
        temp = {[nodes(connectivity(n,1),:)] ; [nodes(connectivity(n,2),:)]};
        [vector(n,:), lengths(n)] = element_length(temp{1}, temp{2});
    end
    
    % Calculate total area of the cross-section
    area = sum(lengths .* t, "all");

    % Calculate angles of elements relative to the x-axis
    angle = zeros(n_element, 1);
    for n = 1:n_element
        angle(n) = acos(dot([1 0]', vector(n,:)') / norm(vector(n,:)));
    end
    
    %% Discretisation

    % Calculate number of discrete parts per element
    n_disc = lengths * (disc_density - 1); 
    increment = vector ./ n_disc; % Incremental distance between discrete points

    % Initialize midpoints and calculate first midpoint for each element
    counter = 1;
    first_points = increment / 2 + nodes(connectivity(:, 1), :);
   
    % Calculate midpoints for all discrete sections of each element
    for n = 1:n_element
        midpoints(counter, :) = first_points(n, :);
        little_n_element(n) = counter;
        counter = counter + 1;
        for i = 1:n_disc(n) - 1
            midpoints(counter, :) = midpoints(counter - 1, :) + increment(n, :);
            counter = counter + 1;
        end
    end
    % Adjust little_n_element for correct counting
    little_n_element(1) = [];
    little_n_element(n_element) = counter - 1;
    temp = [0 little_n_element(1:n_element - 1)];
    little_n_element = little_n_element - temp;

    %% Plot the shape if required
    if plot_yn == 1
        plot(midpoints(:, 1), midpoints(:, 2), '.r') % Plot midpoints in red
        hold on
        plot(nodes(:, 1), nodes(:, 2), 'Xk', 'MarkerSize', 12, 'LineWidth', 3) % Plot nodes
        axis equal
        box off
        set(gcf, 'color', 'w') % Set background to white
        xlabel('X axis')
        ylabel('Y axis')
    end
    
    %% Assign parent elements to discrete sections
    parent_element = [];
    length_disc = zeros(n_element, 1);
    for n = 1:n_element
        length_disc(n) = norm(increment(n, :)); % Length of each discrete section
        parent_element = [parent_element; ones(little_n_element(n), 1) * n];
    end
    % Fix first and last parent elements
    parent_element(1) = [];
    parent_element = [parent_element; parent_element(end, 1)];

    %% Calculate areas of discrete sections
    area_disc = zeros(sum(little_n_element), 1);
    counter = 1;
    for n = 1:n_element
        for i = 1:little_n_element(n)
            area_disc(counter) = t(n) * length_disc(n);
            counter = counter + 1;
        end
    end
    
    %% Calculate centroid of the cross-section
    centroid = sum(midpoints .* area_disc) / sum(area_disc);

    % Plot centroid if required
    if plot_yn == 1
        hold on
        plot(centroid(1), centroid(2), 'Xr') % Plot centroid in red
    end
    
    %% Calculate moments of inertia (Ixx, Iyy) using parallel axis theorem
    sexy_disc_xg = zeros(sum(little_n_element), 1);
    sexy_disc_yg = sexy_disc_xg;
 
    counter = 1;
    for n = 1:n_element
        for i = 1:little_n_element(n)
            % Moments of inertia for discrete sections
            sexy_disc_xg(counter) = t(n) * length_disc(n) / 12 * ((t(n) * cos(angle(n))^2 + (length_disc(n) * sin(angle(n))^2)));
            sexy_disc_yg(counter) = t(n) * length_disc(n) / 12 * ((t(n) * sin(angle(n))^2 + (length_disc(n) * cos(angle(n))^2)));
            counter = counter + 1;
        end
    end

    % Shift midpoints relative to the centroid
    centroid_midpoints = zeros(length(midpoints), 2);
    for n = 1:length(midpoints)
        centroid_midpoints(n, :) = midpoints(n, :) - centroid;
    end
    
    % Apply the parallel axis theorem to compute Ixx and Iyy
    Ixx = 0;
    Iyy = 0;
    for n = 1:length(centroid_midpoints)
        Iyy = Iyy + sexy_disc_xg(n) + area_disc(n) * centroid_midpoints(n, 1)^2;
        Ixx = Ixx + sexy_disc_yg(n) + area_disc(n) * centroid_midpoints(n, 2)^2;
    end

    % Element length function
    function [v, l] = element_length(a, b)
        v = b - a; % Vector between two nodes
        x = v.^2;
        l = sqrt(sum(x)); % Euclidean distance between the two nodes
    end
end
