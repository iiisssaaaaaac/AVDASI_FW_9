function [wingbox_mass, wingbox_volume, total_mass, CS] = wingbox_weight_v1(wingbox, const, CS, nSpan, disc_density)
    % Function to calculate wingbox mass, volume, and total mass
    % Inputs:
    %   wingbox: structure containing all wingbox geometry info (in ratios of chord/span)
    %   const: constants and material properties (e.g., g = 9.81, rho = 2300)
    %   CS: cross-section struct with spanwise information, adapted from Terence AVDASI 3
    %   nSpan: number of spanwise nodes to consider
    %   disc_density: discretization density (nodes per length) for cross-section analysis
    
    % Loop over each spanwise section
    for span_pos = 1:nSpan
        % Prepare the geometry for the cross-section analyzer for this section
        geometry = geom_cs(CS, wingbox, span_pos);
        
        % Perform cross-sectional analysis to calculate area and moments of inertia
        [CS(span_pos).Iyy, CS(span_pos).Izz, CS(span_pos).Centroid, CS(span_pos).area, ~, CS(span_pos).small_nodes, ...
            CS(span_pos).small_node_areas] = ...
            x_section_anal(geometry.nodes(:,2:3), geometry.connectivity, geometry.t, disc_density, 0);
    end
    
    % Calculate wingbox volume by integrating cross-sectional areas from tip to root
    wingbox_volume = flip(cumtrapz(1:nSpan, flip([CS(:).area]', 1)), 1);
    
    % Calculate wingbox mass using material density
    wingbox_mass = wingbox_volume .* const.rho;
    
    % Extract the total mass at the root (first element of wingbox_mass)
    total_mass = wingbox_mass(1);

    %% Nested function for geometric setup of cross-sections

    function out = geom_cs(CS, x, section_no)
        % Function to set up geometric data for cross-section analysis
        % Takes nodes, connectivity, and thicknesses for each spanwise section
        % Inputs:
        %   CS: cross-section struct containing all geometric and material data
        %   x: wingbox input with stringer information
        %   section_no: current spanwise section to analyze
        
        % NODES
        % Initialize nodes array with wingbox corner coordinates
        nodes = [CS(section_no).WingBoxCornerXYZ];
        
        % Append top and bottom stringer coordinates to nodes array
        for ii = 1:x.Stringer
            nodes = [nodes; CS(section_no).TopStringerXYZ(ii)];
        end
        for ii = 1:x.Stringer
            nodes = [nodes; CS(section_no).BotStringerXYZ(ii)];
        end
        nodes = cell2mat(nodes); % Convert cell array to numeric matrix
    
        % THICKNESSES
        % Initialize thickness vector for wingbox components
        t = zeros(4 + (2 * x.Stringer), 1);
        t([1 3], 1) = CS(section_no).tSkin; % Set skin thickness at two locations
        t([2 4], 1) = CS(section_no).tWeb; % Set web thickness at two locations
        t(5:end, 1) = CS(section_no).StringerThickness; % Set thickness for each stringer
        
        % CONNECTIVITY
        % Define connectivity between nodes:
        % - First four nodes are corners of the wingbox
        connectivity = [1 2; 3 2; 4 3; 4 1];
        
        % Add stringer connectivity
        j = (5:2:length(nodes) - 1)'; % Odd indices for stringers
        jj = (j + 1); % Corresponding even indices for stringers
        ii(1:x.Stringer, 1:2) = [jj(1:x.Stringer) j(1:x.Stringer)]; % Top stringers
        ii = [ii; [j(x.Stringer+1:end) jj(x.Stringer+1:end)]]; % Bottom stringers
        connectivity = [connectivity; ii]; % Append stringer connectivity to main connectivity
        
        % Output struct containing geometry data for cross-section analysis
        out.nodes = nodes;
        out.connectivity = connectivity;
        out.t = t;
    end
end
