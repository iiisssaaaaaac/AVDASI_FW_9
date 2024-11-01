function [shear, bending, nodes, new_bending, total_lift] = shear_moment_v1(lift_distribution, point_load_distribution, wingbox_mass, wing_semi_span,...
    nSpan, const, load_case)
    % what is this function

    % Import the lift distribution data
    % lift_distribution = readmatrix("lift_distribution_estimate_v1.xlsx");
    [~,index] = sort(lift_distribution(:,1),1,"ascend"); % find ascending order of position (in case not in order)
    lift_distribution = lift_distribution(index,:); % reorder so in positional order

    % Import the point loads data
    % point_load_distribution = readmatrix("point_loads_v1.xlsx");
    point_load_distribution = point_load_distribution(:,2:4);
    [~,index] = sort(point_load_distribution(:,1),1,"ascend"); % find ascending order of position (in case not in order)
    point_load_distribution = point_load_distribution(index,:); % reorder so in positional order
    
    
    nodes.pos = linspace(0, wing_semi_span,nSpan)'; % find position of nodes
    nodes.lift_per_m = interp1(lift_distribution(:,1), lift_distribution(:,2), nodes.pos, "spline");
    
    % point loads are applied to nearest node, and combined into the original struct
    nodes = point_load_on_node(nodes, point_load_distribution, const.g);
    
    %% Shear force
    % Shear force is sum of all forces from tip to root (so integrate lift, as it is distributed)

    
    shear.lift = flip(cumtrapz(nodes.pos, flip(nodes.lift_per_m, 1)), 1); % finding the sum of lift at each point from tip to root
    nodes.lift = [0; shear.lift] - [shear.lift; 0]; % finding the difference between sum of lift at each node to find lift at each node
    nodes.lift(1) = [];


    shear.wingbox_weight = - wingbox_mass .* const.g;
    nodes.wingbox_weight = [0; shear.wingbox_weight] - [shear.wingbox_weight; 0]; % finding the difference between sum of lift at each node to find lift at each node
    nodes.wingbox_weight(1) = [];
    nodes.min_total = nodes.lift + nodes.wingbox_weight + nodes.min_point_loads;
    nodes.max_total = nodes.lift + nodes.wingbox_weight + nodes.max_point_loads;
    

    shear.lift = - shear.lift .* load_case;
    total_lift = - shear.lift(1,1);

    
    shear.wingbox_weight = load_case .* - shear.wingbox_weight;

    shear.min_point_loads = load_case .* - flip(cumsum(flip(nodes.min_point_loads,1)),1); 
    shear.max_point_loads = load_case .* - flip(cumsum(flip(nodes.max_point_loads,1)),1);

    shear.min_total = shear.lift + shear.min_point_loads + shear.wingbox_weight;
    shear.max_total = shear.lift + shear.max_point_loads + shear.wingbox_weight;
        
        
    
    % load_case_no = 1;
    % for load_case_no = 1:width(load_case)
        % LEF = load_case(load_case_no);

        % shear.min_point_loads(:, load_case_no) = LEF .* - flip(cumsum(flip(nodes.min_point_loads,1)),1); 
        % shear.max_point_loads(:, load_case_no) = LEF .* - flip(cumsum(flip(nodes.max_point_loads,1)),1);
        % 
        
        
        % shear.lift(:, load_case_no) = LEF .* - flip(cumsum(flip(nodes.lift, 1)), 1); % commented out so not going back and forth
        % between and maybe introducing extra rounding errors

        % shear.min_total(:, load_case_no) = shear.lift(:, load_case_no) + shear.min_point_loads(:, load_case_no) + shear.wingbox_weight(:, load_case_no);
        % shear.max_total(:, load_case_no) = shear.lift(:, load_case_no) + shear.max_point_loads(:, load_case_no) + shear.wingbox_weight(:, load_case_no);
        % 
    % end
    %% Bending moment
    
    bending.lift = -flip(cumtrapz(nodes.pos, flip(shear.lift, 1)), 1);
    bending.min_point_loads = -flip(cumtrapz(nodes.pos, flip(shear.min_point_loads, 1)), 1);
    bending.max_point_loads = -flip(cumtrapz(nodes.pos, flip(shear.max_point_loads, 1)), 1);

    bending.max_total = -flip(cumtrapz(nodes.pos, flip(shear.min_total, 1)), 1);
    bending.min_total = -flip(cumtrapz(nodes.pos, flip(shear.max_total, 1)), 1);


    % Reformatting into a better struct shape
    new_bending = repmat(struct('lift', [], 'min_point_loads', [], 'max_point_loads', [], 'max_total', [], 'min_total', []), 1, 50);
    for span_pos = 1:nSpan
        new_bending(span_pos).lift = bending.lift(span_pos, :);
        new_bending(span_pos).min_point_loads = bending.min_point_loads(span_pos, :);
        new_bending(span_pos).max_point_loads = bending.max_point_loads(span_pos, :);
        new_bending(span_pos).max_total = bending.max_total(span_pos, :);
        new_bending(span_pos).min_total = bending.min_total(span_pos, :);
    end

    %% Nested function
    function [nodes] = point_load_on_node(nodes, point_load_distribution, g)
        nodes.max_point_loads = zeros(length(nodes.pos),1);
        nodes.min_point_loads = zeros(length(nodes.pos),1);
        
        % Loop through each point load and distribute it using interpolation
        for i = 1:size(point_load_distribution, 1)
            % first is case where load applied exactly at node
            if any(point_load_distribution(i,1) == nodes.pos)
                node_in = find(point_load_distribution(i,1) == nodes.pos, 1);
    
                nodes.max_point_loads(node_in,1) = nodes.max_point_loads(node_in,1) + point_load_distribution(i,2);
                nodes.min_point_loads(node_in,1) = nodes.min_point_loads(node_in,1) + point_load_distribution(i,3);
            else
                % kinda complex one here, but essentially finds the ratio of weight distribution between two nodes
                node_in = find((nodes.pos - point_load_distribution(i,1))<=0, 1, "last");
                node_out = find((nodes.pos - point_load_distribution(i,1))>=0, 1, "first");
        
                weight_in = -(nodes.pos(node_in) - point_load_distribution(i,1)) / nodes.pos(2,1);
                weight_out = (nodes.pos(node_out) - point_load_distribution(i,1)) / nodes.pos(2,1);
        
                nodes.max_point_loads(node_in:node_out,1)  = nodes.max_point_loads(node_in:node_out,1) + ...
                                                             [weight_in * point_load_distribution(i,2) .* g; 
                                                              weight_out * point_load_distribution(i,2) .* g];
                nodes.min_point_loads(node_in:node_out,1)  = nodes.max_point_loads(node_in:node_out,1) + ...
                                                             [weight_in * point_load_distribution(i,3) .* g; 
                                                              weight_out * point_load_distribution(i,3) .* g];
            end
        end
    end
end

