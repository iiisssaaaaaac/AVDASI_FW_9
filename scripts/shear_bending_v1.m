% %% 
% v1 initial code
clc; clear; close all
% set path
old_path = path;
if ispc
    folder = "C:\Users\cj20794\OneDrive - University of Bristol\AEROSPACE Y4\AVDASI 4\structures\wing";
    addpath(genpath(folder));
    addpath(genpath("C:\Users\cj20794\OneDrive - University of Bristol\AEROSPACE Y3\Colour_Palettes"))
elseif ismac
    folder = "/Users/isaacgardner/Library/CloudStorage/OneDrive-UniversityofBristol/AEROSPACE Y4/AVDASI 4/structures/wing";
    fig_folder = "/Users/isaacgardner/Library/CloudStorage/OneDrive-UniversityofBristol/AEROSPACE Y4/AVDASI 4/structures/wing";
    addpath(genpath(folder));
    addpath("/Users/isaacgardner/Library/CloudStorage/OneDrive-UniversityofBristol/AEROSPACE Y3/Colour_Palettes")
    savepath
end
colours_six
%% Importing 
g = 9.81;

master_datasheet = readtable("DataSheetMASTER_19_10.xlsx",VariableNamingRule = "preserve"); % master datasheet as a table

nSpan = 50; % number of discrete nodes to consider wing as

% Import the lift distribution data
lift_distribution = readmatrix("lift_distribution_estimate_v1.xlsx");
[~,index] = sort(lift_distribution(:,1),1,"ascend"); % find ascending order of position (in case not in order)
lift_distribution = lift_distribution(index,:); % reorder so in positional order

% Import the point loads data
point_load_distribution = readmatrix("point_loads_v1.xlsx");
point_load_distribution = point_load_distribution(:,2:4);
[~,index] = sort(point_load_distribution(:,1),1,"ascend"); % find ascending order of position (in case not in order)
point_load_distribution = point_load_distribution(index,:); % reorder so in positional order


nodes.pos = linspace(0,(master_datasheet{35,2}/2),nSpan)'; % find position of nodes
nodes.lift = interp1(lift_distribution(:,1),lift_distribution(:,2),nodes.pos,"spline");

% point loads are applied to nearest node, and combined into the original struct
nodes = point_load_on_node(nodes, point_load_distribution, g);


%% Shear force
% Shear force is sum of all forces from tip to root (so integrate lift, as it is distributed)

nodes.min_total = nodes.lift + nodes.min_point_loads;
nodes.max_total = nodes.lift + nodes.max_point_loads;

% LEF = 1; % load factor multiple
column = 1;
for LEF = [-1 1 2.5]
    shear.lift(:, column) = LEF .* -flip(cumtrapz(nodes.pos, flip(nodes.lift, 1)), 1);
    shear.min_point_loads(:, column) = LEF .* -flip(cumtrapz(nodes.pos, flip(nodes.min_point_loads, 1)), 1); 
    shear.max_point_loads(:, column) = LEF .* -flip(cumtrapz(nodes.pos, flip(nodes.max_point_loads, 1)), 1);
    
    shear.min_total(:, column) = LEF .* -flip(cumtrapz(nodes.pos, flip(nodes.min_total, 1)), 1);
    shear.max_total(:, column) = LEF .* -flip(cumtrapz(nodes.pos, flip(nodes.max_total, 1)), 1);
    column = column + 1;
end

number_size = 20;
font_size = 25;
plot_name = 'shear_force';
figure(Name=plot_name)
set(gcf,'color','w')
hold on
h = plot(nodes.pos,shear.max_total,LineWidth=2);
h(1).Color = C.b;
h(2).Color = C.c;
h(3).Color = C.d;
box on
fontsize(number_size, "points")
ylabel('Shear Force [$N$]','Interpreter','latex','fontsize', font_size)
xlabel('Spanwise Location [$m$]','Interpreter','latex','fontsize', font_size)
grid
legend box off
legend('1g','-1g','2.5g','Interpreter','latex')
set(gcf, 'Position',  [100, 100, 1200, 250])
exportgraphics(gcf, fullfile(fig_folder, [plot_name '.pdf']), ContentType='vector');


bending.lift = -flip(cumtrapz(nodes.pos, flip(shear.lift, 1)), 1);
bending.min_point_loads = -flip(cumtrapz(nodes.pos, flip(shear.min_point_loads, 1)), 1);
bending.max_point_loads = -flip(cumtrapz(nodes.pos, flip(shear.max_point_loads, 1)), 1);

bending.max_total = -flip(cumtrapz(nodes.pos, flip(shear.min_total, 1)), 1);
bending.min_total = -flip(cumtrapz(nodes.pos, flip(shear.max_total, 1)), 1);

plot_name = 'Bending_moment';
figure(Name=plot_name)
set(gcf,'color','w')
hold on
h = plot(nodes.pos, bending.max_total, LineWidth=2);
h(1).Color = C.b;
h(2).Color = C.c;
h(3).Color = C.d;
box on
fontsize(number_size, "points")
ylabel('Bending Moment [$Nm$]','Interpreter','latex','fontsize', font_size)
xlabel('Spanwise Location [$m$]','Interpreter','latex','fontsize', font_size)
grid
legend box off
legend('1g','-1g','2.5g','Interpreter','latex')
set(gcf, 'Position',  [100, 100, 1200, 250])
exportgraphics(gcf, fullfile(fig_folder, [plot_name '.pdf']), ContentType='vector');

%% functions init
function nodes = point_load_on_node(nodes, point_load_distribution, g)
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


