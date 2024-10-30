%% AVDASI 4 WINGBOX CALCULATIONS (USER DEFINED WINGBOX PARAMETERS)
% functions are run in order here (testing functions work, ready for using optimisation script)
clc; close all
% set path
old_path = path;
if ispc
    folder = "C:\Users\cj20794\OneDrive - University of Bristol\AEROSPACE Y4\AVDASI 4\structures\wing";
    addpath(genpath(folder));
    addpath(genpath("C:\Users\cj20794\OneDrive - University of Bristol\AEROSPACE Y3\Colour_Palettes"))
elseif ismac
    folder = "/Users/isaacgardner/Library/CloudStorage/OneDrive-UniversityofBristol/AEROSPACE Y4/AVDASI 4/structures/wing";
    addpath(genpath(folder));
    addpath("/Users/isaacgardner/Library/CloudStorage/OneDrive-UniversityofBristol/AEROSPACE Y3/Colour_Palettes")
    savepath
    output_folder = "/Users/isaacgardner/Library/CloudStorage/OneDrive-UniversityofBristol/AEROSPACE Y4/AVDASI 4/structures/wing/output_data";
end
colours_six

master_datasheet = readtable("DataSheetMASTER_19_10.xlsx",VariableNamingRule = "preserve"); % master datasheet as a table
aero_datasheet = readtable("Aero Data Sheet(Sheet1).csv",VariableNamingRule = "preserve"); % master datasheet as a table

%% MATERIAL PROPERTIES

E     = 70e9;              % Young's modulus [Pa] 
nu    = 0.34;              % Poisson's ratio
G     = E/(2*(1+nu));       % shear modulus [Pa]
rho   = 2700;              % density [kg/m^3]
YieldStrength = 276*1e6;   % [Pa]
ShearStrength = 207*1e6;   % [Pa]
g = 9.80665;

const.E = E;
const.nu = nu;
const.G = G;
const.rho = rho;
const.YieldStrength = YieldStrength;
const.ShearStrength = ShearStrength;
const.g = g;

%% 
wingbox = reconstruction_geom_v2(optim);


%% WINGBOX SELF WEIGHT

nSpan = 50; % Number of spanwise wing nodes considered
disc_density = 350; % nodes per length for cross section analyser script

wing_semi_span = (master_datasheet{35,2}/2); % wingspan from main datasheet (divided by two for semispan)
wingbox.semi_span = wing_semi_span; % adding semi span here so less inputs to functions later

Chord = aero_datasheet{8:9,4}; % root/tip chord from aero datasheet 

CS = WingParameterisation_HIPPO(wingbox, wing_semi_span, Chord, nSpan, true);

[wingbox_mass, wingbox_volume, total_mass, CS] = wingbox_weight_v1(wingbox, const, CS, nSpan, disc_density);

%% SHEAR FORCE AND BENDING MOMENT

lift_distribution = readmatrix("lift_distribution_estimate_v1.xlsx");
point_load_distribution = readmatrix("point_loads_v1.xlsx");
LEF = [1 -1 2.5]; % load cases 1, 2, 3 etc

[shear, bending, nodes] = shear_moment_v1(lift_distribution, point_load_distribution, wingbox_mass, wing_semi_span,...
    nSpan, const, LEF);

plot_force_moment(shear, bending, nodes, C, true)

%% AXIAL STRESSES
load_case = [1 2 3];

[stresses, Mz] = axial_stress_v1(CS, bending, nSpan, load_case, nodes, false);

% Validation of axial stress by comparison of moment back calculated
% plot([nodes.pos], [bending.max_total(:,1)])
% hold on
% plot([nodes.pos], Mz)

%% RESERVE FACTORS

[RF, allowable, CS, RF_crit] = structural_strength_RF_v1(CS, const, stresses, wingbox, nSpan);

writematrix(RF_crit,fullfile(output_folder,'critical_RFs'),"FileType","spreadsheet")



%% FUNCTIONS
function plot_force_moment(shear, bending, nodes, C, plot_yn)
    if plot_yn
        number_size = 20;
        font_size = 25;
    
        plot_name = 'shear_force';
        figure(Name=plot_name)
        set(gcf,'color','w')
        hold on
        h = plot(nodes.pos, shear.max_total,LineWidth=2);
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
        % exportgraphics(gcf, fullfile(fig_folder, [plot_name '.pdf']), ContentType='vector');
    
        plot_name = 'Bending_moment';
        figure(Name=plot_name)
        set(gcf,'color','w')
        hold on
        h = plot([nodes.pos], [bending.max_total], LineWidth=2);
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
        % exportgraphics(gcf, fullfile(fig_folder, [plot_name '.pdf']), ContentType='vector');
    end
end






