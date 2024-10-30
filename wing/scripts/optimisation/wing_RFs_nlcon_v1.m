function [g, ceq] = wing_RFs_nlcon_v1(x0)
    wingbox = reconstruction_geom_v2(x0);
    master_datasheet = readtable("DataSheetMASTER_19_10.xlsx",VariableNamingRule = "preserve"); % master datasheet as a table
    aero_datasheet = readtable("Aero Data Sheet(Sheet1).csv",VariableNamingRule = "preserve"); % master datasheet as a table
    wingbox.Stringer = 4;
    % %% MATERIAL PROPERTIES

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
    % %% WINGBOX SELF WEIGHT
    
    nSpan = 50; % Number of spanwise wing nodes considered
    disc_density = 350; % nodes per length for cross section analyser script
    
    wing_semi_span = (master_datasheet{35,2}/2); % wingspan from main datasheet (divided by two for semispan)
    wingbox.semi_span = wing_semi_span; % adding semi span here so less inputs to functions later
    
    Chord = aero_datasheet{8:9,4}; % root/tip chord from aero datasheet 
    
    CS = WingParameterisation_HIPPO(wingbox, wing_semi_span, Chord, nSpan, false);
    
    [wingbox_mass, ~, ~, CS] = wingbox_weight_v1(wingbox, const, CS, nSpan, disc_density);
    
    % %% SHEAR FORCE AND BENDING MOMENT
    
    lift_distribution = readmatrix("lift_distribution_estimate.xlsx");
    point_load_distribution = readmatrix("point_loads_v1.xlsx");
    LEF = [1 -1 2.5]; % load cases 1, 2, 3 etc
    
    [~, bending, nodes] = shear_moment_v1(lift_distribution, point_load_distribution, wingbox_mass, wing_semi_span,...
        nSpan, const, LEF);
    % C = [];
    % plot_force_moment(shear, bending, nodes, C, false)
    
    % %% AXIAL STRESSES
    load_case = [1 2 3];
    
    [stresses, ~] = axial_stress_v1(CS, bending, nSpan, load_case, nodes, false);
    
    % Validation of axial stress by comparison of moment back calculated
    % plot([nodes.pos], [bending.max_total(:,1)])
    % hold on
    % plot([nodes.pos], Mz)
    
    % %% RESERVE FACTORS
    
    [~, ~, CS, RF_crit] = structural_strength_RF_v1(CS, const, stresses, wingbox, nSpan);
    
    g =  ones((nSpan - 1)*4,1) - reshape(RF_crit(1:49,:), (nSpan - 1)*4,1);
    ceq = [];



end