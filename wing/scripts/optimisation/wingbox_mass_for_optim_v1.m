function [total_mass] = wingbox_mass_for_optim_v1(x0)
    wingbox = reconstruction_geom_v2(x0);
    master_datasheet = readtable("DataSheetMASTER_19_10.xlsx",VariableNamingRule = "preserve"); % master datasheet as a table
    aero_datasheet = readtable("Aero Data Sheet(Sheet1).csv",VariableNamingRule = "preserve"); % master datasheet as a table
    wingbox.Stringer = 4;
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

    % %% WINGBOX SELF WEIGHT
    
    nSpan = 50; % Number of spanwise wing nodes considered
    disc_density = 350; % nodes per length for cross section analyser script
    
    wing_semi_span = (master_datasheet{35,2}/2); % wingspan from main datasheet (divided by two for semispan)
    wingbox.semi_span = wing_semi_span; % adding semi span here so less inputs to functions later
    
    Chord = aero_datasheet{8:9,4}; % root/tip chord from aero datasheet 
    
    CS = WingParameterisation_HIPPO(wingbox, wing_semi_span, Chord, nSpan, false);
    
    [~, ~, total_mass, ~] = wingbox_weight_v1(wingbox, const, CS, nSpan, disc_density);
end