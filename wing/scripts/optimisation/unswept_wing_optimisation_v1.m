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

%% WINGBOX DIMENSIONS
clc

wingbox_base.BoxGeo = [0.2  ; 0.6; 
            0.2  ; 0.5]; % fraction of chord. Root front and rear

wingbox_base.tSkin = [   0.03     
              0.008]; 

wingbox_base.tWeb  = [   0.01    
              0.005]; 

wingbox_base.Stringer          = 4;         
wingbox_base.StringerHeight    = [   0.1    
                          0.02]; 

wingbox_base.StringerThickness = [   0.008 ;
                          0.003]; 

x0 = cell2mat(struct2cell(wingbox_base));

% LOWER BOUND

wingbox_low.BoxGeo = [0.1  ; 0.4; 
            0.1 ;  0.4]; % fraction of chord. Root front and rear

wingbox_low.tSkin = [   0.01     
              0.001]; 

wingbox_low.tWeb  = [   0.005    
              0.001]; 

wingbox_low.Stringer          = 4;         
wingbox_low.StringerHeight    = [   0.01    
                          0.002]; 

wingbox_low.StringerThickness = [   0.002 ;
                          0.001]; 

xlb = cell2mat(struct2cell(wingbox_low));

% UPPER BOUND

wingbox_up.BoxGeo = [0.6 ;  0.9; 
            0.45 ;  0.9]; % fraction of chord. Root front and rear

wingbox_up.tSkin = [   0.1     
              0.1]; 

wingbox_up.tWeb  = [   0.15    
              0.1]; 

wingbox_up.Stringer          = 4;         
wingbox_up.StringerHeight    = [   0.2    
                          0.05]; 

wingbox_up.StringerThickness = [   0.1 ;
                          0.05];

xub = cell2mat(struct2cell(wingbox_up));

% Constraints:
Aa = [];
b = [];
Aeq = [];
beq = [];
%

options = optimoptions('fmincon','Display','iter','MaxIterations',50,'FunctionTolerance',1e-3);
optim = fmincon(@wingbox_mass_for_optim_v1, x0 ,Aa,b,Aeq,beq,xlb,xub,@wing_RFs_nlcon_v1,options);

%%
optimised_wingbox = reconstruction_geom_v2(optim);
wingbox_mass = wingbox_mass_for_optim_v1(optim);




