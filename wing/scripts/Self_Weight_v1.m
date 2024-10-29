% SELF WEIGHT
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
    addpath(genpath(folder));
    addpath("/Users/isaacgardner/Library/CloudStorage/OneDrive-UniversityofBristol/AEROSPACE Y3/Colour_Palettes")
    savepath
end
colours_six

master_datasheet = readtable("DataSheetMASTER_19_10.xlsx",VariableNamingRule = "preserve"); % master datasheet as a table
aero_datasheet = readtable("Aero Data Sheet(Sheet1).csv",VariableNamingRule = "preserve"); % master datasheet as a table


%% WINGBOX DIMENSIONS

wingbox.BoxGeo = [0  0.2   0.6; 
            1  0.2   0.5]; % fraction of chord. Root front and rear

wingbox.tSkin = [0   0.01     
           1   0.005]; 

wingbox.tWeb  = [0   0.01    
           1   0.005]; 

wingbox.Stringer          = 4;         
wingbox.StringerHeight    = [0   0.05    
                       1   0.02]; 

wingbox.StringerThickness = [0   0.003 ;
                       1   0.002]; 



%% MATERIAL PROPERTIES

% material properties
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

%% find area of cross sections
nSpan = 50;
disc_density = 350;
WingSpan = (master_datasheet{35,2}/2); % wingspan from main datasheet (divided by two for semispan)
Chord = aero_datasheet{8:9,4}; % root/tip chord from aero datasheet 

CS = WingParameterisation_HIPPO(wingbox,WingSpan,Chord,nSpan,false);

for span_pos = 1:nSpan
    geometry = geom_cs(CS,wingbox,span_pos);
    [CS(span_pos).Ixx, CS(span_pos).Iyy, CS(span_pos).Centroid, CS(span_pos).area ]= ...
        x_section_anal(geometry.nodes(:,2:3), geometry.connectivity, geometry.t, disc_density, 0);
end

span_locs = linspace(0, WingSpan, 50);
wingbox_areas = interp1(0:(nSpan-1), [CS(:).area]', span_locs)'; 
wingbox_volume = flip(cumtrapz(span_locs, flip(wingbox_areas, 1)), 1);
wingbox_mass = wingbox_volume .* rho;




%% local functions

function out = geom_cs(CS,x,section_no)
    % func for preparing the geometric data about ribs to be used in a cross section analyser (ie getting nodes, connectivity, etc)
    % NODES
    % adds wingbox corners, then top edge, then bottom edge to the nodes variable
    nodes = [CS(section_no).WingBoxCornerXYZ];
    
    for ii = 1:x.Stringer
        nodes = [nodes; CS(section_no).TopStringerXYZ(ii)];
    end
    for ii = 1:x.Stringer
        nodes = [nodes; CS(section_no).BotStringerXYZ(ii)] ;       
    end
    nodes = cell2mat(nodes);

    % THICKNESSES
    t = zeros(4 + (2*x.Stringer),1);
    t([1 3],1) = CS(section_no).tSkin; % skin
    t([2 4],1) = CS(section_no).tWeb; % web
    t(5:end,1) = CS(section_no).StringerThickness; % stringers
    
    % CONNECTIVITY
    % first 4 are always corners so this is always the start first ascending y, then ascending z (always in x plane)
    
    connectivity = [1 2; 3 2; 4 3; 4 1]; % should always be the same (wing box)
    
    % adding the stringer connectivity
    j = (5:2:length(nodes)-1)';
    jj = (j+1);
    % so that the top stringers connectivity is in +ve direction:
    ii(1:x.Stringer,1:2) = [jj(1:x.Stringer) j(1:x.Stringer)]; 
    ii = [ii; [j(x.Stringer+1:end) jj(x.Stringer+1:end)]];
    connectivity = [connectivity; ii];
    
    out.nodes = nodes;
    out.connectivity = connectivity;
    out.t = t;
    % out.norm_pos = CS(section_no).xSpan./WingSpan;
end
% end




