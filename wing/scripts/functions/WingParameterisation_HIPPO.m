% Model Parameterisation
function [CSred,CS] = WingParameterisation_HIPPO(x,WingSpan,chord,nSpan,PLOT)
% addpath("/Users/isaacgardner/Library/CloudStorage/OneDrive-UniversityofBristol/AEROSPACE Y4/AVDASI 4/structures/wing/input_data")
run NACA_632615  % aerofoil
% run NACA_2412  % tip aerofoil

%% interpolating data
nStringer         = x.Stringer;
% WingSpan    = 10; % meters
xSpan             = linspace(0,WingSpan,nSpan);
Chord             = linspace(chord(1,:),chord(2,:),nSpan); % chord distribution along the span (m)
BoxGeo            = interp1(x.BoxGeo(:,1),x.BoxGeo,xSpan/WingSpan);
tSkin             = interp1(x.tSkin(:,1),x.tSkin,xSpan/WingSpan);
tWeb              = interp1(x.tWeb(:,1),x.tWeb,xSpan/WingSpan);
StringerHeight    = interp1(x.StringerHeight(:,1),x.StringerHeight,xSpan/WingSpan);
StringerThickness = interp1(x.StringerThickness(:,1),x.StringerThickness,xSpan/WingSpan); 


for i=1:nSpan
   CS(i).xSpan             = xSpan(i);             % position of the cross-section along the span (m)
   CS(i).Chord             = Chord(i);
   CS(i).beta              = xSpan(i)/WingSpan;    % interpolation coefficient
   CS(i).ZYProfile_Norm    = NACA632615*(1-CS(i).beta) + NACA632615*CS(i).beta; % for blending between two airfoils (tip and root)


   CS(i).ZProfileOffset    = 0.5*(BoxGeo(i,2)+BoxGeo(i,3))*Chord(i);  % middle of the wing box along the chord
   CS(i).nProfilePoint     = size(CS(i).ZYProfile_Norm,1);
   CS(i).XYZProfile        = [CS(i).xSpan*ones(CS(i).nProfilePoint,1) fliplr(CS(i).ZYProfile_Norm*Chord(i))];
   CS(i).XYZProfile(:,3)   = CS(i).XYZProfile(:,3) - CS(i).ZProfileOffset; 

   %% wing box geometry (rectangular approximation)
   zMid   = 0.5*(BoxGeo(i,2)+BoxGeo(i,3)); % normalised
   zStart = (BoxGeo(i,2)-zMid)*Chord(i);   % normalised
   zEnd   = (BoxGeo(i,3)-zMid)*Chord(i);   % normalised

   UpperProfile = CS(i).XYZProfile(1:21,:);   % upper skin = suction side
   LowerProfile = CS(i).XYZProfile(21:end,:); % lower skin = pressure side

   % 1 ---- 2 
   % 4 ---- 3
   P1 = interp1(UpperProfile(:,3),UpperProfile,zStart);   % 1
   P2 = interp1(UpperProfile(:,3),UpperProfile,zEnd);     % 2
   P3 = interp1(LowerProfile(:,3),LowerProfile,zEnd);     % 3
   P4 = interp1(LowerProfile(:,3),LowerProfile,zStart);   % 4
   P1(2) = 0.5*(P1(2)+P2(2)); % average the height to get a rectangle
   P2(2) = P1(2);
   P3(2) = 0.5*(P3(2)+P3(2)); % average the height to get a rectangle
   P4(2) = P3(2);
   CS(i).WingBoxCornerXYZ = [P1; P2;  P3; P4]; 

   %%% wing box thickness
   CS(i).tSkin = tSkin(i,2);
   CS(i).tWeb  = tWeb(i,2);

   %% offset the profile and wing box to be more or less centre around elastic axis
   Zmid = 0.5*(CS(i).WingBoxCornerXYZ(1,3) + CS(i).WingBoxCornerXYZ(2,3));
   Ymid = 0.5*(CS(i).WingBoxCornerXYZ(1,2) + CS(i).WingBoxCornerXYZ(4,2));
   CS(i).YProfileOffset = Ymid;
   CS(i).XYZProfile(:,2)         = CS(i).XYZProfile(:,2) - CS(i).YProfileOffset; 
   CS(i).WingBoxCornerXYZ(:,2)   = CS(i).WingBoxCornerXYZ(:,2) - CS(i).YProfileOffset; 

   P1 = CS(i).WingBoxCornerXYZ(1,:);
   P2 = CS(i).WingBoxCornerXYZ(2,:);
   P3 = CS(i).WingBoxCornerXYZ(3,:); 
   P4 = CS(i).WingBoxCornerXYZ(4,:); 

   %% Top skin stringer positions and nodes
   CS(i).StringerHeight    = StringerHeight(i,2);
   CS(i).StringerThickness = StringerThickness(i,2);

   alfa       = linspace(0,1,nStringer+2);
   alfa       = alfa(2:end-1);
   PStringers = [P1(1)*(1-alfa) + P2(1)*alfa;
                 P1(2)*(1-alfa) + P2(2)*alfa;
                 P1(3)*(1-alfa) + P2(3)*alfa]; % [X;Y;Z]

   for j = 1:nStringer
      CS(i).TopStringerXYZ{j} = (PStringers(:,j) + [0 0 ;-CS(i).tSkin/2 -CS(i).tSkin/2-CS(i).StringerHeight;0 0])'; % [X Y Z]
   end

   %% Bottom skin stringer positions and nodes
   PStringers = [P4(1)*(1-alfa) + P3(1)*alfa;
                 P4(2)*(1-alfa) + P3(2)*alfa;
                 P4(3)*(1-alfa) + P3(3)*alfa]; % [X;Y;Z]

   for j = 1:nStringer
      CS(i).BotStringerXYZ{j} = (PStringers(:,j) + [0 0 ;+CS(i).tSkin/2 +CS(i).tSkin/2+CS(i).StringerHeight;0 0])'; % [X Y Z]
   end
end


if PLOT
   figure
   tg = uitabgroup; % tabgroup
   MainTab  = uitab(tg,'title','Wing Plot'); % build iith tab
   MainAXES = axes('Parent',MainTab);
   hold on
   grid on
   xlabel('x (m)')
   ylabel('z (m)')
   zlabel('y (m)')
   view(3)
   set(gca, 'Ydir', 'reverse')
   axis equal
   MainAXES.Clipping = 'off';
   plot3(xSpan, zeros(nSpan,1), zeros(nSpan,1),'red','LineWidth',2,'displayName','BeamAxis')

   for ii=1:2 
      for i=1:nSpan
         if ii == 1 
            ax = MainAXES;
         else     
            tab{i} = uitab(tg,'title',['CS ' num2str(i)]); %#ok<AGROW> % build iith tab
            ax     = axes('Parent',tab{i});
            view([90,0])
            hold on
            grid on
            xlabel('x (m)')
            ylabel('z (m)')
            zlabel('y (m)')
            axis equal
            ax.Clipping = 'off';

            CS(i).PlotAxes = ax;
         end

         plot3(ax,CS(i).XYZProfile(:,1), ...
                  CS(i).XYZProfile(:,3), ...
                  CS(i).XYZProfile(:,2),'blue .-')
   
         %%% wing box
         X = CS(i).WingBoxCornerXYZ([1 2 3 4 1],1);
         Y = CS(i).WingBoxCornerXYZ([1 2 3 4 1],2);
         Z = CS(i).WingBoxCornerXYZ([1 2 3 4 1],3);
         plot3(X,Z,Y,'black .-','linewidth',1)
   
         %%% wing box thickness
         % 1 ---- 2 
         % 4 ---- 3
         fill3(ax,[X([1 2]); X([2 1])], ...
               [Z([1 2]); Z([2 1])],...
               [Y([1 2])+CS(i).tSkin/2; Y([2 1])-CS(i).tSkin/2],'red','FaceAlpha', 0.25)
   
         fill3(ax,[X([4 3]); X([3 3])], ...
               [Z([4 3]); Z([3 4])],...
               [Y([4 3])+CS(i).tSkin/2; Y([3 4])-CS(i).tSkin/2],'red','FaceAlpha', 0.25)
   
         fill3(ax,[X([2 3]); X([3 2])], ...
               [Z([2 3])+CS(i).tWeb/2; Z([3 2])-CS(i).tWeb/2],...
               [Y([2 3]); Y([3 2])],'red','FaceAlpha', 0.5)
   
         fill3(ax,[X([1 4]); X([4 1])], ...
               [Z([1 4])+CS(i).tWeb/2; Z([4 1])-CS(i).tWeb/2],...
               [Y([1 4]); Y([4 1])],'red','FaceAlpha', 0.5)


         %%% stringers
         for jj=1:nStringer
            %%% top 
            X = CS(i).TopStringerXYZ{jj}(:,1);
            Y = CS(i).TopStringerXYZ{jj}(:,2);
            Z = CS(i).TopStringerXYZ{jj}(:,3);

            plot3(ax,X,Z,Y,'green .-')
            fill3(ax,[X([1 2]); X([2 1])], ...
                     [Z([1 2])+CS(i).StringerThickness/2; Z([2 1])-CS(i).StringerThickness/2],...
                     [Y([1 2]); Y([2 1])],'green','FaceAlpha', 0.25)

            %%% Bottom 
            X = CS(i).BotStringerXYZ{jj}(:,1);
            Y = CS(i).BotStringerXYZ{jj}(:,2);
            Z = CS(i).BotStringerXYZ{jj}(:,3);

            plot3(ax,X,Z,Y,'green .-')
            fill3(ax,[X([1 2]); X([2 1])], ...
                     [Z([1 2])+CS(i).StringerThickness/2; Z([2 1])-CS(i).StringerThickness/2],...
                     [Y([1 2]); Y([2 1])],'green','FaceAlpha', 0.25)
        
         end
      end
   end

end

%% reduced version of CS to make is clearer to students
CSred = CS;
CSred = rmfield(CSred,'beta');
CSred = rmfield(CSred,'ZYProfile_Norm');
CSred = rmfield(CSred,'ZProfileOffset');
CSred = rmfield(CSred,'YProfileOffset');

end








