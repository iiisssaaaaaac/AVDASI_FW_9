function out = column_buckling_prep_v1(CS, wingbox, section_no)
    % func for preparing the geometric data of top surface of winbox, for buckling calc (ie getting nodes, connectivity, etc)
    % normalised spanwise location
    % CS: wing cross section struct
    % wingbox: input into CS, contains geometric info for wingbox design
    % wing_semi_span: semi span of wing (obviously)
    % section_no: index of discrete wing node (0 at root)

    section_norm_pos = CS(section_no).xSpan/wingbox.semi_span;
     % adds wingbox corners, then top edge, then bottom edge to the nodes variable
    nodes = [CS(section_no).WingBoxCornerXYZ(1:2,:)];
    
    % forming the thickness list. different values for skin and stringers
    t = zeros(1 + (wingbox.Stringer),1);
    t(1,1) = CS(section_no).tSkin; % interp1(x.tSkin(:,1),x.tSkin(:,2),section_norm_pos); % skin
    t(2:end,1) = CS(section_no).StringerThickness; % interp1(x.StringerThickness(:,1),x.StringerThickness(:,2),section_norm_pos); % stringers
    
    for ii = 1:wingbox.Stringer
        nodes = [nodes; CS(section_no).TopStringerXYZ(ii)];
    end
    nodes = cell2mat(nodes);
    
    connectivity = [1 2]; % should always be the same 
    
    % adding the stringer connectivity
    j = (3:2:length(nodes)-1)';
    jj = (j+1);
    % so that the top stringers connectivity is in +ve direction:
    ii(1:wingbox.Stringer,1:2) = [jj(1:wingbox.Stringer) j(1:wingbox.Stringer)]; 
    ii = [ii; [j(wingbox.Stringer+1:end) jj(wingbox.Stringer+1:end)]];
    connectivity = [connectivity; ii];
    
    out.nodes = nodes;
    out.connectivity = connectivity;
    out.t = t;
    out.norm_pos = section_norm_pos;
end