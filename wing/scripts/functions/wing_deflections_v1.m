function [new_nodes, tip_deflection] = wing_deflections_v1(nodes, CS, const, LEF)

    nNodes = length(CS);
    wing_nodes = [[CS.xSpan]', zeros(nNodes,2)];
    connectivity = [[0; (1:nNodes)'], [(1:nNodes)'; 0]];
    connectivity([1,end],:) = [];


    % FORCE INPUT INTO FE
    F = [zeros(nNodes,1), LEF .* nodes.max_total, zeros(nNodes,4)];
    F = reshape(F',6*nNodes,1); % put into single list, thing, as needed...



    new_nodes = bar_anal(wing_nodes, connectivity, [CS.Iyy]', [CS.Izz]', [CS.Iyy]'+[CS.Izz]' , const.E, const.G, [CS.area]', F, false);

    tip_deflection = new_nodes(end,:) - wing_nodes(end,:);
    

 

    function new_nodes = bar_anal(nodes, connectivity, Iyy, Izz, J, E, G, A, F, PLOT)
        % A lot of this is a bit overcomplicated, but it is easier this way as I am copying my 3D solver from before
        % new_nodes = bar_anal(nodes, connectivity, Iyy, Izz, J, E, G, A, F)
        % 
        if PLOT
            plot_elements(nodes,connectivity,'3D', 'blue .-', true, true)
           
        end
        n_nodes = length(nodes);
        n_elements = length(connectivity);
        lengths = zeros(n_elements,1);
        vector = zeros(n_elements,3);
    
        for n = 1:n_elements
            penis = {[nodes(connectivity(n,1),:)] ; [nodes(connectivity(n,2),:)]};
            lengths(n) = element_length(penis{1},penis{2});
            vector(n,:) = penis{2}-penis{1};
        end
        
        % empty rotation matrix. Will remain identity matrix as no elements are angled differently
        for n = 1:n_elements
            r(1:3,1:3,n) = eye(3);
        end
        % convert to big rotation matrix
        for n = 1:n_elements
            r_big(:,:,n) = blkdiag(r(:,:,n),r(:,:,n),r(:,:,n),r(:,:,n));
        end
        % Local stiffies
        k_elem = zeros(12,12,n_elements);
        k_glob = k_elem;
        
        for n = 1:n_elements
            k_elem(:,:,n) = give_me_a_stiffy(E,A(n+1),Iyy(n+1),Izz(n+1),G,J(n+1),lengths(n));
            k_glob(:,:,n) = r_big(:,:,n)*k_elem(:,:,n)*r_big(:,:,n)';
        end
    
        % |Dof = description|
        % |1   = node 1 translation in x|
        % 
        % |2   = node 1 translation in y|
        % 
        % |3   = node 1 translation in z|
        % 
        % |4   = node 1 rotation    in x|
        % 
        % |5   = node 1 rotation    in y|
        % 
        % |6   = node 1 rotation    in z|
        % 
        % |7   = node 2 translation in x|
        % 
        % |8   = node 2 translation in y|
        % 
        % |etc...|
    
        % DoFs linked to each element
        dof_index = zeros(n_elements,12);
        for n = 1:n_elements
            dof_index(n,:) = ([(6*connectivity(n,1)-5):(6*connectivity(n,1)) (6*connectivity(n,2)-5):(6*connectivity(n,2))]);
        end
        % Big stiffy
        k_big = zeros(6*n_nodes);
        for n = 1:n_elements
            k_big(dof_index(n,:),dof_index(n,:)) = k_big(dof_index(n,:),dof_index(n,:)) + k_glob(:,:,n);
        end
    
        % Boundary and forces
        node_BC = [1]; %fully constrained nodes
        node_free = [1:n_nodes]; node_free = setdiff(node_free,node_BC); % any nodes without BCs are considered free
        dof_BC = zeros(1,6*length(node_BC));
        
        for n = 1:length(node_BC)
            dof_BC(1,6*n-5:6*n) = [6*node_BC(n)-5:6*node_BC(n)];
        end
        dof_free = [1:(6*n_nodes)]; dof_free = setdiff(dof_free,dof_BC);
        
        u = zeros(6*n_nodes,1); % initialise the displacement matrix
    
        k_big_free = k_big(dof_free,dof_free); % stiffness matrix containing everything but the BC nodes
        F_free = F(dof_free);
        
        % Runnning the fatty booooooyahhhhhhh
        u(dof_free) = k_big_free^-1*F_free;
    
        u_xy  = reshape(u,6,[]);
        Ux = u_xy(1,:)';
        Uy = u_xy(2,:)';
        Uz = u_xy(3,:)';
        Tx = u_xy(4,:)';
        Ty = u_xy(5,:)';
        Tz = u_xy(6,:)';
    
        new_nodes = nodes + [Ux, Uy, Uz];
    
        if PLOT
            plot_elements(new_nodes, connectivity, '3D', 'red o-', false, false)
            axis equal
            set(gca, 'Ydir', 'reverse')
        end
    
    end
    
    % Stiffness matrix function
    function Ke_Beam = give_me_a_stiffy(E,A,Iyy,Izz,G,J,L)
        % E - Young modularse
        % G - Torsional modularse
        % A - XXXSexional area
        % Iyy - sexy moment of area yy
        % Izz - same as above but other axis
        % J - same as above but torsional
        % L - Element length
        Kbar    = E*A/L*[...
                              1,  -1, ;
                             -1,   1, ];
        
        KShaft  = G*J/L*[...
                           1,  -1, ;
                          -1,   1, ];                
        
        KbeamXZ = E*Iyy/L^3*[...
                            12        -6*L    -12     -6*L;
                            -6*L   4*L^2  6*L  2*L^2 ;
                            -12        6*L    12       6*L;
                            -6*L   2*L^2  6*L  4*L^2];  
                        
        KbeamXY = E*Izz/L^3*[...
                            12         6*L    -12      6*L;
                             6*L   4*L^2  -6*L 2*L^2 ;
                            -12        -6*L   12       -6*L;
                            6*L    2*L^2  -6*L  4*L^2];  
                         
        Ke_Beam = zeros(12,12);
        Ke_Beam([1 7],[1 7])             = Kbar;
        Ke_Beam([4 10],[4 10])           = KShaft;
        Ke_Beam([2 6 8 12],[2 6 8 12])   = KbeamXY;
        Ke_Beam([3 5 9 11],[3 5 9 11])   = KbeamXZ;
    end
    
    function l = element_length(a,b)
        x = b-a;
        x = x.^2;
        l = sqrt(sum(x));
    end

end