function [stresses, Mz] = axial_stress_v1(CS, bending, nSpan, load_case, nodes, plot)
    % AXIAL STRESS IN RIBS
    % INPUTS:
    % moments <- shear_moment function (bending)
    % second moments of area <- wingbox_weight (CS)
    % position from neutral axis <- wingbox_weight (CS) (neutral axis at 0 due to symmetry)


    
    for span_pos = 1:nSpan
        % CS(span_pos).small_nodes(:,1) is the y position of each small node in the rib discretisation
        % bending.max_total(:,load_case) is the bending moments for a given load case (or cases)

        stresses(span_pos).axial = stress_func(CS(span_pos).Izz, bending.max_total(span_pos, load_case), CS(span_pos).small_nodes(:,1));
        stresses(span_pos).max_axial = max(stresses(span_pos).axial,[],1);
    end

    if plot
        figure(Name='Axial stress due to bending')
        hold on
        for span_pos = 1:nSpan
            scatter3(CS(span_pos).small_nodes(:,2), ones(1,length(CS(span_pos).small_nodes(:,1))).* nodes.pos(span_pos),...
            CS(span_pos).small_nodes(:,1), [], stresses(span_pos).axial/1e6) % axial/1e6 THESE ARE NOT THE RIGHT AXES... JUST MAKING IT LOOK RIGHT WHEN PLOTTED
        end
        font_size = 23;
        fontsize(20, "points") % this is actually adjusting the size of the numbers
        axis equal
        axis([-3 3 0 max(nodes.pos) -3 3])
        view([-125 25])
        yticks(0:2:max(nodes.pos))
        set(gcf,'color','w')
        box off
        xlabel('x [m]','Interpreter','latex','fontsize', font_size)
        ylabel('z [m]','Interpreter','latex','fontsize', font_size)
        zlabel('y [m]','Interpreter','latex','fontsize', font_size)
        set(gca, 'Ydir', 'reverse')
        c = colorbar;
        c.Location = 'eastoutside';
        c.Label.String = 'Axial Stress [MPa]';
    end


    %% Validation
        Mz = zeros(nSpan,1);
        for span_pos = 1:nSpan
            for ii = 1:size(stresses(span_pos).axial)
                Mz(span_pos) = Mz(span_pos) - (stresses(span_pos).axial(ii) * CS(span_pos).small_nodes(ii,1) * CS(span_pos).small_node_areas(ii,1));
            end
        end

    
    %% function to calcualte the stress given node positions, moments, moments of area
    function stress = stress_func(I, M, yi)
        stress = - yi.*(M./I);
    end

end