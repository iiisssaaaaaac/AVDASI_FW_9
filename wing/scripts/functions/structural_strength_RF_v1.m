function [RF, allowable, CS, RF_crit] = structural_strength_RF_v1(CS, const, stresses, wingbox, nSpan)
    % finding RFs for different failure modes given the material strength
    no_load_cases               = width(stresses(1).axial);
    disc_density                = 1000;
    
    RF                          = zeros(nSpan, 4, no_load_cases);
    allow_yield                 = zeros(nSpan, 1);
    RF_yield                    = zeros(nSpan, no_load_cases);
    allow_skin_panel_buckling   = zeros(nSpan, 1);
    RF_skin_panel_buckling      = zeros(nSpan, no_load_cases);
    allow_plate_between_string  = zeros(nSpan, 1);
    RF_plate_between_string     = zeros(nSpan, no_load_cases);
    allow_string_plate          = zeros(nSpan, 1);
    RF_string_plate             = zeros(nSpan, no_load_cases);

    b_sk                        = zeros(nSpan,1);

    for span_pos = 1:nSpan
        
    % Yield from axial stress
        allow_yield(span_pos,:) = const.YieldStrength;
        RF_yield(span_pos,:)    = allow_yield(span_pos,:) ./ stresses(span_pos).max_axial;
        
    % Column buckling of top and bottom stiffened panels, finding critical stress value to cause it (finding the allowable)
        % first find radius of gyration of whole sitffened skin (top = bottom)
        gyration_calc           = column_buckling_prep_v1(CS, wingbox, span_pos);
        [Iyy, Izz, ~, area]     = x_section_anal(gyration_calc.nodes(:,2:3), gyration_calc.connectivity, gyration_calc.t, disc_density, 0);
        CS(span_pos).R_g        = sqrt(min(Iyy, Izz)/area);

        % CS(2).xSpan is L, length (span direction) of the stiffened panel for all. think about it and it makes sense :)

        allow_skin_panel_buckling(span_pos,:)   = (pi^2 * const.E / (CS(2).xSpan / CS(span_pos).R_g)^2);
        RF_skin_panel_buckling(span_pos,:)      = allow_skin_panel_buckling(span_pos,:) ./ stresses(span_pos).max_axial;
    
    % Plate buckling in between stringers, like above, but the other direction
        b_sk(span_pos,1)                        = CS(span_pos).TopStringerXYZ{1,2}(2,3) - CS(span_pos).TopStringerXYZ{1,1}(1,3);
        allow_plate_between_string(span_pos,:)  = (CS(span_pos).tSkin/b_sk(span_pos))^2 * (4 * pi^2 * const.E)/(12 * (1 - const.nu^2));
        RF_plate_between_string(span_pos,:)     = allow_plate_between_string(span_pos,:) ./ stresses(span_pos).max_axial;

    % Stringers buckling as plates (stringers themselves buckling between ribs)
        allow_string_plate(span_pos,:)  = (CS(span_pos).StringerThickness/CS(span_pos).StringerHeight)^2 * (0.43 * pi^2 * const.E)/(12 * (1 - const.nu^2));
        RF_string_plate(span_pos,:)     = allow_string_plate(span_pos,:) ./ stresses(span_pos).max_axial;
    end
    
    allowable = [allow_yield, allow_skin_panel_buckling, allow_plate_between_string, allow_string_plate];

    for load_case = 1:no_load_cases
        RF(:,:,load_case) = [RF_yield(:,load_case), RF_skin_panel_buckling(:,load_case), RF_plate_between_string(:,load_case), RF_string_plate(:,load_case)];
    end
    RF_crit = min(RF,[],3);
end