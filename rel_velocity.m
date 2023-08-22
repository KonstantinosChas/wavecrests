function[crest_avg_vel,bulge_avg_vel,crest_vel,crest_vel_idx_top,vel_idx,ixvel,ixvel_1] = rel_velocity(bulge,crit_idx_bot,NUM,NUM_in,im_name)
% max_surf_el_idx_min,max_surf_el_idx_mtx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% moving frame of reference to compare crest/bulge evolution 
% find global (crest) and local (bulge) kinematic parameters
% : apply galilelean tranformations to x-dir to caclculate relative velocities
%% A COMPLICATED SCRIPT WITH MANY CLAUSES AND SUB STATEMENTS THAT ACCOUNT FOR MANY 
%% DIFFERENT CASES OF CRESTS
% Konstantinos Chasaspis 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc; 
    close all;
    dt = 0.002; % time step in sec
    cutoff_vel_idx = 0; % bulge additional space, randomly chosen
%--------------------------------------------------------------------------
    for i = NUM_in:NUM
        time.(im_name{i-NUM_in+1}) = dt.*i;  % in seconds
    end
%--------------------------------------------------------------------------
    for i = NUM_in:NUM
        bulgey.(im_name{i-NUM_in+1}) = bulge.(im_name{i-NUM_in+1})(:,2);
        bulgex.(im_name{i-NUM_in+1}) = bulge.(im_name{i-NUM_in+1})(:,1);
    end
%--------------------------------------------------------------------------
% find the field with the minimum elements 
    for i = NUM_in:NUM
        bulge_length.(im_name{i-NUM_in+1}) = length((bulgex.(im_name{i-NUM_in+1})));
    end
%--------------------------------------------------------------------------
    for i = NUM_in:NUM
        bulge_length_aux(i-NUM_in+1) = bulge_length.(im_name{i-NUM_in+1});  % put lengths in array 
    end
%--------------------------------------------------------------------------
    for i = NUM_in:NUM
        bulgey_mov_in.(im_name{i-NUM_in+1}) = bulgey.(im_name{i-NUM_in+1});
        bulgex_mov_in.(im_name{i-NUM_in+1}) = bulgex.(im_name{i-NUM_in+1});

    end
%--------------------------------------------------------------------------
%% Generate moving frame of reference 
% calculate crest speeds from extracted crests  
% v_crest = crest_vel = dx/dt 
%--------------------------------------------------------------------------
    crest_vel.(im_name{1}) = zeros(length(bulgex_mov_in.(im_name{1})),1);
    for i = NUM_in:NUM
% locate common y points between crests to calculate velocities
% when first crest is too short then all of them are cut short: SO BE
% CAREFUL!
% ixvel : the common indices in bulge.im(i), ixvel_1: the common indices in
% bulge(1) (for each bulge(i))
        [vel_idx.(im_name{i-NUM_in+1}),ixvel.(im_name{i-NUM_in+1}),ixvel_1.(im_name{i-NUM_in+1})] =...
            intersect(bulgey_mov_in.(im_name{i-NUM_in+1})(:), bulgey_mov_in.(im_name{1})(:));   
        
%         [vel_idx.(im_name{i-NUM_in+1}),ixvel.(im_name{i-NUM_in+1}),ixvel_1.(im_name{i-NUM_in+1})] =...
%             intersect(bulgey_mov_in.(im_name{i-NUM_in+1})(:), bulgey_mov_in.(im_name{i-NUM_in})(:));  
    end  
%--------------------------------------------------------------------------
% create arrays of crests with common y points 
% THIS COULD BE TRICKY DEPENDING ON WHAT FIRST EXTRACTED CREST LOOKS LIKE
% but in general is correct
    for i = NUM_in:NUM
        bulgex_mov_in_aux.(im_name{i-NUM_in+1}) = bulgex_mov_in.(im_name{i-NUM_in+1})(ixvel.(im_name{i-NUM_in+1}));
        bulgey_mov_in_aux.(im_name{i-NUM_in+1}) = bulgey_mov_in.(im_name{i-NUM_in+1})(ixvel.(im_name{i-NUM_in+1}));
%--------------------------------------------------------------------------
% calculate velocities in x-dir of each point. vel(x) = (x_n-x_1)/(t_n-t_1)  
        crest_vel.(im_name{i-NUM_in+1}) = ((bulgex_mov_in.(im_name{i-NUM_in+1})(ixvel.(im_name{i-NUM_in+1})) - ...
            bulgex_mov_in.(im_name{1})(ixvel_1.(im_name{i-NUM_in+1})))./(time.(im_name{i-NUM_in+1}) - time.(im_name{1})));  % THIS IS IMPORTANT
        
%         crest_vel.(im_name{i-NUM_in+1}) = ((bulgex_mov_in.(im_name{i-NUM_in+1})(ixvel.(im_name{i-NUM_in+1})) - ...
%             bulgex_mov_in.(im_name{i-NUM_in})(ixvel_1.(im_name{i-NUM_in+1})))./(time.(im_name{i-NUM_in+1}) - time.(im_name{i-NUM_in})));  % THIS IS IMPORTANT
    end
%--------------------------------------------------------------------------    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END OF PART 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Separate crests from bulges to calculate AVERAGE velocities. 
% There are two ways of doing this:
% 1. the latest crest is cut at the inflection point (minimum) found in
% bulge_slope.m 
% 2. the min crest_velocity array length is found and used as cut-off point that separates crests from bulges.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    for i = NUM_in:NUM
        crest_vel_length_aux.(im_name{i-NUM_in+1}) = numel(crest_vel.(im_name{i-NUM_in+1}));
    end
%--------------------------------------------------------------------------
    crest_vel_length = zeros(1,NUM-NUM_in);
    crit_idx_bot_mtx = zeros(1,NUM-NUM_in);
% put values in arrays:
    for i = NUM_in:NUM
        crest_vel_length(i-NUM_in+1) = crest_vel_length_aux.(im_name{i-NUM_in+1}); 
        crit_idx_bot_mtx(i-NUM_in+1) = min(crit_idx_bot.(im_name{i-NUM_in+1}));
        if isempty(crit_idx_bot_mtx(i-NUM_in+1)) == 1
            crit_idx_bot_mtx(i-NUM_in+1) = 0;
        else
            crit_idx_bot_mtx(i-NUM_in+1) = crit_idx_bot_mtx(i-NUM_in+1);
        end
    end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% 1.
    for i = NUM_in:NUM
          crest_vel_idx_top.(im_name{i-NUM_in+1}) = [];  %create structure
    end
    if numel(crit_idx_bot.(im_name{end})) == 1
        if crit_idx_bot.(im_name{end}) ~= 0   % if there is an infection point, do 1.
           
% find which point the index "crit_idx_bot" corresponds to in the the bulgey_mov_in_aux
%% this is an index:
            crest_vel_idx_top.(im_name{end}) = find(bulgey_mov_in_aux.(im_name{end}) == bulgey.(im_name{end})(crit_idx_bot.(im_name{end}) - cutoff_vel_idx));
% if it nothing is found: 
            if isempty(crest_vel_idx_top.(im_name{end}))
% search again within larger range of values:
                crest_vel_idx_top.(im_name{end}) = find(bulgey_mov_in_aux.(im_name{end})<= (bulgey.(im_name{end})(crit_idx_bot.(im_name{end}) - cutoff_vel_idx)+ 0.2) & ...
                bulgey_mov_in_aux.(im_name{end})>= (bulgey.(im_name{end})(crit_idx_bot.(im_name{end}) - cutoff_vel_idx) - 0.2));
            end
            if isempty(crest_vel_idx_top.(im_name{end}))
% seatch again within larger range of values:
                crest_vel_idx_top.(im_name{end}) = find(bulgey_mov_in_aux.(im_name{end})<= (bulgey.(im_name{end})(crit_idx_bot.(im_name{end}) - cutoff_vel_idx)+ 0.5) & ...
                bulgey_mov_in_aux.(im_name{end})>= (bulgey.(im_name{end})(crit_idx_bot.(im_name{end}) - cutoff_vel_idx) - 0.5));
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the value that crest_vel_idx_top(end) corresponds to:
            crest_vel_top_value = bulgey_mov_in_aux.(im_name{end})(crest_vel_idx_top.(im_name{end}));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isempty(crest_vel_top_value)
                disp('top value is empty');
                pause(1)
                crest_vel_top_value = bulgey_mov_in_aux.(im_name{end-1})(crest_vel_idx_top.(im_name{end}));
            end
%check again:
            if isempty(crest_vel_top_value)
                disp('top value is empty');
                pause(1)
                crest_vel_top_value = bulgey_mov_in_aux.(im_name{end-2})(crest_vel_idx_top.(im_name{end}));
            end
%--------------------------------------------------------------------------
% if 1. use crest_vel_top_value to find the indices that this value corresponds to at the other crests:
            for i = NUM_in:NUM-1
%--------------------------------------------------------------------------
                if isempty(crest_vel_top_value) == 0
                    crest_vel_idx_top.(im_name{i-NUM_in+1}) = find(crest_vel_top_value(1)+0.1 >= (bulgey_mov_in_aux.(im_name{i-NUM_in+1})) & ...
                        crest_vel_top_value(1)-0.1 <= (bulgey_mov_in_aux.(im_name{i-NUM_in+1})));
                end
%-------------------------------------------------------------------------- 
% search again if nothing is found
                if isempty(crest_vel_idx_top.(im_name{i-NUM_in+1}))
                 crest_vel_idx_top.(im_name{i-NUM_in+1}) = find(crest_vel_top_value(1)+0.2 >= (bulgey_mov_in_aux.(im_name{i-NUM_in+1})) & ...
                    crest_vel_top_value(1)-0.2 <= (bulgey_mov_in_aux.(im_name{i-NUM_in+1})));
                end
                if isempty(crest_vel_idx_top.(im_name{i-NUM_in+1}))
                crest_vel_idx_top.(im_name{i-NUM_in+1}) = find(crest_vel_top_value(1)+0.5 >= (bulgey_mov_in_aux.(im_name{i-NUM_in+1})) & ...
                    crest_vel_top_value(1)-0.5 <= (bulgey_mov_in_aux.(im_name{i-NUM_in+1})));
                end
% if still nothing found, assign values of index from other crests (not last):
                if isempty(crest_vel_idx_top.(im_name{i-NUM_in+1}))
                    crest_vel_idx_top.(im_name{i-NUM_in+1}) = crest_vel_idx_top.(im_name{i-NUM_in});
                end 
                if isempty(crest_vel_idx_top.(im_name{i-NUM_in+1}))
                    crest_vel_idx_top.(im_name{i-NUM_in+1}) = crest_vel_idx_top.(im_name{i-NUM_in-1});
                end
% a special clause that mostly applies when very short crests are detected (edges of image):  
                if crest_vel_idx_top.(im_name{i-NUM_in+1})>numel(bulgey_mov_in_aux.(im_name{i-NUM_in+1}))
                    crest_vel_idx_top.(im_name{i-NUM_in+1}) = numel(bulgey_mov_in_aux.(im_name{i-NUM_in+1}));
                end
            end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% 2.
        elseif crit_idx_bot.(im_name{end}) == 0  % if there is NOT an inflection point, do 2.
          disp('no critical points');
%--------------------------------------------------------------------------
            for i = NUM_in:NUM
                crest_vel_idx_top.(im_name{i-NUM_in+1}) = (crest_vel_length(i-NUM_in+1)) - cutoff_vel_idx;  % find non-zero minimum number of elements
            end
        end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% the same for numel>1
    elseif numel(crit_idx_bot.(im_name{end})) > 1
        if crit_idx_bot.(im_name{end}) ~= 0   % if there is an infection point, do 1.
% find which point the index "crit_idx_bot" corresponds to in the the bulgey_mov_in_aux
%% same as before:
%% this is an index:
            disp('more than 1 index found. choose one:')
            disp(num2str(crit_idx_bot.(im_name{end})));
            crit_idx_bot.(im_name{end}) = input('');
            crest_vel_idx_top.(im_name{end}) = find(bulgey_mov_in_aux.(im_name{end}) == bulgey.(im_name{end})(crit_idx_bot.(im_name{end}) - cutoff_vel_idx));
% if it nothing is found: 
            if isempty(crest_vel_idx_top.(im_name{end}))
% search again within larger range of values:
                crest_vel_idx_top.(im_name{end}) = find(bulgey_mov_in_aux.(im_name{end})<= (bulgey.(im_name{end})(crit_idx_bot.(im_name{end}) - cutoff_vel_idx)+ 0.4) & ...
                bulgey_mov_in_aux.(im_name{end})>= (bulgey.(im_name{end})(crit_idx_bot.(im_name{end}) - cutoff_vel_idx) - 0.4));
            end
            if isempty(crest_vel_idx_top.(im_name{end}))
% seatch again within larger range of values:
                crest_vel_idx_top.(im_name{end}) = find(bulgey_mov_in_aux.(im_name{end})<= (bulgey.(im_name{end})(crit_idx_bot.(im_name{end}) - cutoff_vel_idx)+ 1.5) & ...
                bulgey_mov_in_aux.(im_name{end})>= (bulgey.(im_name{end})(crit_idx_bot.(im_name{end}) - cutoff_vel_idx) - 1.5));
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the value that crest_vel_idx_top(end) corresponds to:
            crest_vel_top_value = bulgey_mov_in_aux.(im_name{end})(crest_vel_idx_top.(im_name{end}));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isempty(crest_vel_top_value)
                disp('top value is empty');
                pause(1)
                crest_vel_top_value = bulgey_mov_in_aux.(im_name{end-1})(crest_vel_idx_top.(im_name{end}));
            end
            if isempty(crest_vel_top_value)
                disp('top value is still empty');
                pause(1)
                crest_vel_top_value = 40;
%                 bulgey_mov_in_aux.(im_name{end-5})(crest_vel_idx_top.(im_name{end})); % or e.g. 70 arb in mm
            end
%--------------------------------------------------------------------------
% if 1. use crest_vel_top_value to find the indices that this value corresponds to at the other crests:
            for i = NUM_in:NUM-1
%--------------------------------------------------------------------------
                crest_vel_idx_top.(im_name{i-NUM_in+1}) = find(crest_vel_top_value(1)+0.1 >= (bulgey_mov_in_aux.(im_name{i-NUM_in+1})) & ...
                    crest_vel_top_value(1)-0.1 <= (bulgey_mov_in_aux.(im_name{i-NUM_in+1})));
%-------------------------------------------------------------------------- 
% search again if nothing is found
                if isempty(crest_vel_idx_top.(im_name{i-NUM_in+1}))
                 crest_vel_idx_top.(im_name{i-NUM_in+1}) = find(crest_vel_top_value(1)+0.2 >= (bulgey_mov_in_aux.(im_name{i-NUM_in+1})) & ...
                    crest_vel_top_value(1)-0.2 <= (bulgey_mov_in_aux.(im_name{i-NUM_in+1})));
                end
                if isempty(crest_vel_idx_top.(im_name{i-NUM_in+1}))
                crest_vel_idx_top.(im_name{i-NUM_in+1}) = find(crest_vel_top_value(1)+0.5 >= (bulgey_mov_in_aux.(im_name{i-NUM_in+1})) & ...
                    crest_vel_top_value(1)-0.5 <= (bulgey_mov_in_aux.(im_name{i-NUM_in+1})));
                end
% if still nothing found, assing values of index of previous crests (not last):
                if isempty(crest_vel_idx_top.(im_name{i-NUM_in+1}))
                    crest_vel_idx_top.(im_name{i-NUM_in+1}) = crest_vel_idx_top.(im_name{i-NUM_in});
                end 
                if isempty(crest_vel_idx_top.(im_name{i-NUM_in+1}))
                    crest_vel_idx_top.(im_name{i-NUM_in+1}) = crest_vel_idx_top.(im_name{i-NUM_in-1});
                end
            end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%2.
        elseif crit_idx_bot.(im_name{end}) == 0  % if there is NOT an inflection point, do 2.
%--------------------------------------------------------------------------
            for i = NUM_in:NUM
                crest_vel_idx_top.(im_name{i-NUM_in+1}) = (crest_vel_length(i-NUM_in+1)) - cutoff_vel_idx;  % find non-zero minimum number of elements
            end
%--------------------------------------------------------------------------
        end
    end
    if isempty(crest_vel_idx_top) 
        disp('no equal values found, calculating velocities from min crest length');
        crest_vel_idx_top.(im_name{end}) = (crest_vel_length) - cutoff_vel_idx;  % find non-zero minimum number of elements
    else
    end
%% last option: if all crests have critical point
    if all(crit_idx_bot_mtx)==1
% for each crest find index of crit_idx_bot 
        for i = NUM_in:NUM
            vel_idx.(im_name{i-NUM_in+1}) = find((bulgey.(im_name{i-NUM_in+1})(min(crit_idx_bot.(im_name{i-NUM_in+1})) - cutoff_vel_idx)) == (bulgey_mov_in_aux.(im_name{i-NUM_in+1})));
            if isempty(vel_idx.(im_name{i-NUM_in+1}))
                vel_idx.(im_name{i-NUM_in+1}) = find(bulgey_mov_in_aux.(im_name{i-NUM_in+1}) >= (bulgey.(im_name{i-NUM_in+1})(min(crit_idx_bot.(im_name{i-NUM_in+1})) - cutoff_vel_idx)+ 0.1) ...
                    & bulgey_mov_in_aux.(im_name{i-NUM_in+1}) <= (bulgey.(im_name{i-NUM_in+1})(min(crit_idx_bot.(im_name{i-NUM_in+1})) - cutoff_vel_idx)- 0.1));
            end
            if isempty(vel_idx.(im_name{i-NUM_in+1}))
                vel_idx.(im_name{i-NUM_in+1}) = find(bulgey_mov_in_aux.(im_name{i-NUM_in+1}) >= (bulgey.(im_name{i-NUM_in+1})(min(crit_idx_bot.(im_name{i-NUM_in+1})) - cutoff_vel_idx)+ 0.2) ...
                    & bulgey_mov_in_aux.(im_name{i-NUM_in+1}) <= (bulgey.(im_name{i-NUM_in+1})(min(crit_idx_bot.(im_name{i-NUM_in+1})) - cutoff_vel_idx)- 0.2));
            end
            if isempty(vel_idx.(im_name{i-NUM_in+1}))
                vel_idx.(im_name{i-NUM_in+1}) = find(bulgey_mov_in_aux.(im_name{i-NUM_in+1}) >= (bulgey.(im_name{i-NUM_in+1})(min(crit_idx_bot.(im_name{i-NUM_in+1})) - cutoff_vel_idx)+ 0.5) ...
                    & bulgey_mov_in_aux.(im_name{i-NUM_in+1}) <= (bulgey.(im_name{i-NUM_in+1})(min(crit_idx_bot.(im_name{i-NUM_in+1})) - cutoff_vel_idx)- 0.5));
            end
% if still nothing found, assign value of index of previous crest:
            if isempty(vel_idx.(im_name{i-NUM_in+1}))
                vel_idx.(im_name{i-NUM_in+1}) = vel_idx.(im_name{i-NUM_in});
            end
        end
    end
% finally, keep only one value:    
    for i = NUM_in:NUM
        crest_vel_idx_top.(im_name{i-NUM_in+1}) = min(crest_vel_idx_top.(im_name{i-NUM_in+1}));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END OF PART 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% CALCULATE AVERAGE VELOCITIES (1 VEL PER CREST AND 1 PER BULGE) 
%--------------------------------------------------------------------------
%% first find an index in each crest corresponding to cuty=50mm and use as cut-off
    for i = NUM_in:NUM
        cuty_down = 30;  % cutoff for crest vel index 
        cutoff_down.(im_name{i-NUM_in+1}) = find(bulgey_mov_in_aux.(im_name{i-NUM_in+1}) == cuty_down);
        cutoff_down.(im_name{i-NUM_in+1}) = min(cutoff_down.(im_name{i-NUM_in+1}));
        if isempty(cutoff_down.(im_name{i-NUM_in+1}))
            cutoff_down.(im_name{i-NUM_in+1}) = find(bulgey_mov_in_aux.(im_name{i-NUM_in+1}) <= (cuty_down+ 5) & ...
                bulgey_mov_in_aux.(im_name{i-NUM_in+1}) >= (cuty_down- 5));
            cutoff_down.(im_name{i-NUM_in+1}) = min(cutoff_down.(im_name{i-NUM_in+1}));
        end
        if isempty(cutoff_down.(im_name{i-NUM_in+1}))
            disp('no cutoff'); 
            cutoff_down.(im_name{i-NUM_in+1}) = find(min(bulgey_mov_in_aux.(im_name{i-NUM_in+1})));  
        end
        if isempty(cutoff_down.(im_name{i-NUM_in+1}))
            disp('no cutoff'); 
            cutoff_down.(im_name{i-NUM_in+1}) = 1;  
        end        
    end
    
% if all crests have critical points:
    if all(crit_idx_bot_mtx)==1
        for i = NUM_in:NUM
%--------------------------------------------------------------------------
            crest_avg_vel.(im_name{i-NUM_in+1}) = sum(crest_vel.(im_name{i-NUM_in+1})(cutoff_down.(im_name{i-NUM_in+1}):vel_idx.(im_name{i-NUM_in+1})))./...
                length(crest_vel.(im_name{i-NUM_in+1})(cutoff_down.(im_name{i-NUM_in+1}):vel_idx.(im_name{i-NUM_in+1})));
% (min(vel_idx.(im_name{i-NUM_in+1}))-cutoff_down.(im_name{i-NUM_in+1}));
            crest_avg_vel_low.(im_name{i-NUM_in+1}) = sum(crest_vel.(im_name{i-NUM_in+1})(1:cutoff_down.(im_name{i-NUM_in+1})))./cutoff_down.(im_name{i-NUM_in+1});
%--------------------------------------------------------------------------
            bulge_avg_vel.(im_name{i-NUM_in+1}) = sum(crest_vel.(im_name{i-NUM_in+1})((vel_idx.(im_name{i-NUM_in+1})):end))/...
                length(crest_vel.(im_name{i-NUM_in+1})(vel_idx.(im_name{i-NUM_in+1}):end));
%--------------------------------------------------------------------------
        end
% if not all crests have critical points:
    else  
        for i = NUM_in:NUM
            if crit_idx_bot.(im_name{end}) ~= 0 
                disp('not all crests have critical points');
                pause(0.5)
%--------------------------------------------------------------------------
            crest_avg_vel.(im_name{i-NUM_in+1}) = sum(crest_vel.(im_name{i-NUM_in+1})(cutoff_down.(im_name{i-NUM_in+1}):(crest_vel_idx_top.(im_name{i-NUM_in+1}))))./...
                    length(crest_vel.(im_name{i-NUM_in+1})(cutoff_down.(im_name{i-NUM_in+1}):(crest_vel_idx_top.(im_name{i-NUM_in+1}))));
% (min(crest_vel_idx_top.(im_name{i-NUM_in+1}))-cutoff_down.(im_name{i-NUM_in+1}));
                crest_avg_vel_low.(im_name{i-NUM_in+1}) = sum(crest_vel.(im_name{i-NUM_in+1})(1:cutoff_down.(im_name{i-NUM_in+1})))./cutoff_down.(im_name{i-NUM_in+1});
%--------------------------------------------------------------------------
                bulge_avg_vel.(im_name{i-NUM_in+1}) = sum(crest_vel.(im_name{i-NUM_in+1})(min(crest_vel_idx_top.(im_name{i-NUM_in+1})):end))/...
                    length(crest_vel.(im_name{i-NUM_in+1})(min(crest_vel_idx_top.(im_name{i-NUM_in+1})):end));
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
            elseif crit_idx_bot.(im_name{end}) == 0 
                crest_avg_vel.(im_name{i-NUM_in+1}) = sum(crest_vel.(im_name{i-NUM_in+1})(cutoff_down.(im_name{i-NUM_in+1}):(crest_vel_idx_top.(im_name{i-NUM_in+1}))))./...
                    length(crest_vel.(im_name{i-NUM_in+1})(cutoff_down.(im_name{i-NUM_in+1}):min(crest_vel_idx_top.(im_name{i-NUM_in+1}))));
% (min(crest_vel_idx_top.(im_name{end}))-cutoff_down.(im_name{i-NUM_in+1}));
                crest_avg_vel_low.(im_name{i-NUM_in+1}) = sum(crest_vel.(im_name{i-NUM_in+1})(1:cutoff_down.(im_name{i-NUM_in+1})))./(cutoff_down.(im_name{i-NUM_in+1}));
                bulge_avg_vel.(im_name{i-NUM_in+1}) = crest_avg_vel.(im_name{i-NUM_in+1});
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
%--------------------------------------------------------------------------
%     crest_avg_vel_tot = sum(crest_avg_vel_array(1:end))./(length(crest_avg_vel_array)); % not used
    for i = NUM_in:NUM
        crest_avg_vel.(im_name{i-NUM_in+1}) = abs(crest_avg_vel.(im_name{i-NUM_in+1}));
        bulge_avg_vel.(im_name{i-NUM_in+1}) = abs(bulge_avg_vel.(im_name{i-NUM_in+1}));
        crest_avg_vel_low.(im_name{i-NUM_in+1}) = abs(crest_avg_vel_low.(im_name{i-NUM_in+1}));

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% calculate another crest velocity based on more arbirtrary criteria (the
%     %lowest part of the crest, that is far away from breaking
%     
%     %find max surface elevation
%     for i = NUM_in:NUM
%         max_surf_el.(im_name{i-NUM_in+1}) = round(max(bulge.(im_name{i-NUM_in+1})(:,2)),2);
%         % calculate 0.4 (= arbitrary value)* of max surf elevation
%         max_surf_el_idx.(im_name{i-NUM_in+1}) = find(bulge.(im_name{i-NUM_in+1})(:,2) >= 0.7*max_surf_el.(im_name{i-NUM_in+1})(1) - 5 &...
%             bulge.(im_name{i-NUM_in+1})(:,2) <= 0.7*max_surf_el.(im_name{i-NUM_in+1})(1)+ 5);
%         
%         if isempty(max_surf_el_idx.(im_name{i-NUM_in+1}))
%             max_surf_el_idx.(im_name{i-NUM_in+1}) = 0;
%         end
%     end
%     for i = NUM_in:NUM
%         max_surf_el_idx_mtx(i-NUM_in+1) = max_surf_el_idx.(im_name{i-NUM_in+1})(1);
%     end
%      max_surf_el_idx_min = min(max_surf_el_idx_mtx(max_surf_el_idx_mtx>40 & max_surf_el_idx_mtx<300));
%     for i = NUM_in:NUM
%            if max_surf_el_idx_mtx(i-NUM_in+1) == 0
%                 max_surf_el_idx_mtx(i-NUM_in+1) = 40;
%            end
%            if max_surf_el_idx_mtx(i-NUM_in+1) <=10
%                max_surf_el_idx_mtx(i-NUM_in+1) = 40;  % arbitrary value
%            end
%            if max_surf_el_idx_mtx(i-NUM_in+1) > numel(crest_vel.(im_name{i-NUM_in+1}))
%               max_surf_el_idx_mtx(i-NUM_in+1) = round(numel(crest_vel.(im_name{i-NUM_in+1}))./6,0);  
% 
%            end
%     end
%     for i = NUM_in:NUM 

%         crest_avg_vel_alt.(im_name{i-NUM_in+1}) = abs(sum(crest_vel.(im_name{i-NUM_in+1})(1:max_surf_el_idx_mtx(i-NUM_in+1)))./...
%                     (max_surf_el_idx_mtx(i-NUM_in+1)));
%     end
%     for i = NUM_in:NUM - 10
%         crest_avg_vel_alt.(im_name{i-NUM_in+1}) = abs(sum(crest_vel.(im_name{i-NUM_in+1})(1:20))./20);
%     end
%     for i = NUM - 10 : NUM
%              crest_avg_vel_alt.(im_name{i-NUM_in+1}) = abs(sum(crest_vel.(im_name{i-NUM_in+1})(1:end))./numel(crest_vel.(im_name{i-NUM_in+1})));
%     end
%     for i = NUM_in:NUM
%         crest_avg_vel_alt.(im_name{i-NUM_in+1}) = abs(sum(crest_vel.(im_name{i-NUM_in+1}))./numel(crest_vel.(im_name{i-NUM_in+1})));
%         
%     end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%