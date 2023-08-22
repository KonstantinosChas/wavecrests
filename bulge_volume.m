function [y_hypoten,hor_points,vert_points,num_hor_grid,bulge_area,bulge_area_points] = bulge_volume(bulge,crit_idx_bot,ifpx,NUM,NUM_in,im_name)
% function [hor_points,vert_points,num_hor_grid] = bulge_volume(bulge,crit_idx_bot,ifpx,NUM,NUM_in,im_name)

% 
% : Calculates area of bulge for crests with or without critical points
% Creates closed area with fictional bottom and side borders
% using bot criterion or inflection point AND max point of bulge, 
% both of them found in other functions

% find index corresponding to ifpx (which comes from fitted curves...complicated as F&CK and as usual)
   for i = NUM_in:NUM
        if crit_idx_bot.(im_name{i-NUM_in+1}) ~= 0  
            ind.(im_name{i-NUM_in+1}) = 1;
        elseif crit_idx_bot.(im_name{i-NUM_in+1}) == 0  
            ind.(im_name{i-NUM_in+1}) = find(bulge.(im_name{i-NUM_in+1})(:,1) <= ifpx.(im_name{i-NUM_in+1}) + 0.4 & ...
                bulge.(im_name{i-NUM_in+1})(:,1) >= ifpx.(im_name{i-NUM_in+1}) - 0.4);
            ind.(im_name{i-NUM_in+1}) = ind.(im_name{i-NUM_in+1})(1);
        end
   end
%  create matrix of horizontal line points starting from bulge crit point and ending in vertical line
    for i = NUM_in:NUM
        if crit_idx_bot.(im_name{i-NUM_in+1}) ~= 0  
    %check if there are multiple y max points
    %            if bulge.(im_name{i-NUM_in+1})(end,2) == bulge.(im_name{i-NUM_in+1})(end-1,2)
    %                hor_points.(im_name{i-NUM_in+1}) = (bulge.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1}),1):step:(bulge.(im_name{i-NUM_in+1})(end-1,1)))';
    %                num_hor_grid.(im_name{i-NUM_in+1}) = numel(hor_points.(im_name{i-NUM_in+1}));
    %            elseif bulge.(im_name{i-NUM_in+1})(end,2) ~= bulge.(im_name{i-NUM_in+1})(end-1,2)
            if bulge.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(end),1) > bulge.(im_name{i-NUM_in+1})(end,1) % clause when overturning results in bulge being far in front of top of wave (weird clause - 19/04/23)
                step = -0.1;
                hor_points.(im_name{i-NUM_in+1}) = (bulge.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(end),1):step:(bulge.(im_name{i-NUM_in+1})(end,1)))';
            else 
                step = 0.1;
                hor_points.(im_name{i-NUM_in+1}) = (bulge.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(end),1):step:(bulge.(im_name{i-NUM_in+1})(end,1)))';
            end
            num_hor_grid.(im_name{i-NUM_in+1}) = numel(hor_points.(im_name{i-NUM_in+1}));
    %            end
        elseif crit_idx_bot.(im_name{i-NUM_in+1}) == 0  
            step = 0.1;
           hor_points.(im_name{i-NUM_in+1}) = (bulge.(im_name{i-NUM_in+1})(ind.(im_name{i-NUM_in+1}),1):step:(bulge.(im_name{i-NUM_in+1})(end,1)))';
           num_hor_grid.(im_name{i-NUM_in+1}) = numel(hor_points.(im_name{i-NUM_in+1}));
        end
    end
    for i = NUM_in:NUM
        if crit_idx_bot.(im_name{i-NUM_in+1}) ~= 0  
            for j = 1:num_hor_grid.(im_name{i-NUM_in+1}) 
                hor_points.(im_name{i-NUM_in+1})(j,2) = bulge.(im_name{i-NUM_in+1})(max(crit_idx_bot.(im_name{i-NUM_in+1})),2);
            end
        elseif crit_idx_bot.(im_name{i-NUM_in+1}) == 0  
           for j = 1:num_hor_grid.(im_name{i-NUM_in+1}) 
               hor_points.(im_name{i-NUM_in+1})(j,2) = bulge.(im_name{i-NUM_in+1})(ind.(im_name{i-NUM_in+1}),2);
           end            
        end
    end
    %same for vertical    
    for i = NUM_in:NUM
        if crit_idx_bot.(im_name{i-NUM_in+1}) ~= 0  

    %           if bulge.(im_name{i-NUM_in+1})(end,2) == bulge.(im_name{i-NUM_in+1})(end-1,2)
    %                vert_points.(im_name{i-NUM_in+1})(:,2) = (bulge.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1}),2):step:(bulge.(im_name{i-NUM_in+1})(end-1,2)));
    %                num_vert_grid.(im_name{i-NUM_in+1}) = numel(vert_points.(im_name{i-NUM_in+1})(:,2));
    %           elseif bulge.(im_name{i-NUM_in+1})(end,2) ~= bulge.(im_name{i-NUM_in+1})(end-1,2)
            if bulge.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(end),2) > bulge.(im_name{i-NUM_in+1})(end,2) % same as in hor_points
                step = -0.1;
                vert_points.(im_name{i-NUM_in+1})(:,2) = (bulge.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(end),2):step:(bulge.(im_name{i-NUM_in+1})(end,2)));
            else 
                step = 0.1;
                vert_points.(im_name{i-NUM_in+1})(:,2) = (bulge.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(end),2):step:(bulge.(im_name{i-NUM_in+1})(end,2)));
            end
                num_vert_grid.(im_name{i-NUM_in+1}) = numel(vert_points.(im_name{i-NUM_in+1})(:,2));
    %           end    
        elseif crit_idx_bot.(im_name{i-NUM_in+1}) == 0  
            step = 0.1;
            vert_points.(im_name{i-NUM_in+1})(:,2) = (bulge.(im_name{i-NUM_in+1})(ind.(im_name{i-NUM_in+1}),2):step:bulge.(im_name{i-NUM_in+1})(end,2));
            num_vert_grid.(im_name{i-NUM_in+1}) = numel(vert_points.(im_name{i-NUM_in+1})(:,2));            
        end
    end
    for i = NUM_in:NUM
        if crit_idx_bot.(im_name{i-NUM_in+1}) ~= 0  

    %           if bulge.(im_name{i-NUM_in+1})(end,2) == bulge.(im_name{i-NUM_in+1})(end-1,2)
    %                for j = 1:num_vert_grid.(im_name{i-NUM_in+1}) 
    %                     vert_points.(im_name{i-NUM_in+1})(j,1) = (bulge.(im_name{i-NUM_in+1})(end-1,1));
    %                end
    %           elseif bulge.(im_name{i-NUM_in+1})(end,2) ~= bulge.(im_name{i-NUM_in+1})(end-1,2)
            for j = 1:num_vert_grid.(im_name{i-NUM_in+1}) 
                vert_points.(im_name{i-NUM_in+1})(j,1) = (bulge.(im_name{i-NUM_in+1})(end,1));
            end
    %           end
        elseif crit_idx_bot.(im_name{i-NUM_in+1}) == 0  
            for j = 1:num_vert_grid.(im_name{i-NUM_in+1}) 
                vert_points.(im_name{i-NUM_in+1})(j,1) = bulge.(im_name{i-NUM_in+1})(end,1);
            end
        end
    end 
    % find points lying on line y=tan(theta)x between (z_bulge,x_bulge) and
    % (z_max,x_max)  or(z_inf,x_inf) and maxs
    % part of these calculations are repeated in theta.m function
    % this is used instead of hor+vert points to calculate bulge volume
    %     if crit_idx_bot.(im_name{i-NUM_in+1}) ~= 0  
    for i = NUM_in:NUM
        dx.(im_name{i-NUM_in+1}) = hor_points.(im_name{i-NUM_in+1})(end,1) - hor_points.(im_name{i-NUM_in+1})(1,1);
        dy.(im_name{i-NUM_in+1}) = vert_points.(im_name{i-NUM_in+1})(end,2)- vert_points.(im_name{i-NUM_in+1})(1,2);
    end
    for i = NUM_in:NUM
        slopetheta.(im_name{i-NUM_in+1}) = dy.(im_name{i-NUM_in+1})/dx.(im_name{i-NUM_in+1});
    end
    for i = NUM_in:NUM
    %y = ax 
        hypotenuse.(im_name{i-NUM_in+1}) = @(x) slopetheta.(im_name{i-NUM_in+1})*x;
    % b = y1 - ax1;
        beta.(im_name{i-NUM_in+1}) = hor_points.(im_name{i-NUM_in+1})(1,2) - hypotenuse.(im_name{i-NUM_in+1})(hor_points.(im_name{i-NUM_in+1})(1,1));
    % y = ax +b;
        y_hypoten.(im_name{i-NUM_in+1}) = hypotenuse.(im_name{i-NUM_in+1})(hor_points.(im_name{i-NUM_in+1})(1:end,1)) + beta.(im_name{i-NUM_in+1});
    end
    %    elseif crit_idx_bot.(im_name{i-NUM_in+1}) == 0  
    %        for i = NUM_in:NUM
    %            dx.(im_name{i-NUM_in+1}) = 0; slopetheta.(im_name{i-NUM_in+1}) = 0; beta.(im_name{i-NUM_in+1}) = 0;
    %            dy.(im_name{i-NUM_in+1}) = 0; hypotenuse.(im_name{i-NUM_in+1}) = 0;
    %            y_hypoten.(im_name{i-NUM_in+1}) = 0;
    %        end
    %     end
    % concacenate bulge points and hor/vert in matrix   
    for i = NUM_in:NUM
        if crit_idx_bot.(im_name{i-NUM_in+1}) ~= 0  
        %                 bulge_area_points.(im_name{i-NUM_in+1})(:,1) =[bulge.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1}):end,1);hor_points.(im_name{i-NUM_in+1})(:,1);vert_points.(im_name{i-NUM_in+1})(:,1)];
        %                 bulge_area_points.(im_name{i-NUM_in+1})(:,2) =[bulge.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1}):end,2);hor_points.(im_name{i-NUM_in+1})(:,2);vert_points.(im_name{i-NUM_in+1})(:,2)];
            bulge_area_points.(im_name{i-NUM_in+1})(:,1) =[bulge.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(end):end,1);hor_points.(im_name{i-NUM_in+1})(:,1)];
            disp(i);
            pause(0.5)
            bulge_area_points.(im_name{i-NUM_in+1})(:,2) =[bulge.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(end):end,2);y_hypoten.(im_name{i-NUM_in+1})];

        elseif crit_idx_bot.(im_name{i-NUM_in+1}) == 0  
            bulge_area_points.(im_name{i-NUM_in+1})(:,1) =[(bulge.(im_name{i-NUM_in+1})(ind.(im_name{i-NUM_in+1}):end,1));hor_points.(im_name{i-NUM_in+1})(:,1)];
            disp(i);
            pause(0.5)
            bulge_area_points.(im_name{i-NUM_in+1})(:,2) =[(bulge.(im_name{i-NUM_in+1})(ind.(im_name{i-NUM_in+1}):end,2));y_hypoten.(im_name{i-NUM_in+1})];  
        end 
   end

% finally, find area using boundary function
% --> alphashape
    for i = NUM_in:NUM
        [bulge_area_index.(im_name{i-NUM_in+1}),bulge_area.(im_name{i-NUM_in+1})] = ...
                boundary(bulge_area_points.(im_name{i-NUM_in+1})(:,1),bulge_area_points.(im_name{i-NUM_in+1})(:,2),0.75);
    end






end










