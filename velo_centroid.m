function[velo_ycm_dt,velo_xcm_dt,velo_xcm,velo_ycm,velo_crest_dt,velo_crest_top_dt,velo_crest_top_dt_smooth,velo_crest_dt_smooth,velo_crest_tip_dt,velo_crest_tip_dt_smooth] = velo_centroid(bulge,max_idx,crit_idx_top,xcm,ycm,NUM,NUM_in,im_name)
    clc; 
    close all;
    dt = 0.002; % time step in sec
%--------------------------------------------------------------------------
    for i = NUM_in:NUM
        time.(im_name{i-NUM_in+1}) = dt.*i;  % in seconds
    end
%--------------------------------------------------------------------------
    if isnan(xcm.(im_name{1})) == 0
        for i = NUM_in:NUM
            velo_xcm.(im_name{i-NUM_in+1}) = abs((xcm.(im_name{i-NUM_in+1})-xcm.(im_name{1}))./(time.(im_name{i-NUM_in+1}) - time.(im_name{1})));
%--------------------------------------------------------------------------
        % we allow directionality in z for now:
%--------------------------------------------------------------------------
            velo_ycm.(im_name{i-NUM_in+1}) = ((ycm.(im_name{i-NUM_in+1})-ycm.(im_name{1}))./(time.(im_name{i-NUM_in+1}) - time.(im_name{1})));
        end
    elseif isnan(xcm.(im_name{1})) == 1
        if isnan(xcm.(im_name{2})) == 0
            for i = NUM_in:NUM
                velo_xcm.(im_name{i-NUM_in+1}) = abs((xcm.(im_name{i-NUM_in+1})-xcm.(im_name{2}))./(time.(im_name{i-NUM_in+1}) - time.(im_name{2})));
        % velocity between two crests (tn-tn-1)
%                 velo_xcm_dt.(im_name{i-NUM_in+1}) = abs((xcm.(im_name{i-NUM_in+2})-xcm.(im_name{i-NUM_in+1}))./(time.(im_name{i-NUM_in+2}) - time.(im_name{i-NUM_in+1})));
%--------------------------------------------------------------------------
        % we allow directionality in z for now:
%--------------------------------------------------------------------------
                velo_ycm.(im_name{i-NUM_in+1}) = ((ycm.(im_name{i-NUM_in+1})-ycm.(im_name{2}))./(time.(im_name{i-NUM_in+1}) - time.(im_name{2})));  
            end
        elseif isnan(xcm.(im_name{2})) == 1
            for i = NUM_in:NUM
                velo_xcm.(im_name{i-NUM_in+1}) = abs((xcm.(im_name{i-NUM_in+1})-xcm.(im_name{3}))./(time.(im_name{i-NUM_in+1}) - time.(im_name{3})));
        % velocity between two crests (tn-tn-1)
%                 velo_xcm_dt.(im_name{i-NUM_in+1}) = abs((xcm.(im_name{i-NUM_in+2})-xcm.(im_name{i-NUM_in+1}))./(time.(im_name{i-NUM_in+2}) - time.(im_name{i-NUM_in+1})));
%--------------------------------------------------------------------------
        % we allow directionality in z for now:
%--------------------------------------------------------------------------
                velo_ycm.(im_name{i-NUM_in+1}) = ((ycm.(im_name{i-NUM_in+1})-ycm.(im_name{3}))./(time.(im_name{i-NUM_in+1}) - time.(im_name{3})));  
            end
        end
    end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    % instantaneous velocity
    for i = NUM_in:NUM-1
            velo_xcm_dt.(im_name{i-NUM_in+1}) = abs((xcm.(im_name{i-NUM_in+2})-xcm.(im_name{i-NUM_in+1}))./(time.(im_name{i-NUM_in+2}) - time.(im_name{i-NUM_in+1})));
            velo_ycm_dt.(im_name{i-NUM_in+1}) = ((ycm.(im_name{i-NUM_in+2})-ycm.(im_name{i-NUM_in+1}))./(time.(im_name{i-NUM_in+2}) - time.(im_name{i-NUM_in+1})));
    end
%--------------------------------------------------------------------------
    % instantaneous centroid velocity
    for i = NUM_in:NUM
        mean_crest_xcm_idx.(im_name{i-NUM_in+1}) = find(bulge.(im_name{i-NUM_in+1})(:,2) == mean(bulge.(im_name{i-NUM_in+1})(round(1:end,1),2)));
    if isempty(mean_crest_xcm_idx.(im_name{i-NUM_in+1}))
        mean_crest_xcm_idx.(im_name{i-NUM_in+1}) = find(bulge.(im_name{i-NUM_in+1})(:,2) <= mean(bulge.(im_name{i-NUM_in+1})(round(1:end,1),2)) + 0.2 & bulge.(im_name{i-NUM_in+1})(:,2) >= mean(bulge.(im_name{i-NUM_in+1})(round(1:end/10,1),2)) - 0.2 );
    end
    end
    for i = NUM_in:NUM
        if mean_crest_xcm_idx.(im_name{i-NUM_in+1})>0
        mean_crest_xcm_idx_mtx(i-NUM_in+1) = mean_crest_xcm_idx.(im_name{i-NUM_in+1})(end);
        else
          mean_crest_xcm_idx_mtx(i-NUM_in+1) = 0;
        end
    end
    for i = NUM_in:NUM
        if mean_crest_xcm_idx.(im_name{i-NUM_in+1})>0
            mean_crest_xcm(i-NUM_in+1,1) = bulge.(im_name{i-NUM_in+1})(mean_crest_xcm_idx_mtx(i-NUM_in+1),1);
        else
            mean_crest_xcm(i-NUM_in+1,1) = 0;
        end
            
    end
    for i = NUM_in:NUM-1
        if mean_crest_xcm_idx.(im_name{i-NUM_in+1})>0
            velo_crest_dt.(im_name{i-NUM_in+1}) = abs((mean_crest_xcm(i-NUM_in+2,1)- mean_crest_xcm(i-NUM_in+1,1))./(time.(im_name{i-NUM_in+2}) - time.(im_name{i-NUM_in+1})));
        else
            velo_crest_dt.(im_name{i-NUM_in+1}) =[];
        end
    end
    % filter : not tested window and filter choise
    for i = NUM_in:NUM-1
        if mean_crest_xcm_idx.(im_name{i-NUM_in+1})>0
            velo_crest_dt_smooth.(im_name{i-NUM_in+1}) = smooth(velo_crest_dt.(im_name{i-NUM_in+1}),0.25,'rloess');
        else
            velo_crest_dt_smooth.(im_name{i-NUM_in+1}) = [];
            
        end
    end
%--------------------------------------------------------------------------
    % instant velocity of max surf elevation
    for i = NUM_in:NUM-1
        if mean_crest_xcm_idx.(im_name{i-NUM_in+1})>0
            velo_crest_top_dt.(im_name{i-NUM_in+1}) = abs((bulge.(im_name{i-NUM_in+2})(end,1) - bulge.(im_name{i-NUM_in+1})(end,1))./(time.(im_name{i-NUM_in+2}) - time.(im_name{i-NUM_in+1})));
        else
            velo_crest_top_dt.(im_name{i-NUM_in+1}) = 0;
        end
    end
    % filter : not tested window and filter choise
    for i = NUM_in:NUM-1
        if mean_crest_xcm_idx.(im_name{i-NUM_in+1})>0
            velo_crest_top_dt_smooth.(im_name{i-NUM_in+1}) = smooth(velo_crest_top_dt.(im_name{i-NUM_in+1}),'moving',3);
        else
            velo_crest_top_dt_smooth.(im_name{i-NUM_in+1}) = 0;
        end
    end
%--------------------------------------------------------------------------    
    % instant velocity of bulge tip 

    
    for i = NUM_in:NUM-1
        if crit_idx_top.(im_name{i-NUM_in+1})(end)>0 
            if crit_idx_top.(im_name{i-NUM_in+2})(end)>0
                velo_crest_tip_dt.(im_name{i-NUM_in+1}) = abs((bulge.(im_name{i-NUM_in+2})(crit_idx_top.(im_name{i-NUM_in+2})(end),1) - bulge.(im_name{i-NUM_in+1})(crit_idx_top.(im_name{i-NUM_in+1})(end),1))./(time.(im_name{i-NUM_in+2}) - time.(im_name{i-NUM_in+1})));
            elseif crit_idx_top.(im_name{i-NUM_in+2}) == 0
                velo_crest_tip_dt.(im_name{i-NUM_in+1}) = 0;
            end
        else 
            disp('sorry cannot calculate tip velocities');
            velo_crest_tip_dt.(im_name{i-NUM_in+1}) = 0;
        end
    end
    % filter : not tested window and filter choise
    for i = NUM_in:NUM-1
        if crit_idx_top.(im_name{i-NUM_in+1})>0 
            if crit_idx_top.(im_name{i-NUM_in+2})>0
                velo_crest_tip_dt_smooth.(im_name{i-NUM_in+1}) = smooth(velo_crest_tip_dt.(im_name{i-NUM_in+1}),0.25,'rloess');
            elseif crit_idx_top.(im_name{i-NUM_in+2}) == 0
                velo_crest_tip_dt_smooth.(im_name{i-NUM_in+1}) = 0;
            end
        else
            velo_crest_tip_dt_smooth.(im_name{i-NUM_in+1}) = 0;
        end
    end
    
 
    
    
%     for i = NUM_in:NUM-1
%         velo_crest_avg_dt.(im_name{i-NUM_in+1}) = abs((mean_crest_xcm(i-NUM_in+2,1)- mean_crest_xcm(1,1))./(time.(im_name{i-NUM_in+2}) - time.(im_name{1})));
%     end    
%         for i = NUM_in:NUM-1
%         velo_crest_top_avg_dt.(im_name{i-NUM_in+1}) = abs((bulge.(im_name{i-NUM_in+2})(end,1) - bulge.(im_name{1})(end,1))./(time.(im_name{i-NUM_in+2}) - time.(im_name{1})));
%     end
    
end