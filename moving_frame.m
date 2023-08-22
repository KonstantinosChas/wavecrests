function [bulgex_mov,bulgey_mov] = moving_frame(crest_avg_vel,bulgex,bulgey,NUM,NUM_in,im_name)
% calculate moving frame of reference coordinates:     
% X' = X-VT
%% y'=y, x'=x-crest_vel*time

    clc; 
    close all;
    dt = 0.002; % time step in sec
%--------------------------------------------------------------------------
    for i = NUM_in:NUM
        time.(im_name{i-NUM_in+1}) = dt.*i;  % in seconds
    end
    bulgex_mov.(im_name{1}) = bulgex.(im_name{1});
    bulgey_mov.(im_name{1}) = bulgey.(im_name{1});

%--------------------------------------------------------------------------    
    for i = NUM_in+1:NUM
%         crestx_dif.(im_name{i}) = (crest_avg_vel_tot.*time.(im_name{i})); % crest_vel*time
%% MINUS MINUS MINUS MINUS!!!
%--------------------------------------------------------------------------
        crestx_dif.(im_name{i-NUM_in+1}) = (-crest_avg_vel.(im_name{i-NUM_in+1}).*time.(im_name{i-NUM_in+1})); % crest_vel*time
%         crestx_dif_low.(im_name{i-NUM_in+1}) = (crest_avg_vel_low.(im_name{i-NUM_in+1}).*time.(im_name{i-NUM_in+1})); % crest_vel*time
%--------------------------------------------------------------------------
%         bulgex_mov.(im_name{i-NUM_in+1}) = bulgex_mov_in_aux.(im_name{i-NUM_in+1}) - crestx_dif.(im_name{i-NUM_in+1}); % X' = X-VT
%         bulgey_mov.(im_name{i-NUM_in+1}) = bulgey_mov_in_aux.(im_name{i-NUM_in+1}); 
%--------------------------------------------------------------------------
        bulgex_mov.(im_name{i-NUM_in+1}) = bulgex.(im_name{i-NUM_in+1}) - crestx_dif.(im_name{i-NUM_in+1}); % X' = X-VT
        bulgey_mov.(im_name{i-NUM_in+1}) = bulgey.(im_name{i-NUM_in+1}); 
%         bulgex_mov_low.(im_name{i-NUM_in+1}) = bulgex.(im_name{i-NUM_in+1}) - crestx_dif_low.(im_name{i-NUM_in+1}); % X' = X-VT
%--------------------------------------------------------------------------
    end
end

