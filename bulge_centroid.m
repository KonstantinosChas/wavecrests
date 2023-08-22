function[xcm,ycm,ycm_test] = bulge_centroid(bulgex_fin,bulgey_fin,NUM,NUM_in,im_name)
    clc; 
    close all;
%     .(im_name{i-NUM_in+1})
% find centroid coordinates 
    for i = NUM_in:NUM
        ycm.(im_name{i-NUM_in+1}) = mean(bulgey_fin .(im_name{i-NUM_in+1}));
        xcm.(im_name{i-NUM_in+1}) = mean(bulgex_fin .(im_name{i-NUM_in+1}));
        ycm_test.(im_name{i-NUM_in+1}) = sum(bulgey_fin.(im_name{i-NUM_in+1}))./numel(bulgey_fin.(im_name{i-NUM_in+1}));
    end

end