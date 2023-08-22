function [localmaxima,localminima,crit_idx_top,crit_idx_bot] = bulge_slope(bulge,bulgex,bulgey,Amp_mm)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates the gradient of the crest at each instance t = T. 
% find local maxima and minima in order to formulate B, i.e. t=0 sec criterion
% ROTATE CRESTS CLOCK OR ANTI-CLOCKWISE and find CRITICAL POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    theta = -90; % rotate minus theta degrees 
    theta_90 = 90;  % rotate theta degrees
%--------------------------------------------------------------------------
    Rotation_90 = [cosd(theta_90) -sind(theta_90); sind(theta_90) cosd(theta_90)]; % Rotation matrix
    Rotation_minus90 = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];  % 2nd Rotation matrix
%--------------------------------------------------------------------------
    crest_rot = Rotation_minus90*bulge(1:end,:)'; % Rotate points minus 90 deg
    crest_rot = crest_rot';
    cresty_rot = crest_rot(:,2);
    crestx_rot = crest_rot(:,1);
%--------------------------------------------------------------------------  
    crest_rot_90 = Rotation_90*bulge'; % Rotate your point plus 90 deg
    crest_rot_90 = crest_rot_90';
    cresty_rot_90 = crest_rot_90(:,2);
    crestx_rot_90 = crest_rot_90(:,1);
%--------------------------------------------------------------------------
% flip both sides to test the stationary points
    localmaxima_90  = islocalmax(cresty_rot_90,'MinProminence',0.2);
    localminima_90  = islocalmin(cresty_rot_90,'MinProminence',0.2);
%     localmaxima = islocalmax(cresty_rot,1,'MinProminence',0.2);
%     localmaxima = isoutlier(cresty_rot,'mean');
%     localminima = islocalmin(cresty_rot,1,'FlatSelection','all');
%      localmaxima  = ischange(cresty_rot,'Threshold',1000);
    localmaxima = islocalmax(cresty_rot,1,'MinProminence',0.1);
    localminima = islocalmin(cresty_rot,1,'MinProminence',0.1);
%--------------------------------------------------------------------------
% find the indices to which the critical points correspond to
    crit_idx_top = find(localmaxima);
    crit_idx_bot = find(localminima);
    
    % if wave has no bulge, then caclulate change of slope manually, by
    % finding an abrupt change in slope on the curve: 
% num_g = 5; % number of elements between gradient 
% if isempty(crit_idx_top) 
%         for k = 1 : numel(bulgey) - num_g % running gradient of bulge
%             bulge_gradient(k,1) = (bulgey(k+num_g,1) - bulgey(k,1))./...
%                 (bulgex(k+num_g,1) - bulgex(k,1));
%         end
%         for k = 1 : numel(bulgey) - num_g
%         if bulge_gradient(k,1)> 1.7  % default 1.7 **important**
%             crit_idx_bot = k; % this index becomes critical point
%             crit_idx_top = crit_idx_bot + 15;
%         end
%         end
% else
%     for k = 1 : numel(bulgey) - num_g % calculate gradient anyway
%             bulge_gradient(k,1) = (bulgey(k+num_g,1) - bulgey(k,1))./...
%                 (bulgex(k+num_g,1) - bulgex(k,1));
%     end
% end

%--------------------------------------------------------------------------
    if isempty(crit_idx_top)
            crit_idx_top = 0;    % replace empty cells with zeros
            crit_idx_bot = 0;

%     else
%             crit_idx_top =  crit_idx_top;
%             crit_idx_bot =  crit_idx_bot;
    end
    if isempty(crit_idx_bot)
            crit_idx_top = 0;    % replace empty cells with zeros
            crit_idx_bot = 0;
    end
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




