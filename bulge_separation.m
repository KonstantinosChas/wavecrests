function[bulgex_fin,bulgey_fin,bulge_idx_bot,bulge_idx_top,bulge_minima,bulge_maxima] = bulge_separation(bulge,crit_idx_bot,crit_idx_top,NUM,NUM_in,im_name,Amp_pseudo_mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% separates bulge from the rest of the crest according to the critical
% points found in bulge slope.m
% run with wave_crests_main.m
% note: it is a bit complicated for the moment, but it is hard to impose general
% conditions for all types of breakers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc; 
    close all;
%--------------------------------------------------------------------------
%-----------------------------PARAMETERS-----------------------------------
    cutoff = 0;  % add a small cutoff threshold index in 
    amp = Amp_pseudo_mm;     % amplitude 
    amp_add = 0;
    arb_height = amp + amp_add; % amp in mm
%--------------------------------------------------------------------------        
    if  any(structfun(@all,crit_idx_bot))== 0 
%--------------------------------------------------------------------------        
        for i = NUM_in:NUM
            T.(im_name{i-NUM_in+1}) = 1; % clause index for last statement
% if critical points have NOT been found
% arbitary cutoff is chosen, because no actual bulge exists
            bulge_idx_gen.(im_name{i-NUM_in+1}) = find(bulge.(im_name{i-NUM_in+1})(:,2) <= arb_height + 0.1 & bulge.(im_name{i-NUM_in+1})(:,2) >= arb_height - 0.1); 
            if isempty(bulge_idx_gen.(im_name{i-NUM_in+1}))
                bulge_idx_gen.(im_name{i-NUM_in+1}) = find(bulge.(im_name{i-NUM_in+1})(:,2) <= arb_height + 0.2 & bulge.(im_name{i-NUM_in+1})(:,2) >= arb_height - 0.2);
            end
            if isempty(bulge_idx_gen.(im_name{i-NUM_in+1}))
                bulge_idx_gen.(im_name{i-NUM_in+1}) = find(bulge.(im_name{i-NUM_in+1})(:,2) <= arb_height + 2.3 & bulge.(im_name{i-NUM_in+1})(:,2) >= arb_height - 2.3);
            end
            % for crests with no bulge early in the recording
%             if bulge.(im_name{i-NUM_in+1})(bulge_idx_gen.(im_name{i-NUM_in+1}),2) <= arb_height
    %ANCIENT REMAINS
%             end
            % for crests with no bulge late in the recording
%             if bulge.(im_name{i-NUM_in+1})(1,2) >= arb_height - 5
                bulgey_fin.(im_name{i-NUM_in+1}) = bulge.(im_name{i-NUM_in+1})(3*end/4:end,2);
                bulgex_fin.(im_name{i-NUM_in+1}) = bulge.(im_name{i-NUM_in+1})(3*end/4:end,1);  
%             else
%                 bulgey_fin.(im_name{i-NUM_in+1}) = bulge.(im_name{i-NUM_in+1})(bulge_idx_gen.(im_name{i-NUM_in+1})(1):end,2);
%                 bulgex_fin.(im_name{i-NUM_in+1}) = bulge.(im_name{i-NUM_in+1})(bulge_idx_gen.(im_name{i-NUM_in+1})(1):end,1);
%             end
            
            
%--------------------------------------------------------------------------        
            bulge_idx_bot.(im_name{i-NUM_in+1}) =  numel(bulge.(im_name{i-NUM_in+1})(:,2)) - max(bulge_idx_gen.(im_name{i-NUM_in+1}));
            bulge_idx_top.(im_name{i-NUM_in+1}) =  numel(bulgey_fin.(im_name{i-NUM_in+1}));
        end
%--------------------------------------------------------------------------        
    elseif any(structfun(@all,crit_idx_bot))== 1
%--------------------------------------------------------------------------        
% if any critical points have been found
        for i = NUM_in:NUM
%             if crit_idx_bot.(im_name{i-NUM_in+1}) ~= 0 
            if any(crit_idx_bot.(im_name{i-NUM_in+1})) == 1
% find the index of the crest to which the first critical point belongs to
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                nonzero_idx = i;
                disp('found one!');
                pause(0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                break
            elseif any(crit_idx_bot.(im_name{i-NUM_in+1})) == 0
                disp('still searching for critical point...');
                pause(0.5);
            end
        end
%--------------------------------------------------------------------------        
        for i = NUM_in:NUM
%--------------------------------------------------------------------------
%-------------------------------------------------------------------------- 
            if crit_idx_bot.(im_name{i-NUM_in+1}) ~= 0 
%--------------------------------------------------------------------------
%-------------------------------------------------------------------------- 
% for crests with critical points calculate bulge according to those
                T.(im_name{i-NUM_in+1}) = 1; % clause for last statement
                bulgey_fin.(im_name{i-NUM_in+1}) = bulge.(im_name{i-NUM_in+1})(max(crit_idx_bot.(im_name{i-NUM_in+1}))-cutoff:end,2);
                bulgex_fin.(im_name{i-NUM_in+1}) = bulge.(im_name{i-NUM_in+1})(max(crit_idx_bot.(im_name{i-NUM_in+1}))-cutoff:end,1);
%--------------------------------------------------------------------------  
% find the critical indices
                bulge_idx_bot.(im_name{i-NUM_in+1}) = find(bulgey_fin.(im_name{i-NUM_in+1})  == bulge.(im_name{i-NUM_in+1})(max(crit_idx_bot.(im_name{i-NUM_in+1})),2));
                bulge_idx_top.(im_name{i-NUM_in+1}) = find(bulgey_fin.(im_name{i-NUM_in+1})  == bulge.(im_name{i-NUM_in+1})(max(crit_idx_top.(im_name{i-NUM_in+1})),2));
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
            elseif crit_idx_bot.(im_name{i-NUM_in+1}) == 0
%--------------------------------------------------------------------------
%-------------------------------------------------------------------------- 
                T.(im_name{i-NUM_in+1}) = 0; % clause for last statement 
%--------------------------------------------------------------------------
                if numel(bulge.(im_name{i-NUM_in+1})(:,2)) >= numel(bulge.(im_name{nonzero_idx-NUM_in+1})(:,2))
% for fields larger than the (first) critical field - find the cutoff index where those 2 fields are equal
                    bulge_zero_idx.(im_name{i-NUM_in+1}) = find(bulge.(im_name{i-NUM_in+1})(:,2) <= bulge.(im_name{nonzero_idx-NUM_in+1})(crit_idx_bot.(im_name{nonzero_idx-NUM_in+1})(1),2) + 1 ...
                        & bulge.(im_name{i-NUM_in+1})(:,2) >= bulge.(im_name{nonzero_idx-NUM_in+1})(crit_idx_bot.(im_name{nonzero_idx-NUM_in+1})(1),2) - 1);
%-------------------------------------------------------------------------- 
                        disp('Hello');
                else
% for fields smaller than the critical field - the cutoff index is arbitary
% (cannot be anything else)
                    bulge_zero_idx.(im_name{i-NUM_in+1}) = find(bulge.(im_name{i-NUM_in+1})(:,2) <= arb_height + 3 & bulge.(im_name{i-NUM_in+1})(:,2) >= arb_height - 3);
%                         bulge_zero_idx.(im_name{i-NUM_in+1}) = numel(bulge.(im_name{i-NUM_in+1})(:,2));
                        disp('or here')
                end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% these are the final bulges
                bulgey_fin.(im_name{i-NUM_in+1}) = bulge.(im_name{i-NUM_in+1})(max(bulge_zero_idx.(im_name{i-NUM_in+1}))-cutoff:end,2);
                bulgex_fin.(im_name{i-NUM_in+1}) = bulge.(im_name{i-NUM_in+1})(max(bulge_zero_idx.(im_name{i-NUM_in+1}))-cutoff:end,1);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% check number of elements of each field and then reshape according to bulge idx_top 
                    if max(crit_idx_top.(im_name{nonzero_idx-NUM_in+1}))<= numel(bulgey_fin.(im_name{i-NUM_in+1}))
                        bulge_idx_bot.(im_name{i-NUM_in+1}) = find(bulgey_fin.(im_name{i-NUM_in+1}) == bulge.(im_name{i-NUM_in+1})(max(bulge_zero_idx.(im_name{i-NUM_in+1})),2));
                        bulge_idx_top.(im_name{i-NUM_in+1}) = find(bulgey_fin.(im_name{i-NUM_in+1}) == bulge.(im_name{i-NUM_in+1})(max(crit_idx_top.(im_name{nonzero_idx-NUM_in+1})),2));
                    else
                        bulge_idx_bot.(im_name{i-NUM_in+1}) = find(bulgey_fin.(im_name{i-NUM_in+1}) == bulge.(im_name{i-NUM_in+1})(max(bulge_zero_idx.(im_name{i-NUM_in+1})),2));
                        bulge_idx_top.(im_name{i-NUM_in+1}) = numel(bulgey_fin.(im_name{i-NUM_in+1}));
                    end
            end
        end
%--------------------------------------------------------------------------        
    end
% generate new localminima/localmaxima arrays that correspond to new bulge matrix indices 
    for i = NUM_in:NUM
        if T.(im_name{i-NUM_in+1}) == 1 
            bulge_minima.(im_name{i-NUM_in+1}) = zeros(numel(bulgey_fin.(im_name{i-NUM_in+1})),1);
            bulge_maxima.(im_name{i-NUM_in+1}) = zeros(numel(bulgey_fin.(im_name{i-NUM_in+1})),1);
%             bulge_minima.(im_name{i-NUM_in+1})(max(bulge_idx_bot.(im_name{i-NUM_in+1}))) = 1;
            bulge_minima.(im_name{i-NUM_in+1})(1) = 1;
            bulge_maxima.(im_name{i-NUM_in+1})(max(bulge_idx_top.(im_name{i-NUM_in+1}))) = 1;
            bulge_zero_idx = 0;
        elseif T.(im_name{i-NUM_in+1}) == 0
            bulge_minima.(im_name{i-NUM_in+1}) = zeros(numel(bulgey_fin.(im_name{i-NUM_in+1})),1);
            bulge_maxima.(im_name{i-NUM_in+1}) = zeros(numel(bulgey_fin.(im_name{i-NUM_in+1})),1);
%             bulge_minima.(im_name{i-NUM_in+1})(bulge_idx_bot.(im_name{nonzero_idx-NUM_in+1})) = 1;
            bulge_minima.(im_name{i-NUM_in+1})(1) = 1;
            bulge_maxima.(im_name{i-NUM_in+1})(bulge_idx_top.(im_name{nonzero_idx-NUM_in+1})) = 1;
        end
    end
% convert to logical
    for i = NUM_in:NUM
        bulge_minima.(im_name{i-NUM_in+1}) = logical(bulge_minima.(im_name{i-NUM_in+1}));
        bulge_maxima.(im_name{i-NUM_in+1}) = logical(bulge_maxima.(im_name{i-NUM_in+1}));
    end
%--------------------------------------------------------------------------        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%