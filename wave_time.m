function [ROI,NUM_c,frames,time] = wave_time(wave,phase,wave_type)
% assigns ROI, number of images (NUM_c) and time to each wave case whether
% in freshwater or in saltwater.
dt = 0.002; % timestep

if strcmp(wave_type,'freshwater')
    if strcmp(wave,'GW_175') && strcmp(phase,'peak')
        ROI = {'ROI3','ROI3','ROI3'};
        NUM_c = 58;
        frames_in.GW_175.peak.ROI3 = 31855;  
        frames_end.GW_175.peak.ROI3 = 31913;
        frames.GW_175.peak.ROI3 = frames_in.GW_175.peak.ROI3:frames_end.GW_175.peak.ROI3;
        time.GW_175.peak = round(dt.*frames.GW_175.peak.ROI3(1:end-1),3) - 64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'GW_175') && strcmp(phase,'trough')
        ROI = {'ROI4','ROI4','ROI4'};
        NUM_c = 57;
        frames_in.GW_175.trough.ROI4 = 31423;  
        frames_end.GW_175.trough.ROI4 = 31480;
        frames.GW_175.trough.ROI4 = frames_in.GW_175.trough.ROI4:frames_end.GW_175.trough.ROI4;
        time.GW_175.trough = round(dt.*frames.GW_175.trough.ROI4(1:end-1),3) - 64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'GW_175') && strcmp(phase,'minp2')
        ROI = {'ROI4','ROI4','ROI4'};
        NUM_c = 60;
        frames_in.GW_175.minp2.ROI4 = 31627;  
        frames_end.GW_175.minp2.ROI4 = 31687;
        frames.GW_175.minp2.ROI4 = frames_in.GW_175.minp2.ROI4:frames_end.GW_175.minp2.ROI4;
        time.GW_175.minp2 = round(dt.*frames.GW_175.minp2.ROI4(1:end-1),3) -64;
     %-------------------------------------------------------------------------
    elseif strcmp(wave,'GW_175') && strcmp(phase,'posp2')
        ROI = {'ROI3','ROI3','ROI3'};
        NUM_c = 48;
        frames_in.GW_175.posp2.ROI3 = 32120;  
        frames_end.GW_175.posp2.ROI3 = 32168;
        frames.GW_175.posp2.ROI3 = frames_in.GW_175.posp2.ROI3:frames_end.GW_175.posp2.ROI3;  
        time.GW_175.posp2 = round(dt.*frames.GW_175.posp2.ROI3(1:end-1),3) -64;
    elseif strcmp(wave,'GW_18') && strcmp(phase,'peak')
    %----------------------------------------------------------------------
    %% ------------------------------GW_18---------------------------------
    %%---------------------------------------------------------------------
        ROI = {'ROI4','ROI4','ROI4'};
        NUM_c = 53;
        frames_in.GW_18.peak.ROI4 = 31831;  
        frames_end.GW_18.peak.ROI4 = 31890;
        frames.GW_18.peak.ROI4 = frames_in.GW_18.peak.ROI4:frames_end.GW_18.peak.ROI4;
        time.GW_18.peak = round(dt.*frames.GW_18.peak.ROI4(1:end-1),3) - 64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'GW_18') && strcmp(phase,'trough')
        ROI = {'ROI5','ROI5','ROI5'};
        NUM_c = 55;
        frames_in.GW_18.trough.ROI5 = 31400;  
        frames_end.GW_18.trough.ROI5 = 31460;
        frames.GW_18.trough.ROI5 = frames_in.GW_18.trough.ROI5:frames_end.GW_18.trough.ROI5;
        time.GW_18.trough = round(dt.*frames.GW_18.trough.ROI5(1:end-1),3) -64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'GW_18') && strcmp(phase,'minp2')
        ROI = {'ROI5','ROI5','ROI5'};
        NUM_c = 53;
        frames_in.GW_18.minp2.ROI5 = 31600;  
        frames_end.GW_18.minp2.ROI5 = 31643;
        frames.GW_18.minp2.ROI5 = frames_in.GW_18.minp2.ROI5:frames_end.GW_18.minp2.ROI5;
        time.GW_18.minp2 = round(dt.*frames.GW_18.minp2.ROI5(1:end-1),3) -64;
     %-------------------------------------------------------------------------
     elseif strcmp(wave,'GW_18') && strcmp(phase,'posp2')
        ROI = {'ROI3','ROI3','ROI3'};
         NUM_c = 50;
        frames_in.GW_18.posp2.ROI3 = 32073;  
        frames_end.GW_18.posp2.ROI3 = 32123;
        frames.GW_18.posp2.ROI3 = frames_in.GW_18.posp2.ROI3:frames_end.GW_18.posp2.ROI3;   
        time.GW_18.posp2 = round(dt.*frames.GW_18.posp2.ROI3(1:end-1),3) -64;
    elseif strcmp(wave,'JS_15') && strcmp(phase,'peak')
    %--------------------------------------------------------------------------
    %% ------------------------------JS_15-------------------------------------
    %--------------------------------------------------------------------------
        ROI = {'ROI1','ROI1','ROI1'};
        NUM_c = 45;
        frames_in.JS_15.peak.ROI1 = 31890;  
        frames_end.JS_15.peak.ROI1 = 31935;
        frames.JS_15.peak.ROI1 = frames_in.JS_15.peak.ROI1:frames_end.JS_15.peak.ROI1;
        time.JS_15.peak = round(dt.*frames.JS_15.peak.ROI1(1:end-1),3) - 64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'JS_15') && strcmp(phase,'trough')
        ROI = {'ROI1','ROI1','ROI1'};
        NUM_c = 65;
        frames_in.JS_15.trough.ROI1 = 31594;  
        frames_end.JS_15.trough.ROI1 = 31659;
        frames.JS_15.trough.ROI1 = frames_in.JS_15.trough.ROI1:frames_end.JS_15.trough.ROI1;
        time.JS_15.trough = round(dt.*frames.JS_15.trough.ROI1(1:end-1),3) -64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'JS_15') && strcmp(phase,'minp2')
        ROI = {'ROI2','ROI2','ROI2'};
        NUM_c = 47;
        frames_in.JS_15.minp2.ROI2 = 31730;  
        frames_end.JS_15.minp2.ROI2 = 31777;
        frames.JS_15.minp2.ROI2 = frames_in.JS_15.minp2.ROI2:frames_end.JS_15.minp2.ROI2;
        time.JS_15.minp2 = round(dt.*frames.JS_15.minp2.ROI2(1:end-1),3) -64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'JS_15') && strcmp(phase,'posp2')
        ROI = {'ROI2','ROI2','ROI2'};
        NUM_c = 60;
        frames_in.JS_15.posp2.ROI2 = 32053;  
        frames_end.JS_15.posp2.ROI2 = 32113;
        frames.JS_15.posp2.ROI2 = frames_in.JS_15.posp2.ROI2:frames_end.JS_15.posp2.ROI2;  
        time.JS_15.posp2 = round(dt.*frames.JS_15.posp2.ROI2(1:end-1),3) -64;
    elseif strcmp(wave,'JS_155') && strcmp(phase,'peak')
    %-------------------------------------------------------------------------- 
    %% ------------------------------JS_155------------------------------------
    %--------------------------------------------------------------------------
        ROI = {'ROI2','ROI2','ROI2'};
        NUM_c = 55;
        frames_in.JS_155.peak.ROI2 = 31890;  
        frames_end.JS_155.peak.ROI2 = 31945;
        frames.JS_155.peak.ROI2 = frames_in.JS_155.peak.ROI2:frames_end.JS_155.peak.ROI2;
        time.JS_155.peak = round(dt.*frames.JS_155.peak.ROI2(1:end-1),3) - 64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'JS_155') && strcmp(phase,'trough')
        ROI = {'ROI1','ROI1','ROI1'};
        NUM_c = 39;
        frames_in.JS_155.trough.ROI1 = 31618;  
        frames_end.JS_155.trough.ROI1 = 31657;
        frames.JS_155.trough.ROI1 = frames_in.JS_155.trough.ROI1:frames_end.JS_155.trough.ROI1;
        time.JS_155.trough = round(dt.*frames.JS_155.trough.ROI1(1:end-1),3) -64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'JS_155') && strcmp(phase,'minp2')
        ROI = {'ROI2','ROI2','ROI2'};
        NUM_c = 49;
        frames_in.JS_155.minp2.ROI2 = 31744;  
        frames_end.JS_155.minp2.ROI2 = 31794;
        frames.JS_155.minp2.ROI2 = frames_in.JS_155.minp2.ROI2:frames_end.JS_155.minp2.ROI2;
        time.JS_155.minp2 = round(dt.*frames.JS_155.minp2.ROI2(1:end-1),3) -64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'JS_155') && strcmp(phase,'posp2')
        ROI = {'ROI2','ROI2','ROI2'};
        NUM_c = 50;
        frames_in.JS_155.posp2.ROI2 = 32055;  
        frames_end.JS_155.posp2.ROI2 = 32105;
        frames.JS_155.posp2.ROI2 = frames_in.JS_155.posp2.ROI2:frames_end.JS_155.posp2.ROI2;  
        time.JS_155.posp2 = round(dt.*frames.JS_155.posp2.ROI2(1:end-1),3) -64;
    %--------------------------------------------------------------------------
    end

%---------------------------------------------------------------------------
elseif strcmp(wave_type,'saltwater')
%---------------------------------------------------------------------------

    if strcmp(wave,'GW_175') && strcmp(phase,'peak')
        ROI = {'ROI3','ROI3','ROI3'};
        NUM_c = 62;
        frames_in.GW_175.peak.ROI3 = 31860;  
        frames_end.GW_175.peak.ROI3 = 31925;
        frames.GW_175.peak.ROI3 = frames_in.GW_175.peak.ROI3:frames_end.GW_175.peak.ROI3;
        time.GW_175.peak = round(dt.*frames.GW_175.peak.ROI3(1:end-1),3) - 64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'GW_175') && strcmp(phase,'trough')
        ROI = {'ROI4','ROI4','ROI4'};
        NUM_c = 60;
        frames_in.GW_175.trough.ROI4 = 31427;  
        frames_end.GW_175.trough.ROI4 = 31487;
        frames.GW_175.trough.ROI4 = frames_in.GW_175.trough.ROI4:frames_end.GW_175.trough.ROI4;
        time.GW_175.trough = round(dt.*frames.GW_175.trough.ROI4(1:end-1),3) -64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'GW_175') && strcmp(phase,'minp2')
        ROI = {'ROI4','ROI4','ROI4'};
        NUM_c = 60;
        frames_in.GW_175.minp2.ROI4 = 31632;  
        frames_end.GW_175.minp2.ROI4 = 31692;
        frames.GW_175.minp2.ROI4 = frames_in.GW_175.minp2.ROI4:frames_end.GW_175.minp2.ROI4;
        time.GW_175.minp2 = round(dt.*frames.GW_175.minp2.ROI4(1:end-1),3) -64;
     %-------------------------------------------------------------------------
    elseif strcmp(wave,'GW_175') && strcmp(phase,'posp2')
        ROI = {'ROI2','ROI2','ROI2'};
        NUM_c = 60;
        frames_in.GW_175.posp2.ROI3 = 32110;  
        frames_end.GW_175.posp2.ROI3 = 32170;
        frames.GW_175.posp2.ROI3 = frames_in.GW_175.posp2.ROI3:frames_end.GW_175.posp2.ROI3;  
        time.GW_175.posp2 = round(dt.*frames.GW_175.posp2.ROI3(1:end-1),3) -64;
    elseif strcmp(wave,'GW_18') && strcmp(phase,'peak')
    %--------------------------------------------------------------------------
    %% ------------------------------GW_18-------------------------------------
    % %------------------------------------------------------------------------
        ROI = {'ROI3','ROI3','ROI3'};
        NUM_c = 67;
        frames_in.GW_18.peak.ROI4 = 31820;  
        frames_end.GW_18.peak.ROI4 = 31887;
        frames.GW_18.peak.ROI4 = frames_in.GW_18.peak.ROI4:frames_end.GW_18.peak.ROI4;
        time.GW_18.peak = round(dt.*frames.GW_18.peak.ROI4(1:end-1),3) - 64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'GW_18') && strcmp(phase,'trough')
        ROI = {'ROI3','ROI3','ROI3'};
        NUM_c = 62;
        frames_in.GW_18.trough.ROI5 = 31377;  
        frames_end.GW_18.trough.ROI5 = 31445;
        frames.GW_18.trough.ROI5 = frames_in.GW_18.trough.ROI5:frames_end.GW_18.trough.ROI5;
        time.GW_18.trough = round(dt.*frames.GW_18.trough.ROI5(1:end-1),3) -64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'GW_18') && strcmp(phase,'minp2')
        % ROI = {'ROI3','ROI3','ROI3'};
        ROI = {'ROI4','ROI4','ROI4'};
        NUM_c = 67;   % 67 for ROI3 60 for ROI4
        frames_in.GW_18.minp2.ROI5 = 31575;  
        frames_end.GW_18.minp2.ROI5 = 31642;
        frames.GW_18.minp2.ROI5 = frames_in.GW_18.minp2.ROI5:frames_end.GW_18.minp2.ROI5;
        time.GW_18.minp2 = round(dt.*frames.GW_18.minp2.ROI5(1:end-1),3) -64;
     %-------------------------------------------------------------------------
     elseif strcmp(wave,'GW_18') && strcmp(phase,'posp2')
        ROI = {'ROI2','ROI2','ROI2'};
        NUM_c = 69;
        frames_in.GW_18.posp2.ROI3 = 32067;  
        frames_end.GW_18.posp2.ROI3 = 32136;
        frames.GW_18.posp2.ROI3 = frames_in.GW_18.posp2.ROI3:frames_end.GW_18.posp2.ROI3;   
        time.GW_18.posp2 = round(dt.*frames.GW_18.posp2.ROI3(1:end-1),3) -64;
    elseif strcmp(wave,'JS_15') && strcmp(phase,'peak')
    %--------------------------------------------------------------------------
    %% ------------------------------JS_15-------------------------------------
    %--------------------------------------------------------------------------
        ROI = {'ROI2','ROI2','ROI2'};
        NUM_c = 54;
        frames_in.JS_15.peak.ROI2 = 31895;  
        frames_end.JS_15.peak.ROI2 = 31950;
        frames.JS_15.peak.ROI2 = frames_in.JS_15.peak.ROI2:frames_end.JS_15.peak.ROI2;
        time.JS_15.peak = round(dt.*frames.JS_15.peak.ROI2(1:end-1),3) - 64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'JS_15') && strcmp(phase,'trough')
        ROI = {'ROI1','ROI1','ROI1'};
        NUM_c = 50;
        frames_in.JS_15.trough.ROI1 = 31600;  
        frames_end.JS_15.trough.ROI1 = 31655;
        frames.JS_15.trough.ROI1 = frames_in.JS_15.trough.ROI1:frames_end.JS_15.trough.ROI1;
        time.JS_15.trough = round(dt.*frames.JS_15.trough.ROI1(1:end-1),3) -64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'JS_15') && strcmp(phase,'minp2')
        ROI = {'ROI1','ROI1','ROI1'};
        NUM_c = 54;
        frames_in.JS_15.minp2.ROI1 = 31722;  
        frames_end.JS_15.minp2.ROI1 = 31777;
        frames.JS_15.minp2.ROI1 = frames_in.JS_15.minp2.ROI1:frames_end.JS_15.minp2.ROI1;
        time.JS_15.minp2 = round(dt.*frames.JS_15.minp2.ROI1(1:end-1),3) - 64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'JS_15') && strcmp(phase,'posp2')
        ROI = {'ROI1','ROI1','ROI1'};
        NUM_c = 55;
        frames_in.JS_15.posp2.ROI1 = 32053;  
        frames_end.JS_15.posp2.ROI1 = 32113;
        frames.JS_15.posp2.ROI1 = frames_in.JS_15.posp2.ROI1:frames_end.JS_15.posp2.ROI1;  
        time.JS_15.posp2 = round(dt.*frames.JS_15.posp2.ROI1(1:end-1),3) -64;
    elseif strcmp(wave,'JS_155') && strcmp(phase,'peak')
    %-------------------------------------------------------------------------- 
    %% ------------------------------JS_155------------------------------------
    %--------------------------------------------------------------------------
        ROI = {'ROI2','ROI2','ROI2'};
        NUM_c = 50;
        frames_in.JS_155.peak.ROI2 = 31897;  
        frames_end.JS_155.peak.ROI2 = 31950;
        frames.JS_155.peak.ROI2 = frames_in.JS_155.peak.ROI2:frames_end.JS_155.peak.ROI2;
        time.JS_155.peak = round(dt.*frames.JS_155.peak.ROI2(1:end-1),3) - 64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'JS_155') && strcmp(phase,'trough')
        ROI = {'ROI1','ROI1','ROI1'};
        NUM_c = 50;
        frames_in.JS_155.trough.ROI1 = 31615;  
        frames_end.JS_155.trough.ROI1 = 31665;
        frames.JS_155.trough.ROI1 = frames_in.JS_155.trough.ROI1:frames_end.JS_155.trough.ROI1;
        time.JS_155.trough = round(dt.*frames.JS_155.trough.ROI1(1:end-1),3) -64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'JS_155') && strcmp(phase,'minp2')
        ROI = {'ROI2','ROI2','ROI2'};
        NUM_c = 60;
        frames_in.JS_155.minp2.ROI2 = 31740;  
        frames_end.JS_155.minp2.ROI2 = 31805;
        frames.JS_155.minp2.ROI2 = frames_in.JS_155.minp2.ROI2:frames_end.JS_155.minp2.ROI2;
        time.JS_155.minp2 = round(dt.*frames.JS_155.minp2.ROI2(1:end-1),3) -64;
    %--------------------------------------------------------------------------
    elseif strcmp(wave,'JS_155') && strcmp(phase,'posp2')
        ROI = {'ROI2','ROI2','ROI2'};
        NUM_c =50;
        frames_in.JS_155.posp2.ROI2 = 32065;  
        frames_end.JS_155.posp2.ROI2 = 32115;
        frames.JS_155.posp2.ROI2 = frames_in.JS_155.posp2.ROI2:frames_end.JS_155.posp2.ROI2;  
        time.JS_155.posp2 = round(dt.*frames.JS_155.posp2.ROI2(1:end-1),3) -64;
    %--------------------------------------------------------------------------
    end
end
end