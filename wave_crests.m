
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------Wave Crests Main Program---------------------------
% what is does: 
% 1. Load Calibrated images from wave_calibration.m for each run
% 2. Wave Crest/Bulge Extraction
% 3. Calculates Kinematics and Area of deformation
% 4. More Detailed Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%% Konstantinos Chasapis 2019-21
clear;
clc; 
close all;
beep off;
clear dir
clear run
clear NUM_in NUM
%--------------------------------------------------------------------------
%% Part 0 : parameter allocation
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
wave_type = {'saltwater'};
wave = {'JS_15'};
gain = 1.5; 
% phase = {'peak'};
phase = {'trough'};
% phase = {'minp2'};
% phase = {'posp2'};
run = {'run3'};
%--------------------------------------------------------------------------
T = 1;  % : unified coordinate system 1:YES/0:NO  
P = 1;  % : test plot 1:YES/0:NO
%--------------------------------------------------------------------------
%% --------------- ----- time,ROI and number of images---------------------
[ROI,NUM_c,frames,time] = wave_time(wave{1},phase{1},wave_type{1});
%--------------------------------------------------------------------------
NUM_case =  1;
%--------------------------------------------------------------------------
NUM_in = 1;
NUM = NUM_c;                  % number of images; 
%% ---------------------------parameters-----------------------------------
dir_target = '60mm';       % directory of calibration targets
lens = '60mm';
imtype = 'jpg';
num_pix_x = 1280;
num_pix_y = 800;
gain_pseudo = 1.2;
Amp_initial = 40;
Amp_mm = Amp_initial.*gain;   % real linear amp in mm
Amp_pseudo_mm = Amp_initial.*gain_pseudo; % pseudo linear amplitude in mm for cut off purposes..
eta_cor = 10;                 % correction for eta(x) in mm (not used)
lambda_mm = 1450;             % central (f_p) wavelength in mm
phase_speed = 1460;           % linear phase speed U_p or U_c.

dt = 1/500; % interframe time
%--------------------------------------------------------------------------
for case_idx = 1:NUM_case
%% ------------------------------Paths-------------------------------------
    clear crest bulge bulgex bulgey 
    clear bulgex_fin bulgey_fin
    clear buglex_mov bulgey_mov
    clear xcm ycm velo_xcm velo_ycm
    clear crest_vel crest_avg_vel bulge_avg_vel bulgex_mov
    clear bulge_idx_top bulge_idx_bot crit_idx_bot crit_idx_top 
    clear criterion criterion_bulge
%--------------------------------------------------------------------------    
    % C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\wave crests\images\saltwater\60mm\GW_18\ROI2\run1
    dir = ['C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\wave crests\images\',wave_type{case_idx},'\',lens,'\',wave{case_idx},'\',ROI{case_idx},'\',run{case_idx}];  % directory of images to analyze
    dir_save = ['C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\wave crests\data\',wave_type{case_idx},'\',wave{case_idx},'\',ROI{case_idx},'\',run{case_idx}];
%     dir_plot_save = ['C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\wave crests\plots\freshwater\',wave,'\',ROI{case_idx},'\',run{case_idx}];
dir_plot_save = ['C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\wave crests\plots\',wave_type{case_idx}];
%-------------------------------------------------------------------------
    addpath(dir);
    addpath('C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\wave crests\matlab');
    cd('C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\wave crests\matlab');
    cd(dir);
%% ----------------------------Load Images---------------------------------
    im_name = cell(1,NUM-NUM_in+1);
    name =[wave{case_idx},'_',phase{case_idx}];   % names of images (e.g. ...peak_0001)
    name_0 = [wave{case_idx},'_',phase{case_idx},'_000'];
    name_10 = [wave{case_idx},'_',phase{case_idx},'_00'];
    name_100 = [wave{case_idx},'_',phase{case_idx},'_0'];
    name_1000 = [wave{case_idx},'_',phase{case_idx},'_'];
    for i = NUM_in:NUM
            if i>=1000
                im_name{i-NUM_in+1} = ([name_1000,num2str(i)]);
            elseif i>=100 && i<1000
                im_name{i-NUM_in+1} = ([name_100,num2str(i)]);
            elseif i<100 && i>=10
                im_name{i-NUM_in+1} = ([name_10,num2str(i)]);
            elseif i<10
                im_name{i-NUM_in+1} = ([name_0,num2str(i)]);
            end
    end
    for i = NUM_in : NUM
        folderName        =   imread([im_name{i-NUM_in+1},'_cal.',imtype]);
        images.(im_name{i-NUM_in+1}) = folderName;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------------Wave Extraction--------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = NUM_in:NUM
        [crest.(im_name{i-NUM_in+1}), bulge.(im_name{i-NUM_in+1}),bulge_2.(im_name{i-NUM_in+1}),max_idx.(im_name{i-NUM_in+1}),...
            max_bulge.(im_name{i-NUM_in+1}), imbin_g.(im_name{i-NUM_in+1}),imbin_g_aux.(im_name{i-NUM_in+1}),Gmag.(im_name{i-NUM_in+1})] = wave_extraction(images.(im_name{i-NUM_in+1}));  
    end
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------------Coordinate System------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if T ==1
        [bulgex,bulgey,bulge,crestx,cresty,crest,x_cor] = coordinates2D(bulge,crest,ROI{case_idx},NUM_in,NUM,im_name);
    elseif T==0
        for i = NUM_in:NUM
            cresty.(im_name{i-NUM_in+1}) = crest.(im_name{i-NUM_in+1})(:,2);
            crestx.(im_name{i-NUM_in+1}) = crest.(im_name{i-NUM_in+1})(:,1);
            bulgey.(im_name{i-NUM_in+1}) = bulge.(im_name{i-NUM_in+1})(:,2);
            bulgex.(im_name{i-NUM_in+1}) = bulge.(im_name{i-NUM_in+1})(:,1);
            bulge.(im_name{i-NUM_in+1})(:,2) = bulge.(im_name{i-NUM_in+1})(:,2);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------------Bulge Slope Calculation I-------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = NUM_in:NUM
        [localmaxima.(im_name{i-NUM_in+1}),localminima.(im_name{i-NUM_in+1}),...
            crit_idx_top.(im_name{i-NUM_in+1}),crit_idx_bot.(im_name{i-NUM_in+1})] = bulge_slope(bulge.(im_name{i-NUM_in+1}),bulgex.(im_name{i-NUM_in+1}),bulgey.(im_name{i-NUM_in+1}),Amp_mm);
    end   
%     ,bulge_gradient.(im_name{i-NUM_in+1})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------------Bulge Separation-------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [bulgex_fin,bulgey_fin,bulge_idx_bot,bulge_idx_top,bulge_minima,bulge_maxima] = ...
        bulge_separation(bulge,crit_idx_bot,crit_idx_top,NUM,NUM_in,im_name,Amp_pseudo_mm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------Relative Velocities Calculation---------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [crest_avg_vel,bulge_avg_vel,crest_vel,crest_vel_idx_top,vel_idx,ixvel,ixvel_1] = rel_velocity(bulge,crit_idx_bot,NUM,NUM_in,im_name);
%--------------------------------------------------------------------------    
    [bulgex_mov,bulgey_mov] = moving_frame(crest_avg_vel,bulgex,bulgey,NUM,NUM_in,im_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------Centroid Velocities Calculation-----------------------
% centroid coordinates
    [xcm,ycm,ycm_test] = bulge_centroid(bulgex_fin,bulgey_fin,NUM,NUM_in,im_name);
%-------------------------------------------------------------------------- 
% centroid velocities
%     [velo_ycm_dt,velo_xcm_dt,velo_xcm,velo_ycm,velo_crest_dt,velo_crest_top_dt,velo_crest_top_dt_smooth,velo_crest_dt_smooth,velo_crest_tip_dt,velo_crest_tip_dt_smooth] = velo_centroid(bulge,max_idx,crit_idx_top,xcm,ycm,NUM,NUM_in,im_name);
    
%% --------------------Transformed Coordinate System-----------------------
%     [bulge_m_1,bulge_m,bulgey_poly,bulgey_fit,bulgey_fitlimit,bulgex_fitlimit] = crest_polyfit(crit_idx_bot,bulgex,bulgey,im_name,NUM_in,NUM);
% 
%     for i = NUM_in:NUM
%         bulge_angle_ratio.(im_name{i-NUM_in+1}) = bulge_theta.(im_name{i-NUM_in+1})./bulge_m.(im_name{i-NUM_in+1});
%     end
%--------------------------------------------------------------------------
%% ----------------------------Curve fitting-------------------------------
    for i = NUM_in:NUM
        [x_pp.(im_name{i-NUM_in+1})] = crest_polyfit(bulgex.(im_name{i-NUM_in+1}),bulgey.(im_name{i-NUM_in+1}));
    end
%% ---------------------------Inflection Point-----------------------------
%     for i = NUM_in: 26
%         [x_pp.(im_name{i-NUM_in+1}),ifpy.(im_name{i-NUM_in+1}),ifpx.(im_name{i-NUM_in+1}),...
%             allrts.(im_name{i-NUM_in+1}),ifpx_all.(im_name{i-NUM_in+1})] = crvfnob(bulgex.(im_name{i-NUM_in+1}),bulgey.(im_name{i-NUM_in+1}));
%     end
    for i = NUM_in:NUM
        [ifpy.(im_name{i-NUM_in+1}),ifpx.(im_name{i-NUM_in+1}),allrts.(im_name{i-NUM_in+1}),der2fy.(im_name{i-NUM_in+1})] = crvf(bulgex.(im_name{i-NUM_in+1}),bulgey.(im_name{i-NUM_in+1}),wave);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -----------------------Bulge Geometrical Analysis-----------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -----------------------Bulge Volume Calculation-------------------------
    [y_hypoten,hor_points,vert_points,num_hor_grid,bulge_area,bulge_area_points] = bulge_volume(bulge,crit_idx_bot,ifpx,NUM,NUM_in,im_name);
    % [hor_points,vert_points,num_hor_grid] = bulge_volume(bulge,crit_idx_bot,ifpx,NUM,NUM_in,im_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -----------criteria of tip toe,ifp: Dz (h), Dx(l), slope -------------
% the easy way : this is in a function as well. - but is simplified here
    for i = NUM_in:NUM
        if crit_idx_top.(im_name{i-NUM_in+1})(end) > 0
            tip_l.(im_name{i-NUM_in+1}) = bulgex.(im_name{i-NUM_in+1})(end) - bulgex.(im_name{i-NUM_in+1})(crit_idx_top.(im_name{i-NUM_in+1})(end)); 
            tip_l.(im_name{i-NUM_in+1}) = round(tip_l.(im_name{i-NUM_in+1}),2);
            tip_h.(im_name{i-NUM_in+1}) = bulgey.(im_name{i-NUM_in+1})(end) - bulgey.(im_name{i-NUM_in+1})(crit_idx_top.(im_name{i-NUM_in+1})(end));
            slope_tip.(im_name{i-NUM_in+1}) =  tip_h.(im_name{i-NUM_in+1})./tip_l.(im_name{i-NUM_in+1});
            tipel.(im_name{i-NUM_in+1}) = bulgey.(im_name{i-NUM_in+1})(crit_idx_top.(im_name{i-NUM_in+1})(end));
            tiple.(im_name{i-NUM_in+1}) = bulgex.(im_name{i-NUM_in+1})(crit_idx_top.(im_name{i-NUM_in+1})(end));

            toe_l.(im_name{i-NUM_in+1}) = bulgex.(im_name{i-NUM_in+1})(end) - bulgex.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(end)); 
            toe_l.(im_name{i-NUM_in+1}) = round(toe_l.(im_name{i-NUM_in+1}),2);
            toe_h.(im_name{i-NUM_in+1}) = bulgey.(im_name{i-NUM_in+1})(end) - bulgey.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(end));
            slope_toe.(im_name{i-NUM_in+1}) =  toe_h.(im_name{i-NUM_in+1})./toe_l.(im_name{i-NUM_in+1});
            
            toeel.(im_name{i-NUM_in+1}) = bulgey.(im_name{i-NUM_in+1})(crit_idx_top.(im_name{i-NUM_in+1})(end));
            toele.(im_name{i-NUM_in+1}) = bulgex.(im_name{i-NUM_in+1})(crit_idx_top.(im_name{i-NUM_in+1})(end));
        else
            tip_l.(im_name{i-NUM_in+1}) = 0;
            tip_h.(im_name{i-NUM_in+1}) = 0;
            tipel.(im_name{i-NUM_in+1}) = 0;
            tiple.(im_name{i-NUM_in+1}) = 0;
            slope_tip.(im_name{i-NUM_in+1}) = 0 ;
            slope_toe.(im_name{i-NUM_in+1}) = 0 ;

            toe_l.(im_name{i-NUM_in+1}) = 0;
            toe_h.(im_name{i-NUM_in+1}) = 0;
            toeel.(im_name{i-NUM_in+1}) = 0;
            toele.(im_name{i-NUM_in+1}) = 0;
        end
    end
    for i = NUM_in:NUM
% for inflection points and crests with no bulges
        ifp_h.(im_name{i-NUM_in+1}) = bulgey.(im_name{i-NUM_in+1})(end) - ifpy.(im_name{i-NUM_in+1});
        ifp_l.(im_name{i-NUM_in+1}) = bulgex.(im_name{i-NUM_in+1})(end) - ifpx.(im_name{i-NUM_in+1});
        slope_ifp.(im_name{i-NUM_in+1}) = ifp_h.(im_name{i-NUM_in+1})./ifp_l.(im_name{i-NUM_in+1});
    end  
% find mean gamma (slope)
    for i = NUM_in:NUM
        slopetipaux(i-NUM_in+1) = slope_tip.(im_name{i-NUM_in+1});
    end
    slopetipm = mean(slopetipaux(slopetipaux>0));
    slope_tipmm = zeros(1,NUM);slope_tipmm(1,:) = slopetipm;
% find min gamma (slope) and which index it corresponds to
    for i = NUM_in:NUM
        crit_idx_top_aux(i-NUM_in+1) = crit_idx_top.(im_name{i-NUM_in+1})(end);
    end
    if any(crit_idx_top_aux(:)) == 1
        [slope_tipmin,slopetip_minidx] = min(slopetipaux(slopetipaux>0));
    % gamma min corresponds to breaking location - find it
            breaking_loc_run = x_pp.(im_name{1,slopetip_minidx})(end);
    end
    % for crests without bulges
    for i = NUM_in:NUM
        slopeifpaux(i-NUM_in+1) = slope_ifp.(im_name{i-NUM_in+1});
        bulge_areaaux(i-NUM_in+1) = bulge_area.(im_name{i-NUM_in+1});
    end
% find min gamma_ifp (slope) and which index it corresponds to
    [slope_ifpmin,slope_ifpminidx] = min(slopeifpaux(slopeifpaux>0));
    [break_vol_ifp,break_vol_ifpidx] = max(bulge_areaaux(bulge_areaaux>0));
% gamma min corresponds to breaking location and time: 
%         breaking_loc_run_ifp = x_pp.(im_name{1,slope_ifpminidx})(end);
        breaking_loc_run_purespiller = x_pp.(im_name{1,break_vol_ifpidx})(end);
        time_loc_run = time.(wave{case_idx}).(phase{case_idx})(break_vol_ifpidx) + 64; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------Theta Caclulation----------------------------
%     [bulge_theta,bulge_slope,dx,dy] = theta(hor_points,vert_points,NUM,NUM_in,im_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------------------Arrays---------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    velo_xcm_mtx = zeros(1,NUM);
    velo_ycm_mtx = zeros(1,NUM);
    ycm_mtx = zeros(1,NUM);
    xcm_mtx = zeros(1,NUM);
    velo_ycm_dt_mtx = zeros(1,NUM-1);
    velo_xcm_dt_mtx = zeros(1,NUM-1);
    velo_crest_top_dt_smooth_mtx = zeros(1,NUM-1);
    velo_crest_dt_smooth_mtx = zeros(1,NUM-1);
    velo_crest_tip_dt_smooth_mtx = zeros(1,NUM-1);
    B = zeros(1,NUM-1);    

%% -------------------------------Some plots-------------------------------
if P == 1
%% -------------------------------Some plots-------------------------------
    cc = jet(NUM-NUM_in+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % plot bulges for inspection:
%         sizefont = 18;
%         sizefont1 = 18;
%         left = 0.11;
%         width = 0.75;
%         height = 0.75;
%         bottom = 0.19;
%         fig_title = ['stpts_',wave{case_idx},'_',phase{case_idx},'_',ROI{1},'_',run{1}];  
%         a1 = figure('Name',fig_title,'NumberTitle','off'); 
%         axes('Position',[left bottom width height]); hold on;
%         set(a1,'units','normalized','outerposition',[0 0 0.6 0.6]);
%         set(a1,'PaperSize',[8 10]);
%         set(a1,'Color',[1 1 1]); 
% %         title(['wave crests and bulges: ', phase{case_idx}]);
%         for i = NUM_in:NUM
%             windowSize = 35; 
%             b = (1/windowSize)*ones(1,windowSize);
%             a = 1;
% 
% %             bulgeys.(im_name{1,i-NUM_in+1}) = filter(b,a,bulgey.(im_name{1,i-NUM_in+1}));
% 
%             bulgeys.(im_name{1,i-NUM_in+1}) = smooth(bulgey.(im_name{1,i-NUM_in+1}),'moving',55);
% %             bulgexs.(im_name{1,i-NUM_in+1}) = smooth(bulgex.(im_name{1,i-NUM_in+1}),'moving',55);
%             plot(-bulgex.(im_name{1,i-NUM_in+1})./lambda_mm,bulgey.(im_name{1,i-NUM_in+1})./Amp_mm,'linewidth',1,'color','black');
% %             plot(-bulgex.(im_name{4})./lambda_mm,bulgey.(im_name{1,4})./Amp_mm,'linewidth',2,'color','black');
% %             plot(-bulgex.(im_name{34})./lambda_mm,bulgey.(im_name{1,34})./Amp_mm,'linewidth',2,'color','black');
% %             plot(-bulgex.(im_name{45})./lambda_mm,bulgey.(im_name{1,45})./Amp_mm,'linewidth',2,'color','black');
%             hold on;
% %             plot(-bulgex.(im_name{1,i-NUM_in+1})./lambda_mm,bulgeys.(im_name{1,i-NUM_in+1})./Amp_mm,'linewidth',1,'color','blue');
% %             plot(-x_pp.(im_name{4})./lambda_mm,bulgey.(im_name{1,4})./Amp_mm,'linewidth',1,'color','red');
% %             plot(-x_pp.(im_name{34})./lambda_mm,bulgey.(im_name{1,34})./Amp_mm,'linewidth',1,'color','red');
% %             plot(-x_pp.(im_name{45})./lambda_mm,bulgey.(im_name{1,45})./Amp_mm,'linewidth',1,'color','red');
% 
% %             plot(-x_pp.(im_name{i-NUM_in+1})./lambda_mm,bulgey.(im_name{1,i-NUM_in+1})./Amp_mm,'linewidth',1,'color','black');
%             hold on;
%             %             plot(-bulgex.(im_name{i-NUM_in+1})(localmaxima.(im_name{i-NUM_in+1}))./lambda_mm,bulgey.(im_name{i-NUM_in+1})(localmaxima.(im_name{i-NUM_in+1}))./Amp_mm,'r*'); 
% %             plot(-bulgex.(im_name{i-NUM_in+1})(localminima.(im_name{i-NUM_in+1}))./lambda_mm,bulgey.(im_name{i-NUM_in+1})(localminima.(im_name{i-NUM_in+1}))./Amp_mm,'b*');
%             plot(-x_pp.(im_name{i-NUM_in+1})(end)./lambda_mm,bulgey.(im_name{i-NUM_in+1})(end)./Amp_mm,'*','color','black'); 
% 
%             if crit_idx_bot.(im_name{i-NUM_in+1})>0
%             plot(-x_pp.(im_name{i-NUM_in+1})(crit_idx_top.(im_name{i-NUM_in+1})(1))./lambda_mm,bulgey.(im_name{i-NUM_in+1})(crit_idx_top.(im_name{i-NUM_in+1})(1))./Amp_mm,'r*');
%             plot(-x_pp.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(1))./lambda_mm,bulgey.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(1))./Amp_mm,'b*');
%             end
%         end
%         for i = NUM_in:NUM
%             for j = 1:numel(allrts.(im_name{i-NUM_in+1}))
%                 if bulgey.(im_name{i-NUM_in+1})(localmaxima.(im_name{i-NUM_in+1}))> 0 
% % 
% %                     if allrts.(im_name{i-NUM_in+1})(j) <= bulgey.(im_name{i-NUM_in+1})(localmaxima.(im_name{i-NUM_in+1})) && allrts.(im_name{i-NUM_in+1})(j)>= bulgey.(im_name{i-NUM_in+1})(localminima.(im_name{i-NUM_in+1}))
% %         %             scatter(-ifpx.(im_name{i-NUM_in+1})./lambda_mm,ifpy.(im_name{i-NUM_in+1})./Amp_mm,100,'green','+','LineWidth',2);
% %                     scatter(-ifpx_all.(im_name{i-NUM_in+1})(j)./lambda_mm,allrts.(im_name{i-NUM_in+1})(j)./Amp_mm,100,'green','+','LineWidth',2);
% %                     end
%                 else
%                     scatter(-ifpx.(im_name{i-NUM_in+1})./lambda_mm,ifpy.(im_name{i-NUM_in+1})./Amp_mm,100,'green','*','LineWidth',1);
% %                     scatter(-ifpx_all.(im_name{i-NUM_in+1})./lambda_mm,allrts.(im_name{i-NUM_in+1})./Amp_mm,100,'green','+','LineWidth',2);
%                 end
%             end
%         end
% %         for i = 58
% %             scatter(-bulge_area_xpoints.(im_name{i-NUM_in+1})./lambda_mm,bulge_area_ypoints.(im_name{i-NUM_in+1})./Amp_mm,'LineWidth',.5,'MarkerEdgeColor','green');
% %             plot(-bulgex_mov.(im_name{1,i-NUM_in+1})./lambda_mm,bulgey_mov.(im_name{1,i-NUM_in+1})./Amp_mm,'linewidth',3,'color','blue');
% %             plot(-bulgex_sp.(im_name{1,i-NUM_in+1})./lambda_mm,bulgey_sp.(im_name{1,i-NUM_in+1})./Amp_mm,'linewidth',3,'color','magenta');
% %         end
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Xdir','reverse');
%             set(gca,'Color',[1 1 1]);
%             xlabel('$\mathrm{x''/\lambda_c}$','Fontsize',sizefont,'interpreter','latex');
%             ylabel('$\eta/A_{linear}$','Fontsize', sizefont,'interpreter','latex');
% %             xlabel('$\mathrm{x''(mm)}$','Fontsize',sizefont,'interpreter','latex');
% %             ylabel('$\mathrm{\eta (mm)}$','Fontsize', sizefont,'interpreter','latex');
% %             ylim([0.8 1.4]);
% %             xlim([-.185 -0.04]);
% %             xlim([-.15 -0.02]);
%             grid on;
%             box on;
% %--------------------------------------------------------------------------
%         cd(dir_plot_save);
%         % export_fig([dir_plot_save,'\',fig_title],'-png','-nocrop');
% 
% %%-------------------------------------------------------------------------
% % %% time evolution of ratios 
% % %     
%         sizefont = 18;
%         sizefont1 = 18;
%         bottom = 0.33;
%         left = 0.14; % 0.14
%         width = 0.8; % 0.8
%         height = 0.62;
%         fig_title = ['parameters',wave{case_idx},'_',phase{case_idx},'_',ROI{1},'_',run{1}];  
%         a1 = figure('Name',[fig_title],'NumberTitle','off'); 
%         axes('Position',[left bottom width height]); hold on;
%         set(a1,'units','normalized','outerposition',[0 0 0.85 0.45]);
%         set(a1,'PaperSize',[8 10]);
%         set(a1,'Color',[1 1 1]); 
%         for i = NUM_in:NUM
%             if crit_idx_bot.(im_name{i-NUM_in+1}) > 0
% %                 scatter(time(i-NUM_in+1),tip_l.(im_name{i-NUM_in+1}),90,'^','MarkerEdgeColor','red','LineWidth',1);
% %                 scatter(time(i-NUM_in+1),tip_h.(im_name{i-NUM_in+1}),40,'*','MarkerEdgeColor','red','LineWidth',1);
% %                 
% %%this one
%                 scatter(time.(wave{case_idx}).(phase{case_idx})(i-NUM_in+1),bulgey.(im_name{i-NUM_in+1})(crit_idx_top.(im_name{i-NUM_in+1})(1))./Amp_mm,40,'*','MarkerEdgeColor','red','LineWidth',1);
% 
% %                 scatter(time(i-NUM_in+1),slope_tip.(im_name{i-NUM_in+1}),40,'+','MarkerEdgeColor','red','LineWidth',1);
% 
% 
% %                 scatter(time(i-NUM_in+1),toe_l.(im_name{i-NUM_in+1}),50,'^','MarkerEdgeColor','blue','LineWidth',1);
% %                 scatter(time(i-NUM_in+1),toe_h.(im_name{i-NUM_in+1}),20,'*','MarkerEdgeColor','blue','LineWidth',1);
%                 scatter(time.(wave{case_idx}).(phase{case_idx})(i-NUM_in+1),bulgey.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(1))./Amp_mm,40,'*','MarkerEdgeColor','blue');
% 
% %                 scatter(time(i-NUM_in+1),slope_toe.(im_name{i-NUM_in+1}),40,'+','MarkerEdgeColor','blue','LineWidth',1);
% 
%             end
%         end
%             hold on;
%             for i = NUM_in:NUM 
% %%this one
%                 scatter(time.(wave{case_idx}).(phase{case_idx})(i-NUM_in+1),bulgey.(im_name{i-NUM_in+1})(end)./Amp_mm,40,'*','MarkerEdgeColor','black','LineWidth',1);
%                 for j = 1:numel(allrts.(im_name{i-NUM_in+1}))
%                     if bulgey.(im_name{i-NUM_in+1})(localmaxima.(im_name{i-NUM_in+1}))> 0 
% % %                         if allrts.(im_name{i-NUM_in+1})(j) <= bulgey.(im_name{i-NUM_in+1})(localmaxima.(im_name{i-NUM_in+1})) && allrts.(im_name{i-NUM_in+1})(j)>= bulgey.(im_name{i-NUM_in+1})(localminima.(im_name{i-NUM_in+1}))
% % %                             scatter(time(i-NUM_in+1),ifp_l.(im_name{i-NUM_in+1}),40,'^','green','LineWidth',1);
% % %                             scatter(time(i-NUM_in+1),ifp_h.(im_name{i-NUM_in+1}),40,'*','green','LineWidth',1);
% % %                             scatter(time(i-NUM_in+1),ifpy.(im_name{i-NUM_in+1}),40,'d','green','LineWidth',1);                   
% % %                         else
% % %                         end
%                     else
% %                         scatter(time(i-NUM_in+1),ifp_l.(im_name{i-NUM_in+1}),40,'^','green','LineWidth',1);
% %                         scatter(time(i-NUM_in+1),ifp_h.(im_name{i-NUM_in+1}),40,'*','green','LineWidth',1);
% %%this one
%                         scatter(time.(wave{case_idx}).(phase{case_idx})(i-NUM_in+1),ifpy.(im_name{i-NUM_in+1})./Amp_mm,40,'*','green','LineWidth',1); 
% %                         
% % %                         scatter(time(i-NUM_in+1),slope_ifp.(im_name{i-NUM_in+1}),40,'+','green','LineWidth',1);
% % 
%                     end
%                 end
%             end
% %             scatter(time(i-NUM_in+1),ifp_l.(im_name{i-NUM_in+1}),40,'^','green','LineWidth',1);
% %             scatter(time(i-NUM_in+1),ifp_h.(im_name{i-NUM_in+1}),40,'*','green','LineWidth',1);
% %             scatter(time(i-NUM_in+1),ifpy.(im_name{i-NUM_in+1}),40,'d','green','LineWidth',1);
% 
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Fontsize',sizefont1);
%             set (gca,'Xdir','reverse')
%             set(gca,'Color',[1 1 1]);
% %             ylabel('\Gamma','Fontsize',sizefont);
% %             ylabel('elevation/A_{linear}','Fontsize',sizefont);
%             ylabel('$\eta/A_{linear}$','Fontsize', sizefont,'interpreter','latex');
%             xlabel('$t*f_c$','Fontsize',sizefont,'interpreter','latex');
%             box on;
%             grid on;
%             xlim([time.(wave{case_idx}).(phase{case_idx})(1)-0.002 time.(wave{case_idx}).(phase{case_idx})(end)]);
%             ylim([0.8 1.4]);
%             yticks([0.8 1 1.2 1.4]);
% %             xlim([62.5 64]);
% %             
% %             
% % % --------------------------------------------------------------------------
% %         cd(dir_plot_save);
% %         export_fig([dir_plot_save,'\',fig_title],'-png','-nocrop');
% % % --------------------------------------------------------------------------  
% % %% time evolution of gammas
% % 
%         sizefont = 18;
%         sizefont1 = 18;
%         bottom = 0.33;
%         left = 0.14; % 0.14
%         width = 0.8; % 0.8
%         height = 0.62;
%         fig_title = ['gamma_',wave{case_idx},'_',phase{case_idx},'_',ROI{1},'_',run{1}];  
%         a1 = figure('Name',[fig_title],'NumberTitle','off'); 
%         axes('Position',[left bottom width height]); hold on;
%         set(a1,'units','normalized','outerposition',[0 0 0.85 0.45]);
%         set(a1,'PaperSize',[8 10]);
%         set(a1,'Color',[1 1 1]); 
%         for i = NUM_in:NUM
%             if crit_idx_bot.(im_name{i-NUM_in+1}) > 0
%                 scatter(time.(wave{case_idx}).(phase{case_idx})(i-NUM_in+1),slope_tip.(im_name{i-NUM_in+1}),40,'+','MarkerEdgeColor','red','LineWidth',1);
%                 scatter(time.(wave{case_idx}).(phase{case_idx})(i-NUM_in+1),slope_toe.(im_name{i-NUM_in+1}),40,'+','MarkerEdgeColor','blue','LineWidth',1);
% 
%             end
%         end
%         hold on;
%         if crit_idx_bot.(im_name{i-NUM_in+1}) > 0
%             plot(time.(wave{case_idx}).(phase{case_idx})(1:numel(slope_tipmm)),slope_tipmm,'-','Linewidth',0.5,'color','black');
%         end
%         for i = NUM_in:NUM 
%             for j = 1:numel(allrts.(im_name{i-NUM_in+1}))
%                 if bulgey.(im_name{i-NUM_in+1})(localmaxima.(im_name{i-NUM_in+1}))> 0 
%                 else
%                     scatter(time.(wave{case_idx}).(phase{case_idx})(i-NUM_in+1),slope_ifp.(im_name{i-NUM_in+1}),40,'+','green','LineWidth',1);
% 
%                 end
%             end
%         end
%         set(gca,'Fontsize',sizefont1);
%         set(gca,'Fontsize',sizefont1);
%         set (gca,'Xdir','reverse')
%         set(gca,'Color',[1 1 1]);
%         ylabel('$\Gamma$','Fontsize',sizefont,'interpreter','latex');
%         xlabel('$t*f_c$','Fontsize',sizefont,'interpreter','latex');
%         box on;
%         grid on;
%         xlim([time.(wave{case_idx}).(phase{case_idx})(1)-0.002 time.(wave{case_idx}).(phase{case_idx})(end)]);
%         ylim([0.3 2.2]);
%         yticks([0.4 0.6 0.8 1 1.6 2.2]);
%         set(gca,'YScale','log','Fontsize',sizefont);
% %             xlim([62.5 64]);
% 
% %             
% % % --------------------------------------------------------------------------
% %         cd(dir_plot_save);
% %         export_fig([dir_plot_save,'\',fig_title],'-png','-nocrop');
% % % -------------------------------------------------------------------------- 
% % %% plot bulge area 
%         sizefont = 15;
%         sizefont1 = 15;
%         bottom = 0.24;
%         left = 0.125;
%         width = 0.85;
%         height = 0.65;
%         bottom1 = 0.27;
%         height1 = 0.6;
%         fig_title = ['bulge area_',wave,'_',phase{case_idx},'_',ROI{1},'_',run{1}];  
%         a1 = figure(10);
%         axes('Position',[left bottom width height]); hold on;
%         set(a1,'units','normalized','outerposition',[0 0 0.9 0.9]);
%         set(a1,'PaperSize',[8 10]);
%         set(a1,'Color',[1 1 1]); 
% %         title(['wave crests and bulges: ', phase{case_idx}]);
%         for i = NUM_in:NUM
%             plot(-bulgex.(im_name{i-NUM_in+1}),bulgey.(im_name{i-NUM_in+1}),'linewidth',0.5,'color','black');
% %             plot(-bulgex_fitlimit.(im_name{i-NUM_in+1}),bulgey_poly.(im_name{i-NUM_in+1}),'linewidth',2,'color','green');
%            hold on;
% %             if bulge_area_points.(im_name{i-NUM_in+1})(:,:) ~=0
%                 scatter(-bulge_area_points.(im_name{i-NUM_in+1})(:,1),bulge_area_points.(im_name{i-NUM_in+1})(:,2),'LineWidth',.5);
% %                 scatter(hor_points.(im_name{i-NUM_in+1})(:,1),y_hypoten.(im_name{i-NUM_in+1}));
% %             end
% %             'MarkerEdgeColor','red','MarkerFaceColor','red'
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Color',[1 1 1]);
%             set(gca,'Xdir','reverse');
%             xlabel('x (mm)','Fontsize',sizefont);
%             ylabel('$\eta$(mm)','Fontsize', sizefont,'interpreter','latex');
% %             ylim([Amp_pseudo_mm-20  max(bulgey_fin.(im_name{end}))+25]);
% %             xlim([-270 -150]);
%             box on;
%             grid on;
%         end
%         pause(2);
% % % %--------------------------------------------------------------------------
% % %         cd(dir_plot_save);
% % %         export_fig([dir_plot_save,'\',fig_title],'-png','-nocrop');
% % % %--------------------------------------------------------------------------
%% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % moving frame
%         sizefont = 22;
%         sizefont1 = 22;
%         bottom = 0.24;
%         left = 0.14;
%         width = 0.79;
%         height = 0.65;
%         bottom1 = 0.27;
%         height1 = 0.6;
%         % fig_title = ['movingframe_plunger',wave,'_',phase{1},'_',ROI{1},'_',run{1}];  
%         fig_title = 'movingframe_purespiller';  
% 
%         a3 = figure(12);
%         axes('Position',[left bottom width height]); hold on;
%         set(a3,'units','normalized','outerposition',[0 0 0.35 0.9]);
%         set(a3,'PaperSize',[8 10]);
%         set(a3,'Color',[1 1 1]); 
% %         title(['wave crests and bulges: ', phase{case_idx}]);
%         hold on;
% %         plot(-x_pp.(im_name{1,NUM_in+25})./lambda_mm,bulgey_mov.(im_name{1,NUM_in+25})./Amp_mm,'linewidth',6,'color','green');
% 
%         plot(-bulgex_mov.(im_name{1,NUM_in+1})./lambda_mm,bulgey_mov.(im_name{1,NUM_in+1})./Amp_mm,'linewidth',6,'color','green');
%         for i = NUM_in+1:3:NUM-4
% %         plot(-x_pp.(im_name{1,i-NUM_in+1})./lambda_mm,bulgey_mov.(im_name{1,i-NUM_in+1})./Amp_mm,'--','linewidth',1,'color','black');
% 
%             plot(-bulgex_mov.(im_name{1,i-NUM_in+1})./lambda_mm,bulgey_mov.(im_name{1,i-NUM_in+1})./Amp_mm,'--','linewidth',1,'color','black');
%         end
%         plot(-bulgex_mov.(im_name{1,end-4})./lambda_mm,bulgey_mov.(im_name{1,end-4})./Amp_mm,'-','linewidth',3,'color','blue');
%         box on;
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Xdir','reverse');
%             set(gca,'Color',[1 1 1]);
%             % xlabel('${x''/\lambda_c}$','Fontsize',sizefont,'interpreter','latex');
%             % ylabel('$\eta/A_{linear}$','Fontsize', sizefont,'interpreter','latex');
%             % set(gca,'TickLabels',[]);
%             ylim([0.7 1.3]);
%             yticks([0.7 0.9 1.1 1.3]);
% %             xlim([-0.23 -0.20]);
%                 grid on; box on;
% % % %--------------------------------------------------------------------------
%         cd(dir_plot_save);
%         exportgraphics(a3,[fig_title,'.emf'],'ContentType','vector');  % other formats pdf, eps, tiff
% %         export_fig([dir_plot_save,'\',fig_title],'-png','-nocrop');
% %--------------------------------------------------------------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    % %--------------------------------------------------------------------------
% % vertical centroid and max elevation trajectory with time
%         sizefont = 15;
%         sizefont1 = 15;
%         bottom = 0.24;
%         left = 0.125;
%         width = 0.85;
%         height = 0.65;
%         fig_title = ['freefall_t_',wave,'_',phase{case_idx},'_',ROI{1},'_',run{1}];  
%         a1 = figure('Name',[fig_title],'NumberTitle','off'); 
%         axes('Position',[left bottom width height]); hold on;
%         set(a1,'units','normalized','outerposition',[0 0 0.5 0.85]);
%         set(a1,'PaperSize',[8 10]);
%         set(a1,'Color',[1 1 1]); 
%         for i = NUM_in+4:NUM-4
%             if crit_idx_bot.(im_name{i-NUM_in+1}) > 0
%             scatter(time(i-NUM_in+1),bulge.(im_name{i-NUM_in+1})(crit_idx_top.(im_name{i-NUM_in+1})(end),2),'MarkerEdgeColor',cc(i-NUM_in+1,:,:),'LineWidth',5)
%             scatter(time(i-NUM_in+1),bulge.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(end),2),'MarkerEdgeColor',cc(i-NUM_in+1,:,:),'LineWidth',5)
% %             scatter(time(i-NUM_in+1),ycm_mtx(i-NUM_in+1),'MarkerEdgeColor',cc(i-NUM_in+1,:,:),'LineWidth',5)
%             end
%             scatter(time(i-NUM_in+1),bulge.(im_name{i-NUM_in+1})(end,2),'MarkerEdgeColor','blue','LineWidth',5)
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Color',[1 1 1]);
%             xlabel('time (sec)','Fontsize',sizefont);
%             ylabel('elevation (mm)','Fontsize', sizefont,'interpreter','latex');
% %             ylim([70 100]);
%             box on;
%             grid on;
%             hold on;
% 
%         end
% %--------------------------------------------------------------------------
%         cd(dir_plot_save);
%         export_fig([dir_plot_save,'\',fig_title],'-png','-nocrop');
% %--------------------------------------------------------------------------
%         sizefont = 15;
%         sizefont1 = 15;
%         bottom = 0.24;
%         left = 0.125;
%         width = 0.85;
%         height = 0.65;
%         time = linspace(63,64,NUM-NUM_in+1);
%         fig_title = ['norm_dt_velx_t_',wave,'_',phase{case_idx},'_',ROI{1},'_',run{1}];  
%         a1 = figure('Name',[fig_title],'NumberTitle','off'); 
%         axes('Position',[left bottom width height]); hold on;
%         set(a1,'units','normalized','outerposition',[0 0 0.5 0.85]);
%         set(a1,'PaperSize',[8 10]);
%         set(a1,'Color',[1 1 1]);
%         for i = NUM_in:NUM-1
% %             if crit_idx_bot.(im_name{i-NUM_in+1}) > 0
% %             scatter(time(i-NUM_in+1),velo_crest_dt.(im_name{i-NUM_in+1})./phase_speed,'MarkerEdgeColor','blue','LineWidth',4)
% %             hold on;
% %             scatter(time(i-NUM_in+1),velo_crest_top_dt.(im_name{i-NUM_in+1})./phase_speed,'MarkerEdgeColor','red','LineWidth',4);
% %             scatter(time(i-NUM_in+1),velo_crest_tip_dt.(im_name{i-NUM_in+1})./phase_speed,'MarkerEdgeColor','black','LineWidth',4);
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Color',[1 1 1]);
%             xlabel('time (sec)','Fontsize',sizefont);
% %             ylabel('instant crest speed/linear phase speed','Fontsize', sizefont,'interpreter','latex');
%             ylabel('U/C','Fontsize', sizefont,'interpreter','latex');
%             ylim([0 3]);
%             box on;
%             grid on;
% %             end
%         end
%         hold on;
%             plot(time(1:end-1),B,'blue','LineWidth',1);
% %         plot(time(1:end-1),velo_crest_dt_smooth_mtx./phase_speed,'blue','LineWidth',1);
% %         plot(time(1:end-1),velo_crest_top_dt_smooth_mtx./phase_speed,'red','LineWidth',1);
% %         plot(time(1:end-1),velo_crest_tip_dt_smooth_mtx./phase_speed,'black','LineWidth',1)
% %--------------------------------------------------------------------------
%         cd(dir_plot_save);
%         export_fig([dir_plot_save,'\',fig_title],'-png','-nocrop');
% % --------------------------------------------------------------------------

% %--------------------------------------------------------------------------
%     fig_title  = ('gradient_edges');
%     a1 = figure('Name',[fig_title],'NumberTitle','off');
%     set(a1,'units','normalized','outerposition',[0. 0. 1. 1.]);
%     set(a1,'PaperSize',[8 10]);
%     set(a1,'Color',[1 1 1]);
% %     for i = NUM_in:NUM-1
% %         imshow(flipud(imbin_g.(im_name{i-NUM_in+1})));
% %         imshow(imbin_g_aux.(im_name{i-NUM_in+1}));
%         imshow(Gmag.(im_name{i-NUM_in+1}));
% 
% %     imshow(images);
% %     end
% %--------------------------------------------------------------------------
%         cd(dir_plot_save);
%         export_fig([dir_plot_save,'\',fig_title],'-png','-nocrop');
% %%-------------------------------------------------------------------------

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------------SAVE RESULTS-----------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd(dir_save);
%--------------------------------------------------------------------------
    save(['bulgex_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'bulgex','-v7.3');
    save(['bulgey_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'bulgey','-v7.3');
    save(['crestx_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'crestx','-v7.3');
    save(['cresty_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'cresty','-v7.3');
%--------------------------------------------------------------------------
    save(['x_pp_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'x_pp','-v7.3');
    save(['bulge_area_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'bulge_area','-v7.3');
if exist('tip_l','var') == 1    
    save(['tip_l_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'tip_l','-v7.3');
    save(['tip_h_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'tip_h','-v7.3');
    save(['tipel_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'tipel','-v7.3');
    save(['tiple_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'tiple','-v7.3');
    save(['gamme_tip_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'slope_tip','-v7.3');

    save(['toe_l_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'toe_l','-v7.3');
    save(['toe_h_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'toe_h','-v7.3');
    save(['toeel_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'toeel','-v7.3');
    save(['toele_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'toele','-v7.3');
    save(['gamme_toe_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'slope_toe','-v7.3');
end

    save(['ifp_h',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'ifp_h','-v7.3');
    save(['ifp_l',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'ifp_l','-v7.3');
    save(['gamme_ifp_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'slope_ifp','-v7.3');

    save(['break_vol_idx_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'break_vol_ifpidx','-v7.3');
    save(['timeorig_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'time_loc_run','-v7.3');
    

end





%  %% lost project: fitting a curve to non-breaking crests and comparing with crests with bulges   
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
%     for i = NUM_in:NUM
%         bulge_mov.(im_name{i-NUM_in+1})(:,1) =  round(bulgex_mov.(im_name{i-NUM_in+1}),2);
%         bulge_mov.(im_name{i-NUM_in+1})(:,2) =  round(bulgey_mov.(im_name{i-NUM_in+1}),2);
%     
%     end
% 
%     
%     % find which index of crest 2 corresponds to max crest height of other
%     % crests
%     for i = NUM_in+1:NUM
%         bulge_mov_top_idx.(im_name{i-NUM_in+1}) = find(bulge_mov.(im_name{2})(:,2) == bulge_mov.(im_name{i-NUM_in+1})(end,2));
%     end
%     for i = NUM_in+1:NUM
%         if isempty(bulge_mov_top_idx.(im_name{i-NUM_in+1}))
%             bulge_mov_top_idx.(im_name{i-NUM_in+1}) = 0;
%         else
%              bulge_mov_top_idx.(im_name{i-NUM_in+1}) =  bulge_mov_top_idx.(im_name{i-NUM_in+1})(1);
%         end
%     end
%     for i = NUM_in+1:NUM
%         if bulge_mov_top_idx.(im_name{i-NUM_in+1})(1) == 0
%             bulge_mov_top_idx.(im_name{i-NUM_in+1})(1) = bulge_mov_top_idx.(im_name{i-NUM_in})(1);
%         end
%     end
%     
%     for i = NUM_in+1:NUM
%         idx_down.(im_name{i-NUM_in+1}) = bulge_mov_top_idx.(im_name{i-NUM_in+1}) - numel(bulge_mov.(im_name{i-NUM_in+1})(:,2));
%     end
% 
% 
%   % rearrange second bulge to compare with each crest  
%     for i = NUM_in+1:NUM
%         if idx_down.(im_name{i-NUM_in+1})>0
%             bulge_mov_tr.(im_name{i-NUM_in+1}) = bulge_mov.(im_name{2})(idx_down.(im_name{i-NUM_in+1}):bulge_mov_top_idx.(im_name{i-NUM_in+1}),:);
%         else 
%             bulge_mov_tr.(im_name{i-NUM_in+1}) = bulge_mov.(im_name{2});
%         end
%     end
%     % points of intersection
%     for i = NUM_in+1:NUM
%         [x0.(im_name{i-NUM_in+1}),y0.(im_name{i-NUM_in+1})] = intersections(bulge_mov.(im_name{i-NUM_in+1})(:,1),...
%             bulge_mov.(im_name{i-NUM_in+1})(:,2),bulge_mov_tr.(im_name{i-NUM_in+1})(:,1),bulge_mov_tr.(im_name{i-NUM_in+1})(:,2));
%         
%         x0.(im_name{i-NUM_in+1}) = round(x0.(im_name{i-NUM_in+1}),1);
%         y0.(im_name{i-NUM_in+1}) = round(y0.(im_name{i-NUM_in+1}),1);
%         bulge_mov.(im_name{i-NUM_in+1}) = round(bulge_mov.(im_name{i-NUM_in+1}),1);
%     end    
%     
%     for i = NUM_in+1:NUM
%         if numel(x0.(im_name{i-NUM_in+1})) == 1
%             idx_com.(im_name{i-NUM_in+1}) = find(bulge_mov.(im_name{i-NUM_in+1})(:,1) <= x0.(im_name{i-NUM_in+1}) + 0.2 ...
%             & bulge_mov.(im_name{i-NUM_in+1})(:,1) >= x0.(im_name{i-NUM_in+1}) - 0.2);            
%         else
%             idx_com.(im_name{i-NUM_in+1}) = find(bulge_mov.(im_name{i-NUM_in+1})(:,1) <= x0.(im_name{i-NUM_in+1})(end-1) + 0.2 ...
%             & bulge_mov.(im_name{i-NUM_in+1})(:,1) >= x0.(im_name{i-NUM_in+1})(end-1) - 0.2); 
%         end
%     end
%     % find indices of points of intersection in second crest
%     for i = NUM_in+1:NUM
%         if numel(x0.(im_name{i-NUM_in+1})) == 1
%             idx_com_first.(im_name{i-NUM_in+1}) = find(bulge_mov_tr.(im_name{i-NUM_in+1})(:,1) <= x0.(im_name{i-NUM_in+1}) + 0.2 ...
%             & bulge_mov_tr.(im_name{i-NUM_in+1})(:,1) >= x0.(im_name{i-NUM_in+1}) - 0.2);            
%         else
%         idx_com_first.(im_name{i-NUM_in+1}) = find(bulge_mov_tr.(im_name{i-NUM_in+1})(:,1) <= x0.(im_name{i-NUM_in+1})(end-1) + 0.2 ...
%             & bulge_mov_tr.(im_name{i-NUM_in+1})(:,1) >= x0.(im_name{i-NUM_in+1})(end-1) - 0.2);
%         end
%     end
%     for i = NUM_in+1:NUM
%         %fill empties with zeros
%         if isempty(idx_com.(im_name{i-NUM_in+1})) 
%             idx_com.(im_name{i-NUM_in+1}) = 0;
%         end
%          if isempty(idx_com_first.(im_name{i-NUM_in+1})) 
%             idx_com_first.(im_name{i-NUM_in+1}) = 0;
%          end   
%          idx_com.(im_name{i-NUM_in+1}) = idx_com.(im_name{i-NUM_in+1})(1);
%          idx_com_first.(im_name{i-NUM_in+1}) = idx_com_first.(im_name{i-NUM_in+1})(1);
%         %put values in matrix
%         idx_com_mtx(i-NUM_in+1) = idx_com.(im_name{i-NUM_in+1})(1);
%         idx_com_first_mtx(i-NUM_in+1) = idx_com_first.(im_name{i-NUM_in+1})(1);
% 
%     end
%     % these are the final intersection indices
%     for i = NUM_in+1:NUM
%         if idx_com_mtx(i-NUM_in+1) == 0
%             idx_com_mtx(i-NUM_in+1) = idx_com_mtx(i-NUM_in) - 10;
%         end
%         if idx_com_first_mtx(i-NUM_in+1) == 0
%             idx_com_first_mtx(i-NUM_in+1) = idx_com_first_mtx(i-NUM_in) - 10;
%         end
%     end
% 
%     
%     % choose part of wave crest above first crest
%   for i = NUM_in+1:NUM
%     bulge_mov_sel_crt.(im_name{i-NUM_in+1}) = bulge_mov.(im_name{i-NUM_in+1})(idx_com_mtx(i-NUM_in+1):end,:);
%   end
%   % choose part of first crest below each crest
%    for i = NUM_in+1:NUM
%     bulge_mov_sel_fi.(im_name{i-NUM_in+1}) = bulge_mov_tr.(im_name{i-NUM_in+1})(idx_com_first_mtx(i-NUM_in+1):end,:);
%    end
%  
% % concacenate points found just before 
%     for i = NUM_in+1:NUM
%         % x
%     bulge_area_xpoints.(im_name{i-NUM_in+1}) =[bulge_mov_sel_crt.(im_name{i-NUM_in+1})(:,1);bulge_mov_sel_fi.(im_name{i-NUM_in+1})(:,1)];
%         % y
%     bulge_area_ypoints.(im_name{i-NUM_in+1}) =[bulge_mov_sel_crt.(im_name{i-NUM_in+1})(:,2);bulge_mov_sel_fi.(im_name{i-NUM_in+1})(:,2)];
% 
%     end
%     for i = NUM_in+1:NUM
%             [bulge_area_index.(im_name{i-NUM_in+1}),bulge_vol_spi.(im_name{i-NUM_in+1}).(im_name{i-NUM_in+1})] = ...
%                 boundary(bulge_area_xpoints.(im_name{i-NUM_in+1}),bulge_area_ypoints.(im_name{i-NUM_in+1}),0);
%         
%     end









% %%--------------------------------------------------------------------------
%         sizefont = 15;
%         sizefont1 = 15;
%         bottom = 0.24;
%         left = 0.125;
%         width = 0.85;
%         height = 0.65;
%         time = linspace(63,64,NUM-NUM_in+1);
%         fig_title = ['surface_elevation_t_',wave,'_',phase{case_idx},'_',ROI{1},'_',run{1}];  
%         a1 = figure('Name',[fig_title],'NumberTitle','off'); 
%         axes('Position',[left bottom width height]); hold on;
%         set(a1,'units','normalized','outerposition',[0 0 0.5 0.85]);
%         set(a1,'PaperSize',[8 10]);
%         set(a1,'Color',[1 1 1]);
%         for i = NUM_in:NUM
%             scatter(time(i-NUM_in+1),bulge.(im_name{i-NUM_in+1})(end,2),'MarkerEdgeColor','blue','LineWidth',5)
%             hold on;
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Color',[1 1 1]);
%             xlabel('time (sec)','Fontsize',sizefont);
%             ylabel('\eta','Fontsize', sizefont,'interpreter','latex');
% %             ylim([0 3]);
%             box on;
%             grid on;
%         end
% %--------------------------------------------------------------------------
%         cd(dir_plot_save);
%         export_fig([dir_plot_save,'\',fig_title],'-png','-nocrop');
% horizontal crest velocity to linear phase velocity ratio
%         sizefont = 15;
%         sizefont1 = 15;
%         bottom = 0.24;
%         left = 0.125;
%         width = 0.85;
%         height = 0.65;
%         time = linspace(63,64,NUM-NUM_in+1);
%         fig_title = ['norm_velx_t_comparisons_',wave,'_',phase{case_idx},'_',ROI{1},'_',run{1}];  
%         a1 = figure('Name',[fig_title],'NumberTitle','off'); 
%         axes('Position',[left bottom width height]); hold on;
%         set(a1,'units','normalized','outerposition',[0 0 0.5 0.85]);
%         set(a1,'PaperSize',[8 10]);
%         set(a1,'Color',[1 1 1]);
%         for i = NUM_in:NUM-1
% %             if crit_idx_bot.(im_name{i-NUM_in+1}) > 0
% %             scatter(time(i-NUM_in+1),crest_avg_vel_alt_mtx(i-NUM_in+1)./phase_speed,'MarkerEdgeColor','black','LineWidth',5);
%             scatter(time(i-NUM_in+1),velo_xcm_mtx(i-NUM_in+1)./phase_speed,'MarkerEdgeColor','black','LineWidth',5)
% 
%             hold on;
%             scatter(time(i-NUM_in+1),velo_crest_avg_dt.(im_name{i-NUM_in+1})./phase_speed,'MarkerEdgeColor','blue','LineWidth',5)
%             scatter(time(i-NUM_in+1),velo_crest_top_avg_dt.(im_name{i-NUM_in+1})./phase_speed,'MarkerEdgeColor','red','LineWidth',5)
%             scatter(time(i-NUM_in+1),crest_avg_vel.(im_name{i-NUM_in+1})./phase_speed,'MarkerEdgeColor','green','LineWidth',5)
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Color',[1 1 1]);
%             xlabel('time (sec)','Fontsize',sizefont);
%             ylabel('norm centroid crest - entire crest speed','Fontsize', sizefont,'interpreter','latex');
%             ylim([0. 1.5]);
%             box on;
%             grid on;
% %             end
%         end
% %--------------------------------------------------------------------------
%         cd(dir_plot_save);
%         export_fig([dir_plot_save,'\',fig_title],'-png','-nocrop');
% %%--------------------------------------------------------------------------
% % correlation between bulge velocity (x) and bulge area
%         sizefont = 15;
%         sizefont1 = 15;
%         bottom = 0.24;
%         left = 0.125;
%         width = 0.85;
%         height = 0.65;
% 
%         fig_title = ['correlation_b_velx_',wave,'_',phase{case_idx},'_',ROI{1},'_',run{1}];  
%         a1 = figure('Name',[fig_title],'NumberTitle','off'); 
%         axes('Position',[left bottom width height]); hold on;
%         set(a1,'units','normalized','outerposition',[0 0 0.5 0.85]);
%         set(a1,'PaperSize',[8 10]);
%         set(a1,'Color',[1 1 1]); 
%         for i = NUM_in:NUM
%             if crit_idx_bot.(im_name{i-NUM_in+1}) > 0
%             scatter(velo_xcm_mtx(i-NUM_in+1),bulge_area_mtx(i-NUM_in+1),'MarkerEdgeColor',cc(i-NUM_in+1,:,:),'LineWidth',5)
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Fontsize',sizefont1);
%             set(gca,'Color',[1 1 1]);
%             xlabel('horizontal velocity (mm/sec)','Fontsize',sizefont);
%             ylabel('$\mathrm{area}$ ($\mathrm{mm^{2})}$','Fontsize', sizefont,'interpreter','latex');
% %             xlim([1200 2000]);
%             box on;
%             grid on;
%             hold on;
%             end
%         end
% %--------------------------------------------------------------------------
%         cd(dir_plot_save);
%         export_fig([dir_plot_save,'\',fig_title],'-png','-nocrop');
% %-------------------------------------------------------------------------- 





















% %% save text files
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     % frames_in.peak = 31855;  % initial frame peak ROI3
%     % frames_end.peak = 31913;
%     frames_in = 31423;  %first frame
%     frames_end = 31480;
%     time = round((frames_in:frames_end)./500 - 64,3);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % x-z coordinates
% cd('C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\cloud data\paper1\wave crest profiles');
% for i = 1:NUM
%     fid = fopen(['crest.',wave{1},'.',phase{1},'.',num2str(time(i),'%05.3f'),'sec.txt'],'wt');
%     header1 = 'x (mm)';
%     header2 = 'z (mm)';
%     fprintf(fid, [ ' ' header1 '   ' header2 '\r\n']);
%     fprintf(fid, '%4.1f     %4.1f\r\n',bulge.(im_name{i-NUM_in+1})');
%     fclose(fid);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % time with tip toe inf and gamma
% for i = NUM_in:NUM
%     if crit_idx_top.(im_name{i-NUM_in+1})(1) == 0
%         tip(i,1) = NaN;
%         tip(i,2) = NaN;
% 
% 
%         gamma_tip(i,1) = NaN;
%         gamma_toe(i,1) = NaN; 
%     elseif crit_idx_top.(im_name{i-NUM_in+1})(1) > 0
%         tip(i,1) = -bulgex.(im_name{i-NUM_in+1})(crit_idx_top.(im_name{i-NUM_in+1})(1)); % tip
%         tip(i,2) = bulgey.(im_name{i-NUM_in+1})(crit_idx_top.(im_name{i-NUM_in+1})(1)); % tip
% 
% 
%         gamma_tip(i,1) = slope_tip.(im_name{i-NUM_in+1});
%         gamma_toe(i,1) = slope_toe.(im_name{i-NUM_in+1});
% 
%     end
%     if crit_idx_bot.(im_name{i-NUM_in+1})(1) == 0
%         toe(i,1) = NaN;
%         toe(i,2) = NaN;
%     elseif crit_idx_bot.(im_name{i-NUM_in+1})(1) > 0
%         toe(i,1) = -bulgex.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(1)); % toe
%         toe(i,2) = bulgey.(im_name{i-NUM_in+1})(crit_idx_bot.(im_name{i-NUM_in+1})(1)); % toe
%     end
% 
%      if ifpy.(im_name{i-NUM_in+1}) > 0
%         inf(i,1) = -ifpx.(im_name{i-NUM_in+1});
%         inf(i,2) = ifpy.(im_name{i-NUM_in+1});
% 
%         gamma_ifp(i,1) = slope_ifp.(im_name{i-NUM_in+1});
%     elseif ifpy.(im_name{i-NUM_in+1}) == 0
%         inf(i,1) = NaN;
%         inf(i,2) = NaN;
% 
%         gamma_ifp(i,1) = NaN;
%      end  
% end
% % put all in one matrix
% for i = NUM_in:NUM
%     crestparam = [time(1:end-1)', tip(:,1), tip(:,2), toe(:,1), toe(:,2), inf(:,1), inf(:,2), gamma_tip(:,1), gamma_toe(:,1),gamma_ifp(:,1)];
% end
% 
% 
% % for i = 1:NUM
% fid = fopen(['crestparam.',wave{1},'.',phase{1},'.txt'],'wt');
% header1 = 'time (s)';
% header2 = 'tip x-co (mm)'; header3 = 'tip z-co (mm)'; header4 = 'toe x-co (mm)'; header5 = 'toe z-co'; header6 = 'inf x-co'; header7 = 'inf z-co'; 
% header8 = 'gamma tip'; header9 = 'gamma toe'; header10 = 'gamma inf';
% fprintf(fid, [ ' ' header1 '   ' header2 '   ' header3 '   ' header4 '   ' header5 '   ' header6 '  ' header7 '  ' header8  '  ' header9 '  ' header10  '\r\n']);
% fprintf(fid,        '%4.3f       %4.1f        %4.1f            %4.1f       %4.1f         %4.1f        %4.1f         %4.1f        %4.1f        %4.1f \r\n',crestparam');
% fclose(fid);
% % end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %--------------------------------------------------------------------------
%     % save(['criterion_aspect_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'criterion_aspect','-v7.3');
% 
% %     save(['criterion_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['criterion'],'-v7.3');
% %     save(['criterion_bulge_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['criterion_bulge'],'-v7.3');
% %     save(['criterion_no_dim',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['criterion_no_dim'],'-v7.3');
% %--------------------------------------------------------------------------
% %     save(['localmaxima_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['localmaxima'],'-v7.3');
% %     save(['localminima_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['localminima'],'-v7.3');
%     save(['crit_idx_top_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'crit_idx_top','-v7.3');
%     save(['crit_idx_bot_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'crit_idx_bot','-v7.3');
%     save(['slope_tip_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'slope_tip','-v7.3');
%     %--------------------------------------------------------------------------    
% %     save(['bulgex_mov_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['bulgex_mov'],'-v7.3');
% %     save(['bulgey_mov_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['bulgey_mov'],'-v7.3');
% %--------------------------------------------------------------------------    
% %     save(['crest_avg_vel_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['crest_avg_vel'],'-v7.3');
% %     save(['bulge_avg_vel_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['bulge_avg_vel'],'-v7.3');
% %--------------------------------------------------------------------------
% %     save(['bulgex_fin_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['bulgex_fin'],'-v7.3');
% %     save(['bulgey_fin_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['bulgey_fin'],'-v7.3');
% %     save(['bulge_maxima_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['bulge_maxima'],'-v7.3');
% %     save(['bulge_minima_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['bulge_minima'],'-v7.3');
% % %--------------------------------------------------------------------------
% save(['xcm_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'xcm','-v7.3');
%     save(['ycm_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'ycm','-v7.3');
% %     save(['velo_xcm_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['velo_xcm'],'-v7.3');
% %     save(['velo_ycm_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['velo_ycm'],'-v7.3');
%     % save(['velo_ycm_dt_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'velo_ycm_dt','-v7.3');
%     % save(['velo_xcm_dt_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'velo_xcm_dt','-v7.3');
%     % save(['velo_crest_top_dt_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'velo_crest_top_dt_smooth','-v7.3');
%     % save(['velo_crest_tip_dt_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'velo_crest_tip_dt_smooth','-v7.3');
% 
% %--------------------------------------------------------------------------
%     save(['bulge_area_points_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'bulge_area_points','-v7.3');
%     save(['bulge_area_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'bulge_area','-v7.3');
%     % save(['bulge_theta_',wave{case_idx},'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],'bulge_theta','-v7.3');
% %     save(['bulge_m_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['bulge_m'],'-v7.3');
% %     save(['bulge_angle_ratio_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'_.mat'],['bulge_angle_ratio'],'-v7.3');



    








