%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------Wave Crests Calibration---------------------------
% what is does: 
% 1. Calibration for lens distortion
% 2. Saves cal. images
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
%--------------------------------------------------------------------------
%% Part 0 : parameter allocation
%--------------------------------------------------------------------------
wave_type = 'saltwater';
wave = {'GW_18'};
% phase = {'peak','peak','peak'};
% phase = {'trough','trough','trough'};
% phase = {'minp2','minp2','minp2'};
phase = {'posp2','posp2','posp2'};

run = {'run1','run2','run3'};
% run = {'run1'};
[ROI,NUM_c,frames,time] = wave_time(wave{1},phase{1},wave_type);
%--------------------------------------------------------------------------
NUM_case = numel(run);
%--------------------------------------------------------------------------
NUM_in = 1;
NUM = NUM_c;
%-----------------------------Parameters-----------------------------------
% date = ['191020'];  % not used: freshwater/saltwater replaced it
lens = '60mm';
imtype = 'jpg';
num_pix_x = 1280;
num_pix_y = 800;
% --------------------------------------------------------------------------
for case_idx = 1:NUM_case
%% ------------------------------Paths-------------------------------------
%--------------------------------------------------------------------------    
    dir_target = '60mm';       % u save calibrated images in same folder as raw
    dir = ['C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\wave crests\images\',wave_type,'\',lens,'\',wave{1},'\',ROI{case_idx},'\',run{case_idx}];  % directory of images to analyze
%-------------------------------------------------------------------------
    addpath(dir);
    addpath('C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\wave crests\matlab');
    cd('C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\wave crests\matlab');
    cd(dir);
%% ----------------------------Load Images---------------------------------
    im_name = cell(1,NUM-NUM_in+1);
    name =[wave{1},'_',phase{case_idx}];   % names of images (e.g. ...peak_0001)
    name_0 = [wave{1},'_',phase{case_idx},'_000'];
    name_10 = [wave{1},'_',phase{case_idx},'_00'];
    name_100 = [wave{1},'_',phase{case_idx},'_0'];
    name_1000 = [wave{1},'_',phase{case_idx},'_'];
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
    clear images
    for i = NUM_in : NUM
        folderName        =   imread([im_name{i-NUM_in+1},'.',imtype]);
        images.(im_name{i-NUM_in+1}) = folderName;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------------Calibration----------------------------------
%% --------------------------Save images-----------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = NUM_in : NUM
        images.(im_name{i-NUM_in+1}) = camcal(20,[im_name{i-NUM_in+1},'.',imtype],dir,dir_target);
        cd(dir)
        images_aux = images.(im_name{i-NUM_in+1});
%         save([im_name{i-NUM_in+1},'_cal.jpg'],images_aux);
        imwrite(images_aux, [im_name{i-NUM_in+1},'_cal.jpg']);
    end
end




