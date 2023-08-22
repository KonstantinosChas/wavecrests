function [crest,bulge,bulge_2,max_idx,max_bulge,imbin_g,imbin_g_aux,Gmag] = wave_extraction(images)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extracts wave profile from image
%input: image
% output: x-z wave crest coordinates.
% parameters
    scaling_factor = 1./5.5;   % -->px to mm  1/5.5 for freshwater
    cut = 1;             % cut crests: 0:NO 1:YES cuts them at max. surf. el.
    ordinate_min = 0;    % specify the cut-off ordinate 
    abscissa_max = 300;  % specify the cut-off abscissa - old parameter was 220 - now is obsolete - set large value so it doesn't affect anything
    cutofflim = 0;       % cut-off lim in number of elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%% Part 1: Image Pre-processing
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fsize = 35;  % pixel size of filter, default: 21, 9-11 for JS, 35-55 for GW and smoother results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rgb_gray = rgb2gray(images);
    % imbin = im2bw(images,graythresh(images)); 
    fd = im2double(images);  % convert to double
    h = fspecial('gaussian',fsize);
    g = sqrt(imfilter(fd,h,'replicate').^2 + imfilter(fd,h','replicate').^2);  % smooth
    % h1 = fspecial('prewitt');
    % g1  = sqrt(imfilter(fd,h1,'replicate').^2 + imfilter(fd,h1','replicate').^2);   
    h2 = fspecial('average',fsize);
    g2 = sqrt(imfilter(fd,h2,'replicate').^2 + imfilter(fd,h2','replicate').^2);  %  smooth
%-----------------image intensity gradients:-------------------------------
    imbin_fil = im2bw(g2,graythresh(g2));   % convert to double 
    [Gmag,Gdir] = imgradient(imbin_fil, 'central'); %intensity gradient
    % [Gx,Gy] = imgradientxy(imbin_fil, 'central');
    imbin_g_aux = im2bw(g2,graythresh(Gmag));
    % imbin_g = im2bw(imbin_g,graythresh(Gmag-1));
    % [Gmag_2,Gdir_2] = imgradient(imbin_g_aux, 'central');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%--------------------------------------------------------------------------
%% Part 2 : Main Image Processing (crest extraction) 
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imbin_g = edge(imbin_g_aux,'Canny',0.85); %0.85: default ,  0.95 alternative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imbin_g = flipud(imbin_g);  % flip y-coordinates, because in images (0,0) is at top left corner
    [crest,L] = bwboundaries(imbin_g,'holes');   % trace wave crest boundaries in bin image
    % stats =  bwconncomp(imbin_g);
    % trace only in the vicinity of crest pixels (i.e. zoomed image)
    [crest_aux,crest_init_idx] = max(cellfun(@length,crest));
    boundary = crest{crest_init_idx}; %(:,:);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    cresty = boundary(:,1).*scaling_factor;   
    crestx = boundary(:,2).*scaling_factor;
    crest = [crestx, cresty];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% APPLY Additional filters on extracted crests to smooth result
    windowsize = 35; % 15-35  default
    % windowend = 800; % 500
    b = (1/windowsize)*ones(1,windowsize); % function for smoothing filter
    a = 1;
    cresty_2 = filter(b,a,cresty); 
    crestx_2 = filter(b,a,crestx);
    crest_2  = [crestx_2(windowsize:end),cresty_2(windowsize:end)];
    % cresty_2 = cresty_2(windowsize:end); 
    % crestx_2 = crestx_2(windowsize:end);
%        clear crest;
   crest = crest_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3 : Post-Processing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Initial crest Analysis: cut to elliminate non physical shapes 
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [maximum_x ,ix] = min(abs(abscissa_max  - crestx));  % find the min distance from the abscissa and its index 
%--------------------------------------------------------------------------
    if cut == 1 
        bulgey = cresty(1:ix);    % cut the extracted crests: x-dir  
        bulgex = crestx(1:ix);    
    elseif cut == 0
        bulgey = cresty(1:end);     % don't cut the crests  
        bulgex = crestx(1:end); 
    end
%--------------------------------------------------------------------------
    minimum_ordinate = (abs(ordinate_min  + bulgey(1,1)));
%--------------------------------------------------------------------------
%     [minimum_y.(folderfieldname{i}),iy.(folderfieldname{i})] = min(abs(ordinate_min  - bulgey.(folderfieldname{i})));  % find the min distance from the ordinate and its index 
   iy = find(bulgey == minimum_ordinate);
%--------------------------------------------------------------------------
    bulgey_2 = bulgey(iy:end);    % reshape again the crests: z-dir
    bulgex_2 = bulgex(iy:end);     
%--------------------------------------------------------------------------
    if isempty(bulgey_2) == 1
%     disp('field is empty');
        bulgey_2 = 0;
    else
        bulgey_2 =  bulgey_2;
    end
%--------------------------------------------------------------------------
    if isempty(bulgex_2) == 1
%         disp('field is empty');
        bulgex_2 = 0;
    else
        bulgex_2 =  bulgex_2;
    end
%--------------------------------------------------------------------------
    bulge_2  = [bulgex_2,bulgey_2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eliminate erroneous crests with very few pixels  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ny = numel(bulgey_2);
    nx = numel(bulgex_2);
%--------------------------------------------------------------------------
    if nx < 100
        bulgex_2 =  zeros(max(nx),1);
    else
    end
%--------------------------------------------------------------------------
    if ny < 100
        bulgey_2 =  zeros(max(ny),1);
    else
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the cutting again but depending on where the top of the bulge is.
% perform method twice for better results (find and the re-find the max of maxima)
% then cut the crest after the maxima and eliminate non-physical
% distrurbances of the wave crests in the images
    localmaxima_crest = islocalmax(bulgey_2,'FlatSelection', 'first');   
%--------------------------------------------------------------------------
    idx = find(localmaxima_crest);
%--------------------------------------------------------------------------
    max_bulge = max(bulgey_2);
%--------------------------------------------------------------------------
    max_idx = find(bulgey_2 == max_bulge(end));
    max_idx = max_idx(end);
%--------------------------------------------------------------------------
%     clear bulge bulgex bulgey;
%--------------------------------------------------------------------------
    bulgey_3 = bulgey_2(1:max_idx(1) + cutofflim);     % re-shape the crest data in z-dir: cut non-physical shapes
    bulgex_3 = bulgex_2(1:max_idx(1) + cutofflim);     % the same for x-dir
% apply filters on bulge_3
    bulgey_3 = filter(b,a,bulgey_3); 
    bulgex_3 = filter(b,a,bulgex_3);
    bulgey_3 = bulgey_3(windowsize:end); 
    bulgex_3 = bulgex_3(windowsize:end);
% keep significant figures based on image resolution
    bulgey =  round(bulgey_3,1);    
    bulgex =  round(bulgex_3,1);    
%--------------------------------------------------------------------------
    % bulge_3  = [bulgex_3,bulgey_3];
    bulge  = [bulgex,bulgey];
%     max_idx = find(bulgey == max(bulgey));

%--------------------------------------------------------------------------
%     fig_a1  = ('filters+gradients+edge detection');
%     a1 = figure('Name',[fig_a1],'NumberTitle','off');
%     set(a1,'units','normalized','outerposition',[0. 0. 1. 1.]);
%     set(a1,'PaperSize',[8 10]);
%     set(a1,'Color',[1 1 1]);
%     imshow(flipud(imbin_g));
% %     imshow(images);
%     pause(0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


