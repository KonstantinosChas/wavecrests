function [ifpy,ifpx,allrts,der2fy] = crvf(bulgex,bulgey,wave)
% finds inflection points on crest profiles
% allrts,ifpx_all,der2fy
%% 1. not an output: curve fitting x = f(y) (y = f(x) is not a single valued function)
% smooth;
window = 1;
bulgey_smooth = smooth(bulgey,'loess',window);
pp = polyfit(bulgey_smooth,bulgex,14); 
syms x

% 18th order polynomial
% fit_x = @(x) pp(1).*x.^18 + pp(2).*x.^17 + pp(3).*x.^16+ pp(4).*x.^15 + pp(5).*x.^14 + pp(6).*x.^13 + ...
%     pp(7).*x.^12 + pp(8).*x.^11 + pp(9).*x.^10 + pp(10).*x.^9 + pp(11).*x.^8 + pp(12).*x.^7 + ...
%     pp(13).*x.^6 + pp(14).*x.^5 + pp(15).*x.^4 + pp(16).*x.^3 + pp(17).*x.^2 + pp(18).*x + pp(19);

% 14th order polynomial ---
fit_x = @(x) pp(1).*x.^14 + pp(2).*x.^13 + ...
    pp(3).*x.^12 + pp(4).*x.^11 + pp(5).*x.^10 + pp(6).*x.^9 + pp(7).*x.^8 + pp(8).*x.^7 + ...
    pp(9).*x.^6 + pp(10).*x.^5 + pp(11).*x.^4 + pp(12).*x.^3 + pp(13).*x.^2 + pp(14).*x + pp(15);

% 12th order polynomial
% fit_x = @(x) pp(1).*x.^12 + pp(2).*x.^11 + pp(3).*x.^10+ pp(4).*x.^9 + pp(5).*x.^8 + pp(6).*x.^7 + ...
%     pp(7).*x.^6 + pp(8).*x.^5 + pp(9).*x.^4 + pp(10).*x.^3 + pp(11).*x.^2 + pp(12).*x + pp(13);

% 10th order polynomial
% fit_x = @(x) pp(1).*x.^10 + pp(2).*x.^9 + pp(3).*x.^8+ pp(4).*x.^7 + pp(5).*x.^6 + pp(6).*x.^5 + ...
%     pp(7).*x.^4 + pp(8).*x.^3 + pp(9).*x.^2 + pp(10).*x + pp(11);

% 8th order polynomial
% fit_x = @(x) pp(1).*x.^8+ pp(2).*x.^7 + pp(3).*x.^6 + pp(4).*x.^5 + ...
%     pp(5).*x.^4 + pp(6).*x.^3 + pp(7).*x.^2 + pp(8).*x + pp(9);

% 6th order polynomial
% fit_x = @(x) pp(1).*x.^6 + pp(2).*x.^5 + ...
%     pp(3).*x.^4 + pp(4).*x.^3 + pp(5).*x.^2 + pp(6).*x + pp(7);

%% 2. inflection points 
% 1st derivative
der1f = (diff(fit_x(x))); 
% 2nd derivative
der2f = diff(der1f);

% find 2nd der. roots
der2fy = double(solve(der2f == 0));

% set im. = 0 
rootsdyidx = imag(der2fy)==0; 
% keep real values
der2fy = der2fy(rootsdyidx);
% keep values inside data range;
der2fy = der2fy(der2fy<=max(bulgey) & der2fy>35); 
% keep all real roots for checking
allrts = der2fy; % z
% ifpx_all = fit_x(allrts); %
% 
if strcmp(wave,'JS_15') || strcmp(wave,'JS_155') % some special clauses for pure spillers
    ifpy = max(der2fy(der2fy>=0.60 & der2fy<=68));
    if ifpy < 50
      ifpy = max(der2fy(der2fy>=0.50 & der2fy<65)); 
    end
else
    % keep largest real value (realistic)
    ifpy = max(der2fy); % inflection point z-coordinate
end

ifpx = fit_x(ifpy); % inflection point x - coordinate
ifpx = round(ifpx,1);
ifpy = round(ifpy,1);

end

