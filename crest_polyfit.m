function [x_pp] = crest_polyfit(bulgex,bulgey)
% fit polunomial curve to original discrete crest profile
% smooth
window = 41;
bulgey_smooth = smooth(bulgey,'loess',window);
[pp,~,mu] = polyfit(bulgey_smooth,bulgex,14); 
x_pp = polyval(pp,bulgey,[],mu);
% smooth
window = 31;
x_pp = smooth(x_pp,'moving',window);
end