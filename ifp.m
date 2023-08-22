function [outputArg1,outputArg2] = ifp(inputArg1,inputArg2)


der1f = (diff(fit_x(x))); 
der2f = diff(der1f);

% der1fy = roots(polyder(polyder(pp)));

der2fy = double(solve(der2f == 0));
rootsdyidx = imag(der2fy)==0; 
der2fy = der2fy(rootsdyidx);% keep reals 
der2fy = der2fy(der2fy<=max(cr_1_smu(:,3))); % keep values inside discrere points;
ifp = max(der2fy);% keep largest value (realistic)


% rootsdy = max(rootsdy);  % find max;

der1fy = smooth(der1fy,'moving',window);


der2f = smooth(der2f,'moving',window);

for j = 1:numel(ind)-3
    crvf(j) = der2f(j).*der2f(j+1);
end
a = 1;
for j = 1:numel(ind)-3
    if crvf(j) < 0
        cp(a) = j;
        a = a + 1;
    end
end










end

