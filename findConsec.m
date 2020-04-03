function tau = findConsec(t,N,first)


% Code to find first or all indices of t with followed by N consecutive numbers from http://www.mathworks.com/matlabcentral/answers/86420-find-a-series-of-consecutive-numbers-in-a-vector
x = diff(t)==1;
b = size(x);
if b(1) == 1 % row vec
    f = find([false x]~=[x false]);
else % col vec
    f = find([false;x]~=[x;false]);
end
if first == 1 % Just find first time this happens
    g = find(f(2:2:end)-f(1:2:end-1)>=N,1,'first');
else
    g = find(f(2:2:end)-f(1:2:end-1)>=N); % Find all instances where N in a row in range
end

tau = t(f(2*g-1)); % t followed by >=N consecutive numbers