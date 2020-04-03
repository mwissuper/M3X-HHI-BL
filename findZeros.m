function x = findZeros(x,n)

% Look for long periods of zeros in marker data since Nexus fills gaps with
% zeros. Replace with nans before processing. Find periods of at least n
% length

ind = find(x==0);
a = diff(ind);
b = findConsec(a,n,0); % indices of x where there are 20 or more 0's in a row