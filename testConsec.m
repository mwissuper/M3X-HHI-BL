x = [34 35 36 78 79 80 81 82 84 85 86 102 103 104 105 106 107 201 202 203 204];
[b, n, idx] = RunLength(x - (1:length(x)));
match  = (n > 5);
result = x(idx(match));