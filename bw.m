function [n, m, sc] = bw(j)
% sc = 's' = sin(),  sc = 'c' = cos(),  sc = ' ' for neither.
% sc = 1 = sin(),  sc = 2 = cos().

% calculate the n value
n1 = (-1 + sqrt(1 + 8 * j)) / 2;
n = floor(n1);
if n1 == n
    n = n - 1;
end % if statement

% calculate the m value
k = (n+1)*(n+2)/2;
d = k - j;
m1 = n - 2*d;
m = abs(m1);

% calculate the sc value
sc = ' ';
if m1 > 0
    sc = 's';
end % if statement
if m1 < 0
    sc = 'c';
end % if statement

end % function bw