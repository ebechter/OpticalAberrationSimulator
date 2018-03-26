function [n, m, sc] = mahajan(j)
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
m = n - 2 * floor((k - j)/2);

% calculate the sc value
sc = ' ';
if (m ~= 0) && (mod(j,2) ~= 0)
    % m ~= 0 and j odd
    sc = 's';
end % if statement
if (m ~= 0) && (mod(j,2) == 0)
    % m ~= 0 and j even
    sc = 'c';
end % if statement

end % function mahajan