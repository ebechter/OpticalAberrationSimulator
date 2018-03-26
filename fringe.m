function [n, m, sc] = fringe(j)
%
% Note that 'FRINGE', 'ISO', 'WYANT' use this function
% to assign (n, m) pairs to the j values.
%
% sc = 's' = sin(),  sc = 'c' = cos(),  sc = ' ' for neither.
%

switch j
    case 1
        n = 0; m = 0; sc = ' ';
    case 2
        n = 1; m = 1; sc = 'c';
    case 3
        n = 1; m = 1; sc = 's';
    case 4
        n = 2; m = 0; sc = ' ';
    case 5
        n = 2; m = 2; sc = 'c';
    case 6
        n = 2; m = 2; sc = 's';
    case 7
        n = 3; m = 1; sc = 'c';
    case 8
        n = 3; m = 1; sc = 's';
    case 9
        n = 4; m = 0; sc = ' ';
    case 10
        n = 3; m = 3; sc = 'c';
    case 11
        n = 3; m = 3; sc = 's';
    case 12
        n = 4; m = 2; sc = 'c';
    case 13
        n = 4; m = 2; sc = 's';
    case 14
        n = 5; m = 1; sc = 'c';
    case 15
        n = 5; m = 1; sc = 's';
    case 16
        n = 6; m = 0; sc = ' ';
    case 17
        n = 4; m = 4; sc = 'c';
    case 18
        n = 4; m = 4; sc = 's';
    case 19
        n = 5; m = 3; sc = 'c';
    case 20
        n = 5; m = 3; sc = 's';
    case 21
        n = 6; m = 2; sc = 'c';
    case 22
        n = 6; m = 2; sc = 's';
    case 23
        n = 7; m = 1; sc = 'c';
    case 24
        n = 7; m = 1; sc = 's';
    case 25
        n = 8; m = 0; sc = ' ';
    case 26
        n = 5; m = 5; sc = 'c';
    case 27
        n = 5; m = 5; sc = 's';
    case 28
        n = 6; m = 4; sc = 'c';
    case 29
        n = 6; m = 4; sc = 's';
    case 30
        n = 7; m = 3; sc = 'c';
    case 31
        n = 7; m = 3; sc = 's';
    case 32
        n = 8; m = 2; sc = 'c';
    case 33
        n = 8; m = 2; sc = 's';
    case 34
        n = 9; m = 1; sc = 'c';
    case 35
        n = 9; m = 1; sc = 's';
    case 36
        n = 10; m = 0; sc = ' ';
    case 37
        n = 12; m = 0; sc = ' ';
        
end % switch j

end % function fringe