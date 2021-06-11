function m = minmod(A)

%% MINMOD The minmod function
%    m = minmod(A) applies the so-called minmod function to an array.
%
%     » If A is a vector of length n, then minmod(A) returns s*min(abs(A)),
%       if s = sign(A(1)) = ... = sign(A(n)), otherwise it returns 0.
%       
%       Examples: A = [  1  2  3 ]; minmod(A) =  1
%                 B = [ -2 -3 -4 ]; minmod(B) = -2
%                 C = [ -1  2  3 ]; minmod(C) =  0
%
%     » If A is a matrix, then minmod(A) returns a row vector m where 
%       m(i) = minmod(A(:,i)), i.e., it applies the minmod function as
%       described above to each column of A.
%   
%       Example: A = [ 1 -2 -1; 2 -3 2; 3 -4 3 ]; minmod(A) = [ 1 -2 0 ]
%

%% Validate input

ip = inputParser;
vA = @(x)validateattributes(x,{'numeric'},{'2d'});
ip.addRequired('A',vA); 
ip.parse(A);
ip.Results;

%% Apply the minmod function

switch min(size(A))
    case 1
        m = isequal(abs(sum(sign(A))),length(A))*sign(A(1))*min(abs(A));
    otherwise
        m = (abs(sum(sign(A))) == size(A,1)).*sign(A(1,:)).*min(abs(A));
end

end