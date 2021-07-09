function DoFs = slopeLimiter(DoFs,p,N,m)

%% SLOPELIMITER Applies the minmod slope limiter
%    DoFs = slopeLimiter(DoFs,p,N) applies the minmod slope limiter to the
%    degrees of freedom, DoFs, of a DG solution of degree p for N elements.
%    See pages 193 and 194 of [1] for details.
%
%    Note: This implementation assumes the Legendre polynomials as a basis
%    for the DG solution and periodic boundary conditions.
%
%--------------------------------------------------------------------------
%    Reference:
%
%    [1] Bernardo Cockburn and Chi-Wang Shu, "Runge--Kutta Discontinuous
%        Galerkin methods for convection-dominated problems," Journal of
%        Scientific Computing, 16(3), 173–261, 2001.
%

%% Validate input

vDoFs = @(x)validateattributes(x,{'numeric'},{'column','size',[N*(p+1) 1]});
vp = @(x)validateattributes(x,{'numeric'},{'scalar','integer','>',0});
vN = @(x)validateattributes(x,{'numeric'},{'scalar','integer','>',0});
ip = inputParser;
ip.addRequired('DoFs',vDoFs);
ip.addRequired('p',vp); ip.addRequired('N',vN);
ip.parse(DoFs,p,N); ip.Results;

%% Apply the slope limiter

% Reshape the DoFs into a (p+1) by N array
u = reshape(DoFs,p+1,N);
% Store the mean value of the DG solution over each element
uBar = u(1,:);
% Compute the DG solution at the end points of the element
uPlus  = (-ones(1,p+1)).^(0:p)*u;
uMinus = ones(1,p+1)*u;
% Compute the arguments used in the minmod function
u1 = [uMinus - uBar; uBar - uPlus];
u2 = uBar - [uBar(end), uBar(1:end-1)];
u3 = [uBar(2:end), uBar(1)] - uBar;
% Apply the minmod function
[minmodMinus,iMinus] = minmod([u1(1,:); u2; u3]);
[minmodPlus,iPlus] = minmod([u1(2,:); u2; u3]);
% Identify elements where |u1| <= m (see page 195 of [1])
iMinus(abs(u1(1,:))<=m) = 1;
iPlus(abs(u1(2,:))<=m) = 1;
% Find elements to which slope limiter is applied
iMod = or(iMinus~=1,iPlus~=1);
% Compute the modified endpoints (see Eqs 2.10 and 2.11 on page 193 of [1])
uMinusMod(iMod) = uBar(iMod) + minmodMinus(iMod);
uPlusMod(iMod)  = uBar(iMod) - minmodPlus(iMod);
% Compute the new slope of the element
u(2,iMod) = (uMinusMod(iMod) - uPlusMod(iMod))/2;
% If p > 1, then...
if p > 1
    % Zero out higher-order DoFs
    u(3:end,iMod) = 0;
end
% Save the modified DoFs
DoFs = u(:);

end
