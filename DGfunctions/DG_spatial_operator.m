function RHS = DG_spatial_operator(Uh,DG)


% Compute U at integration points
%------------------------------------------------------------------
U.n = DG.PHI.points*Uh;
% Compute U at boundary points
%------------------------------------------------------------------
U.minus = DG.PHI.minus'*Uh; U.plus = DG.PHI.plus'*Uh;
% Enforce periodic boundary conditions
%------------------------------------------------------------------
if strcmp(DG.ProbDef.BCtype, 'periodic')
  U.minus = [ U.minus(end); U.minus ];
  U.plus  = [ U.plus; U.plus(1) ];
end
% Compute numerical flux
%------------------------------------------------------------------
C = max([abs(DG.ProbDef.c(U.minus)) abs(DG.ProbDef.c(U.plus))],[],2);
fhat = 1/2*(DG.ProbDef.f(U.plus) + DG.ProbDef.f(U.minus) - C.*(U.plus-U.minus));
% Compue RHS
%------------------------------------------------------------------
RHS = DG.A*DG.ProbDef.f(U.n) + DG.B*fhat;
