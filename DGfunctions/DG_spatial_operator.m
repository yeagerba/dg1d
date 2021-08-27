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
elseif strcmp(DG.ProbDef.BCtype, 'dirichlet')
  U.minus = [DG.ProbDef.BCL; U.minus];
  U.plus  = [U.plus; DG.ProbDef.BCR];
elseif strcmp(DG.ProbDef.BCtype, 'dirichletfnt')
  U.minus = [DG.ProbDef.BCL(DG.t); U.minus];
  U.plus  = [U.plus; DG.ProbDef.BCR(DG.t)];
end
% Compute Lax-Friedrichs numerical flux
%----------------------------------------------------------------------
% C = max([abs(DG.ProbDef.c(U.minus)) abs(DG.ProbDef.c(U.plus))],[],2);
C = max([abs(DG.ProbDef.dfdu(U.minus)) abs(DG.ProbDef.dfdu(U.plus))],[],2);
fhat = 1/2*(DG.ProbDef.f(U.plus) + DG.ProbDef.f(U.minus) - C.*(U.plus-U.minus));
% % Compute Roe flux with entropy fix
% % ---------------------------------------------------------------------
% C = max([abs(DG.ProbDef.dfdu(U.minus)) abs(DG.ProbDef.dfdu(U.plus))],[],2);
% fhat = 1/2*(DG.ProbDef.f(U.plus) + DG.ProbDef.f(U.minus) - C.*(U.plus-U.minus));
% fL = DG.ProbDef.f(U.minus);
% fR = DG.ProbDef.f(U.plus);
% fLswitch = logical( (DG.ProbDef.dfdu(U.minus) >= 0) ...
%                  .* (DG.ProbDef.dfdu(U.plus) >= 0) );
% fRswitch = logical( (DG.ProbDef.dfdu(U.minus) < 0) ...
%                  .* (DG.ProbDef.dfdu(U.plus) < 0) );
% fhat(fLswitch) = fL(fLswitch);
% fhat(fRswitch) = fR(fRswitch);
% Compue RHS
%------------------------------------------------------------------
RHS = DG.A*DG.ProbDef.f(U.n) + DG.B*fhat;
