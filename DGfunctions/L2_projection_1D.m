function [U_L2] = L2_projection_1D(DG)

% Retrieve a sufficient number of Gauss points
%--------------------------------------------------------------------------
QRule = quadGaussJacobi(2*DG.p+3,0,0);
% Evaluate the basis functions at the Gauss points
%--------------------------------------------------------------------------
for i = 0:DG.p
    phi_l2(i+1).f = polyval(DG.phi{i+1},QRule.Points);
end
% Compute the element L2 projection matrix
%--------------------------------------------------------------------------
m = [phi_l2.f]'*diag(QRule.Weights)*[phi_l2.f]; q = [phi_l2.f]'*diag(QRule.Weights); l2 = m\q;
% Compute the element psi vector
%--------------------------------------------------------------------------
psi(1).l2 = polyval([-1/2 1/2],QRule.Points); psi(2).l2 = polyval([1/2 1/2],QRule.Points);
% Create a global L2 matrix
%--------------------------------------------------------------------------
% L2 = cell(1,DG.Mesh.Nelems); [L2{:}] = deal(l2); L2 = spblkdiag(L2);
L2 = kron(speye(DG.Mesh.Nelems),l2);
% Create the global X vector
%--------------------------------------------------------------------------
PSI.l2 = cell(1,DG.Mesh.Nelems); [PSI.l2{:}] = deal(([psi.l2]));
X.elem = num2cell([DG.Mesh.Points(1:end-1),DG.Mesh.Points(2:end)]',1);
X.l2 = cellfun(@mtimes,PSI.l2,[X.elem],'UniformOutput',0);
X.l2 = reshape([X.l2{:}],numel([X.l2{:}]),1);
% Take the L2 projection of U
%--------------------------------------------------------------------------
if DG.t == 0
  U_L2 = L2*DG.ProbDef.U0(X.l2);
else
  U_L2 = L2*DG.ProbDef.Ue(X.l2, DG.t);
end

end
