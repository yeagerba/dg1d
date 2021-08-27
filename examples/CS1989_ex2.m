clear all; close all;
% Directory containing DG discretization code
addpath('../')
% Add path to included timesteppers so we can use them
addpath('../timesteppers/')

% ==============================================================
% INPUT
% ==============================================================
M = 20; % tolerance for slope limiting (NaN for no limiting)
timestepper = 'RK'; % 'RK' or 'LM'
% ==============================================================

% --------------------------------------------------------------
% Problem definition structure
% --------------------------------------------------------------
% Cockburn, Shu (1989) Example 2
% ------------------------------
ProbDef.name = 'CS_example_2';
ProbDef.xL = -1;
ProbDef.xR = 1;
ProbDef.BCtype = 'dirichletfnt';
ProbDef.BCL = @(t) Ue_example2(ProbDef.xL,t);
ProbDef.BCR = @(t) Ue_example2(ProbDef.xR,t);
ProbDef.c = @(u) u/2;
ProbDef.f = @(u) ProbDef.c(u).*u;
ProbDef.dfdu = @(u) u
ProbDef.U0 = @(x) 1/4 + 1/2 * sin(pi*x);
ProbDef.t0 = 0;
ProbDef.T = 2/pi * 1.0;

% --------------------------
% Construct the mesh object
% --------------------------
Nelems = 80;
points(:,1) = linspace(ProbDef.xL,ProbDef.xR,Nelems+1);
connectivity(:,1) = [1:Nelems];
connectivity(:,2) = [2:Nelems+1];

% Initialize mesh object
Mesh = MeshClass(points, connectivity);

% ------------------------------------
% Initialize the solution in DG space
% ------------------------------------
p = 2; % DG polynomial degree
DG = DG(p, Mesh, ProbDef);

% Calculate initial total variation
TVsnorm1 = DG.TV_seminorm()

% -------------------------------------------------
% Step the solution forward in time
% -------------------------------------------------

if strcmp(timestepper, 'RK')
  % Define RK method
  RK.q = 3;      % Time stepper order
  RK.s = 4;      % Number of stages
  [RK.alpha, RK.beta, RK.cfl] = sspRKold(RK.q, RK.s);
  dt = DG.Mesh.dx_min*RK.cfl
  NT = DG.ProbDef.T/dt;

  % Time stepping
  for j = 1:NT
    DG.RKstep(RK, dt, M);
  end

elseif strcmp(timestepper, 'LM')
  % Define LM method
  LM.q = 3;      % Time stepper order
  LM.r = 4;      % Number of steps
  [LM.C_ssp, LM.alpha, LM.beta] = Rkp(LM.r, LM.q); % SSP-optimized LM methods
  LM.cfl = 0.052; % Theoretical cfl from Google Drive

  dt = DG.Mesh.dx_min*LM.cfl
  NT = DG.ProbDef.T/dt;

  % Initialize necessary DG solution variables
  for i = 1:LM.r
    DG.t = (1-i)*dt;
    DG.Uh(:,i) = L2_projection_1D(DG);
    DG.RHS(:,i) = DG_spatial_operator(DG.Uh(:,i), DG);
  end
  DG.t = 0; % Reset t to 0 before time stepping

  % Time stepping
  for j = 1:NT
    DG.LMstep(LM, dt, M);
  end
end

% ---------------
% Postprocessing
% ---------------
% L2error = DG.L2_error()
TVsnorm2 = DG.TV_seminorm()
DG.plot()


% ---------------------------------------------------------------
% Nested Function for Burgers' equation solution (for BCs)
% ---------------------------------------------------------------
options = optimoptions('fsolve', ...
                       'TolFun',1e-13,'TolX',1e-13, ...
                       'MaxIter', 5000, 'MaxFunEvals', 1e6, ...
                       'SpecifyObjectiveGradient', true, ...
                       'CheckGradients', true, ...
                       'Display', 'iter');

function UE = Ue_example2(x,t)
  Uefunc = @(u) u - ( 1/4 + 1/2 * sin(pi*(x-u*t)) );
  [UE, fval, exitflag, fsolveout] = fzero(Uefunc, sin(pi*(x)));
  if fval > 1e-12
    disp('WARNING: fval > 1e-12')
  end
  if exitflag ~= 1 & exitflag ~= 4
    disp(['WARNING: exitflag = ', string(exitflag)])
    disp(fsolveout)
    error('fsolve failed')
  end
end
