clear all; close all;
% Directory containing DG discretization code
addpath('../')
% Add path to included timesteppers so we can use them
addpath('../timesteppers/')

% ==============================================================
% INPUT
% ==============================================================
M = 20; % tolerance for slope limiting (NaN for no limiting)
timestepper = 'LM'; %'RK' or 'LM'
% ==============================================================

% --------------------------------------------------------------
% Problem definition structure
% --------------------------------------------------------------
% Cockburn, Shu (1989) Example 3
% ------------------------------
ProbDef.name = 'CS_example_3';
ProbDef.xL = -1;
ProbDef.xR = 1;
ProbDef.BCtype = 'dirichlet';
% Flux 1
ProbDef.BCL = -3;
ProbDef.BCR =  3;
ProbDef.f = @(u) (1/4) .* (u.^2 - 1) .* (u.^2 - 4);
ProbDef.U0 = @(x) ProbDef.BCL * (x<0) + ProbDef.BCR * (x>0);
ProbDef.t0 = 0;
ProbDef.T = 0.03; %*0.0077*1e-2;

% Flux 2 - Buckley-Leverett flux
ProbDef.BCL = 0;
ProbDef.BCR = 0;
ProbDef.f = @(u) ( 4*u.^2 ) ./ ( 4*u.^2 + (1-u).^2 );
ProbDef.U0 = @(x) (x > -0.5) .* (x < 0.0);
ProbDef.t0 = 0;
ProbDef.T = 0.4;

% --------------------------
% Construct the mesh object
% --------------------------
Nelems = 40;
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
  LM.cfl = 0.052*0.5; % Theoretical cfl from Google Drive

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
