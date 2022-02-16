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
% Cockburn, Shu (1989) Example 1
% ------------------------------
ProbDef.name = 'CS_example_1';
ProbDef.xL = -1;
ProbDef.xR = 1;
ProbDef.BCtype = 'periodic';
ProbDef.f = @(u) (1/2) * u.^2;
ProbDef.U0 = @(x) 1/4 + 1/2 * sin(pi*x);
ProbDef.t0 = 0;
ProbDef.T = 2/pi;

% --------------------------
% Construct the mesh object
% --------------------------
Nelems = 3;
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


% -------------------------------------------------
% Step the solution forward in time
% -------------------------------------------------

% Calculate initial L2 norm
L2norm1 = DG.L2_norm()
% Calculate initial total variation in the means
TVsnorm1 = DG.TV_seminorm()

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
L2norm2 = DG.L2_norm()
TVsnorm2 = DG.TV_seminorm()

LogicStr = {'False', 'True'};
disp(['L2 stable? ',LogicStr{(L2norm2<L2norm1)+1}, ...
      ', L2 norm % growth = ',num2str( (L2norm2-L2norm1)/L2norm1*100 )])
disp(['TVDM? ',LogicStr{(TVsnorm2<TVsnorm1)+1}, ...
      ', TV seminorm % growth = ',num2str( (TVsnorm2-TVsnorm1)/TVsnorm1*100 )])
DG.plot()
