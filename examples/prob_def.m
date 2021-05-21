% Directory containing DG discretization code
addpath('../')
% Add path to included timesteppers so we can use them
addpath('../timesteppers/')

% --------------------------------------------------------------
% Problem definition structure
% --------------------------------------------------------------

% Gaussian distribution propogating to the right
% ----------------------------------------------
% ProbDef.name = 'linear_gaussian';  % Problem name
% ProbDef.xL = -1;                   % Left endpoint
% ProbDef.xR = 1;                    % Right endpoint
% ProbDef.BCtype = 'periodic';       % Boundary conditions
% ProbDef.c  = @(u) 1;               % Wave propagation speed (linear advection)
% ProbDef.f  = @(u) ProbDef.c(u).*u; % Flux function (linear advection)
% ProbDef.Ue = @(x,t) 1*exp( -((mod(x+1-t,2)-1)/0.25).^2 ); % Exact solution (if known)
% ProbDef.U0 = @(x) ProbDef.Ue(x,0); % Exact solution (if known)
% ProbDef.t0 = 0;                    % Initial time (set to 0 if not specified)
% ProbDef.T = 9;                     % Final time

% Cockburn, Shu (1989) Example 1
% ------------------------------
ProbDef.name = 'Burgers''';
ProbDef.xL = -1;
ProbDef.xR = 1;
ProbDef.BCtype = 'periodic';
ProbDef.c = @(u) u/2;
ProbDef.f = @(u) ProbDef.c(u).*u;
ProbDef.U0 = @(x) 1/4 + 1/2 * sin(pi*x);
ProbDef.t0 = 0;
ProbDef.T = 2/pi * 0.5;


% --------------------------
% Construct the mesh object
% --------------------------
Nelems = 10;
points(:,1) = linspace(ProbDef.xL,ProbDef.xR,Nelems+1);
connectivity(:,1) = [1:Nelems];
connectivity(:,2) = [2:Nelems+1];

% Initialize mesh object
Mesh = MeshClass(points, connectivity);


% ------------------------------------
% Initialize the solution in DG space
% ------------------------------------

p = 2; % DG polynomial degree
DG = DGClass(p, Mesh, ProbDef);

TVsnorm1 = DG.TV_seminorm()


% -------------------------------------------------
% Step the solution forward in time
% -------------------------------------------------

% RK example
% ------------
% % Define RK method
% RK.q = 3;      % Time stepper order
% RK.s = 4;      % Number of stages
% [RK.alpha, RK.beta, RK.cfl] = sspRKold(RK.q, RK.s);
% dt = DG.Mesh.dx_min*RK.cfl;
% NT = DG.ProbDef.T/dt;
%
% % Time stepping
% for j = 1:NT
%   DG.RKstep(RK, dt);
% end

% LM example
% -----------
% Define LM method
LM.q = 3;      % Time stepper order
LM.r = 4;      % Number of steps
[LM.C_ssp, LM.alpha, LM.beta] = Rkp(LM.r, LM.q); % SSP-optimized LM methods

dt = 0.01; % Choose a small number since we don't know the cfl
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
  DG.LMstep(LM, dt);
end

% ---------------
% Postprocessing
% ---------------
% L2error = DG.L2_error()
TVsnorm2 = DG.TV_seminorm()
DG.plot()
