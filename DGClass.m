% Define a DG class
classdef DGClass < handle

  properties
    p       % DG polynomial degree
    Mesh    % Stored copy of mesh object
    ProbDef % Stored copy of ProbDef structure
    NDoFs   % Total number of degrees of freedom in discretized soln
    Uh      % Solution(s) in DG space (Can store multiple step values in (:,j>1))
    RHS     % RHS DG spatial evaluation(s) (Can store multiple step values in (:,j>1))
    M       % Global mass matrix
    A       % Global A matrix
    B       % Global B matrix
    phi     % Local basis function evaluations matrix
    PHI     % Global basis function evaluations matrix
    t       % Current time
  end

  methods
    function DGClass = DGClass(p, Mesh, ProbDef)
      % This function initializes the DG solution
      % p       : DG polynomial degree
      % Mesh    : a Mesh object
      % ProbDef : a ProbDef structure defining the problem to solve

      % Path to additional DG functions
      classpath = mfilename('fullpath')
      addpath([classpath(1:end-7), 'DGfunctions/'])

      DGClass.p = p;
      DGClass.Mesh = Mesh;
      DGClass.ProbDef = ProbDef;
      DGClass.NDoFs = Mesh.Nelems*(DGClass.p+1);

      % Get the quadrature points and weigths
      %--------------------------------------------------------------------------
      QRule = quadGaussJacobi(DGClass.p+1,0,0);
      W = sparse(1:length(QRule.Weights),1:length(QRule.Weights),QRule.Weights);

      for i = 0:DGClass.p
        % Construct the i-th basis function
        %----------------------------------------------------------------------
        phi{i+1} = polyJacobi(i,0,0);
        % Evaluate basis function at the boundary points of the master element
        %----------------------------------------------------------------------
        Phi.plus(i+1,1)  = polyval(phi{i+1},-1);
        Phi.minus(i+1,1) = polyval(phi{i+1}, 1);
        % Evaluate the basis function and derivative at quadrature points
        %----------------------------------------------------------------------
        Phi.n(:,i+1)  = polyval(phi{i+1},QRule.Points);
        dPhi.n(:,i+1) = polyval(polyder(phi{i+1}),QRule.Points);
      end
      Phi.n    = sparse(Phi.n);    dPhi.n    = sparse(dPhi.n);
      Phi.plus = sparse(Phi.plus); Phi.minus = sparse(Phi.minus);

      DGClass.phi = phi; % (This is only used for plotting)

      % Compute the mass matrix
      %--------------------------------------------------------------------------
      M = Phi.n'*W*Phi.n; M(abs(M)<=(eps*100)) = 0;
      DGClass.M = double(sym(M));

      % Compute the element matrics a, b+ and b-
      %--------------------------------------------------------------------------
      a       = inv(M)*dPhi.n'*W; a(abs(a)<=(eps*100)) = 0;
      b.plus  = inv(M)*Phi.plus;  b.plus  = double(sym(b.plus));
      b.minus = inv(M)*Phi.minus; b.minus = double(sym(b.minus));

      % Assemble global matrices
      %--------------------------------------------------------------------------
      DGClass.A = kron(diag(2./Mesh.dx),a);
      DGClass.B = [ kron(diag(2./Mesh.dx),b.plus) zeros([DGClass.NDoFs,1])   ] + ...
          [ zeros([DGClass.NDoFs,1]) kron(diag(2./Mesh.dx),-b.minus) ];

      DGClass.PHI.points = kron(eye(Mesh.Nelems),Phi.n);
      DGClass.PHI.plus   = kron(eye(Mesh.Nelems),Phi.plus);
      DGClass.PHI.minus  = kron(eye(Mesh.Nelems),Phi.minus);

      % Compute the L2 projection of the initial condition
      %------------------------------------------------------------------------
      DGClass.t = ProbDef.t0;
      DGClass.Uh = L2_projection_1D(DGClass);

    end

    function RKstep(DGClass, RK, dt)
      % ===================================================================
      % This function steps the DG solution forward in time one step using
      % the given Runge--Kutta method
      % ===================================================================
      % RK : A structure defining the RK method in Shu--Osher form
      %      RK.alpha - alpha coefficients vector
      %      RK.beta  - beta coefficients vector
      % dt : timestep - Can this live within DGClass or RK?
      % ===================================================================

      % Initialize matrices to store stage values
      y = zeros(DGClass.NDoFs, RK.s+1);
      RHS = zeros(DGClass.NDoFs, RK.s);

      % The first stage solution value is the current DG solution
      y(:,1) = DGClass.Uh;

      % Loop over stages
      for i = 2:RK.s+1
        % Compute right hand side of current stage (i-1)
        RHS(:,i-1) = DG_spatial_operator(y(:,i-1), DGClass);
        % Step to the next stage
        for j = 1:i-1
          y(:,i) = y(:,i) + ...
            RK.alpha(i-1,j)*y(:,j) + dt*RK.beta(i-1,j)*RHS(:,j);
          % Apply slope limiter
          if DGClass.ProbDef.limit_slope
            y(:,i) = slopeLimiter(y(:,i), DGClass.p, DGClass.Mesh.Nelems);
          end
        end
      end

      % Update DG Solution
      DGClass.Uh = y(:,end);
      % Update time
      DGClass.t = DGClass.t + dt;
    end

    function LMstep(DGClass, LM, dt)
      % ===================================================================
      % This function steps the DG solution forward in time one step using
      % the given linear multistep method
      % ===================================================================
      % LM : A structure defining the LM method
      %      LM.alpha - alpha coefficients vector
      %      LM.beta  - beta coefficients vector
      % dt : timestep - Can this live within DGClass or LM?
      % ===================================================================

      % % Store current DG solution in matrix of previous steps
      % DGClass.Uh_prev(:,i) = DGClass.Uh;
      % DGClass.Uh = zeros(NDoFs, 1)

      % Initialize new solution vector
      U = zeros(DGClass.NDoFs, 1);

      % Step to the next stage
      for i = 1:LM.r
        U = U + LM.alpha(i)*DGClass.Uh(:,i) + dt*LM.beta(i)*DGClass.RHS(:,i);
      end

      % Apply slope limiter
      if DGClass.ProbDef.limit_slope
        U = slopeLimiter(U, DGClass.p, DGClass.Mesh.Nelems);
      end

      % Update DG Solution, storing new solution in DGClass.Uh(:,1)
      for i = LM.r:-1:2
        DGClass.Uh(:,i) = DGClass.Uh(:,i-1);
        DGClass.RHS(:,i) = DGClass.RHS(:,i-1);
      end
      DGClass.Uh(:,1) = U;
      DGClass.RHS(:,1) = DG_spatial_operator(U, DGClass);

      % Update time
      DGClass.t = DGClass.t + dt;
    end

    function L2_Error = L2_error(DGClass)

      if ~isfield(DGClass.ProbDef, 'Ue')
        disp('Can''t compute L2 error - no exact solution ProbDef.Ue')
        return
      end

      % Retrieve a sufficient number of Gauss points
      %--------------------------------------------------------------------------
      QRule = quadGaussJacobi(2*DGClass.p+3,0,0);
      % Evaluate the basis functions at the Gauss points
      %--------------------------------------------------------------------------
      for i = 0:DGClass.p
          phi_l2(i+1,:) = polyval(DGClass.phi{i+1},QRule.Points);
      end
      PHI = kron(speye(DGClass.Mesh.Nelems),phi_l2');
      % Compute the DG solution at Gauss points
      %--------------------------------------------------------------------------
      Uh = PHI*DGClass.Uh(:,1);
      % Compute the quadrature matrix
      %--------------------------------------------------------------------------
      q = diag(QRule.Weights); Q = cell(1,DGClass.Mesh.Nelems); [Q{:}] = deal(q);
      Q = cellfun(@times,Q,num2cell(DGClass.Mesh.dx/2).','UniformOutput',0);
      Q = kron(speye(DGClass.Mesh.Nelems),q);
      % Compute the element psi vector
      %--------------------------------------------------------------------------
      psi(1).l2 = polyval([-1/2 1/2],QRule.Points); psi(2).l2 = polyval([1/2 1/2],QRule.Points);
      % Create the global X vector
      %--------------------------------------------------------------------------
      PSI.l2 = cell(1,DGClass.Mesh.Nelems); [PSI.l2{:}] = deal(([psi.l2]));
      X.elem = num2cell([DGClass.Mesh.Points(1:end-1),DGClass.Mesh.Points(2:end)]',1);
      X.l2 = cellfun(@mtimes,PSI.l2,[X.elem],'UniformOutput',0);
      X.l2 = reshape([X.l2{:}],numel([X.l2{:}]),1);
      % Determine U at all points X.l2
      %--------------------------------------------------------------------------
      UE = DGClass.ProbDef.Ue(X.l2,DGClass.t);

      % Compute the L2 error
      %--------------------------------------------------------------------------
      L2_Error = Q*( (UE - Uh).^2 );
      L2_Error = sqrt(sum(L2_Error));

    end

    function TV = TV_seminorm(DGClass)
      %--------------------------------------------------------------------------
      % Input: DGClass.Uh    = Vector of DG degrees of freedom
      %        DGClass.phi  = cell of basis functions
      %        DGClass.Mesh = mesh structure with fields Connectivity and Points
      %--------------------------------------------------------------------------

      % Determine the number of int points
      %--------------------------------------------------------------------------
      nq = length(DGClass.phi);
      % Determine the mesh sizes
      %--------------------------------------------------------------------------
      h = [DGClass.Mesh.Points(2:end)-DGClass.Mesh.Points(1:end-1)]';
      % Retrieve Gauss points
      %--------------------------------------------------------------------------
      QRule = quadGaussJacobi(nq,0,0);
      % Evaluate the basis functions at the Gauss points
      %--------------------------------------------------------------------------
      for i = 0:DGClass.p
          phi_tv(i+1,:) = polyval(DGClass.phi{i+1},QRule.Points);
      end
      PHI = kron(speye(DGClass.Mesh.Nelems),phi_tv);
      % Compute the DG solution at Gauss points
      %--------------------------------------------------------------------------
      Uh = PHI*DGClass.Uh(:,1);
      Uh = reshape(Uh, [length(QRule.Weights), DGClass.Mesh.Nelems]).';
      % Compute mean DG solution over each element
      % -------------------------------------------------------------------------
      Uint = Uh*QRule.Weights.*h.'./2;
      Umean = Uint./h.';

      % Compute total variation
      TV = sum(abs([Umean;Umean(1)] - [Umean(end); Umean]));

    end

    function plot(DGClass)

      % Construct the global basis for plotting
      %----------------------------------------
      xi = [ 2/DGClass.Mesh.dx(1) -1 ];
      for i = 0:DGClass.p
          m = length(DGClass.phi{i+1});
          phiplot(i+1,:) = [ zeros([1 DGClass.p+1-m]), polycompose(DGClass.phi{i+1},xi) ];
      end
      PHI.Plot = kron(speye(DGClass.Mesh.Nelems),phiplot);

      % Setup the figure and axes
      % -------------------------
      figure('Renderer','zbuffer');
      set(gca,'NextPlot','replaceChildren');
      set(gcf, 'Position', [100 100 995 452]);
      set(gca,'FontSize',13,'FontWeight','bold');
      axis([DGClass.Mesh.Points(1), DGClass.Mesh.Points(end), -1.7, 1.7])
      box on
      xplot = linspace(DGClass.ProbDef.xL,DGClass.ProbDef.xR,DGClass.Mesh.Nelems*10);
      hold on

      % Plot the exact solution if it exists
      % ------------------------------------
      if isfield(DGClass.ProbDef, 'Ue')
        exactline = plot(xplot,DGClass.ProbDef.Ue(xplot,DGClass.t),...
                      'Color',0.7*ones(1,3),...
                      'LineWidth',3,...
                      'DisplayName','Exact Solution');
        plot(xplot,DGClass.ProbDef.Ue(xplot,DGClass.t),...
                      'Color',0.5*ones(1,3),...
                      'LineWidth',1,...
                      'DisplayName','Exact Solution');
      end

      h1 = gca;

      % Plot the DG solution in MATLAB piecewise polynomial (pp) form
      % -------------------------------------------------------------
      breaks = DGClass.Mesh.Points(:,1)';
      coeffs = full(DGClass.Uh(:,1)'*PHI.Plot);
      DGpp = ppmak(breaks,coeffs);
      fnplt(DGpp,'jumps') % plot DG solution within elements
      ax = gca; DGline = ax.Children(1);
      set(DGline,'Color',[0/255, 24/255, 168/255],...
          'DisplayName',['DG Solution, \it p \rm = ',num2str(DGClass.p)]')

      xb = repmat(breaks,3,1);
      yb = repmat([200;-200;NaN],1,length(breaks));
      % Mark element boundaries
      plot(xb(:),yb(:),':','Color',0.8*ones(1,3))
      % Plot DG solution on right side of element boundary
      plot(breaks,fnval(DGpp,breaks),'bo',...
          'LineWidth',2,...
          'MarkerEdgeColor',[0/255, 24/255, 168/255],...
          'MarkerFaceColor',[.49 1 .63],...
          'MarkerSize',5);
      % Plot DG solution on left side of element boundary
      plot(breaks,fnval(DGpp,breaks,'l'),'bo',...
          'LineWidth',2,...
          'MarkerEdgeColor',[0/255, 24/255, 168/255],...
          'MarkerFaceColor',[.49 1 .63],...
          'MarkerSize',5)

      % Label axes, title, legend, and draw
      % -----------------------------------
      title(h1,['t = ',num2str(DGClass.t),', p = ',num2str(DGClass.p),')']);
      ylabel('u, u_h','FontSize',13,'FontWeight','bold')
      xlabel('x','FontSize',13,'FontWeight','bold')
      if isfield(DGClass.ProbDef, 'Ue')
        legend([exactline, DGline])
      else
        legend(DGline)
      end
      box on
      drawnow
    end

  end
end
