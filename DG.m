% Define a DG class
classdef DG < handle

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
    function DG = DG(p, Mesh, ProbDef)
      % This function initializes the DG solution
      % p       : DG polynomial degree
      % Mesh    : a Mesh object
      % ProbDef : a ProbDef structure defining the problem to solve

      % Path to additional DG functions
      classpath = mfilename('fullpath');
      addpath([classpath(1:end-2), 'DGfunctions/']);

      DG.p = p;
      DG.Mesh = Mesh;
      DG.NDoFs = Mesh.Nelems*(DG.p+1);
      DG.ProbDef = ProbDef;
      % Automatically compute the derivative of the flux function
      syms u
      DG.ProbDef.dfdu = matlabFunction(diff(DG.ProbDef.f(u)), 'vars', [u]);

      % Get the quadrature points and weigths
      %--------------------------------------------------------------------------
      QRule = quadGaussJacobi(DG.p+1,0,0);
      W = sparse(1:length(QRule.Weights),1:length(QRule.Weights),QRule.Weights);

      for i = 0:DG.p
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

      DG.phi = phi; % (This is only used for plotting)

      % Compute the mass matrix
      %--------------------------------------------------------------------------
      M = Phi.n'*W*Phi.n; M(abs(M)<=(eps*100)) = 0;
      DG.M = double(sym(M));

      % Compute the element matrics a, b+ and b-
      %--------------------------------------------------------------------------
      a       = inv(M)*dPhi.n'*W; a(abs(a)<=(eps*100)) = 0;
      b.plus  = inv(M)*Phi.plus;  b.plus  = double(sym(b.plus));
      b.minus = inv(M)*Phi.minus; b.minus = double(sym(b.minus));

      % Assemble global matrices
      %--------------------------------------------------------------------------
      DG.A = kron(diag(2./Mesh.dx),a);
      DG.B = [ kron(diag(2./Mesh.dx),b.plus) zeros([DG.NDoFs,1])   ] + ...
          [ zeros([DG.NDoFs,1]) kron(diag(2./Mesh.dx),-b.minus) ];

      DG.PHI.points = kron(eye(Mesh.Nelems),Phi.n);
      DG.PHI.plus   = kron(eye(Mesh.Nelems),Phi.plus);
      DG.PHI.minus  = kron(eye(Mesh.Nelems),Phi.minus);

      % Compute the L2 projection of the initial condition
      %------------------------------------------------------------------------
      DG.t = ProbDef.t0;
      DG.Uh = L2_projection_1D(DG);

    end

    function RKstep(DG, RK, dt, M)
      arguments
        DG
        RK
        dt (1,1) double
        M  (1,1) double  = 20
      end
      % ===================================================================
      % This function steps the DG solution forward in time one step using
      % the given Runge--Kutta method
      % ===================================================================
      % RK : A structure defining the RK method in Shu--Osher form
      %      RK.alpha - alpha coefficients vector
      %      RK.beta  - beta coefficients vector
      % dt : timestep - Can this live within DG or RK?
      % ===================================================================

      % Initialize matrices to store stage values
      y = zeros(DG.NDoFs, RK.s+1);
      RHS = zeros(DG.NDoFs, RK.s);

      % The first stage solution value is the current DG solution
      y(:,1) = DG.Uh;

      % Loop over stages
      for i = 2:RK.s+1
        % Compute right hand side of current stage (i-1)
        RHS(:,i-1) = DG_spatial_operator(y(:,i-1), DG);
        % Step to the next stage
        for j = 1:i-1
          y(:,i) = y(:,i) + ...
            RK.alpha(i-1,j)*y(:,j) + dt*RK.beta(i-1,j)*RHS(:,j);
          % Apply slope limiter
          if ~isnan(M)
            m = M * DG.Mesh.dx.^2; % Should be row vector!
            % U = slopeLimiter(U, DG.p, DG.Mesh.Nelems, m.');
            y(:,i) = DG.slopeLimiter(y(:,i), m.');
          end
        end
      end

      % Update DG Solution
      DG.Uh = y(:,end);
      % Update time
      DG.t = DG.t + dt;
    end

    function LMstep(DG, LM, dt, M)
      arguments
        DG
        LM
        dt (1,1) double
        M  (1,1) double  = 20
      end
      % ===================================================================
      % This function steps the DG solution forward in time one step using
      % the given linear multistep method
      % ===================================================================
      % LM : A structure defining the LM method
      %      LM.alpha - alpha coefficients vector
      %      LM.beta  - beta coefficients vector
      % dt : timestep - Can this live within DG or LM?
      % M  : slope limiter tolerance (NaN for no limiting - default)
      % ===================================================================

      % Initialize new solution vector
      U = zeros(DG.NDoFs, 1);

      % Step to the next stage
      for i = 1:LM.r
        U = U + LM.alpha(i)*DG.Uh(:,i) + dt*LM.beta(i)*DG.RHS(:,i);
      end

      % Apply slope limiter
      if ~isnan(M)
        m = M * DG.Mesh.dx.^2; % Should be row vector!
        U = DG.slopeLimiter(U, m.');
      end

      % Update DG Solution, storing new solution in DG.Uh(:,1)
      for i = LM.r:-1:2
        DG.Uh(:,i) = DG.Uh(:,i-1);
        DG.RHS(:,i) = DG.RHS(:,i-1);
      end
      DG.Uh(:,1) = U;
      DG.RHS(:,1) = DG_spatial_operator(U, DG);

      % Update time
      DG.t = DG.t + dt;
    end

    function ulim = slopeLimiter(DG, u, m)
      %% SLOPELIMITER Applies the minmod slope limiter
      %    DoFs = slopeLimiter(u,N) applies the minmod slope limiter to the
      %    degrees of freedom, u, of a DG solution of degree p for N elements.
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

      %% Apply the slope limiter

      % Reshape the DoFs into a (p+1) by N array
      u = reshape(u(:,1),DG.p+1,DG.Mesh.Nelems);
      % Store the mean value of the DG solution over each element
      uBar = u(1,:);
      % Compute the DG solution at the end points of the element
      uPlus  = (-ones(1,DG.p+1)).^(0:DG.p)*u;
      uMinus = ones(1,DG.p+1)*u;
      % Compute deviation from the mean at element end points
      enddev = [uMinus - uBar; uBar - uPlus];
      if strcmp(DG.ProbDef.BCtype, 'periodic')
        % Compute mean difference between neighboring elements
        meandiffL = uBar - [uBar(end), uBar(1:end-1)];
        meandiffR = [uBar(2:end), uBar(1)] - uBar;
        % Apply the minmod function
        [minmodMinus,iMinus] = minmod([enddev(1,:); meandiffL; meandiffR]);
        [minmodPlus,iPlus]   = minmod([enddev(2,:); meandiffL; meandiffR]);
      elseif strcmp(DG.ProbDef.BCtype, 'dirichlet')
        % Compute mean difference between neighboring elements
        meandiffL = uBar - [DG.ProbDef.BCL, uBar(1:end-1)];
        meandiffL(1) = 2*meandiffL(1);
        meandiffR = [uBar(2:end), DG.ProbDef.BCR] - uBar;
        meandiffR(end) = sign(enddev(1,1));    % meandiffR(end) = 2*meandiffR(end);
        % Apply the minmod function
        [minmodMinus,iMinus] = minmod([enddev(1,:); meandiffL; meandiffR]);

        meandiffL = uBar - [DG.ProbDef.BCL, uBar(1:end-1)];
        meandiffL(1) = sign(enddev(2,1));
        meandiffR = [uBar(2:end), DG.ProbDef.BCR] - uBar;
        meandiffR(end) = 2*meandiffR(end);
        % Apply the minmod function
        [minmodPlus,iPlus] = minmod([enddev(2,:); meandiffL; meandiffR]);
      elseif strcmp(DG.ProbDef.BCtype, 'dirichletfnt')
        % Compute mean difference between neighboring elements
        meandiffL = uBar - [DG.ProbDef.BCL(DG.t), uBar(1:end-1)];
        meandiffL(1) = 2*meandiffL(1);
        meandiffR = [uBar(2:end), DG.ProbDef.BCR(DG.t)] - uBar;
        meandiffR(end) = sign(enddev(1,1));    % meandiffR(end) = 2*meandiffR(end);
        % Apply the minmod function
        [minmodMinus,iMinus] = minmod([enddev(1,:); meandiffL; meandiffR]);

        meandiffL = uBar - [DG.ProbDef.BCL(DG.t), uBar(1:end-1)];
        meandiffL(1) = sign(enddev(2,1));
        meandiffR = [uBar(2:end), DG.ProbDef.BCR(DG.t)] - uBar;
        meandiffR(end) = 2*meandiffR(end);
        % Apply the minmod function
        [minmodPlus,iPlus] = minmod([enddev(2,:); meandiffL; meandiffR]);
      end

      % Identify elements where |enddev| <= m (see page 195 of [1])
      iMinus(abs(enddev(1,:))<=m) = 1;
      iPlus(abs(enddev(2,:))<=m) = 1;
      % Find elements to which slope limiter is applied
      iMod = or(iMinus~=1,iPlus~=1);
      % Compute the modified endpoints (see Eqs 2.10 and 2.11 on page 193 of [1])
      uMinusMod(iMod) = uBar(iMod) + minmodMinus(iMod);
      uPlusMod(iMod)  = uBar(iMod) - minmodPlus(iMod);
      % Compute the new slope of the element
      u(2,iMod) = (uMinusMod(iMod) - uPlusMod(iMod))/2;
      % If p > 1, then...
      if DG.p > 1
          % Zero out higher-order DoFs
          u(3:end,iMod) = 0;
      end

      % Save the modified DoFs
      ulim(:,1) = u(:);

      % ====================================
      % NESTED FUNCTIONS: minmod
      % ====================================
      function [m,I] = minmod(A)

        %% MINMOD The minmod function
        %    m = minmod(A) applies the so-called minmod function to an array.
        %
        %     � If A is a vector of length n, then minmod(A) returns s*min(abs(A)),
        %       if s = sign(A(1)) = ... = sign(A(n)), otherwise it returns 0.
        %
        %       Examples: A = [  1  2  3 ]; minmod(A) =  1
        %                 B = [ -2 -3 -4 ]; minmod(B) = -2
        %                 C = [ -1  2  3 ]; minmod(C) =  0
        %
        %     � If A is a matrix, then minmod(A) returns a row vector m where
        %       m(i) = minmod(A(:,i)), i.e., it applies the minmod function as
        %       described above to each column of A.
        %
        %       Example: A = [ 1 -2 -1; 2 -1 2; 3 -4 3 ]; minmod(A) = [ 1 -1 0 ]
        %
        %    [m,I] = minmod(A) applies the minmod function as described above and
        %    also returns the the index into the operating dimension that corre-
        %    sponds to the minmod value of A or 0.
        %
        %       Example: A = [ 1 -2 -1; 2 -1 2; 3 -4 3 ]; [m,I] = minmod(A) returns
        %       m = [ 1 -1 0 ], I = [ 1 2 0 ]
        %

        %% Validate input

        ip = inputParser;
        vA = @(x)validateattributes(x,{'numeric'},{'2d'});
        ip.addRequired('A',vA);
        ip.parse(A);
        ip.Results;

        %% Apply the minmod function

        [M,I] = min(abs(A));
        switch min(size(A))
            case 1
                m = isequal(abs(sum(sign(A))),length(A))*sign(A(1))*M;
            otherwise
                m = (abs(sum(sign(A))) == size(A,1)).*sign(A(1,:)).*M;
        end
        I = I.*logical(m);

      end

    end

    function L2_Error = L2_error(DG)

      if ~isfield(DG.ProbDef, 'Ue')
        disp('Can''t compute L2 error - no exact solution ProbDef.Ue')
        return
      end

      % Retrieve a sufficient number of Gauss points
      %--------------------------------------------------------------------------
      QRule = quadGaussJacobi(2*DG.p+3,0,0);
      % Evaluate the basis functions at the Gauss points
      %--------------------------------------------------------------------------
      for i = 0:DG.p
          phi_l2(i+1,:) = polyval(DG.phi{i+1},QRule.Points);
      end
      PHI = kron(speye(DG.Mesh.Nelems),phi_l2');
      % Compute the DG solution at Gauss points
      %--------------------------------------------------------------------------
      Uh = PHI*DG.Uh(:,1);
      % Compute the quadrature matrix
      %--------------------------------------------------------------------------
      q = diag(QRule.Weights); Q = cell(1,DG.Mesh.Nelems); [Q{:}] = deal(q);
      Q = cellfun(@times,Q,num2cell(DG.Mesh.dx/2).','UniformOutput',0);
      Q = kron(speye(DG.Mesh.Nelems),q);
      % Compute the element psi vector
      %--------------------------------------------------------------------------
      psi(1).l2 = polyval([-1/2 1/2],QRule.Points); psi(2).l2 = polyval([1/2 1/2],QRule.Points);
      % Create the global X vector
      %--------------------------------------------------------------------------
      PSI.l2 = cell(1,DG.Mesh.Nelems); [PSI.l2{:}] = deal(([psi.l2]));
      X.elem = num2cell([DG.Mesh.Points(1:end-1),DG.Mesh.Points(2:end)]',1);
      X.l2 = cellfun(@mtimes,PSI.l2,[X.elem],'UniformOutput',0);
      X.l2 = reshape([X.l2{:}],numel([X.l2{:}]),1);
      % Determine U exact at all points X.l2
      %--------------------------------------------------------------------------
      UE = DG.ProbDef.Ue(X.l2,DG.t);

      % Compute the L2 error
      %--------------------------------------------------------------------------
      L2_Error = Q*( (UE - Uh).^2 );
      L2_Error = sqrt(sum(L2_Error));

    end

    function L2_Norm = L2_norm(DG)

      % Retrieve a sufficient number of Gauss points
      %--------------------------------------------------------------------------
      QRule = quadGaussJacobi(2*DG.p+3,0,0);
      % Evaluate the basis functions at the Gauss points
      %--------------------------------------------------------------------------
      for i = 0:DG.p
          phi_l2(i+1,:) = polyval(DG.phi{i+1},QRule.Points);
      end
      PHI = kron(speye(DG.Mesh.Nelems),phi_l2');
      % Compute the DG solution at Gauss points
      %--------------------------------------------------------------------------
      Uh = PHI*DG.Uh(:,1);
      % Compute the quadrature matrix
      %--------------------------------------------------------------------------
      q = diag(QRule.Weights); Q = cell(1,DG.Mesh.Nelems); [Q{:}] = deal(q);
      Q = cellfun(@times,Q,num2cell(DG.Mesh.dx/2).','UniformOutput',0);
      Q = kron(speye(DG.Mesh.Nelems),q);
      % Compute the element psi vector
      %--------------------------------------------------------------------------
      psi(1).l2 = polyval([-1/2 1/2],QRule.Points); psi(2).l2 = polyval([1/2 1/2],QRule.Points);
      % Create the global X vector
      %--------------------------------------------------------------------------
      PSI.l2 = cell(1,DG.Mesh.Nelems); [PSI.l2{:}] = deal(([psi.l2]));
      X.elem = num2cell([DG.Mesh.Points(1:end-1),DG.Mesh.Points(2:end)]',1);
      X.l2 = cellfun(@mtimes,PSI.l2,[X.elem],'UniformOutput',0);
      X.l2 = reshape([X.l2{:}],numel([X.l2{:}]),1);

      % Compute the L2 error
      %--------------------------------------------------------------------------
      L2_Norm = sqrt(sum(Q*Uh.^2));

    end

    function TV = TV_seminorm(DG)
      %--------------------------------------------------------------------------
      % Input: DG.Uh = Vector of DG degrees of freedom
      %        DG.p  = basis function polynomial degree
      %--------------------------------------------------------------------------

      % Compute mean DG solution over each element (just the first DOF)
      Umean = DG.Uh(1:DG.p+1:end,1);

      % Compute total variation
      TV = sum(abs([Umean;Umean(1)] - [Umean(end); Umean]));

    end

    function plot(DG)

      % Construct the global basis for plotting
      %----------------------------------------
      xi = [ 2/DG.Mesh.dx(1) -1 ];
      for i = 0:DG.p
          m = length(DG.phi{i+1});
          phiplot(i+1,:) = [ zeros([1 DG.p+1-m]), polycompose(DG.phi{i+1},xi) ];
      end
      PHI.Plot = kron(speye(DG.Mesh.Nelems),phiplot);

      % Setup the figure and axes
      % -------------------------
      figure('Renderer','zbuffer');
      set(gca,'NextPlot','replaceChildren');
      set(gcf, 'Position', [100 100 995 452]);
      set(gca,'FontSize',13,'FontWeight','bold');
      axis([DG.Mesh.Points(1), DG.Mesh.Points(end), -2.1, 2.1])
      box on
      xplot = linspace(DG.ProbDef.xL,DG.ProbDef.xR,DG.Mesh.Nelems*10);
      hold on

      % Plot the exact solution if it exists
      % ------------------------------------
      if isfield(DG.ProbDef, 'Ue')
        exactline = plot(xplot,DG.ProbDef.Ue(xplot,DG.t),...
                      'Color',[1, 117/255, 24/255],...
                      'LineWidth',4,...
                      'DisplayName','Exact Solution');
        plot(xplot,DG.ProbDef.Ue(xplot,DG.t),...
                      'Color',[1, 117/255, 24/255],...
                      'LineWidth',4,...
                      'DisplayName','Exact Solution');
      end

      h1 = gca;

      % Plot the DG solution in MATLAB piecewise polynomial (pp) form
      % -------------------------------------------------------------
      breaks = DG.Mesh.Points(:,1)';
      coeffs = full(DG.Uh(:,1)'*PHI.Plot);
      DGpp = ppmak(breaks,coeffs);
      fnplt(DGpp,'jumps') % plot DG solution within elements
      ax = gca; DGline = ax.Children(1);
      set(DGline,'Color',[0/255, 24/255, 168/255],...
          'DisplayName',['DG Solution, \it p \rm = ',num2str(DG.p)]')

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
      title(h1,['t = ',num2str(DG.t),', p = ',num2str(DG.p),')']);
      ylabel('u, u_h','FontSize',13,'FontWeight','bold')
      xlabel('x','FontSize',13,'FontWeight','bold')
      if isfield(DG.ProbDef, 'Ue')
        legend([exactline, DGline])
      else
        legend(DGline)
      end
      box on
      drawnow
    end

  end
end
