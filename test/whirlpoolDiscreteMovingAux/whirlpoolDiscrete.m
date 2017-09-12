% Mode can be either 'S', 'P' or 'GB'
% 'S' = load and draw state space
% 'P' = load and draw obstacles
% 'GB' = Generalized Buchi
% numGoals = number of goals
% numAbs = number of abstraction layers
% controllers = array of length=numGoals, where controllers[i] = number of
% controllers for i-th goal
function whirlpoolDiscrete (mode, numGoals, numAbs, controllers)
  w = [0.05 0.05];
  bx = [0 0];
  a1x = [2.5 4];
  a2x = [4 -0.75];
  a3x = [-2.75 -4];
  bInd = [1 2];
  aInd = [3 4; 5 6; 7 8];
  x = [bx a1x a2x a3x];
   
  addpath(genpath('../..'));
  
  % colors
  colors=get(groot,'DefaultAxesColorOrder');
  
  if (strcmp(mode, 'S'))
    figure
    hold on
    box on    
    drawnow
  
    % load and draw state space
    X = SymbolicSet('plotting/X.bdd');
    lb = X.first();
    ub = X.last();   
    axis([lb(1)-1 ub(1)+1 lb(2)-1 ub(2)+1])
    plotCells(X, 'facecolor', 'none', 'edgec', [0.8 0.8 0.8], 'linew', 0.1);
    drawnow
    disp('Done plotting state space')
   
    savefig('system');
  end
  
  if (strcmp(mode,'P'))
    openfig('system');
    hold on
    drawnow
    
    % load and draw obstacles
    try
      O = SymbolicSet('plotting/O.bdd');
      plotCells(O, 'facecolor', colors(1,:)*0.5+0.5, 'edgec', colors(1,:), 'linew', 0.1)
      drawnow
      disp('Done plotting obstacles')
    catch
      warning('No obstacles');
    end
    savefig('problem');
  end
 
  if (strcmp(mode, 'GB')) % generalized Buchi
    numLoops = 50000;
  
    openfig('problem');
    hold on
    drawnow
    
    v = [];
    
    iStep = 1;
    
    N = SymbolicSet(['G/G1' int2str(numAbs) '.bdd']);
    
    for iLoop = 1:numLoops
      disp(['loop: ' int2str(iLoop)])
      for iGoal = 1:numGoals
	disp(['goal: ' int2str(iGoal)])
	numControllers = controllers(iGoal);
	
	baxdim = size(bx,2) + size(aInd(iGoal,:),2);
	
	
	for iController = numControllers:-1:1
	  disp(['controller: ' int2str(iController)])
	  C = SymbolicSet(['C/C' int2str(iGoal) int2str(iController) '.bdd'], 'projection', 1:baxdim);
	  Z = SymbolicSet(['Z/Z' int2str(iGoal) int2str(iController) '.bdd']);
	  if (iController == 1)
	    G = SymbolicSet(['G/G' int2str(iGoal) int2str(numAbs) '.bdd']);
	  else
	    G = SymbolicSet(['Z/Z' int2str(iGoal) int2str(iController-1) '.bdd']);	   
	  end
	  p = G.points();
	  D = gobjects(1,1);
	  D(1) = plot(p(:,1), p(:,2), 'x', 'MarkerFaceColor', colors(mod(iStep,7)+1,:), 'MarkerSize', 1.5);
	  eta = Z.eta();
	  eta = eta';
	  tau = getTau(eta);
	  disp(['eta: ' num2str(eta)])
	  disp(['tau: ' num2str(tau)])
	 
	  while(1)
	    disp('loop')
	    disp(iLoop)
	    disp('goal')
	    disp(iGoal)
	  
	    boolplot = 0;
	    boolbreak = 0;
	  
	    bx = x(end,bInd);
	    ax = x(end,aInd(iGoal,:));
	    bax = [bx ax];
	  
	    disp(['step: ' int2str(iStep)])
	    disp(['base x: ' num2str(bx)])
	    disp(['aux' int2str(iGoal) ' x: ' num2str(ax)])
	  
	    if (mod(iStep,1) == 0)
	      boolplot = 1;
	    end
	
	    if (G.isElement(bax))
	      boolplot = 1;
	      boolbreak = 1;

	    end
	    
	    if (boolplot)
	      P = plot(x(end,bInd(1)),x(end,bInd(2)), 'ko', 'MarkerSize', 5, 'LineWidth', 2);
	      H = gobjects(numGoals, 1);
	      for jGoal = 1:numGoals
		if (jGoal == iGoal)
		  color = colors(2,:);
		else
		  color = colors(6,:);
		end
		H(jGoal) = plot(x(end,aInd(jGoal,1)), x(end,aInd(jGoal,2)), 'o', 'Color', color, 'MarkerSize', 30, 'LineWidth', 3);
	      end
	      drawnow
	      pause(tau/2);
	      
	      if (~(iLoop == numLoops && iGoal == numGoals && iController == 1 && boolbreak))
		delete(P);
		for jGoal = 1:numGoals
		  delete(H(jGoal));
		end
	      end
	      
	    end
	    if (boolbreak)
	      break;
	    end
	  
	    u = C.getInputs(bax);

	    ran = randi([1 size(u,1)], 1, 1);	
	    v = [v; u(ran,:)];
	    d = disturbance(w);
	    
	    disp(['u: ' num2str(u(ran,:))])
	    disp(['d: ' num2str(d)])
	    
	    
	    [t bphi] = ode45(@bODE, [0 tau], bx, [], u(ran,:), d);
	    phi = bphi;	    
	    
	    for jGoal = 1:numGoals
%  	      [t aphi] = ode45(str2func(['a' int2str(jGoal) 'ODE']), [0 tau], x(end,aInd(jGoal,:)), []);
	      
	      if (jGoal == 1)
		aphi = a1next(x(end,aInd(jGoal,:)), tau);
	      elseif (jGoal == 2)
		aphi = a2next(x(end,aInd(jGoal,:)), tau);
	      elseif (jGoal == 3)
		aphi = a3next(x(end,aInd(jGoal,:)), tau);
	      end	      
	      phi = [phi repmat(aphi,size(phi,1),1)];
	    end 
	    x = [x; phi];

	    iStep = iStep + 1;
	  
	  end
	  delete(D(1));
	end
      end
    end

    savefig('gBuchi');
  end
end

function d = disturbance(w)
  d = -w + (2 * w .* rand(size(w)));
end

function dbxdt = bODE(t,bx,u,d)
  dbxdt = zeros(size(bx));
  dbxdt(1) = (-0.5 * bx(1) + 1 * bx(2)) * u(2);
  dbxdt(2) = (-0.5 * bx(1) + 0.5 * bx(2) + u(1)) * u(2);
  dbxdt = dbxdt + d';
end

function a1x = a1next(a1x, tau)
  a1x(1) = a1x(1) + 1 * tau/0.6;
  if (a1x(1) > 2.75)
    a1x(1) = -2.75;
  end
end

function a2x = a2next(a2x, tau)
  a2x(2) = a2x(2) - 1 * tau/0.6;
  if (a2x(2) < -2.75)
    a2x(2) = 2.75;
  end
end

function a3x = a3next(a3x, tau)
  a3x(1) = a3x(1) - 1 * tau/0.6;
  if (a3x(1) < -2.75)
    a3x(1) = 2.75;
  end
end

function tau = getTau(eta)
  tau = eta(1);
end