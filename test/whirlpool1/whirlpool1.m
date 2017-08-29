
function whirlpool1 (mode, numGoals, numAbs, controllers)
  w = [0 0];
  bx = [3.45 3.45];
  a1x = [3.5 -0.9];
  bInd = [1 2];
  aInd = [3 4];
  x = [bx a1x];
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
      O = SymbolicSet('plotting/O.bdd', 'projection', [1 2]);
      plotCells(O, 'facecolor', colors(1,:)*0.5+0.5, 'edgec', colors(1,:), 'linew', 0.1)
      drawnow
      disp('Done plotting obstacles')
    catch
      warning('No obstacles');
    end
    savefig('problem');
  end
 
  if (strcmp(mode, 'GB')) % generalized Buchi
    numLoops = 100;
  
    openfig('problem');
    hold on
    drawnow
    pause;
    
    v = [];
    
    iStep = 1;
    
    N = SymbolicSet(['G/G1' int2str(numAbs) '.bdd']);
    
    for iLoop = 1:numLoops
      disp(['loop: ' int2str(iLoop)])
      for iGoal = 1:numGoals
	disp(['goal: ' int2str(iGoal)])
	numControllers = controllers(iGoal);
	
	baxdim = size(bx,2) + size(aInd(iGoal,:),2);
	
	hitGoal = 0;
	
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
	  plot(p(:,1), p(:,2), 'x', 'MarkerFaceColor', colors(mod(iStep,7)+1,:), 'MarkerSize', 1.5);
	  eta = Z.eta();
	  eta = eta';
	  tau = getTau(eta);
	  disp(['eta: ' num2str(eta)])
	  disp(['tau: ' num2str(tau)]) 
	 
	  while(1)
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
	      if (iController == numControllers)
		if (hitGoal == 0)
		  hitGoal = 1;
		else
		  boolplot = 1;
		  boolbreak = 1;
		end
	      else
		boolplot = 1;
		boolbreak = 1;
	      end
	    end
	    
	    if (boolplot)
	      plot(x(:,bInd(1)),x(:,bInd(2)),'k.-')
	      H = gobjects(numGoals, 1);
	      for jGoal = 1:numGoals
		H(jGoal) = plot(0, x(end,aInd(jGoal,1)), 'o', 'Color', colors(mod(iStep,7)+1,:), 'MarkerSize', 20, 'LineWidth', 3);
	      end
	      drawnow
	      pause(tau/3)
	      
	      if (~(iLoop == numLoops && iGoal == numGoals && iController == 1))
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
	      [t aphi] = ode45(str2func(['a' int2str(jGoal) 'ODE']), [0 tau], x(end,aInd(jGoal,:)), []);
	      aphi = aphi(end,:);	      
	      phi = [phi repmat(aphi,size(phi,1),1)];
	    end 
	    x = [x; phi];

	    iStep = iStep + 1;
	  
	  end   
	end
      end
    end

    savefig('gBuchi');
  end
end

function tau = getTau(eta)
  tau = eta(1) * 3 / 2;
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

function da1xdt = a1ODE(t,a1x)
  da1xdt = zeros(size(a1x));
end