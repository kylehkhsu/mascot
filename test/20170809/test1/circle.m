
function circle (mode, controllers)

  addpath(genpath('/home/kylehsu/control/SCOTS+Adaptive'));
  
  % colors
  colors=get(groot,'DefaultAxesColorOrder');  
 
  if (mode == 'S')
    figure
    hold on
    box on    
    drawnow
  
    % load and draw state space
    X = SymbolicSet('Xgrid.bdd');
    lb = X.first();
    ub = X.last();   
    axis([lb(1)-1 ub(1)+1 lb(2)-1 ub(2)+1])
    plotCells(X, 'facecolor', 'none', 'edgec', [0.8 0.8 0.8], 'linew', 0.1);
    drawnow
    disp('Done plotting state space')
   
    savefig('system');
  end
  
  if (mode == 'P')
    openfig('system');
    hold on
    drawnow
    
    % load and draw obstacles
    try
      O = SymbolicSet('O.bdd');
      plotCells(O, 'facecolor', colors(1,:)*0.5+0.5, 'edgec', colors(1,:), 'linew', 0.1)
      drawnow
      disp('Done plotting obstacles')
    catch
      warning('No obstacles');
    end
        
    % load and draw goal
    G = SymbolicSet('G.bdd');
    plotCells(G, 'facecolor', colors(2,:)*0.5+0.5, 'edgec', colors(2,:), 'linew', 0.1)
    drawnow
    disp('Done plotting goal')
    savefig('problem');

  end
  
  if (mode == 'T')
    openfig('system');
    hold on
    drawnow
    
    Zf = SymbolicSet('test/Zf.bdd');
    Zc = SymbolicSet('test/Zc.bdd');
    plotCells(Zc, 'facecolor', colors(2,:)*0.5+0.5, 'edgec', colors(2,:), 'linew', 0.1)
    drawnow
    pause
    plotCells(Zf, 'facecolor', colors(3,:)*0.5+0.5, 'edgec', colors(3,:), 'linew', 0.1)

  end
  
  if (mode == 'A')
    loops = 5;
  
    openfig('problem');
    hold on
    drawnow
    
    I = SymbolicSet('I.bdd');
    plotCells(I, 'facecolor', colors(3,:)*0.5+0.5, 'edgec', colors(3,:), 'linew', 0.1)
    drawnow
    disp('Done plotting initial')
    
    x = I.points();
    x = x(1,:);
    
    v = [];
    
    j = 1;
    k = 1
    
    while(1)
      if (k > loops)
	break;
      end
      for i = controllers-1:-1:0
	if (i == 0)
	  i = controllers;
	end      
      
	disp(['controller: ' int2str(i)])
	
	C = SymbolicSet(['C/C' int2str(i) '.bdd'], 'projection', [1 2]);
	if (i == 1)
	  G = SymbolicSet(['G.bdd']);
	elseif (i == controllers)
	  G = SymbolicSet(['E.bdd']);
	else
	  G = SymbolicSet(['Z/Z' int2str(i-1) '.bdd']);
	end
	
	p = G.points();
	plot(p(:,1), p(:,2), 'x', 'MarkerFaceColor', colors(mod(i,7)+1,:), 'MarkerSize', 1.5);
	
	Z = SymbolicSet(['Z/Z' int2str(i) '.bdd']);
	eta = Z.eta();
	eta = eta';
	tau = eta(1) * 3 / 2;
	
	disp('eta')
	disp(eta)
	disp('tau')
	disp(tau)
	
	while(1)
	  disp(j)
	  disp('x')
	  disp(x(end,:))
	  
	  if (mod(j,1) == 0)
	    plot(x(:,1),x(:,2),'k.-')
	    drawnow
	    pause
	  end
	
	  if (G.isElement(x(end,:)))
	    plot(x(:,1),x(:,2),'k.-')
	    drawnow
	    break
	  end
	  
	  u = C.getInputs(x(end,:));
	  ran = randi([1 size(u,1)], 1, 1);	
	  v = [v; u(ran,:)];
	  r = radNext(eta, u(ran,:), tau);
	  [t phi] = ode45(@ODE, [0 tau], x(end,:), [], u(ran,:));
	  x = [x; phi];
	  xNext = disturbance(phi(end,:), r);
	  x = [x; xNext];
	
	  disp('u')
	  disp(u(ran,:))
	
	  j = j + 1;
	end
      end
      k = k + 1;
    end
    savefig('alwaysEventually');
  end
      
      
	  
	  
    

  if (mode == 'R')
    openfig('problem');
    hold on
    drawnow
    
    I = SymbolicSet('I.bdd');
    plotCells(I, 'facecolor', colors(3,:)*0.5+0.5, 'edgec', colors(3,:), 'linew', 0.1)
    drawnow
    disp('Done plotting initial')
   
    
    x = I.points();
    x = x(1,:);
    
    v = [];
    
    j = 1;
    start = 1;
    for i = controllers:-1:1
      disp(['iteration: ' int2str(i)])
      
      C = SymbolicSet(['C/C' int2str(i) '.bdd'], 'projection', [1 2]);
      if (i == 1)
	G = SymbolicSet(['G.bdd']);
      else
	G = SymbolicSet(['Z/Z' int2str(i-1) '.bdd']);
      end
      
      points = G.points();
      plot(points(:,1), points(:,2), 'x', 'MarkerFaceColor', colors(mod(i,7)+1,:)*0.3+0.3, 'MarkerSize', 1.5);
      
      Z = SymbolicSet(['Z/Z' int2str(i) '.bdd']);      
      eta = Z.eta();
      eta = eta';
      tau = eta(1) * 3 / 2;
      
      disp('eta')
      disp(eta)
      disp('tau')
      disp(tau)
      
      while (1)
	disp(j)
	disp('x')
	disp(x(end,:))
	
	if (mod(j,1) == 0)
	  plot(x(:,1),x(:,2),'k.-')
	  drawnow
	  pause
	end
	
	if (G.isElement(x(end,:)))
	  plot(x(:,1),x(:,2),'k.-')
	  drawnow
	  break
	end
	
	u = C.getInputs(x(end,:));
	ran = randi([1 size(u,1)], 1, 1);	
	v = [v; u(ran,:)];
	r = radNext(eta, u(ran,:), tau);
	[t phi] = ode45(@ODE, [0 tau], x(end,:), [], u(ran,:));
	x = [x; phi];
	xNext = disturbance(phi(end,:), r);
	x = [x; xNext];
	
	disp('u')
	disp(u(ran,:))
	
	j = j + 1;
      end
    end
%      plotCells(G,'facecolor',colors(2,:)*0.5+0.5,'edgec',colors(2,:),'linew',.1)
    savefig('reach');
  end
end

function r = radNext(eta, u, tau)
  r = eta/2;

end

function xNext = disturbance(x, r)
  w = -r + (2 * r .* rand(size(x)));
  xNext = x + w; 
end

function dxdt = ODE(t,x,u)
  dxdt = zeros(size(x));
  dxdt(1) = (-1 * x(1) + 2 * x(2)) * u(2);
  dxdt(2) = (-1 * x(1) + 1 * x(2) + u(1)) * u(2);
end