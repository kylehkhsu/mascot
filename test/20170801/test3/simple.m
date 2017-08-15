%
% simple.m
%
%

function simple (mode, controllers)



  addpath(genpath('/home/kylehsu/control/SCOTS+Adaptive'));
  
  % colors
  colors=get(groot,'DefaultAxesColorOrder');
  
  if (mode == 'T')
    for i = 1:5
      X = SymbolicSet('X.bdd');
    end
    p = X.points();
    plot(p(:,1),p(:,2),'kx');
  end

  
  if (mode == 'S')
    figure
    hold on
    box on    
    drawnow
  
    % load and draw state space
    X = SymbolicSet('X.bdd');
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
   
    % load and draw initial
    I = SymbolicSet('I.bdd');
    plotCells(I, 'facecolor', colors(3,:)*0.5+0.5, 'edgec', colors(3,:), 'linew', 0.1)
    drawnow
    disp('Done plotting initial')
   
    savefig('problem');

  end
  
  if (mode == 'T')
    openfig('system');
    hold on
    drawnow
    
    Ff = SymbolicSet('T/Ff.bdd');
    Fc = SymbolicSet('T/Fc.bdd');
    plotCells(Ff, 'facecolor', colors(3,:)*0.5+0.5, 'edgec', colors(3,:), 'linew', 0.1)
    drawnow
    pause
    plotCells(Fc, 'facecolor', colors(2,:)*0.5+0.5, 'edgec', colors(2,:), 'linew', 0.1)
  end
    

  if (mode == 'C')
    openfig('problem');
    hold on
    drawnow
    
    I = SymbolicSet('I.bdd');
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
	G = SymbolicSet(['F/F' int2str(i-1) '.bdd']);
      end
      
      points = G.points();
      plot(points(:,1), points(:,2), 'x', 'MarkerFaceColor', colors(mod(i,7)+1,:)*0.3+0.3, 'MarkerSize', 1.5);
      
      F = SymbolicSet(['F/F' int2str(i) '.bdd']);      
      eta = F.eta();
      eta = eta';
      tau = eta(1);
      
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
	  pause
	  break
	end
	
	u = C.getInputs(x(end,:));
	ran = randi([1 size(u,1)], 1, 1);	
	v = [v; u(ran,:)];
	r = radNext(eta, u(ran,:), tau);
	phi = sysNext(x(end,:), u(ran,:), tau);
	xNext = disturbance(phi(end,:), r);
	x = [x; xNext];
	
	disp('u')
	disp(u(ran,:))
	
	j = j + 1;
      end
    end
%      plotCells(G,'facecolor',colors(2,:)*0.5+0.5,'edgec',colors(2,:),'linew',.1)
    savefig('simulation');
  end
  
  if (mode == 'Q')
    openfig('problem');
    hold on
    drawnow
    
    I = SymbolicSet('I.bdd');
    x = I.points();
    x = x(1,:);
    v = [];
    
    G = SymbolicSet('G.bdd');
    eta = G.eta();
    tau = 0.2;
    
    C = SymbolicSet('CSCOTS.bdd', 'projection', [1 2]);
    
    j = 1;
    
    while(1)
      disp(j)
      disp(x(end-1,:))
      disp(x(end,:))
      
      if (mod(j,3) == 0)
	plot(x(:,1),x(:,2),'k.-')
	drawnow
	pause
      end
      
      if (G.isElement(x(end,:)))
	plot(x(:,1),x(:,2),'k.-')
	drawnow
	pause
	break
      end
	
      u = C.getInputs(x(end,:));
      ran = randi([1 size(u,1)], 1, 1);        
      v = [v; u(ran,:)];
      r = radNext(eta, u(ran,:), tau);
      phi = sysNext(x(end,:), u(ran,:), tau);
      xNext = disturbance(phi(end,:), r);
      x = [x; xNext];
      
      j = j + 1;
    end       
  end    
end

function r = radNext(eta, u, tau)
  r = eta + abs(u) .* tau ./ 2;
  r = r ./ 2;  
end

function xNext = disturbance(x, r)
  w = -r + (2 * r .* rand(size(x)));
  xNext = x + w; 
end

function phi = sysNext(x, u, tau)
  phi = x + u * tau;
end
