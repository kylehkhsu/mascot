%
% simple.m
%
%

function unicycle (mode, controllers)

%    taus = tau;
%    for i = 1:10
%      tau = taus(i);
%      taus = [taus tau/3];
%    end

  addpath(genpath('/home/kylehsu/control/SCOTS+Adaptive'));
  
  % colors
  colors=get(groot,'DefaultAxesColorOrder');
  
  if (mode == 'S')
    figure
    hold on
    box on    
    drawnow
  
    % load and draw state space
    X = SymbolicSet('X.bdd', 'projection', [1 2]);
    lb = X.first();
    ub = X.last();   
    axis([lb(1)-1 ub(1)+1 lb(2)-1 ub(2)+1])
    plotCells(X, 'facecolor', 'none', 'edgec', [0.8 0.8 0.8], 'linew', 0.1);
    drawnow
    disp('Done plotting state space')
%      try
%        O = SymbolicSet('O.bdd', 'projection', [1 2]);
%        plotCells(O, 'facecolor', 'none', 'edgec', [0.8 0.8 0.8], 'linew', 0.1);
%        drawnow
%      catch
%        warning('No obstacles');
%      end
    
    savefig('system');
  end
  
  if (mode == 'P')
    openfig('system');
    hold on
    drawnow
    
    % load and draw obstacles
    try
      O = SymbolicSet('O.bdd', 'projection', [1 2]);
      plotCells(O, 'facecolor', colors(1,:)*0.5+0.5, 'edgec', colors(1,:), 'linew', 0.1)
      drawnow
      disp('Done plotting obstacles')
    catch
      warning('No obstacles');
    end
        
    % load and draw goal
    G = SymbolicSet('G.bdd', 'projection', [1 2]);
    plotCells(G, 'facecolor', colors(2,:)*0.5+0.5, 'edgec', colors(2,:), 'linew', 0.1)
    drawnow
    disp('Done plotting goal')
   
    % load and draw initial
    I = SymbolicSet('I.bdd', 'projection', [1 2]);
    plotCells(I, 'facecolor', colors(3,:)*0.5+0.5, 'edgec', colors(3,:), 'linew', 0.1)
    drawnow
    disp('Done plotting initial')
   
    savefig('problem');

  end
  
  if (mode == 'C')
    openfig('problem');
    hold on
    drawnow
    
    I = SymbolicSet('I.bdd', 'projection', [1 2]);
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
      
      F = SymbolicSet(['F/F' int2str(i) '.bdd']);      
      eta = F.eta();
      eta = eta';
      tau = eta(1);
      
      disp(eta)
      disp(tau)
      
      while (1)
	disp(j)
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
	
	disp(u(ran,:))
	
	v = [v; u(ran,:)];
	r = radNext(eta, u(ran,:), tau);
	phi = sysNext(x(end,:), u(ran,:), tau);
	x = [x;phi];
	
	xNext = disturbance(phi(end,:), r);
	x = [x; xNext];
	
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
    x = [x; x];
    v = [];
    
    G = SymbolicSet('G.bdd', 'projection', [1 2]);
    eta = [0.6 0.6 0.3];
    tau = 0.9;
    
    C = SymbolicSet('CSCOTS.bdd', 'projection', [1 2 3]);
    
    j = 1;
    
    while(1)
      disp(j)
      disp(x(end-1,:))
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
      [t phi] = ode45(@unicycle_ode, [0 tau], x(end,:), [], u(ran,:));
      x = [x; phi];
      xNext = disturbance(phi(end,:), r);
      xNext(3) = mod(xNext(3), pi);
      x = [x; xNext];
      
      j = j + 1;
    end       
  end    
end

function r = radNext(eta, u, tau)
  r = eta;
  r(1) = eta(1) + eta(3)*abs(u(1)) * tau;
  r(2) = eta(2) + eta(3)*abs(u(1)) * tau;
  r = r / 1.2;
end

function xNext = disturbance(x, r)
  w = -r + (2 * r .* rand(size(x)));
  xNext = x + w;
end

function dxdt = unicycle_ode(t,x,u)
  dxdt = zeros(3,1);
  dxdt(1)=u(1)*cos(x(3));
  dxdt(2)=u(1)*sin(x(3));
  dxdt(3)=u(2);
end
