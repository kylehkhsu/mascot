function whirlpoolSafe (mode, numAbs)
	w = [0.05 0.05];
  addpath(genpath('../..'));
  colors=get(groot,'DefaultAxesColorOrder');
  
  if (mode == 'S')
    figure
    hold on
    box on    
    drawnow
  
    % load and draw state space
    X = SymbolicSet('plotting/X.bdd');
    lb = X.first();
    ub = X.last();
    axis([lb(1)-(ub(1)-lb(1))/10 ub(1)+(ub(1)-lb(1))/10 lb(2)-(ub(2)-lb(2))/10 ub(2)+(ub(2)-lb(2))/10])
    plotCells(X, 'facecolor', 'none', 'edgec', [0.8 0.8 0.8], 'linew', 0.1);
%    x = X.points;
%    plot(x(:,1),x(:,2),'.','color',[0.8 0.8 0.8]);
    drawnow
    disp('Done plotting state space')
    savefig('system');
  end
  
  if (mode == 'P')
    openfig('system');
    hold on
    drawnow

    S = SymbolicSet(['S/S' num2str(numAbs) '.bdd']);
    plotCells(S, 'facecolor', colors(3,:)*0.5+0.5, 'edgec', colors(3,:), 'linew', 0.1)

    O = SymbolicSet('plotting/O.bdd');
    try
      plotCells(O, 'facecolor', colors(1,:)*0.5+0.5, 'edgec', colors(1,:), 'linew', 0.1)
    end

    hold on
    drawnow
    disp('Done plotting possible safe zone')
    
    savefig('problem');
  end
  
  if (mode == 'V') % visualization
    openfig('problem');
    hold on
    drawnow
    
    D = SymbolicSet('plotting/D.bdd');
    d = D.points;
    plot(d(:,1),d(:,2),'.','color', [0.75 0.85 0.95], 'MarkerSize', 20);
    
    for i = numAbs:-1:1
      try
	Z = SymbolicSet(['Z/Z' int2str(i) '.bdd']);
	z = Z.points;
	plot(z(:,1),z(:,2),'.','color', colors(i,:)*0.5+0.5, 'MarkerSize', 20);
	drawnow
	pause
      end
    end
    
    savefig('visualization');
    
    openfig('problem');
    hold on
    drawnow
    
    D = SymbolicSet('plotting/D.bdd');
    d = D.points;
    plot(d(:,1),d(:,2),'.','color', [0.75 0.85 0.95], 'MarkerSize', 5);
  end    
    
  
  if (strcmp(mode, 'safe')) % adaptive safe
    openfig('visualization');
    hold on
    drawnow
    
    x = [5.7 -1];
    v = [];
    N = 50;        
    
    for j = 1:N
      disp(j)
      disp(x(end,:))
      
      if (mod(j,1) == 0)
	plot(x(:,1),x(:,2),'k.-')
	drawnow
	pause
      end
      
      foundController = 0;      
      for i = 1:numAbs
	Z = SymbolicSet(['Z/Z' int2str(i) '.bdd']);
	if (Z.isElement(x(end,:)))
	  eta = Z.eta;
	  eta = eta';
	  tau = getTau(eta);
	  C = SymbolicSet(['C/C' int2str(i) '.bdd']);
	  foundController = 1;
	  disp(['Controller: ' int2str(i)])
	  break;
	end
      end
      if (~foundController)
	disp('No controller for state found. Breaking.')
	break;
      end
      
      u = C.getInputs(x(end,:));
      ran = randi([1 size(u,1)], 1, 1);
      v = [v; u(ran,:)];
      d = disturbance(w);
      [t phi] = ode45(@sysODE, [0 tau], x(end,:), [], u(ran,:), d);
      x = [x; phi];
      
     disp('u')
     disp(u(ran,:))
     disp('d')
     disp(d)

    end
    savefig('safe'); 
  end
end

function d = disturbance(w)
  d = -w + (2 * w .* rand(size(w)));
end

function dxdt = sysODE(t,x,u,d)
  dxdt = zeros(size(x));
  dxdt(1) = (-0.5 * x(1) + 1 * x(2)) * u(2);
  dxdt(2) = (-0.5 * x(1) + 0.5 * x(2) + u(1)) * u(2);
  dxdt = dxdt + d';
end
  
function tau = getTau(eta)
  tau = eta(1)*3/2;
end
