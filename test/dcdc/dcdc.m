function dcdc (mode, abstractions)
  addpath(genpath('/home/kylehsu/control/SCOTS+Adaptive'));
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
    (ub(1)-lb(1))/10
    axis([lb(1)-(ub(1)-lb(1))/10 ub(1)+(ub(1)-lb(1))/10 lb(2)-(ub(2)-lb(2))/10 ub(2)+(ub(2)-lb(2))/10])
%      plotCells(X, 'facecolor', 'none', 'edgec', [0.8 0.8 0.8], 'linew', 0.1);
    x = X.points;
    plot(x(:,1),x(:,2),'.','color',[0.8 0.8 0.8]);
    drawnow
    disp('Done plotting state space')
    savefig('system');
  end
  
  if (mode == 'P')
    openfig('system');
    hold on
    drawnow
    
    % load and draw possible safe zone
%      S = SymbolicSet('S.bdd');
%      plotCells(S, 'facecolor', colors(2,:)*0.5+0.5, 'edgec', colors(2,:), 'linew', 0.1)

    v = [1.15 5.45; 1.55 5.45; 1.15 5.85; 1.55 5.85 ];
    patch('vertices',v,'faces',[1 2 4 3],'facecolor','none','edgec',colors(2,:),'linew',1)
    hold on
    drawnow
    disp('Done plotting possible safe zone')
    
    savefig('problem');
  end
  
  if (mode == 'V') % visualization
    openfig('problem');
    hold on
    drawnow
    
    D = SymbolicSet('D.bdd');
    d = D.points;
    plot(d(:,1),d(:,2),'.','color', [0.75 0.85 0.95], 'MarkerSize', 5);
    
    for i = abstractions:-1:1
    
      Z = SymbolicSet(['Z/Z' int2str(i) '.bdd']);
      z = Z.points;
      plot(z(:,1),z(:,2),'.','color', colors(i,:), 'MarkerSize', 5);
      drawnow
      pause
    end
    
    savefig('visualization');
    
    openfig('problem');
    hold on
    drawnow
    
    D = SymbolicSet('D.bdd');
    d = D.points;
    plot(d(:,1),d(:,2),'.','color', [0.75 0.85 0.95], 'MarkerSize', 5);
  end    
    
  
  if (mode == 'A') % adaptive safe
    openfig('visualization');
    hold on
    drawnow
    
    tau = 0.5;
    x = [1.2 5.52];
    v = [];
    T = 100;
    
    Z = SymbolicSet('Z/Z3.bdd');
    
    
    for j = 1:T/tau
      disp(j)
      disp(x(end,:))
      
      if (mod(j,5) == 0)
	plot(x(:,1),x(:,2),'k.-')
	drawnow
	pause
      end
      
      foundController = 0;      
      for i = 1:abstractions
	Z = SymbolicSet(['Z/Z' int2str(i) '.bdd']);
	if (Z.isElement(x(end,:)))
	  eta = Z.eta;
	  eta = eta';
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
      r = eta / 2;
      [t r] = ode45(@radODE, [0 tau], r, [], u(ran,:));
      [t phi] = ode45(@sysODE, [0 tau], x(end,:), [], u(ran,:));
      x = [x; phi];
      xNext = disturbance(phi(end,:), r);
      x = [x; xNext];    
    end
    savefig('safe');
 
  end
  
  if (mode == 'Q') % scots safe
    openfig('problem');
    hold on
    drawnow
    
    C = SymbolicSet('CSCOTS.bdd', 'projection', [1 2]);
    c = C.points;
    plot(c(:,1),c(:,2),'.', 'color', [0.75 0.85 0.95], 'MarkerSize', 5);
    
    
    S = SymbolicSet('S.bdd');
    eta = S.eta();
    eta = eta';
    tau = 0.5;
    x = [1.2 5.6];
    v = [];
    
    T = 100;
    
    for j = 1:T/tau
      disp(j)
      disp(x(end,:))
      
      if (mod(j,5) == 0)
	plot(x(:,1),x(:,2),'k.-')
	drawnow
	pause
      end    
    
      u = C.getInputs(x(end,:));
      ran = randi([1 size(u,1)], 1, 1);
      v = [v; u(ran,:)];
      r = eta / 2;
      [t r] = ode45(@radODE, [0 tau], r, [], u(ran,:));
      [t phi] = ode45(@sysODE, [0 tau], x(end,:), [], u(ran,:));
      x = [x; phi];
      xNext = disturbance(phi(end,:), r);
      x = [x; xNext];
    end
    savefig('safeSCOTS');
  end
  
  
end 

function dxdt = sysODE(t, x, u)
  r0 = 1;
  vs = 1;
  rl = 0.05;
  rc = rl / 10;
  xl = 3;
  xc = 70;
  
  b = [(vs / xl); 0];
  
  if (u(1) == 1)
    A = [ -rl / xl  0 ;  0  (-1 / xc) * (1 / (r0 + rc)) ] ;
  else
    A = [ (-1 / xl) * (rl + ((r0 * rc) / (r0 + rc)))  ((-1 / xl) * (r0 / (r0 + rc))) / 5 ;
          5 * (r0 / (r0 + rc)) * (1 / xc)  (-1 / xc) * (1 / (r0 + rc)) ];
  end
  
  dxdt = A*x + b;
end
  
function xNext = disturbance(x, r)
  w = -r + (2 * r .* rand(size(x)));
  xNext = x + w;
end
  
function drdt = radODE(t, r, u)
  r0 = 1;
  rl = 0.05;
  rc = rl / 10;
  xl = 3;
  xc = 70;
  
  if (u(1) == 1)
    A = [ -rl / xl  0 ;  0  (-1 / xc) * (1 / (r0 + rc)) ] ;
  else
    A = [ (-1 / xl) * (rl + ((r0 * rc) / (r0 + rc)))  ((-1 / xl) * (r0 / (r0 + rc))) / 5 ;
          5 * (r0 / (r0 + rc)) * (1 / xc)  (-1 / xc) * (1 / (r0 + rc)) ];
  end
  
  drdt = A * r;

end
  
