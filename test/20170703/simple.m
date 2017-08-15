%
% simple.m
%
%

function simple (mode, paws, iter)

  %% configuration

  % initial state
  x0 = [7 8];
  g0 = [-2, -1];
   
  x = x0;
  v = [];
 
  addpath(genpath('/home/kylehsu/control/SCOTS+Adaptive'));
  clear set
  
  % colors
  colors=get(groot,'DefaultAxesColorOrder');

  openfig('system');
  hold on
  box on
  axis([-11 11 -11 11])
  
  if (mode == 0)
    filename = 'before';
    
    % load and draw goal
    set = SymbolicSet(['G1.bdd']);
    plotCells(set,'facecolor',colors(2,:)*0.5+0.5,'edgec',colors(2,:),'linew',.1)
    drawnow
    disp('Done plotting target set')
    if (paws)
      pause
    end
  end

    if (mode >= 1)
    filename = 'controller';
    lastIter = 1;    
    while (lastIter <= iter)
      R = SymbolicSet(['G' int2str(lastIter) '.bdd']);
      try
	p = R.points();
	break
      catch
	lastIter = lastIter + 1;
      end
    end
    j = 1;
    start = 1;
    for i = iter:-1:lastIter
      clear G C eta
      disp(['iteration: ' int2str(i)])
      try
	% load and draw goal
	G = SymbolicSet(['G' int2str(i) '.bdd']);
	eta = G.eta();
	colorInd = mod(i, size(colors,1) - 2) + 2;	
	plotCells(G,'facecolor',colors(colorInd,:)*0.5+0.5,'edgec',colors(colorInd,:),'linew',.1)
	drawnow
	disp('Done plotting target set')
	if (paws)
	  pause
	end
      catch
	disp('Goal is empty.')
	continue
      end
      C = SymbolicSet(['C' int2str(i) '.bdd'], 'projection', [1 2]);

      while (1)
	disp(j)
	disp(x(end,:))
	
	if (mod(j,1) == 0)
	  stop = j;
	  plot(x(start:stop,1),x(start:stop,2),'k.-')
	  drawnow
	  start = j;
	  pause
        end
	
	if (G.isElement(x(end,:)))
	  break
	end
	
	u = C.getInputs(x(end,:));
	ran = randi([1 size(u,1)], 1, 1);
	v = [v; u(ran,:)];
	r = radNext(eta, u(ran,:));
	phi = sysNext(x(end,:), u(ran,:));
	xNext = disturbance(phi(end,:), r);
	x = [x; xNext];
	
	j = j + 1;
      end
      if (R.isElement(x(end,:)))
	break
      end
    end
    plotCells(R,'facecolor',colors(2,:)*0.5+0.5,'edgec',colors(2,:),'linew',.1)
    drawnow
    plot(x(:,1),x(:,2),'k.-')
    drawnow
  end
  
  % draw starting state
  plot(x(1,1),x(1,2),'.','color',colors(5,:),'markersize',20)
  drawnow
  disp('Done plotting starting state')

  saveas(gcf, filename);
end

function r = radNext(eta, u)
  r = [0.6 0.6];
end

function xNext = disturbance(x, r)
  w = -r + (2 * r .* rand(size(x)));
  xNext = x + w;
  
end

function phi = sysNext(x, u)
  phi = x + u;
end
