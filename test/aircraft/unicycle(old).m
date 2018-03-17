
function unicycle (mode, numAbs, etaAbs, controllers)
  w = [0.05 0.05 0];
  addpath(genpath('../..'));
  addpath(genpath('~/ownCloud/C++/SCOTS_modified/mfiles/'));
  
  % colors
  colors=get(groot,'DefaultAxesColorOrder');
  
  if (strcmp(mode, 'S'))
    figure
    hold on
    box on    
    drawnow
  
    % load and draw state space
    X = SymbolicSet('plotting/X.bdd', 'projection', [1 2]);
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
        
    % load and draw goal
    G = SymbolicSet('plotting/G.bdd', 'projection', [1 2]);
    plotCells(G, 'facecolor', colors(2,:)*0.5+0.5, 'edgec', colors(2,:), 'linew', 0.1)
    drawnow
    disp('Done plotting goal')
   
    % load and draw initial
    I = SymbolicSet('plotting/I.bdd', 'projection', [1 2]);
    plotCells(I, 'facecolor', colors(3,:)*0.5+0.5, 'edgec', colors(3,:), 'linew', 0.1)
    drawnow
    disp('Done plotting initial')
   
    savefig('problem');

  end
  
  if (strcmp(mode, 'C'))
      openfig('problem')
      
      for ii=1:controllers
        C = SymbolicSet(['C/C' int2str(ii) '.bdd']);
        p = C.points;
        plotColor = colors(mod(ii,7)+1,:)*0.3+0.3;
        scatter(p(:,1),p(:,2),'Marker','o','MarkerFaceColor',plotColor);
        pause
      end
  end
  
  if (strcmp(mode,'Q')) % plot controller domains
      openfig('problem');
      hold on;
      cmap = [hot(numAbs) 0.5*ones(numAbs,1)];
      
%       if numAbs~=size(etaAbs,1) % Sanity check
%           error();
%       end
      
      for ii=controllers:-1:1
          C = SymbolicSet(['C/C' int2str(ii) '.bdd'], 'projection', [1 2]);
          eta = C.eta(1:2);
          for jj=1:numAbs
              if isequal(eta',etaAbs(jj,1:2))
                  layer = jj;
                  break;
              end
          end
          
          
          if (ii == 1)
            G = SymbolicSet(['G/G' int2str(numAbs) '.bdd'], 'projection', [1 2]);
          else
            G = SymbolicSet(['Z/Z' int2str(ii-1) '.bdd'], 'projection', [1 2]);
          end
          
          points = setdiff(C.points,G.points,'rows');
          for jj=1:size(points,1)
              rectangle('Position',[points(jj,:)-eta' 2*eta'],'FaceColor',cmap(layer,:));
              %alpha(s,0.5);
          end
      end
  end
           
              
  
  if (strcmp(mode,'R'))
    openfig('problem');
    hold on
    drawnow
    
    I = SymbolicSet('plotting/I.bdd');
    x = I.points();
    x = x(1,:);
    
    v = [];
    
    j = 1;
    start = 1;
    for i = controllers:-1:1
      disp(['iteration: ' int2str(i)])
      
      C = SymbolicSet(['C/C' int2str(i) '.bdd'], 'projection', [1 2 3]);
      if (i == 1)
	G = SymbolicSet(['G/G' int2str(numAbs) '.bdd']);
      else
	G = SymbolicSet(['Z/Z' int2str(i-1) '.bdd']);
      end
      
      points = G.points();
      plot(points(:,1), points(:,2), 'x', 'MarkerFaceColor', colors(mod(i,7)+1,:)*0.3+0.3, 'MarkerSize', 2);       
      
      Z = SymbolicSet(['Z/Z' int2str(i) '.bdd']);      
      eta = Z.eta();
      eta = eta';
      tau = eta(1)*3/2;
      
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
	d = disturbance(w);
	[t phi] = ode45(@sysODE, [0 tau], x(end,:), [], u(ran,:), d);
	x = [x; phi];
	  
	disp('u')
	disp(u(ran,:))
	disp('d')
	disp(d)
	
	j = j + 1;
      end
    end
%      plotCells(G,'facecolor',colors(2,:)*0.5+0.5,'edgec',colors(2,:),'linew',.1)
    savefig('simulation');
  end
  
  if (strcmp(mode,'scots'))
    openfig('problem');
    hold on
    drawnow
    
    I = SymbolicSet('scots/I.bdd');
    x = I.points();
    x = x(1,:);
    x = [x; x];
    v = [];
    
    G = SymbolicSet('scots/G.bdd');
    eta = G.eta;
    tau = eta(1)*3/2;
    
    C = SymbolicSet('scots/C.bdd', 'projection', [1 2 3]);
    
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
      d = disturbance(w);
      [t phi] = ode45(@unicycle_ode, [0 tau], x(end,:), [], u(ran,:), d);
      x = [x; phi];
      
      j = j + 1;
    end
    savefig('scots/simulation')
  end    
end

function d = disturbance(w)
  d = -w + (2 * w .* rand(size(w)));
end

function dxdt = sysODE(t,x,u, d)
  dxdt = zeros(3,1);
  dxdt(1)=u(1)*cos(x(3));
  dxdt(2)=u(1)*sin(x(3));
  dxdt(3)=u(2);
  dxdt = dxdt + d';
end
