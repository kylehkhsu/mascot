%
% simple.m
%
%

function simple (mode, paws)

  %% configuration

  %% target set
  L = [0.65 0; 0 .65];
  c = [-7; -7];

  % initial state
  x0 = [7 8];
  g0 = [-2, -1];
   
  x = x0;
  v = [];
 
  addpath(genpath('/home/kylehsu/control/SCOTS+Adaptive'));
  clear set
  % close all
  
  % colors
  colors=get(groot,'DefaultAxesColorOrder');

  openfig('system');
  hold on
  
  if (mode == 0)
    filename = 'before';
    
    % load and draw goal
    set = SymbolicSet('G.bdd');
    plotCells(set,'facecolor',colors(2,:)*0.5+0.5,'edgec',colors(2,:),'linew',.1)
    drawnow
    disp('Done plotting target set')
    if (paws)
      pause
    end
  end

  if (mode >= 1)
    filename = 'controller';
    % load and draw goal
    set = SymbolicSet('G.bdd');
    plotCells(set,'facecolor',colors(2,:)*0.5+0.5,'edgec',colors(2,:),'linew',.1)
    drawnow
    disp('Done plotting goal set')
    if (paws)
      pause
    end

    
    % load controller
    controller = SymbolicSet('C.bdd','projection',[1 2]);

    i = 1;
    start = 1;

    while(1)
      disp(i)
      disp(x(end,:))
  
      if ( (x(end,:)-c')*L'*L*(x(end,:)'-c)<=1 )
	break;
      end

%        if (x(end,:) == g0)
%  	break;
%        end
      
      u = controller.getInputs(x(end,:));
      ran = randi([1 size(u,1)], 1, 1);
      v = [v; u(ran,:)];
      
      r = radNext();

      
      phi = sysNext(x(end,:), u(ran,:));
      xNext = disturbance(phi(end,:), r)

      x = [x; xNext];
      
      if (1)
	if (mod(i,3) == 0)
	  stop = i;
	  plot(x(start:stop,1),x(start:stop,2),'k.-')
	  drawnow
	  start = i;
	  pause
        end
      end
      
      i = i + 1;      
    end
    disp('Done simulating trajectory')
    
    % draw trajectory
    plot(x(:,1),x(:,2),'k.-')
    drawnow
    disp('Done plotting trajectory')
    if (paws)
      pause
    end
  end
  
  % draw starting state
  plot(x(1,1),x(1,2),'.','color',colors(5,:),'markersize',20)
  drawnow
  disp('Done plotting starting state')
  if (paws)
    pause
  end  

  box on
  axis([-11 11 -11 11])

  saveas(gcf, filename);
end

function r = radNext()
  r = [0.6 0.6];
end

function xNext = disturbance(x, r)
  w = -r + (2 * r .* rand(size(x)));
  xNext = x + w;
  
end

function xNext = sysNext(x, u)
  xNext = x + u;
end
