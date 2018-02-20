
function unicycle (mode, numAbs, tauSet)
% Plot results for the unicycle example (for safety)
% Inputs:
%       - mode (char):  'S' = plot state space
%                       'P' = plot safe set
%                       'T' = plot transition domains in various layers
%                       'Safe' = simulate the safety controller, you can
%                       change the initial state in the function body
%           Note: 'S' and 'P' must always be run first (in the respective
%           order) before 'T' or 'Safe'
%       - numAbs (int): number of abstraction as was used during controller
%       synthesis
%       - tauSet (double array): an array (same length as numAbs)
%       containing the time discretization parameters in different layers,
%       with tauSet(1) being for the coarsest and tauSet(end) being the
%       finest

w = [0.05 0.05 0]; % disturbance, *must* be same as was used during controller synthesis
addpath(genpath('../..')); % adding the mfiles/ in MATLAB path

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
    
    % load and draw safe set
    try
        SafeSet = SymbolicSet('plotting/SafeInner.bdd', 'projection', [1 2]);
        plotCells(SafeSet, 'facecolor', colors(1,:)*0.5+0.5, 'edgec', colors(1,:), 'linew', 0.1)
        drawnow
        disp('Done plotting safe set')
    catch
        warning('No safe set');
    end
   
    savefig('problem');
    
end

if (strcmp(mode, 'T'))
    openfig('problem')
    T1 = SymbolicSet('T/T1.bdd', 'projection', [1 2]);
    p1 = T1.points;
    x = plot(p1(:,1),p1(:,2),'ko');
    pause
    delete(x)
    
    T2 = SymbolicSet('T/T2.bdd', 'projection', [1 2]);
    p2 = T2.points;
    x = plot(p2(:,1),p2(:,2),'bo');
    pause    
    delete(x)
    
    T3 = SymbolicSet('T/T3.bdd', 'projection', [1 2]);
    p3 = T3.points;
    plot(p3(:,1),p3(:,2),'ro');
    title('T')
end

% debug purpose
if (strcmp(mode, 'Zs'))
    fig0 = openfig('problem');
    figure(fig0)
    axis0 = gca;
    fig1 = openfig('problem');
    figure(fig1)
    axis1 = gca;
    fig2 = openfig('problem');
    figure(fig2)
    axis2 = gca;
    
    
    for ii=1:2
        T1 = SymbolicSet(['validZsInit/Zs_1.bdd']);
        p1 = T1.points;
        z = plot(axis0,p1(:,1),p1(:,2),'ko');
        
        T1 = SymbolicSet(['Zs' int2str(ii) '/Zs_1.bdd']);
        p1 = T1.points;
        x = plot(axis1,p1(:,1),p1(:,2),'ko');
        
        T1 = SymbolicSet(['validZs' int2str(ii) '/Zs_1.bdd']);
        p1 = T1.points;
        y = plot(axis2,p1(:,1),p1(:,2),'ko');
        pause
        delete(x)
        delete(y)
        delete(z)

        T2 = SymbolicSet(['validZsInit/Zs_2.bdd']);
        p2 = T2.points;
        z = plot(axis0,p2(:,1),p2(:,2),'ko');
        
        T2 = SymbolicSet(['Zs' int2str(ii) '/Zs_2.bdd']);
        p2 = T2.points;
        x = plot(axis1,p2(:,1),p2(:,2),'bo');
        
        T2 = SymbolicSet(['validZs' int2str(ii) '/Zs_2.bdd']);
        p2 = T2.points;
        y = plot(axis2,p2(:,1),p2(:,2),'bo');
        pause    
        delete(x)
        delete(y)
        delete(z)
        
        T3 = SymbolicSet(['validZsInit/Zs_3.bdd']);
        p3 = T3.points;
        z = plot(axis0,p3(:,1),p3(:,2),'ko');

        T3 = SymbolicSet(['Zs' int2str(ii) '/Zs_3.bdd']);
        p3 = T3.points;
        x = plot(axis1,p3(:,1),p3(:,2),'ro');
        
        T3 = SymbolicSet(['validZs' int2str(ii) '/Zs_3.bdd']);
        p3 = T3.points;
        y = plot(axis2,p3(:,1),p3(:,2),'ro');
        pause
        delete(x)
        delete(y)
        delete(z)
    end
    title('Zs')
end


if (strcmp(mode,'Safe'))
    openfig('problem');
    hold on
    drawnow
   
    x = [1 0.5 0];   % Initial state  
    v = [];    
    j = 1;
    
    % Load controllers (C) and controller domains (Z)
    C = cell(numAbs,1); 
    Z = cell(numAbs,1);
    for ii=1:numAbs
        C{ii} = SymbolicSet(['C/C' int2str(ii) '.bdd']);
        Z{ii} = SymbolicSet(['Z/Z' int2str(ii) '.bdd']);
    end
    
    for rounds=1:400 % time steps (variable length, depending on the abstraction layer) for which the simulation will run
        % find the controller whose domain include the present state
        for ii=1:numAbs
            if (Z{ii}.isElement(x(end,:)))
                tau = tauSet(ii);
                break;
            end
            if (ii==numAbs)
                error('Current state is outside the controller domains');
            end
        end
        u = C{ii}.getInputs(x(end,:));
        
        % choose the control input randomly from many possible options
        ran = randi([1 size(u,1)], 1, 1);
        v = [v; u(ran,:)];
        
        % random disturbance which respects the given bound
        d = disturbance(w);
        
        % solve the ODE
        [t phi] = ode45(@sysODE, [0 tau], x(end,:), [], u(ran,:), d);
        % contains all the points along the trajectory
        x = [x; phi];

        plot(x(:,1),x(:,2),'k.-')
        drawnow
%       pause

        disp('u')
        disp(u(ran,:))
        disp('d')
        disp(d)
       
    end

    %      plotCells(G,'facecolor',colors(2,:)*0.5+0.5,'edgec',colors(2,:),'linew',.1)
    savefig('simulation');
end

end

function d = disturbance(w)
d = -w + (2 * w .* rand(size(w)));
end

% system ODE
function dxdt = sysODE(t,x,u, d)
dxdt = zeros(3,1);
dxdt(1)=u(1)*cos(x(3));
dxdt(2)=u(1)*sin(x(3));
dxdt(3)=u(2);
dxdt = dxdt + d';
end
