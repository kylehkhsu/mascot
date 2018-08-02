
function simple (mode, numAbs, controllers, progression)
w = [0.05 0.05];
addpath(genpath('../..'));
addpath(genpath('~/ownCloud/C++/SCOTS_modified/mfiles_v0.2/'));

% colors
colors=get(groot,'DefaultAxesColorOrder');
if (strcmp(mode, 'Z'))
    openfig('problem')
    C1 = StaticController('C/C1');
    p1 = C1.domain;
    plot(p1(:,1),p1(:,2),'ko');
    pause
    
    Z2 = StaticController('Z/Z2');
    p2 = Z2.domain;
    plot(p2(:,1),p2(:,2),'bo');
    pause
    
    %       Z22 = StaticController('D/Z2.scs');
    %       p2 = Z22.domain;
    %       plot(p2(:,1),p2(:,2),'ro');
    
    Z3 = StaticController('Z/Z3');
    p3 = Z3.domain;
    plot(p3(:,1),p3(:,2),'ro');
    pause
    
    Z4 = StaticController('Z/Z4');
    p4 = Z4.domain;
    plot(p4(:,1),p4(:,2),'go');
    title('Z')
    pause
    
    Z5 = StaticController('Z/Z5');
    p5 = Z5.domain;
    plot(p5(:,1),p5(:,2),'go');
end
if (strcmp(mode, 'D'))
    openfig('problem')
    C = StaticController('C/C1.scs');
    p = C.domain;
    plot(p(:,1),p(:,2),'ko');
    pause
    D1 = StaticController('D1/1.scs');
    p = D1.domain;
    plot(p(:,1),p(:,2),'bo');
    pause
    D2 = StaticController('D2/1.scs');
    p = D2.domain;
    plot(p(:,1),p(:,2),'ro');
    pause
end

if (strcmp(mode, 'CZ'))
    % visualize progression. black = sub-controller domain, red = sub-target
    openfig('problem')
    
    if strcmp(progression, 'b')
        for ii=1:controllers
            disp(num2str(ii))
            if ii ~= 1
                Z = Goal(['Z/Z' int2str(ii-1)]);
                p = Z.points;
                Zset = plot(p(:,1),p(:,2),'ro');
                pause
            end
            
            C = StaticController(['C/C' int2str(ii)]);
            p = C.domain;
            %         plotColor = colors(mod(ii,7)+1,:)*0.3+0.3;
            Cset = plot(p(:,1),p(:,2),'ko');
            pause
            
            if ii ~= 1
                delete(Zset)
            end
            delete(Cset)
        end
    elseif strcmp(progression, 'f')
        for ii = controllers:-1:1
            disp(num2str(ii))
            C = StaticController(['C/C' int2str(ii)]);
            p = C.domain;
            %         plotColor = colors(mod(ii,7)+1,:)*0.3+0.3;
            Cset = plot(p(:,1),p(:,2),'ko');
            pause
            
            if ii ~= 1
                Z = Goal(['Z/Z' int2str(ii-1)]);
                p = Z.points;
                Zset = plot(p(:,1),p(:,2),'ro');
                pause
            end       
            
            if ii ~= 1
                delete(Zset)
            end
            delete(Cset)
        end
    end
end

if (strcmp(mode, 'T'))
    openfig('problem')
    T1 = StaticController('T/T1.scs', 'projection', [1 2]);
    p1 = T1.domain;
    plot(p1(:,1),p1(:,2),'ko');
    pause
    
    T2 = StaticController('T/T2.scs', 'projection', [1 2]);
    p2 = T2.domain;
    plot(p2(:,1),p2(:,2),'bo');
    pause
    
    T3 = StaticController('T/T3.scs', 'projection', [1 2]);
    p3 = T3.domain;
    plot(p3(:,1),p3(:,2),'ro');
    title('T')
end

if (strcmp(mode, 'uTs'))
    openfig('problem')
    T1 = StaticController('uTs/T1.scs', 'projection', [1 2]);
    p1 = T1.domain;
    plot(p1(:,1),p1(:,2),'ko');
    pause
    
    T2 = StaticController('uTs/T2.scs', 'projection', [1 2]);
    p2 = T2.domain;
    plot(p2(:,1),p2(:,2),'bo');
    pause
    
    T3 = StaticController('uTs/T3.scs', 'projection', [1 2]);
    p3 = T3.domain;
    plot(p3(:,1),p3(:,2),'ro');
    title('T')
end

if (strcmp(mode, 'X'))
    openfig('problem')
    X1 = StaticController('X/X1.scs');
    p1 = X1.domain;
    plot(p1(:,1),p1(:,2),'ko');
    pause
    
    X2 = StaticController('X/X2.scs');
    p2 = X2.domain;
    plot(p2(:,1),p2(:,2),'bo');
    pause
    
    X3 = StaticController('X/X3.scs');
    p3 = X3.domain;
    plot(p3(:,1),p3(:,2),'ro');
    title('X')
end

if (strcmp(mode, 'G'))
    openfig('problem')
    X1 = Goal('G/G1','init');
    p1 = X1.points;
    plot(p1(:,1),p1(:,2),'ko');
    pause
    
    X2 = Goal('G/G2','init');
    p2 = X2.points;
    plot(p2(:,1),p2(:,2),'bo');
    pause
    
    X3 = Goal('G/G3','init');
    p3 = X3.points;
    plot(p3(:,1),p3(:,2),'ro');
    title('X')
end

if (strcmp(mode, 'S'))
    figure
    hold on
    box on
    drawnow
    
    % load and draw state space
    X = UniformGrid('plotting/X');
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
        O = StaticController('plotting/O.scs');
        plotCells(O, 'facecolor', colors(1,:)*0.5+0.5, 'edgec', colors(1,:), 'linew', 0.1)
        drawnow
        disp('Done plotting obstacles')
    catch
        warning('No obstacles');
    end
    
    % load and draw goal
    G = StaticController('plotting/G.scs');
    plotCells(G, 'facecolor', colors(2,:)*0.5+0.5, 'edgec', colors(2,:), 'linew', 0.1)
    drawnow
    disp('Done plotting goal')
    
    %     % load and draw initial
    %     I = StaticController('plotting/I.scs');
    %     plotCells(I, 'facecolor', colors(3,:)*0.5+0.5, 'edgec', colors(3,:), 'linew', 0.1)
    %     drawnow
    %     disp('Done plotting initial')
    
    savefig('problem');
    
end

if (strcmp(mode,'R'))
    
    disp('w')
    disp(w)
    
    openfig('problem');
    hold on
    drawnow
    
    %     I = StaticController('plotting/I.scs');
    %     x = I.domain();
    %     x = x(1,:);
    x = [0.5 1];
    
    v = [];
    
    j = 1;
    start = 1;
    for i = controllers:-1:1
        disp(['iteration: ' int2str(i)])
        
        C = StaticController(['C/C' int2str(i)]);
        if (i == 1)
%             G = StaticController(['G/G' int2str(numAbs) '.scs']);
            G = Goal(['G/G' int2str(numAbs)]);
        else
%             G = StaticController(['Z/Z' int2str(i-1) '.scs']);
            G = Goal(['Z/Z' int2str(i-1)]);
        end
        
%         points = G.domain();
        points = G.points;
        set = plot(points(:,1),points(:,2),'ko');
%         Z = StaticController(['Z/Z' int2str(i) '.scs']);
        Z = Goal(['Z/Z' int2str(i)]);
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
%                 pause
            end
            
            if (G.isElement(x(end,:)))
%             if (ismember(x(end,:),G.points))
                plot(x(:,1),x(:,2),'k.-')
                drawnow
                delete(set)
                break
            end
            
            u = C.control(x(end,:));
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
if strcmp(mode,'plot')
    openfig('system');
    for i = controllers:-1:1
        Z = StaticController(['G/G' int2str(numAbs) '.scs']);
        plotCells(Z, 'facecolor', colors(1,:), 'edgec', colors(1,:), 'linew', 0.1)
        drawnow
    end
    % load and draw obstacles
    try
        O = StaticController('plotting/O.scs', 'projection', [1 2]);
        plotCells(O, 'facecolor', colors(1,:)*0.5+0.5, 'edgec', colors(1,:), 'linew', 0.1)
        drawnow
        disp('Done plotting obstacles')
    catch
        warning('No obstacles');
    end
    
    % load and draw goal
    G = StaticController('plotting/G.scs', 'projection', [1 2]);
    plotCells(G, 'facecolor', colors(2,:)*0.5+0.5, 'edgec', colors(2,:), 'linew', 0.1)
    drawnow
    disp('Done plotting goal')
    
    % load and draw initial
    I = StaticController('plotting/I.scs', 'projection', [1 2]);
    plotCells(I, 'facecolor', colors(3,:)*0.5+0.5, 'edgec', colors(3,:), 'linew', 0.1)
    drawnow
    disp('Done plotting initial')
    
    savefig('figure')
end
if (strcmp(mode,'scots'))
    openfig('problem');
    hold on
    drawnow
    
    I = StaticController('scots/I.scs');
    x = I.domain();
    x = x(1,:);
    x = [x; x];
    v = [];
    
    G = StaticController('scots/G.scs');
    eta = G.eta;
    tau = eta(1)*3/2;
    
    C = StaticController('scots/C.scs', 'projection', [1 2]);
    
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
        
        u = C.control(x(end,:));
        ran = randi([1 size(u,1)], 1, 1);
        v = [v; u(ran,:)];
        d = disturbance(w);
        [t phi] = ode45(@sysODE, [0 tau], x(end,:), [], u(ran,:), d);
        x = [x; phi];
        
        j = j + 1;
    end
    savefig('scots/simulation')
end

end

function d = disturbance(w)
d = -w + (2 * w .* rand(size(w)));
end

function dxdt = sysODE(t,x,u,d)
dxdt = zeros(2,1);
dxdt(1) = u(1);
dxdt(2) = u(2);
dxdt = dxdt + d';
end
