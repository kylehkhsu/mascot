
function unicycle (mode, numAbs, controllers, progression)
w = [0.05 0.05 0];
addpath(genpath('../..'));
addpath(genpath('~/ownCloud/C++/SCOTS_modified/mfiles/'));

% colors
colors=get(groot,'DefaultAxesColorOrder');
if (strcmp(mode, 'CZ'))
    % visualize progression. black = sub-controller domain, red = sub-target
    openfig('problem')
    
    if strcmp(progression, 'b')
        for ii=1:controllers
            disp(num2str(ii))
            if ii ~= 1
                Z = SymbolicSet(['Z/Z' int2str(ii-1) '.bdd']);
                p = Z.points;
                Zset = plot(p(:,1),p(:,2),'ro');
                pause
            end
            
            C = SymbolicSet(['C/C' int2str(ii) '.bdd']);
            p = C.points;
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
            C = SymbolicSet(['C/C' int2str(ii) '.bdd']);
            p = C.points;
            %         plotColor = colors(mod(ii,7)+1,:)*0.3+0.3;
            Cset = plot(p(:,1),p(:,2),'ko');
            pause
            
            if ii ~= 1
                Z = SymbolicSet(['Z/Z' int2str(ii-1) '.bdd']);
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


if (strcmp(mode,'Q')) % plot controller domains
    openfig('problem');
    hold on;
    cmap = [hot(numAbs) 0.5*ones(numAbs,1)];
    
    %       if numAbs~=size(etaAbs,1) % Sanity check
    %           error();
    %       end
    etaAbs = [0.6 0.6; 0.3 0.3; 0.15 0.15];
    
    for ii=controllers:-1:1
        ii
        C = SymbolicSet(['C/C' int2str(ii) '.bdd'], 'projection', [1 2]);
        eta = C.eta(1:2);
        for jj=1:numAbs
            if isequal(eta',etaAbs(jj,1:2))
                layer = jj;
                break;
            end
        end
        
        
        if (ii == 1)
            Z = SymbolicSet(['SafeSets/S' int2str(numAbs) '.bdd'], 'projection', [1 2]);
        else
            Z = SymbolicSet(['Z/Z' int2str(ii-1) '.bdd'], 'projection', [1 2]);
        end
        
        points = setdiff(C.points,Z.points,'rows');
        for jj=1:size(points,1)
            rectangle('Position',[points(jj,:)-eta' 2*eta'],'FaceColor',cmap(layer,:));
            %alpha(s,0.5);
        end
        pause;
    end
end

if (strcmp(mode,'Safe'))
    openfig('problem');
    hold on
    drawnow
    
    %     I = SymbolicSet('plotting/I.bdd');
    %     x = I.points();
    %     x = x(1,:);
    x = [1 0.5 0];    
    v = [];    
    j = 1;
    tauSet = [0.9;0.9/2;0.9/2/2];
    
    C = cell(numAbs,1);
    Z = cell(numAbs,1);
    for ii=1:numAbs
        C{ii} = SymbolicSet(['C/C' int2str(ii) '.bdd']);
        Z{ii} = SymbolicSet(['Z/Z' int2str(ii) '.bdd']);
    end
    
    for rounds=1:400
        for ii=1:numAbs
            if (Z{ii}.isElement(x(end,:)))
                tau = tauSet(ii);
                if (ii==3)
                    a = 1;
                end
                break;
            end
        end
        try
            u = C{ii}.getInputs(x(end,:));
        catch
            debug = 1;
        end
        ran = randi([1 size(u,1)], 1, 1);
        v = [v; u(ran,:)];
        d = disturbance(w);
        [t phi] = ode45(@sysODE, [0 tau], x(end,:), [], u(ran,:), d);
        x = [x; phi];

%         if (mod(j,1) == 0)
            plot(x(:,1),x(:,2),'k.-')
            drawnow
%             pause
%         end
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

function dxdt = sysODE(t,x,u, d)
dxdt = zeros(3,1);
dxdt(1)=u(1)*cos(x(3));
dxdt(2)=u(1)*sin(x(3));
dxdt(3)=u(2);
dxdt = dxdt + d';
end
