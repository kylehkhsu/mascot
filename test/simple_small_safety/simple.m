
function simple (mode, numAbs, controllers, progression)
w = [0];
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
    X = SymbolicSet('plotting/X.bdd');
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
        SafeSet = SymbolicSet('plotting/SafeInner.bdd');
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
    openfig('problem')
    for ii=2:2
        T1 = SymbolicSet(['Zs' int2str(ii) '/Zs_1.bdd']);
        p1 = T1.points;
        x = plot(p1(:,1),p1(:,2),'ko');
        pause
        delete(x)

        T2 = SymbolicSet(['Zs' int2str(ii) '/Zs_2.bdd']);
        p2 = T2.points;
        x = plot(p2(:,1),p2(:,2),'bo');
        pause    
        delete(x)

        T3 = SymbolicSet(['Zs' int2str(ii) '/Zs_3.bdd']);
        p3 = T3.points;
        plot(p3(:,1),p3(:,2),'ro');
        delete(x)
    end
    title('Zs')
end

if (strcmp(mode, 'D'))
    openfig('problem')
    for ii=1:4
        try
            D1 = SymbolicSet(['D1/' int2str(ii) '.bdd']);
            p1 = D1.points;
            x = plot(p1(:,1),p1(:,2),'ko');
            pause
            delete(x)
        catch
            warning(['No points in D1/' int2str(ii) '.bdd']);
        end
    
        try
            D2 = SymbolicSet(['D2/' int2str(ii) '.bdd']);
            p2 = D2.points;
            x = plot(p2(:,1),p2(:,2),'ko');
            pause
            delete(x)
        catch
            warning(['No points in D2/' int2str(ii) '.bdd']);
        end
    end
    title('D')
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
        C = SymbolicSet(['C/C' int2str(ii) '.bdd'],'projection',[1 2]);
        eta = C.eta(1:2);
        for jj=1:numAbs
            if isequal(eta',etaAbs(jj,1:2))
                layer = jj;
                break;
            end
        end
        
        
        if (ii == 1)
            Z = SymbolicSet(['SafeSets/S' int2str(numAbs) '.bdd']);
        else
            Z = SymbolicSet(['Z/Z' int2str(ii-1) '.bdd']);
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