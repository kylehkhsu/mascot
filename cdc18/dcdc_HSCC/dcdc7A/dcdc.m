
function dcdc (mode, numAbs, controllers, progression)
w = [0.001 0.001];
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
    axis([0.95*lb(1) 1.05*ub(1) 0.99*lb(2) 1.01*ub(2)])
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
    cmap = [jet(numAbs) 0.5*ones(numAbs,1)];
    openfig('system')
    title(['T' num2str(numAbs)])
    for i = 1:numAbs
        T = SymbolicSet(['T/T' num2str(i) '.bdd'], 'projection', [1 2]);
        p = T.points;
        plot(p(:,1),p(:,2),'x','color',cmap(i,:))
        disp(['layer ' num2str(i)])
        pause
    end    
end

if (strcmp(mode,'Z'))
    openfig('problem');
    for ii=1:numAbs
        Z = SymbolicSet(['Z/Z' int2str(ii) '.bdd']);
        p = Z.points;
        x = plot(p(:,1),p(:,2),'ko');
        pause
        delete(x)
    end
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
    openfig('system');
    title(['C' num2str(numAbs)])
    hold on;
    cmap = [jet(numAbs) 0.5*ones(numAbs,1)];  
 
    for i=1:numAbs
        disp(['layer ' num2str(i)])
        C = SymbolicSet(['C/C' int2str(i) '.bdd'],'projection',[1 2]);
        try
            p = C.points;
        catch
            warning(['No points in C' int2str(i) '.bdd']);
            continue
        end        
        plot(p(:,1),p(:,2),'x','color',cmap(i,:))
        pause
    end
end

if (strcmp(mode,'Safe'))
    openfig('system');
    hold on
    drawnow
    
    x = [1.175 5.475];
%     x = [1.40700000000000,5.48300000000000];
    v = [];    
    j = 1;
    tauSet = [6.4 3.2 1.6 0.8 0.4 0.2 0.1];
    
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
%         try
            u = C{ii}.getInputs(x(end,:));
%         catch
%             debug = 1;
%         end
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
        disp(ii)
%         disp('u')
%         disp(u(ran,:))
%         disp('d')
%         disp(d)
%         disp('x')
        disp(x(end,:))
%         pause
       
    end

    %      plotCells(G,'facecolor',colors(2,:)*0.5+0.5,'edgec',colors(2,:),'linew',.1)
    savefig('simulation');
end

end

function d = disturbance(w)
d = -w + (2 * w .* rand(size(w)));
end

function dxdt = sysODE(t,x,u, d)
    r0=1.0;
  vs = 1.0 ;
  rl = 0.05 ;
  rc = rl / 10 ;
  xl = 3.0 ;
  xc = 70.0 ;

  if(u==0.5)
    A = [ -rl / xl  0 ;  0  (-1 / xc) * (1 / (r0 + rc)) ] ;
  else
    A = [ (-1 / xl) * (rl + ((r0 * rc) / (r0 + rc)))  ((-1 / xl) * (r0 / (r0 + rc))) / 5 ;
          5 * (r0 / (r0 + rc)) * (1 / xc)  (-1 / xc) * (1 / (r0 + rc)) ];
  end
  b = [(vs / xl) ; 0 ] ;
    dxdt = A*x + b + d';
end
