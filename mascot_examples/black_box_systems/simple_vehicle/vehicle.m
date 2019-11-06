function vehicle (mode, numAbs, controllers, spec)

% mascot path
addpath(genpath('../../../'));


% parameters
sDIM = 2;
W_lb = [0 0];
W_ub = [0 0];
lb = [0 0];
ub = [5 5];
tau = 0.4;

% colors
green = [0.7569    0.8667    0.7765];
purple = [0.8196    0.6549    0.8471];
orange = [0.9137    0.8275    0.3804];
c = [green;purple;orange];
cmap = [];
for ii=1:ceil(max(controllers,numAbs)/3)
    cmap = [cmap;c];
end

switch mode
    case 'S'
        % plot state space and obstacles
%         figure('Units','inches',...
%             'OuterPosition',[0 0 3 3.1],...
%             'PaperPositionMode','auto');
        figure;
        hold on
        box on
        drawnow

        % load and draw state space
        X = SymbolicSet('plotting/X.bdd');
        lb = X.first();
        ub = X.last();
        axis([0.99*lb(1) 1.01*ub(1) 0.99*lb(2) 1.01*ub(2)]);
        plotCells(X, 'facecolor', 'none', 'edgec', [0.8 0.8 0.8], 'linew', 0.1);
        drawnow;
        disp('Done plotting state space')

        savefig('Figures/system');
        
    case 'P'
        openfig('Figures/system')
        axis([0.99*lb(1) 1.01*ub(1) 0.99*lb(2) 1.01*ub(2)]);
        try
            O = SymbolicSet('plotting/O.bdd');
            plotCells(O,'fast','facecolor','b');
            disp('Done plotting obstacles')
        catch
            warning('Something wrong in obstacle set (possibly no points).\n');
        end
        
        try
            E = SymbolicSet('plotting/E.bdd');
            plotCells(E,'fast','facecolor','b');
        catch
            warning('Something wrong in exclusion set (possibly no points).\n');
        end
        disp('Done plotting exclusion zones')
        G = SymbolicSet('plotting/G.bdd');
        if (spec==0)
            plotCells(G,'fast','facecolor','g');
            disp('Done plotting goals')
        else
            lb1=min(G.points(:,1))-spec;
            lb2=min(G.points(:,2))-spec;
            ub1=max(G.points(:,1))+spec;
            ub2=max(G.points(:,2))+spec;
            rectangle('Position',[lb1 lb2 (ub1-lb1) (ub2-lb2)],'Facecolor',[0 .5 .5]);
        end
        X0 = SymbolicSet('plotting/X0.bdd');
        plotCells(X0,'fast','facecolor','y');
        disp('Done plotting initial states')
        savefig('Figures/problem');
        
    case 'O'
        openfig('Figures/system')
        axis([0.99*lb(1) 1.01*ub(1) 0.99*lb(2) 1.01*ub(2)]);
        for ii=2:numAbs
            Z = SymbolicSet(['O/O' num2str(ii) '.bdd']);
            plotCells(Z,'fast','facecolor',cmap(ii,:));
            pause
        end
    case 'E'
        openfig('Figures/system')
        axis([0.99*lb(1) 1.01*ub(1) 0.99*lb(2) 1.01*ub(2)]);
        for ii=2:numAbs
            Z = SymbolicSet(['E/E' num2str(ii) '.bdd']);
            plotCells(Z,'fast','facecolor',cmap(ii,:));
            pause
        end
        
    case 'Z'
        openfig('Figures/problem')
        axis([0.99*lb(1) 1.01*ub(1) 0.99*lb(2) 1.01*ub(2)]);
        for ii=2:controllers
            Z = SymbolicSet(['Z/Z' num2str(ii) '.bdd']);
            plotCells(Z,'fast','facecolor',cmap(ii,:));
            pause
        end
        
    case 'Cdom'
        openfig('Figures/problem')
        axis([0.99*lb(1) 1.01*ub(1) 0.99*lb(2) 1.01*ub(2)]);
        for ii=1:controllers
            Cdom = SymbolicSet(['C/C' num2str(ii) '.bdd'],'projection',[1 2]);
            plotCells(Cdom,'fast','facecolor',cmap(ii,:));
            pause
        end
    case 'T'
%         openfig('Figures/system')
        axis([0.99*lb(1) 1.01*ub(1) 0.99*lb(2) 1.01*ub(2)]);
        for ii=1:numAbs
            try
                T = SymbolicSet(['T_backup_2000_spec=0_3/T' num2str(ii) '.bdd'],'projection',[1 2]);
                plotCells(T,'fast','facecolor',cmap(ii,:));
                pause
            catch
                disp('Something wrong with the SymbolicSet (possibly empty)\n');
            end
        end
    case 'D'
        openfig('Figures/system')
        axis([0.99*lb(1) 1.01*ub(1) 0.99*lb(2) 1.01*ub(2)]);
        for ii=1:numAbs
            D = SymbolicSet(['D/D' num2str(ii) '.bdd']);
            plotCells(D,'facecolor',cmap(ii,:));
            pause
        end
    case 'R'
        openfig('problem');
        hold on
        drawnow

        %     I = SymbolicSet('plotting/I.bdd');
        %     x = I.points();
        %     x = x(1,:);
        x = [0.5 10.9 -pi/2];    
        v = [];    
        j = 1;

        for i = controllers:-1:1
            disp(['iteration: ' int2str(i)])

            C = SymbolicSet(['C/C' int2str(i) '.bdd'], 'projection', [1 2 3]);
            if (i == 1)
                G = SymbolicSet(['G/G' int2str(numAbs) '.bdd']);
            else
                G = SymbolicSet(['Z/Z' int2str(i-1) '.bdd']);
            end

            points = G.points();
            Gset = plot(points(:,1), points(:,2), 'o', 'MarkerFaceColor', colors(mod(i,7)+1,:)*0.3+0.3, 'MarkerSize', 2);

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
                    delete(Gset);
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
end
