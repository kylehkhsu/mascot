function spiral(mode, numAbs) % numAbs = num of controllers
    addpath(genpath('../..'));
    colors=get(groot,'DefaultAxesColorOrder');
    tauset = [0.05*2*2; 0.05*2; 0.05];
    
    if (strcmp(mode, 'S'))
%         figure
        figure('Units','inches',...
            'OuterPosition',[0 0 3 3.1],...
            'PaperPositionMode','auto');
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

        savefig('system');
        
        UnSafe = SymbolicSet('plotting/O.bdd');
        plotCells(UnSafe,'facecolor','b');
        savefig('problem'); % for dcdc, they are same
    end
    
    if (strcmp(mode,'T'))
        figure
        cmap = colormap('gray');
        color1 = cmap(1,:);
        color2 = cmap(floor(size(cmap,1)/2),:);
        for ii=1:numAbs
            T = SymbolicSet(['T/T' num2str(ii) '.bdd'],'projection',[1 2]);
            try
                subplot(1,numAbs,ii)
                plotCells(T,'facecolor','b','EdgeColor','b');
            catch
                warning(['There are no points in T' num2str(ii)]);
            end
        end
    end
    
    if (strcmp(mode,'T&C'))
        figure
        cmap = colormap('gray');
        color1 = cmap(1,:);
        color2 = cmap(floor(size(cmap,1)/2),:);
        for ii=1:numAbs
            T = SymbolicSet(['T/T' num2str(ii) '.bdd'],'projection',[1 2]);
            C = SymbolicSet(['C/C' num2str(ii) '.bdd'],'projection',[1 2]);
            try
                subplot(1,numAbs,ii)
                plotCells(T,'facecolor','b','EdgeColor','b');
                plotCells(C,'facecolor','r','EdgeColor','r');
            catch
                warning(['There are no points in T or C' num2str(ii)]);
            end
        end
    end
    
    if (strcmp(mode,'C_dom'))
        cmap = [jet(numAbs) 0.5*ones(numAbs,1)];
        openfig('problem');
        for ii=1:numAbs
            hold on
            C = SymbolicSet(['C/C' num2str(ii) '.bdd'],'projection',[1 2]);
            try
%                 plotCells(C,'facecolor','k');
                plot(C.points(:,1),C.points(:,2),'x','color',cmap(ii,:));
            catch
                warning(['There are no points in C' num2str(ii)]);
            end
            drawnow
            pause
        end
        set(gca,...
            'Units','normalized',...
            'OuterPosition',[0 0 1 1],...
            'FontUnits','points',...
            'FontWeight','normal',...
            'FontSize',9,...
            'FontName','Times')
        ylabel({'$x_2$'},...
                'FontUnits','points',...
                'interpreter','latex',...
                'FontSize',9,...
                'FontName','Times');%,...
%                 'Units','Normalized',...
%                 'Position',[-0.1,0.5,0])
        xlabel({'$x_1$'},...
                'FontUnits','points',...
                'interpreter','latex',...
                'FontSize',9,...
                'FontName','Times');
        savefig('spiral_cdc');
        print -depsc spiral_cdc.eps
    end
    
    if (strcmp(mode,'Safe'))
        openfig('problem');
        hold on
        drawnow

        x = [1.5,1.5];
    %     x = [1.40700000000000,5.48300000000000];
        v = [];    
        j = 1;
%         tauSet = [0.5*2*2*2;0.5*2*2;0.5*2;0.5];

        C = cell(numAbs,1);
        Z = cell(numAbs,1);
        for ii=1:numAbs
            C{ii} = SymbolicSet(['C/C' int2str(ii) '.bdd']);
            Z{ii} = SymbolicSet(['Z/Z' int2str(ii) '.bdd']);
        end

    %     for rounds=1:400
        while 1
            noController = 1;
            for ii=1:numAbs
                if (Z{ii}.isElement(x(end,:)))
%                     tau = Z{ii}.tau;
                    tau = tauset(ii);
                    noController = 0;
                    break;
                end
            end
            if noController==1
                error('The state is outside the multilayered controller domain.')
            end
    %         try
                u = C{ii}.getInputs(x(end,:));
    %         catch
    %             debug = 1;
    %         end
            ran = randi([1 size(u,1)], 1, 1);
%             u1 = u(ran,:);
            u1 = min(abs(u(:)));
            v = [v; u1];
%             d = disturbance(w);
            [th,r] = cart2pol(x(end,1),x(end,2));
          [t z]=ode45(@sysODE,[0 tau], [r th], odeset('abstol',1e-4,'reltol',1e-4),u1);
          [x1,x2] = pol2cart(z(:,2),z(:,1));

          x=[x; x1(end) x2(end)];
          plot(x1,x2);
          drawnow
    %             pause
    %         end
%             disp(ii);
            disp('u')
            disp(v(end,:))
    %         disp(u(ran,:))
    %         disp('d')
    %         disp(d)
    %         disp('x')
%             disp(x(end,:));
    %         pause

        end

        %      plotCells(G,'facecolor',colors(2,:)*0.5+0.5,'edgec',colors(2,:),'linew',.1)
        savefig('simulation');
    end
end

function dxdt = sysODE(t,x,u)
    % ODE in polar coordinate
    omega = 0.1;
    a = 0.5;
    dxdt = zeros(size(x));
    dxdt(1) = -a*x(1) + u;
    dxdt(2) = omega;
end