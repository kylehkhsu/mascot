function dcdcPlotNice(mode, numAbs) % numAbs = num of controllers
    w = [0.001 0.001];
    addpath(genpath('../../..'));
    colors=get(groot,'DefaultAxesColorOrder');
    
    if (strcmp(mode, 'S'))
        figure
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
        savefig('problem'); % for dcdc, they are same
    end
    
    if (strcmp(mode,'T&C'))
        figure('Units','inches',...
            'OuterPosition',[0 0 8 3.1],...
            'PaperPositionMode','auto');
        [ha, pos] = tight_subplot(3, 2, 0.07, 0.06, 0.06);
        set(ha, 'XTickLabelMode', 'auto') 
        set(ha, 'YTickLabelMode', 'auto')
        cmap = colormap('gray');
        color1 = cmap(1,:);
        color2 = cmap(floor(size(cmap,1)/2),:);
        red = [0.9373    0.2980    0.2980];
		gray = [0.2471    0.2471    0.2471];
        for ii=1:numAbs
            disp(ii)
            axes(ha(ii));
            hold on;
            T = SymbolicSet(['T/T' num2str(ii) '.bdd'],'projection',[1 2]);
            C = SymbolicSet(['C/C' num2str(ii) '.bdd'],'projection',[1 2]);
            try
%                 s = subplot(3,numAbs/3,ii); hold on;
%                 plotCells(T,'facecolor',red,'EdgeColor',red);
%                 plotCells(C,'facecolor',gray,'EdgeColor',gray);
                
                pT = T.points;
                pC = C.points;
                plot(pT(:,1),pT(:,2),'rs', 'MarkerSize', 10/ii);
                plot(pC(:,1),pC(:,2),'ks', 'MarkerSize', 10/ii);
%                 set(s,...
%                     'Units','normalized',...
%                     'OuterPosition',[0 0 1 1],...
%                     'YTick',0:2:12,...
%                     'XTick',0:2:12,...
%                     'FontUnits','points',...
%                     'FontWeight','normal',...
%                     'FontSize',9,...
%                     'FontName','Times')
                axis([1.1 1.6 5.4 5.9]);
                title({['$l_' num2str(numAbs + 1 - ii) '$']},...
                    'FontUnits','points',...
                        'interpreter','latex',...
                        'FontSize',14,...
                        'FontName','Times');
                ylabel({'$x_2$'},...
                        'FontUnits','points',...
                        'interpreter','latex',...
                        'FontSize',10,...
                        'FontName','Times');%,...
        %                 'Units','Normalized',...
        %                 'Position',[-0.1,0.5,0])
                xlabel({'$x_1$'},...
                        'FontUnits','points',...
                        'interpreter','latex',...
                        'FontSize',10,...
                        'FontName','Times');                                              
                box on;
            catch
                warning(['There are no points in T or C' num2str(ii)]);
            end
        end
        
        savefig('dcdc6A_cdc');
        print -depsc dcdc6A_cdc.eps
    end
    
    if (strcmp(mode,'Controller domain'))
        for ii=1:numAbs
            figure
            C = SymbolicSet(['C/C' num2str(ii) '.bdd'],'projection',[1 2]);
            try
                plotCells(C,'facecolor','k');
            catch
                warning(['There are no points in C' num2str(ii)]);
            end
        end
    end
    
    if (strcmp(mode,'Safe'))
        openfig('problem');
        hold on
        drawnow

        x = [1.18250000000000,5.81350000000000];
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
                    tau = Z{ii}.tau;
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