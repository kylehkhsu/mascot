addpath(genpath('../../../mfiles'));
% plot different controller domains
name = 'vehicle';
% [sDIM,~,W_lb,W_ub,lb,ub,tau] = readParams('input.hh');
sDIM = 2;
W_lb = [0 0];
W_ub = [0 0];
lb = [0 0];
ub = [5 5];
tau = 0.4;
numAbs = 4;
ab = 1;

addpath(genpath('../../../'));
C = SymbolicSet(['C/C' int2str(ab) '.bdd']);
O = SymbolicSet(['O/O' int2str(ab) '.bdd']);
G = SymbolicSet(['G/G' int2str(ab) '.bdd']);

figure;
hold on;
axis([lb(1) ub(1) lb(2) ub(2)]);
% try
%     pO = O.points;
%     nO = size(pO,1);
%     plot(pO(:,1),pO(:,2),'rx');
% catch ME
%     warning('Something wrong while reading points in the under-approximation controller domain');
% end
% 
% try
%     pC = C.points;
%     no = size(pC,1);
%     plot(pC(:,1),pC(:,2),'bx');
% catch ME
%     warning('Something wrong while reading points in the over-approximation controller domain');
% end
% 
% try
%     pG = G.points;
%     nw = size(pG,1)
%     plot(pG(:,1),pG(:,2),'gx');
% catch ME 
%     warning('Something wrong while reading points in the worst-case controller domain');
% end

% if (strcmp(mode,'Int_dom'))
        iter = 4; % from the directory InterimGoal
        green = [0.7569    0.8667    0.7765];
		purple = [0.8196    0.6549    0.8471];
		orange = [0.9137    0.8275    0.3804];
        cmap = [green;purple;orange];
        cmap = [cmap;cmap;cmap;cmap];
%         openfig('problem');
        hold on;
        h = cell(3,1);
        for ii=numAbs:-1:1
            disp(['Abs = ' int2str(ii)]);
%             for jj=1:iter
                disp(['Iter = ' int2str(jj)]);
                try
                    Z = SymbolicSet(['T/T' num2str(ii) '.bdd'],'projection',[1 2]);
                    h{jj+1} = plotCells(Z,'facecolor',cmap(ii+1,:),'edgec',cmap(ii+1,:),'linew', 0.1);
                catch
                    warning(['There are probably no points in Z' num2str(ii)]);
                end
                drawnow
                pause
%             end 
%             for jj=1:numAbs
%                 delete(h{jj});
%             end
%             pause(0.01)
        end
%     end


% tu = SymbolicSet('transitions_sure.bdd');
% to = SymbolicSet('transitions_maybe.bdd');
% 
% intsets = cell(2,1);
% for ii=1:size(intsets,1)
%     intsets{ii}=SymbolicSet(['IntSets/S' num2str(ii) '.bdd']);
% end