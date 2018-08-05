function [] = PlotCells(X, eta, varargin)
%     plotCells(X, 'facecolor', 'none', 'edgec', [0.8 0.8 0.8], 'linew', 0.1);
    n=length(X(:,1));
    eh=eta./2;
    dim=size(X,2);
    switch dim
        case 2
            for i=1:n
              x=X(i,1);
              xdata=x+[-1 1 1 -1]*eh(1);
              y=X(i,2);
              ydata=y+[-1 -1 1 1]*eh(2);
              v=[xdata(:) ydata(:)];
              f=[1 2 3 4];
              patch('faces',f,'vertices',v,varargin{:});
            end
            drawnow
        case 3
            n=length(X(:,1));
            eh=eta./2;
            for i=1:n
              x=X(i,1);
              xdata=x+[-1 1 1 -1]*eh(1);
              xdata=[xdata xdata];
              y=X(i,2);
              ydata=y+[-1 -1 1 1]*eh(2);
              ydata=[ydata ydata];
              z=X(i,3);
              zdata=z+[-1 -1 -1 -1]*eh(3);
              zdata=[zdata zdata+2*eh(3)];
              v=[xdata(:) ydata(:) zdata(:)];
              f=[1 2 3 4;
                 1 2 6 5;
                 2 3 7 6;
                 3 4 8 7;
                 4 1 5 8;
                 5 6 7 8];
              patch('faces',f,'vertices',v,varargin{:})
            end
    end
end