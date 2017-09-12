function sandbox(mode)
  addpath(genpath('/home/kylehsu/control/SCOTS+Adaptive'));
  colors=get(groot,'DefaultAxesColorOrder'); 
  figure, hold on
    X = SymbolicSet('X.bdd');
    plotCells(X, 'facecolor', 'none', 'edgec', [0.8 0.8 0.8], 'linew', 0.1);
    hold on
    drawnow
    
    Zf = SymbolicSet('Zf.bdd');
    plotCells(Zf, 'facecolor', colors(2,:)*0.5+0.5, 'edgec', colors(2,:), 'linew', 0.1)
	drawnow
	pause

	Zc = SymbolicSet('Zc.bdd');
    plotCells(Zc, 'facecolor', colors(4,:)*0.5+0.5, 'edgec', colors(4,:), 'linew', 0.1)
end
