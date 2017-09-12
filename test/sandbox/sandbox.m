function sandbox(mode)
  addpath(genpath('/home/kylehsu/control/SCOTS+Adaptive'));
  colors=get(groot,'DefaultAxesColorOrder'); 
  figure, hold on
    Xf = SymbolicSet('Xf.bdd');
    plotCells(Xf, 'facecolor', 'none', 'edgec', [0.8 0.8 0.8], 'linew', 0.1);
    hold on
    drawnow
    
    Zf = SymbolicSet('Zf.bdd');
    plotCells(Zf, 'facecolor', colors(2,:)*0.5+0.5, 'edgec', colors(2,:), 'linew', 0.1)
	drawnow
	pause
    
    Zf2 = SymbolicSet('Zf2.bdd');
    plotCells(Zf2, 'facecolor', colors(3,:)*0.5+0.5, 'edgec', colors(3,:), 'linew', 0.1)
	drawnow
	pause

	Zc = SymbolicSet('Zc.bdd');
    plotCells(Zc, 'facecolor', colors(4,:)*0.5+0.5, 'edgec', colors(4,:), 'linew', 0.1)
end
