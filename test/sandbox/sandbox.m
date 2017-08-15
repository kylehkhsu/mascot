function sandbox(mode)
  addpath(genpath('/home/kylehsu/control/SCOTS+Adaptive'));
  colors=get(groot,'DefaultAxesColorOrder'); 
  
  if (mode == 'T')
    openfig('system');
    hold on
    drawnow
    
    Zc = SymbolicSet('test/Zc.bdd');
    Zf = SymbolicSet('test/Zf.bdd');
    
    plotCells(Zf, 'facecolor', colors(3,:)*0.5+0.5, 'edgec', colors(3,:), 'linew', 0.1)
    drawnow
    pause
    plotCells(Zc, 'facecolor', colors(2,:)*0.5+0.5, 'edgec', colors(2,:), 'linew', 0.1)
    
    
end