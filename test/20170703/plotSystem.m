function plotSystem()
  addpath(genpath('/home/kylehsu/control/SCOTS+Adaptive'));
  figure
  hold on
  drawnow
  colors=get(groot,'DefaultAxesColorOrder');
  
  % load and draw state space
  set = SymbolicSet('X.bdd');
  plotCells(set,'facecolor','none','edgec',[0.8 0.8 0.8],'linew',.1)
  drawnow
  disp('Done plotting state space')
  
  % load and draw the symbolic set containing obstacles
  try
    set = SymbolicSet('X_obst.bdd');
    plotCells(set,'facecolor',colors(1,:)*0.5+0.5,'edgec',colors(1,:),'linew',.1)
    drawnow
    disp('Done plotting obstacles')
  catch
    warning('No obstacles');
  end


  savefig('system');
end




