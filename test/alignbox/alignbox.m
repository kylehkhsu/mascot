function alignbox
  figure
  hold on
  box
  
  addpath(genpath('../..'));
  colors=get(groot,'DefaultAxesColorOrder'); 
  
  X = SymbolicSet('X.bdd');
  lb = X.first();
  ub = X.last(); 
  axis([lb(1)-1 ub(1)+1 lb(2)-1 ub(2)+1])
  plotCells(X, 'facecolor', 'none', 'edgec', [0.8 0.8 0.8], 'linew', 0.1);
  
  R = SymbolicSet('P.bdd');
  plotCells(R, 'facecolor', colors(1,:)*0.5+0.5, 'edgec', colors(1,:), 'linew', 0.1)
  
  E = SymbolicSet('E.bdd');
  plotCells(E, 'facecolor', colors(2,:)*0.5+0.5, 'edgec', colors(2,:), 'linew', 0.1)
    
    
end