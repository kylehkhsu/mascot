function unicycle()

  figure
  hold on
  colors=get(groot,'DefaultAxesColorOrder');
  
  % load and draw state space
  set = SymbolicSet('X.bdd','projection',[1 2]);
  plotCells(set,'facecolor','none','edgec',[0.8 0.8 0.8],'linew',.1)
  drawnow
  disp('Done plotting state space')

  % load and draw the symbolic set containing obstacles
  set = SymbolicSet('X_obst.bdd','projection',[1 2]);
  plotCells(set,'facecolor',colors(1,:)*0.5+0.5,'edgec',colors(1,:),'linew',.1)
  drawnow
  disp('Done plotting obstacles')

  % plot the real obstacles 
  plot_domain
  drawnow
  disp('Done plotting domain')

  savefig('unicycle');
end

function plot_domain

colors=get(groot,'DefaultAxesColorOrder');


v=[1     0  ;1.2  0   ; 1     9    ; 1.2 9   ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[2.2   0  ;2.4  0   ; 2.2   5    ; 2.4 5   ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[2.2   6  ;2.4  6   ; 2.2   10   ; 2.4 10  ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[3.4   0  ;3.6  0   ; 3.4   9    ; 3.6 9   ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[4.6   1  ;4.8  1   ; 4.6   10   ; 4.8 10  ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[5.8   0  ;6    0   ; 5.8   6    ; 6   6   ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[5.8   7  ;6    7   ; 5.8   10   ; 6   10  ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[7     1  ;7.2  1   ; 7     10   ; 7.2 10  ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[8.2   0  ;8.4  0   ; 8.2   8.5  ; 8.4 8.5 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[8.4   8.3;9.3  8.3 ; 8.4   8.5  ; 9.3 8.5 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[9.3   7.1;10   7.1 ; 9.3   7.3  ; 10  7.3 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[8.4   5.9;9.3  5.9 ; 8.4   6.1  ; 9.3 6.1 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[9.3   4.7;10   4.7 ; 9.3   4.9  ; 10  4.9 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[8.4   3.5;9.3  3.5 ; 8.4   3.7  ; 9.3 3.7 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[9.3   2.3;10   2.3 ; 9.3   2.5  ; 10  2.5 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));


end