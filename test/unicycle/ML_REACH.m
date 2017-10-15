h = openfig('ML_REACH');
set(gcf,'position',[0 0 1000 1000])
set(gca,'fontsize',24)
set(gca,'linewidth',0.7)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];