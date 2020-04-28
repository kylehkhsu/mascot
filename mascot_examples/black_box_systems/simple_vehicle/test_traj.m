lb = [0 0];
ub = [5 5];
sys_traj = dlmread("Figures/sys_traj.txt");
abs_traj = dlmread("Figures/abs_traj.txt");
% openfig('problem')
vehicle('P',8,3,0.3)
hold on;
% vehicle('Cdom',6,1)
axis([0.99*lb(1) 1.01*ub(1) 0.99*lb(2) 1.01*ub(2)]);
plot(sys_traj(:,1),sys_traj(:,2),'LineWidth',4);
% plot(abs_traj(:,1),abs_traj(:,2),'ko','Markersize',4);
% hold on;
% traj2 = dlmread("Figures/sys_traj.txt");
% plot(traj2(:,1),traj2(:,2),'r','LineWidth',4);
