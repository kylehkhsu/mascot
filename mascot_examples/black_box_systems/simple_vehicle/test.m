% tests
% openfig('Figures/system')
T0 = SymbolicSet('T/T1.bdd');
% T1 = SymbolicSet('testing/Tc0.bdd');
C1 = T0.points;
% C1 = T1.points;
C0 = T0_init;

range0 = [2 2.5; 2.5 3];
D0 = zeros(size(C0));
D1 = zeros(size(C1));

index = 1;
for ii=1:size(D0,1)
    if (C0(ii,1)>=range0(1,1) && C0(ii,1)<=range0(1,2) ...
            && C0(ii,2)>=range0(2,1) && C0(ii,2)<=range0(2,2))
        D0(index,:) = C0(ii,:);
        index = index+1;
    end
end
D0 = unique(D0(1:index-1,5:6),'rows');
plot(D0(:,1),D0(:,2),'bx');
pause
hold on

% range1 = [1.5 3; 1.5 3];
range1 = range0;
index = 1;
for ii=1:size(D1,1)
    if (C1(ii,1)>=range0(1,1) && C1(ii,1)<=range0(1,2) ...
            && C1(ii,2)>=range0(2,1) && C1(ii,2)<=range0(2,2))
        D1(index,:) = C1(ii,:);
        index = index+1;
    end
end
D1 = unique(D1(1:index-1,5:6),'rows');

plot(D1(:,1),D1(:,2),'rx');