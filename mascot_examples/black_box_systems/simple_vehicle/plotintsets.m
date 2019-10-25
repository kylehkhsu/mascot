figure
hold on
for ii=1:3
    set = SymbolicSet(['IntSets/S' num2str(ii) '.bdd']);
    p = set.points;
    plot(p(:,1),p(:,2),'.');
end