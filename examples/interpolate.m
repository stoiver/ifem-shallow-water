[node,elem] = squaremesh([0 1 0 1],1/2);
p = 1:size(elem,1);
figure(1);
subplot(1,6,1); showsolution(node,elem,p,'EdgeColor','k'); view(2);
[node,elem,~,~,tree] = bisect(node,elem,[1 2]);
p = eleminterpolate(p,tree);
subplot(1,6,2); showsolution(node,elem,p,'EdgeColor','k'); view(2);
[node,elem,~,~,tree] = bisect(node,elem,[1 2]);
p = eleminterpolate(p,tree);
subplot(1,6,3); showsolution(node,elem,p,'EdgeColor','k'); view(2);
[node,elem,~,~,tree] = coarsen(node,elem,'all');
size(p)
p = eleminterpolate(p,tree);
size(p)
subplot(1,6,4); showsolution(node,elem,p,'EdgeColor','k'); view(2);
[node,elem,~,~,tree] = coarsen(node,elem,'all');
size(p)
p = eleminterpolate(p,tree);
size(p)
subplot(1,6,5); showsolution(node,elem,p,'EdgeColor','k'); view(2);
[node,elem,~,~,tree] = coarsen(node,elem,'all');
size(p)
p = eleminterpolate(p,tree);
size(p)
subplot(1,6,6); showsolution(node,elem,p,'EdgeColor','k'); view(2);