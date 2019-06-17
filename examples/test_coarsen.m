load ../../ifem/mesh/meshdata/Lshapeadaptivemesh;
set(gcf,'Units','normal'); set(gcf,'Position',[0.25,0.25,0.7,0.4]);

subplot(1,3,1); showmesh(node,elem);
findelem(node,elem);  % plot indices of all triangles
findnode(node);       % plot indices of all vertices

[node,elem] = coarsen(node,elem,'all');
subplot(1,3,2); showmesh(node,elem);
findelem(node,elem);  % plot indices of all triangles
findnode(node);       % plot indices of all vertices

[node,elem] = coarsen(node,elem,'all');
subplot(1,3,3); showmesh(node,elem);
findelem(node,elem);  % plot indices of all triangles
findnode(node);       % plot indices of all vertices