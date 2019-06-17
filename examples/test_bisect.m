node = [0,0; 1,0; 1,1; 0,1];
elem = [2,3,1; 4,1,3];
figure(1); subplot(1,3,1); showmesh(node,elem);
findelem(node,elem);  % plot indices of all triangles
findnode(node);       % plot indices of all vertices

[node,elem,bdFlag,HB,tree] = bisect(node,elem,'all');
figure(1); subplot(1,3,2); showmesh(node,elem);
findelem(node,elem);  % plot indices of all triangles
findnode(node);       % plot indices of all vertices

bdFlag = setboundary(node,elem,'Dirichlet','all','Neumann','y==1');
[node,elem,bdFlag,HB,tree] = bisect(node,elem,[1 4],bdFlag);
figure(1); subplot(1,3,3); showmesh(node,elem);
findelem(node,elem);  % plot indices of all triangles
findnode(node);       % plot indices of all vertices