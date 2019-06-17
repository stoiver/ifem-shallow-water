function planar_dam
%% Shallow water Planar Dam Break
%--------------------------------------------------------------------------
% Sudi Mungkasi 2012, Australian National University
% Steve Roberts 2019, Australian National University
%--------------------------------------------------------------------------

%%
close all; clear all;

%---------------------- Parameters ----------------------------------------
figure(1); set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.8,0.4]);

t = 0;
g = 9.81;
%---------------------- Initial Grid --------------------------------------
node = [-1 -1; 1 -1; 1 1; -1 1];
elem = [2 3 1; 4 1 3];
%bdFlag = setboundary(node,elem,'Dirichlet','all');

for i = 1:1 % Try different level here. 7 is good.
    [node,elem] = uniformbisect(node,elem);
end

x = [node(elem(:,1),1), node(elem(:,2),1), node(elem(:,3),1)];
y = [node(elem(:,1),2), node(elem(:,2),2), node(elem(:,3),2)];


%-------------------------- Initial Condition -----------------------------
nt = size(elem,1);
C  = find_centroid(node,elem);
h  = zeros(1,nt);
u  = zeros(1,nt);
v  = zeros(1,nt);
for i=1:nt
    if C(i,1) < 0.0
        h(i) = 0.5;
        u(i) = 0.0;
        v(i) = 0.0;
    else
        h(i) = 0.2;
        u(i) = 0.0;
        v(i) = 0.0;
    end
end

Q = zeros(4,nt); %storage initiation for quantity
Q(1,:) = h;      %height
Q(2,:) = u.*h;   %xMomentum
Q(3,:) = v.*h;   %yMomentum
Q(4,:) = 0.5*h.*(u.^2+v.^2) + 0.5*g*h.^2; %entropy

Area      = find_area(node,elem);
tol_Area  = max(Area);

fprintf('The initial number of elements is %d triangles.\n', nt);
fprintf('The initial time is %d seconds.\n', t);


% Plot the initial conditions
showmesh(node,elem);
axis on;
findelem(node,elem);  % plot indices of all triangles
findnode(node);       % plot indices of all vertices


%% MAIN PROGRAM STARTS HERE %%%%%%%%%%%

maxRefineCoarsen = 2;
dt = 0.002;

for iterate=1:maxRefineCoarsen
    nt = size(elem,1); %number of triangles
    
    % Test evolve; Used to decide on new mesh
    [U, Area] = evolve_step(node,elem,Q,dt,g);

    % Begin Refinement ----------------------------------------------------
    h  = U(1,:);
    uh = U(2,:);
    u  = uh./h;
    vh = U(3,:);
    v  = vh./h;    
    Ent_after = 0.5*h.*(u.^2+v.^2) + 0.5*g*h.^2; %entropy
    
    NEP =(1/dt)*(Ent_after-U(4,:));%.*Area';

    tol_NEP = 0.25*max(abs(NEP));
    refineElem  = find(abs(NEP) > tol_NEP);%0.000001*tol_NEP);
    
    Elems_to_be_refined = refineElem
    
    for j = 1: size(refineElem,2)
        if Area(refineElem(1,j))/tol_Area < 0.000001+1/(2^maxRefineCoarsen) %0.126 %1.1 %0.6 %0.26
            Elems_to_be_refined = setdiff(Elems_to_be_refined, refineElem(1,j));
        end
    end
    
    % Do the bisection of the domain
    [node,elem, bdFlag, HB, tree] = bisect(node,elem, Elems_to_be_refined);
    
    % Update quantities after bisection
  
    if isempty(tree) == 0   % 0 means that brother is not empty
        
        nb = size(tree,1);     %number of original triangles refined
        ntr = size(elem,1);    %number of triangles after refinement
        QR = zeros(4,ntr);     %storage initiation for quantities on refined grid
        QR(:,1:nt) = Q;        %Copy over old quantity values
        NEPR = zeros(1,ntr);   %storage initiation for NEP on refined grids
        NEPR(1,1:nt) = NEP;    %Copy over old NEP values
        
        for ii=1:nb
            t1 = tree(ii,1);
            t2 = tree(ii,2);
            t3 = tree(ii,3);
            QR(:,t2) = Q(:,t1);
            QR(:,t3) = Q(:,t1);
            NEPR(1,t2) = NEP(:,t1);
            NEPR(1,t3) = NEP(:,t1);
        end
    else
        QR = Q;
        NEPR = NEP;
    end
    % End Refinement ------------------------------------------------------
    
    
    % Begin Coarsening ----------------------------------------------------
    Area = find_area(node,elem);  %areas of triangles
    coarsenElem = find(abs(NEPR) < tol_NEP);
    Elems_to_be_coarsened = coarsenElem;
    for j = 1: size(coarsenElem,2)
        Area_j = Area(coarsenElem(1,j));
        if Area_j/tol_Area >(2^maxRefineCoarsen)- 0.000001 %2.1%1.1 %0.9
            Elems_to_be_coarsened = setdiff(Elems_to_be_coarsened, Area_j);
        end
    end
    
    markedElem = Elems_to_be_coarsened;
    % Do coarsening and update quatities at the same time
    bdEdge = []; % This empty set bdEdge will return nothing.
    
    [node,elem,bdEdge,brother,QC,NEPC] = coarsen_sudi(node,elem, markedElem, bdEdge, QR, NEPR);
    
    %Collect quantities for next iteration from the refine grids
    Q = QC;
    NEP = NEPC;
    % End Coarsening ------------------------------------------------------
    
    
    % Actual Evolve -------------------------------------------------------
    
    [Q, Area] = evolve_step(node,elem,Q,dt);
    t = t+dt;
    
    
    
    %Plot in colour the results
    %clf(figure);
    clf;
    figure(1);
    h  = Q(1,:);
    uh = Q(2,:);
    u  = uh./h;
    vh = Q(3,:);
    v  = vh./h;
    Ent_after = 0.5*h.*(u.^2+v.^2) + 0.5*g*h.^2; %entropy
    NEP =(1/dt)*(Ent_after-U(4,:));%.*Area';
    %maxNEP = max(abs(NEP));
    %NEP = NEP/maxNEP;
    
    clf(figure(1));
    figure(1);
    subplot(1,2,1);
    showsolution(node,elem,h); view(3); zlim([0. 1.]);
    xlabel('x'); ylabel('y'); %zlabel('Stage')
    subplot(1,2,2);
    showsolution(node,elem,-NEP/max(abs(NEP))); view(3);
    xlabel('x'); ylabel('y'); pause(0.001);% zlabel('CK/max(|CK|)');
    
    clf(figure(2));
    figure(2);
    subplot(1,2,1);
    showsolution(node,elem,uh/max(abs(uh))); view(3); zlim([0. 1.]);
    xlabel('x'); ylabel('y'); %zlabel('uh/max(|uh|)')
    subplot(1,2,2);
    showmesh(node,elem); view(3); xlabel('x'); ylabel('y'); pause(0.001);
end
end 

%% Sub Functions

function [U, Area] = evolve_step(node,elem,Q,dt,g)

% Extract all quantities
nt = size(elem,1);
h  = Q(1,:);
uh = Q(2,:);
u  = uh./h;
vh = Q(3,:);
v  = vh./h;
Q(4,:) = 0.5*h.*(u.^2+v.^2) + 0.5*g*h.^2;

% Extract domain properties
T = auxstructure(elem);
Neighbor = T.neighbor;             %find_neighbor(elem); %neighbor
Normal   = find_normal(node,elem); %normals
Edgelength= find_edgelength(node,elem); %edgelength
Area      = find_area(node,elem);  %areas of triangles
U = Q;                             %initiation for temporary storage of the quantity

% Evolve
for j = 1:nt
    Ul = Q(:,j);
    flux = 0.0;
    for k = 1:3
        neighbor = Neighbor(j,k);
        normal   = Normal(:,k,j);
        edgelength = Edgelength(1,k,j);
        
        Ul_rot = rotate(Ul, normal(1), normal(2));
        Ur = Q(:,neighbor);
        % Use reflective boundary at the moving water
        if neighbor==j
            if Ur(3) ~= 0.0
                Ur(3) = -Ur(3);
            end
        end
        Ur_rot = rotate(Ur, normal(1), normal(2));
        
        edgeflux  = FluxKNP(Ul_rot,Ur_rot,g);
        edgeflux  = rotate(edgeflux, normal(1), -normal(2));
        flux      = flux - edgeflux*edgelength;
    end
    U(:,j) = Ul + dt*flux/Area(j);
end
end %% -------------------- End of Evolve Step ----------------------------


%%
function F = FluxKNP(Ul,Ur,g)
%dt=1000.0 must return F and dt
F=zeros(4,1);
h_l  = Ul(1,1);     %height
u_l  = Ul(2,1)/h_l; %x velocity
ct_l = Ul(3,1)/h_l; %y velocity, that is, Concentration of Tracer
e_l  = Ul(4,1);     %entropy
h_r  = Ur(1,1);
u_r  = Ur(2,1)/h_r;
ct_r = Ur(3,1)/h_r;
e_r  = Ur(4,1);
        
% Compute speeds in x-direction.
soundspeed_l = sqrt(g*h_l);
soundspeed_r = sqrt(g*h_r);
s_max = max([u_l+soundspeed_l, u_r+soundspeed_r, 0.0]);
s_min = min([u_l-soundspeed_l, u_r-soundspeed_r, 0.0]);
      
% Flux formulas
flux_l_h       = h_l*u_l;
flux_l_p       = h_l*u_l^2 + 0.5*g*h_l^2;
flux_l_voltrac = u_l*h_l*ct_l;
flux_r_h       = h_r*u_r;
flux_r_p       = h_r*u_r^2 + 0.5*g*h_r^2;
flux_r_voltrac = u_r*h_r*ct_r;

entropy_flux_l  = (0.5*h_l*(u_l^2 + ct_l^2) + g*h_l^2)*u_l;
entropy_flux_r  = (0.5*h_r*(u_r^2 + ct_r^2) + g*h_r^2)*u_r;
        
% Flux computation
denom = s_max-s_min;
edgeflux_h = (s_max*flux_l_h - s_min*flux_r_h + s_max*s_min*(h_r-h_l))/ denom;
edgeflux_p = (s_max*flux_l_p - s_min*flux_r_p + s_max*s_min*(flux_r_h-flux_l_h))/ denom;
edgeflux_voltrac = (s_max*flux_l_voltrac - s_min*flux_r_voltrac + s_max*s_min*(h_r*ct_r-h_l*ct_l))/ denom;
entropy_edgeflux = (s_max*entropy_flux_l - s_min*entropy_flux_r + s_max*s_min*(e_r-e_l))/denom;
max_speed = max(abs(s_max), abs(s_min));
        
% Update timestep
%dt = min(dt, CFL*0.5*dx/max_speed);

F(1,1) = edgeflux_h;
F(2,1) = edgeflux_p;
F(3,1) = edgeflux_voltrac;
F(4,1) = entropy_edgeflux;
%return F, dt
end

%%
function Q = rotate(q, n1, n2)
% Rotate the momentum component q, that is, q(2) and q(3)
% from x,y coordinates to coordinates based on normal vector (n1, n2).
% Result is returned in array 3x1 r
% To rotate in opposite direction, call rotate with (q, n1, -n2)
% Contents of q are changed by this function */

Q = q;
% Shorthands
q2 = q(2);  % uh momentum
q3 = q(3);  % vh momentum

% Rotate
Q(2) =  n1*q2 + n2*q3;
Q(3) = -n2*q2 + n1*q3;
end

function C = find_centroid(node,elem)
x = (node(elem(:,1),1) + node(elem(:,2),1) + node(elem(:,3),1))/3;
y = (node(elem(:,1),2) + node(elem(:,2),2) + node(elem(:,3),2))/3;
C = [x,y]; %space initiation
end

function A = find_area(node,elem)
x1 = node(elem(:,1),1);
x2 = node(elem(:,2),1);
x3 = node(elem(:,3),1);
y1 = node(elem(:,1),2);
y2 = node(elem(:,2),2);
y3 = node(elem(:,3),2);
A = 0.5* ( x1.*y2 + x2.*y3 + x3.*y1 - x1.*y3 - x2.*y1 - x3.*y2 );
end

%%
function n = find_normal(node,elem)

M = [ 0 1 ; -1 0 ];
nt = size(elem,1);

n = zeros(2,3,nt);
el = zeros(1,3,nt);
for i = 1:nt
n(:,1,i) = M*(node(elem(i,3),:) - node(elem(i,2),:))';
n(:,2,i) = M*(node(elem(i,1),:) - node(elem(i,3),:))';
n(:,3,i) = M*(node(elem(i,2),:) - node(elem(i,1),:))';

el(1,1,i) = sqrt(n(1,1,i)^2 + n(2,1,i)^2);
el(1,2,i) = sqrt(n(1,2,i)^2 + n(2,2,i)^2);
el(1,3,i) = sqrt(n(1,3,i)^2 + n(2,3,i)^2);

n(:,1,i) = n(:,1,i)/el(1,1,i) ;
n(:,2,i) = n(:,2,i)/el(1,2,i) ;
n(:,3,i) = n(:,3,i)/el(1,3,i) ;    
end
end

%%
function el = find_edgelength(node,elem)

M = [ 0 1 ; -1 0 ];
nt = size(elem,1);

n = zeros(2,3,nt);
el = zeros(1,3,nt);
for i = 1:nt
n(:,1,i) = M*(node(elem(i,3),:) - node(elem(i,2),:))';
n(:,2,i) = M*(node(elem(i,1),:) - node(elem(i,3),:))';
n(:,3,i) = M*(node(elem(i,2),:) - node(elem(i,1),:))';

el(1,1,i) = sqrt(n(1,1,i)^2 + n(2,1,i)^2);
el(1,2,i) = sqrt(n(1,2,i)^2 + n(2,2,i)^2);
el(1,3,i) = sqrt(n(1,3,i)^2 + n(2,3,i)^2);

n(:,1,i) = n(:,1,i)/el(1,1,i) ;
n(:,2,i) = n(:,2,i)/el(1,2,i) ;
n(:,3,i) = n(:,3,i)/el(1,3,i) ;    
end
end


%%
function [node,elem,bdEdge,brother,UC,NEPC] = coarsen_sudi(node,elem, markedElem, bdEdge, UR, NEPR)
%---------------------------------------------------------sudi
node_store = node;
elem_store = elem;
UC = UR;
NEPC = NEPR;
%---------------------------------------------------------sudi
%% COARSEN coarsen a 2-D triangulation.
%
% [node,elem] = coarsen(node,elem,markedElem) removes good-to-coarsen
% nodes whose star are marked for coarsening
%
% [node,elem,bdEdge] = coarsen(node,elem,markedElem,bdEdge) updates
% boundary conditions represented by bdEdge.
%
% [node,elem,bdEdge,brother,indexMap] = coarsen(node,elem,markedElem,bdEdge)
% outputs two additional information: brother and indexMap. 
%
% - brother(:,1:2) stores two triangles sharing the same parent with
%   ordering brother(t,1) < brother(t,2). It is useful for the
%   interpolation of elementwise functions; see also eleminterpolate.
%
% - indexMap is the map between nodes in the fine mesh (node in the input)
%   to that in the coarse mesh (node in the output). For example,
%   indexMap(10) = 6 means the 10-th node in the fine grid  is now the 6-th
%   node in the coarse one. indexMap is useful for the interpolation of
%   nodewise function; see also nodeinterpolate
%
% Example
%
%     load data/Lshapemesh;
%     set(gcf,'Units','normal'); set(gcf,'Position',[0.25,0.25,0.7,0.4]);
%     subplot(1,3,1); showmesh(node,elem);
%     [node,elem] = coarsen(node,elem,'all');
%     subplot(1,3,2); showmesh(node,elem);
%     [node,elem] = coarsen(node,elem,'all');
%     subplot(1,3,3); showmesh(node,elem);
% 
% <a href="matlab:ifem coarsendoc">Algorithm/coarsen</a>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

brother = []; indexMap  = (1:size(node,1))';
if (nargin<=3) || isempty(bdEdge), bdEdge = []; end
if isempty(markedElem), return; end
if strcmp(markedElem,'all'), markedElem = (1:size(elem,1))'; end

%% Find good-to-coarsen nodes
N = size(node,1); NT = size(elem,1);
valence = accumarray(elem(:),ones(3*NT,1),[N 1]);
markedVal = accumarray(elem(markedElem,1),ones(length(markedElem),1),[N 1]);
isIntGoodNode = ((markedVal==valence) & (valence==4));
isBdGoodNode = ((markedVal==valence) & (valence==2));
% isIntGoodNode(1:N0) = false;
% isBdGoodNode(1:N0) = false;
NTdead = 2*sum(isIntGoodNode) + sum(isBdGoodNode); 
if (NTdead == 0), return; end

%% Remove interiori good-to-coarsen nodes
t2v = sparse([1:NT,1:NT,1:NT], elem(1:NT,:), 1, NT, N);
% Find stars for good-to-coarsen nodes
[ii,jj] = find(t2v(:,isIntGoodNode));
if length(jj)<0
    error('number of good nodes are not correct')
end
nodeStar = reshape(ii,4,sum(isIntGoodNode));
isIntNode = false(size(nodeStar,2),1);
% isIntNode is used to exclude the bd node whose val = 4
% case: 1 2 3 4
idx = (elem(nodeStar(1,:),3) == elem(nodeStar(4,:),2)) & ...
      (elem(nodeStar(1,:),2) == elem(nodeStar(2,:),3));  
nodeStar(:,idx)  = nodeStar([1 2 3 4],idx);
isIntNode(idx) = true;
% case: 1 2 4 3
idx = (elem(nodeStar(1,:),3) == elem(nodeStar(3,:),2)) & ...
      (elem(nodeStar(1,:),2) == elem(nodeStar(2,:),3));  
nodeStar(:,idx)  = nodeStar([1 2 4 3],idx); 
isIntNode(idx) = true;
% case: 1 3 2 4
idx = (elem(nodeStar(1,:),3) == elem(nodeStar(4,:),2)) & ...
      (elem(nodeStar(1,:),2) == elem(nodeStar(3,:),3));  
nodeStar(:,idx)  = nodeStar([1 3 2 4],idx); 
isIntNode(idx) = true;
% case: 1 3 4 2
idx = (elem(nodeStar(1,:),3) == elem(nodeStar(2,:),2)) & ...
      (elem(nodeStar(1,:),2) == elem(nodeStar(3,:),3));  
nodeStar(:,idx)  = nodeStar([1 3 4 2],idx); 
isIntNode(idx) = true;
% case: 1 4 2 3
idx = (elem(nodeStar(1,:),3) == elem(nodeStar(3,:),2)) & ...
      (elem(nodeStar(1,:),2) == elem(nodeStar(4,:),3));  
nodeStar(:,idx)  = nodeStar([1 4 2 3],idx); 
isIntNode(idx) = true;
% case: 1 4 3 2
idx = (elem(nodeStar(1,:),3) == elem(nodeStar(2,:),2)) & ...
      (elem(nodeStar(1,:),2) == elem(nodeStar(4,:),3));  
nodeStar(:,idx)  = nodeStar([1 4 3 2],idx); 
isIntNode(idx) = true;
% merge t1 with t2, and t3 with t4
t1 = nodeStar(1,isIntNode); 
t2 = nodeStar(2,isIntNode); 
t3 = nodeStar(3,isIntNode);
t4 = nodeStar(4,isIntNode);
p2 = elem(t1,3); 
p3 = elem(t2,2); 
p4 = elem(t1,2); 
p5 = elem(t3,2);
elem(t1,:) = [p4 p2 p3]; 
elem(t2,1) = 0;
elem(t3,:) = [p5 p3 p2]; 
elem(t4,1) = 0;
% %---------------------------------------------------------sudi
% if isempty(p4) == 0
%     posIntGoodNode = find(isIntGoodNode==1);
%     nIGN = numel(posIntGoodNode);
%     pos_A = zeros(1,nIGN);
%     pos_B = zeros(1,nIGN);
%     pos_C = zeros(1,nIGN);
%     pos_D = zeros(1,nIGN);    
%     
%     for i=1:nIGN
%         rr = find(ismember(elem_store,[posIntGoodNode(i,1) p4(i,1) p2(i,1)],'rows')==1);
%         ss = find(ismember(elem_store,[posIntGoodNode(i,1) p3(i,1) p4(i,1)],'rows')==1);
%         tt = find(ismember(elem_store,[posIntGoodNode(i,1) p5(i,1) p3(i,1)],'rows')==1);
%         uu = find(ismember(elem_store,[posIntGoodNode(i,1) p2(i,1) p5(i,1)],'rows')==1);
%         pos_A(1,i) = elem_store(rr);
%         pos_B(1,i) = elem_store(ss);
%         pos_C(1,i) = elem_store(tt);
%         pos_D(1,i) = elem_store(uu);
%     end
% 
%     UC(:,t1') = 0.5*(UR(:,pos_A) + UR(:,pos_B));
%     UC(1,t2') = 1000000; % one million. It is only for marking!
%     UC(:,t3') = 0.5*(UR(:,pos_C) + UR(:,pos_D));
%     UC(1,t4') = 1000000; % one million. It is only for marking!
% 
%     NEPC(:,t1') = 0.5*(NEPR(:,pos_A) + NEPR(:,pos_B));
%     NEPC(1,t2') = 1000000; % one million. It is only for marking!
%     NEPC(:,t3') = 0.5*(NEPR(:,pos_C) + NEPR(:,pos_D));
%     NEPC(1,t4') = 1000000; % one million. It is only for marking!
% end
% %---------------------------------------------------------sudi

% update isIntGoodNode
intGoodNode = find(isIntGoodNode);
isIntGoodNode(intGoodNode(~isIntNode)) = false;
%% Remove boundary good-to-coarsen nodes
% Find stars for good-to-coarsen nodes
[ii,jj] = find(t2v(:,isBdGoodNode));
if length(jj)<0
    error('number of  good nodes are not correct')
end
nodeStar = reshape(ii,2,sum(isBdGoodNode));
idx = (elem(nodeStar(1,:),3) == elem(nodeStar(2,:),2));
nodeStar(:,idx)  = nodeStar([2 1],idx); 
% Corner/feature points could be good nodes. Remove these points will
% change the shape of the domain.
t5 = nodeStar(1,:); 
t6 = nodeStar(2,:);
p1 = elem(t5,1);
p2 = elem(t5,3); 
p3 = elem(t6,2); 
p4 = elem(t5,2);
v13 = node(p1,:) - node(p3,:);
v12 = node(p1,:) - node(p2,:);
v23 = node(p2,:) - node(p3,:);
lengthDiff = (sum(v12.^2,2) + sum(v13.^2,2) - sum(v23.^2,2))./sum(v23.^2,2);
idx = (sqrt(lengthDiff) < 1e-3);
elem(t5(idx),:) = [p4 p2 p3]; 
elem(t6(idx),1) = 0;
% %---------------------------------------------------------sudi
% if isempty(p4) == 0
%     posBdGoodNode = find(isBdGoodNode==1);
%     nBGN = numel(posBdGoodNode);
%     pos_A = zeros(1,nBGN);
%     pos_B = zeros(1,nBGN);
%     
%     for i=1:nBGN
%         j = find(ismember(elem_store, [posBdGoodNode(i,1) p4(i,1) p2(i,1)],'rows')==1);
%         k = find(ismember(elem_store, [posBdGoodNode(i,1) p3(i,1) p4(i,1)],'rows')==1);
%         pos_A(1,i) = elem_store(j);
%         pos_B(1,i) = elem_store(k);
%     end
% 
%     UC(:,t5(idx)') = 0.5*(UR(:,pos_A) + UR(:,pos_B));
%     UC(1,t6(idx)') = 1000000; % one million. It is only for marking!
% 
%     NEPC(:,t5(idx)') = 0.5*(NEPR(:,pos_A) + NEPR(:,pos_B));
%     NEPC(1,t6(idx)') = 1000000; % one million. It is only for marking!
% end
% %---------------------------------------------------------sudi

bdGoodNode = find(isBdGoodNode);
isBdGoodNode(bdGoodNode(~idx)) = false;

%% Update boundary edges
if (nargin==4) && (~isempty(bdEdge))
	bdEdge(t1,:) = [bdEdge(t1,2) bdEdge(t2,1) bdEdge(t1,1)];
	bdEdge(t3,:) = [bdEdge(t3,2) bdEdge(t4,1) bdEdge(t3,1)];	
	bdEdge(t5,:) = [bdEdge(t5,2) bdEdge(t6,1) bdEdge(t5,1)];
    bdEdge((elem(:,1) == 0),:) = [];
else
	bdEdge=[];
end

%% record brother structure
NTdead = 2*sum(isIntGoodNode) + sum(isBdGoodNode);
brother = zeros(NTdead,2,'uint32');
brother(1:NTdead,1) = [t1'; t3'; t5(idx)'];
brother(1:NTdead,2) = [t2'; t4'; t6(idx)'];

%% Clean node and elem matrices
elem((elem(:,1) == 0),:) = [];
node(isIntGoodNode | isBdGoodNode,:) = [];
indexMap = zeros(N,1);
indexMap(~(isIntGoodNode | isBdGoodNode))= 1:size(node,1);
elem = indexMap(elem);

%---------------------------------------------------------sudi
% Reorder UC and NEPC with respect to the new order of elements
% n_NodeFine = size(node_store,1);
% n_ElemFine = size(elem_store,1);
% n_NodeCoarse = size(node,1);
n_ElemCoarse = size(elem,1);
U_temp = zeros(4,n_ElemCoarse);
NEP_temp = zeros(1,n_ElemCoarse);

Elem_Col23 = [elem_store(:,2) elem_store(:,3)];

for i=1:n_ElemCoarse
    elem_i = elem(i,:);
    hist1 = find(indexMap==elem_i(1,1));
    hist2 = find(indexMap==elem_i(1,2));
    hist3 = find(indexMap==elem_i(1,3));
    history = [hist1 hist2 hist3];
    pos_history = find(ismember(elem_store,history,'rows')==1);
    if isempty(pos_history) == 0
        U_temp(:,i) = UR(:,pos_history);
        NEP_temp(1,i) = NEPR(1,pos_history);
    else
        hist_posA = find(ismember(Elem_Col23,[history(1,1) history(1,2)],'rows')==1);
        hist_posB = find(ismember(Elem_Col23,[history(1,3) history(1,1)],'rows')==1);
        U_temp(:,i) = 0.5*(UR(:,hist_posA) + UR(:,hist_posB));
        NEP_temp(1,i) = 0.5*(NEPR(1,hist_posA) + NEPR(1,hist_posB));
    end    
end
UC = U_temp;
NEPC = NEP_temp;
%---------------------------------------------------------sudi
end