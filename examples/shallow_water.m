function shallow_water
%% Shallow water Planar Dam Break
%--------------------------------------------------------------------------
% Sudi Mungkasi 2012, Australian National University
% Steve Roberts 2019, Australian National University
%--------------------------------------------------------------------------

%%
close all; clear all;

%---------------------- Parameters ----------------------------------------
%figure(1); set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.8,0.4]);

t = 0;
g = 9.81;
%---------------------- Initial Grid --------------------------------------
node = [-1 -1; 1 -1; 1 1; -1 1];
elem = [2 3 1; 4 1 3];
%bdFlag = setboundary(node,elem,'Dirichlet','all');

for i = 1:4 % Try different level here. 7 is good.
    [node,elem] = uniformbisect(node,elem);
end

x = [node(elem(:,1),1), node(elem(:,2),1), node(elem(:,3),1)];
y = [node(elem(:,1),2), node(elem(:,2),2), node(elem(:,3),2)];


%-------------------------- Initial Condition -----------------------------
nt = size(elem,1);
C  = find_centroid(node,elem);
h  = zeros(nt,1);
u  = zeros(nt,1);
v  = zeros(nt,1);
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

Q = zeros(nt,4); %storage initiation for quantity
Q(:,1) = h;      %height
Q(:,2) = u.*h;   %xMomentum
Q(:,3) = v.*h;   %yMomentum
Q(:,4) = 0.5*h.*(u.^2+v.^2) + 0.5*g*h.^2; %entropy

Area      = find_area(node,elem);
tol_Area  = max(Area);
init_res = tol_Area;

fprintf('The initial number of elements is %d triangles.\n', nt);
fprintf('The initial time is %d seconds.\n', t);


% Plot the initial conditions
% clf(figure(3));
% figure(3);
% showmesh(node,elem);
% axis on;
% findelem(node,elem);  % plot indices of all triangles
% findnode(node);       % plot indices of all vertices


clf(figure(1));
figure(1);
set(gcf,'Units','normal'); set(gcf,'Position',[0.8,0.05,0.15,0.15]);
showmesh(node,elem); view(2); xlabel('x'); ylabel('y');


clf(figure(2));
figure(2);
set(gcf,'Units','normal'); set(gcf,'Position',[0.8,0.3,0.15,0.15]);
subplot(1,2,1);
showsolution(node,elem,Q(:,1)); view(3); zlim([-1. 1.]);
xlabel('x'); ylabel('y'); zlabel('h')
subplot(1,2,2);
showsolution(node,elem,Q(:,2)); view(3); zlim([-1. 1.]);
xlabel('x'); ylabel('y'); zlabel('uh')


%% MAIN PROGRAM STARTS HERE %%%%%%%%%%%

maxRefineCoarsen = 2;
maxIt = 0;
dt = 0.001;


while (t < 0.1) 
    
    nt = size(elem,1); %number of triangles
    %fprintf('nt = %g \n',nt);
    
    for iter = 1:1
        % Test evolve; Used to decide on new mesh
        
        h  = Q(:,1);
        uh = Q(:,2);
        u  = uh./h;
        vh = Q(:,3);
        v  = vh./h;
        entropy_before = 0.5*h.*(u.^2+v.^2) + 0.5*g*h.^2;
        
        Q(:,4) = entropy_before;
        
        [U, Area, dt] = evolve_step(node,elem,Q,dt,g);
        
        
        % Begin Refinement ----------------------------------------------------
        h  = U(:,1);
        uh = U(:,2);
        u  = uh./h;
        vh = U(:,3);
        v  = vh./h;
        entropy_evolved = U(:,4);
        entropy_formula = 0.5*h.*(u.^2+v.^2) + 0.5*g*h.^2; %entropy
        NEP = abs((1/dt)*(entropy_formula-entropy_evolved));%.*Area';
        
        tol_NEP = 0.1*max(abs(NEP));
        refineElem  = find((abs(NEP) > tol_NEP) .* (Area >= init_res/64) );
        
        %     elems_to_be_refined = refineElem;
        %
        %     for j = 1: size(refineElem,2)
        %         if Area(refineElem(1,j))/tol_Area < 0.000001+1/(2^maxRefineCoarsen) %0.126 %1.1 %0.6 %0.26
        %             elems_to_be_refined = setdiff(elems_to_be_refined, refineElem(1,j));
        %         end
        %     end
        %
        %     refineElem = elems_to_be_refined;
        
        % Do the bisection of the domain
        
        %     fprintf('|refineElem| = %g \n',size(refineElem,1));
        [node,elem,~,~,tree] = bisect(node,elem, refineElem);
        
        
        %     clf(figure(6));
        %     figure(6);
        %     set(gcf,'Units','normal'); set(gcf,'Position',[0.6,0.25,0.15,0.15]);
        %     showmesh(node,elem); view(2); xlabel('x'); ylabel('y');
        
        ntn = size(elem,1);  %number of triangles after refinement
        %fprintf('refine ntn = %g \n',ntn);
        QR = zeros(ntn,4);
        for ii = 1:4
            QR(:,ii) = eleminterpolate(Q(:,ii),tree);
        end
        NEPR = eleminterpolate(NEP,tree);
        
        % End Refinement ------------------------------------------------------
        
        
        % Begin Coarsening ----------------------------------------------------
        Area = find_area(node,elem);  %areas of triangles
        tol_Area  = max(Area);
        coarsenElem = find(abs(NEPR) < 0.5*tol_NEP);
        
        
        %fprintf('|coarsenElem| = %g \n',size(coarsenElem,1));
        [node,elem,~,~,tree] = coarsen(node,elem, coarsenElem, []);
        
        %     clf(figure(5));
        %     figure(5);
        %     set(gcf,'Units','normal'); set(gcf,'Position',[0.6,0.05,0.15,0.15]);
        %     showmesh(node,elem); view(2); xlabel('x'); ylabel('y');
        
        
        
        ntn = size(elem,1);  %number of triangles after refinement
        QC = zeros(ntn,4);
        %     fprintf('coarsen ntn = %g \n',ntn);
        %     fprintf('size(QR) = %g \n',size(QR));
        %     fprintf('size(QC) = %g \n',size(QC));
        %     fprintf('size(tree) = %g \n',size(tree));
        
        QC = zeros(ntn,4);
        for ii = 1:4
            QC(:,ii) = eleminterpolate(QR(:,ii),tree);
            %         fprintf('size(QC(:,ii)) = %g \n',size(QC(:,ii)));
        end
        NEPC = eleminterpolate(NEPR,tree);
        
        % End Coarsen -----------------------------------------------------
        
        
        
        
        %Collect quantities for next iteration from the refine grids
        Q = QC;
        NEP = NEPC;
        % End Coarsening --------------------------------------------------
        
    end
    
    % Actual Evolve -------------------------------------------------------
    
    [Q, Area, dt] = evolve_step(node,elem,Q,dt,g);
    
    
    tol_Area  = max(Area);
    t = t+dt;
    
    fprintf('TIME %g dt %g \n',t, dt);
    
    
    
    
    %Plot in colour the results
    h  = Q(:,1);
    uh = Q(:,2);
    u  = uh./h;
    vh = Q(:,3);
    v  = vh./h;
    entropy = Q(:,4);
    entropy_formula = 0.5*h.*(u.^2+v.^2) + 0.5*g*h.^2; %entropy
    
    NEP =(1/dt)*(entropy_formula-entropy);
    maxNEP = max(abs(NEP));
    NEP = NEP/maxNEP;
    
    fprintf('max h  = %g, min h  = %g \n',max(h),min(h));
    fprintf('max uh = %g, min uh = %g \n',max(uh),min(uh));
    fprintf('max ent = %g, min ent = %g \n',max(entropy),min(entropy));
end


clf(figure(3));
figure(3);
set(gcf,'Units','normal'); set(gcf,'Position',[0.8,0.55,0.15,0.15]);
subplot(1,2,1);
showsolution(node,elem,h); view(3); zlim([0. 1.]);
xlabel('x'); ylabel('y'); zlabel('Stage');

subplot(1,2,2);
showsolution(node,elem,abs(NEP)); view(3);
xlabel('x'); ylabel('y'); zlabel('|NEP|');

clf(figure(4));
figure(4);
set(gcf,'Units','normal'); set(gcf,'Position',[0.8,0.75,0.15,0.15]);
subplot(1,2,1);
showsolution(node,elem,uh/max(abs(uh))); view(3); zlim([0. 1.]);
xlabel('x'); ylabel('y'); zlabel('uh/max(|uh|)')
subplot(1,2,2);
showmesh(node,elem); view(3); xlabel('x'); ylabel('y');



end



%% Sub Functions

function [U, Area, dt] = evolve_step(node,elem,Q,dt,g)

CFL = 0.5;

% Extract all quantities
nt = size(elem,1);
% h  = Q(:,1);
% uh = Q(:,2);
% u  = uh./h;
% vh = Q(:,3);
% v  = vh./h;
% entropy = Q(:,4);
%Q(:,4) = 0.5*h.*(u.^2+v.^2) + 0.5*g*h.^2;

% Extract domain properties
T = auxstructure(elem);
Neighbor = T.neighbor;             %find_neighbor(elem); %neighbor
[Normal, Edgelength] = find_normals_edgelengths(node,elem); %normals
Area      = find_area(node,elem);  %areas of triangles
U = Q; 

%initiation for temporary storage of the quantity
Flux = zeros(size(Q));

% Evolve
for j = 1:nt
    Ul = Q(j,:);
    flux = 0.0;
    max_speed = 0.0;
    for k = 1:3
        neighbor = Neighbor(j,k);
        normal   = Normal(:,k,j);
        edgelength = Edgelength(1,k,j);
        
        Ul_rot = rotate(Ul, normal(1), normal(2));
        Ur = Q(neighbor,:);
        Ur_rot = rotate(Ur, normal(1), normal(2));
        
        % Use reflective boundary at the moving water
        if neighbor==j
            Ur_rot(2) = -Ur_rot(2);      
        end
       
        
        [edgeflux, ms] = FluxKNP(Ul_rot,Ur_rot,g);
        
        max_speed = max(max_speed, ms);
        %dt = min(dt, CFL*0.5*edgelength/max_speed);
        
        edgeflux  = rotate(edgeflux, normal(1), -normal(2));
        flux      = flux - edgeflux*edgelength;
    end
    
    
    Flux(j,:) = flux'/Area(j);
end

U = Q + dt*Flux;
end %% -------------------- End of Evolve Step ----------------------------


%%
function [F, max_speed] = FluxKNP(Ul,Ur,g)
%dt=1000.0 must return F and dt
F=zeros(4,1);
h_l  = Ul(1);     %height
u_l  = Ul(2)/h_l; %x velocity
ct_l = Ul(3)/h_l; %y velocity, that is, Concentration of Tracer
e_l  = Ul(4);     %entropy

h_r  = Ur(1);
u_r  = Ur(2)/h_r;
ct_r = Ur(3)/h_r;
e_r  = Ur(4);
        
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
edgeflux_entropy = (s_max*entropy_flux_l - s_min*entropy_flux_r + s_max*s_min*(e_r-e_l))/denom;
max_speed = max(abs(s_max), abs(s_min));
        
% Update timestep
%dt = min(dt, CFL*0.5*dx/max_speed);

F(1) = edgeflux_h;
F(2) = edgeflux_p;
F(3) = edgeflux_voltrac;
F(4) = edgeflux_entropy;


%return F, dt
end

%%
function Q = rotate(q, n1, n2)
% Rotate the momentum component q, that is, q(2) and q(3)
% from x,y coordinates to coordinates based on normal vector (n1, n2).
% Result is returned in array 4x1 r
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

%%
function C = find_centroid(node,elem)
x = (node(elem(:,1),1) + node(elem(:,2),1) + node(elem(:,3),1))/3;
y = (node(elem(:,1),2) + node(elem(:,2),2) + node(elem(:,3),2))/3;
C = [x,y]; %space initiation
end


%%
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
function [n,el] = find_normals_edgelengths(node,elem)

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



