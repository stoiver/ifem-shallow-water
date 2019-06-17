function planardamstandard_triangular
    % Planar Dam Break
    %--------------------------------------------------------------------------
    % Sudi Mungkasi 2012, Australian National University
    %--------------------------------------------------------------------------

    close all; clear all;
    %---------------------- Parameters ----------------------------------------
    %figure(1); set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.8,0.4]);
    t = 0; 
    g = 9.81;
    %---------------------- Initial Grid --------------------------------------
    node = [-1 -1; 1 -1; 1 1; -1 1];
    elem = [2 3 1; 4 1 3];
    %bdFlag = setboundary(node,elem,'Dirichlet','all');

    for i = 1:6    % Try different level here. 7 is good.
        [node,elem] = uniformbisect(node,elem);
    end
    x = [node(elem(:,1),1), node(elem(:,2),1), node(elem(:,3),1)];
    y = [node(elem(:,1),2), node(elem(:,2),2), node(elem(:,3),2)];
    %r = x.^2+y.^2 - 1.0/16.0;
    nt = size(elem,1);    
    C = find_centroid(node,elem);

    h  = zeros(1,nt);    
    for i=1:nt
        if C(i,1) < 0.0
            h(i) = 0.5;
        else
            h(i) = 0.2;
        end
    end  

    uh = zeros(1,nt);
    u  = zeros(1,nt);
    vh = zeros(1,nt);
    v  = zeros(1,nt);
    str = sprintf('The initial number of elements is %d triangles.', nt);
    disp(str);
    str2 = sprintf('The initial time is %d seconds.', t);
    disp(str2);  

    %Plot the initial conditions
%      subplot(1,2,1); showmesh(node,elem);
%      for i=1:size(elem,1)
%          c = h(i)/max(h);
%          tcolor(1,i,1:3) = [c c c];
%      end
%      subplot(1,2,2); 
%      patch(x',y',[h;h;h],tcolor); %pause(0.001)%; view(3) %colorbar;
%      bla

    %%%%%%% MAIN PROGRAM STARTS HERE %%%%%%%%%%%
    Q = zeros(4,nt); %storage initiation for quantity
    Q(1,:) = h;
    Q(2,:) = uh; %xMomentum
    Q(3,:) = vh; %yMomentum  
    v      = vh./h;
    Q(4,:) = 0.5*h.*(u.^2+v.^2) + 0.5*g*h.^2;
    dt = 0.002;
    maxIt = 100;

    for i = 1:maxIt
        i
        nt = size(elem,1); %number of triangles
        % Extract all quantities
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
    Q = U;
    %Plot in colour the results
    %clf(figure);
%     clf;
%     figure(1);
%     h  = Q(1,:);
%     uh = Q(2,:);
%     u  = uh./h;
%     vh = Q(3,:);
%     v  = vh./h;
%     Ent_after = 0.5*h.*(u.^2+v.^2) + 0.5*g*h.^2; %entropy
%     NEP = (1.0/dt)*abs(Ent_after-U(4,:));
%     subplot(1,2,1);
%     showsolution(node,elem,h); view(3); zlim([0. 1.]); %pause(0.01); 
%     subplot(1,2,2);
%     showsolution(node,elem,NEP); view(3); pause(0.01);      
    end
    
    %Plot the results
    strg1 = sprintf('The final number of elements is %d triangles.', nt);
    disp(strg1);    
    strg2 = sprintf('The final time is %d seconds.', t);
    disp(strg2);   
    
    clf(figure);
    x = [node(elem(:,1),1), node(elem(:,2),1), node(elem(:,3),1)];
    y = [node(elem(:,1),2), node(elem(:,2),2), node(elem(:,3),2)];
    h  = Q(1,:);
    uh = Q(2,:);
    u  = uh./h;
    vh = Q(3,:);
    v  = vh./h;
    Ent_after = 0.5*h.*(u.^2+v.^2) + 0.5*g*h.^2; %entropy
    NEP = (1.0/dt)*abs(Ent_after-U(4,:));
    for i=1:size(elem,1)
        c = h(i)/max(h);
        tcolor(1,i,1:3) = [c c c];
    end
    subplot(1,2,1);% showmesh(node,elem);
    patch(x',y',[h;h;h],tcolor); view(3); zlim([0. 5.]);
    subplot(1,2,2);
    patch(x',y',[NEP;NEP;NEP],tcolor); pause(0.001); view(3)%; zlim([0. 1000.]);%colorbar;
	saveas(gcf,'classic_dam_standard.fig')    
end
%-------------------- End of CIRCULAR DAM --------------------------------


%-------------------- Sub functions ---------------
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
if edgeflux_h >= 0.0
    edgeflux_voltrac = edgeflux_h*ct_l;
else
    edgeflux_voltrac = edgeflux_h*ct_r;
end
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
A = 0.5*( x1.*y2 + x2.*y3 + x3.*y1 - x1.*y3 - x2.*y1 - x3.*y2 );
end

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

function neighbor = find_neighbor(elem)
totalEdge = uint32(sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2));
[edge, i2, j] = unique(totalEdge,'rows');
NT = size(elem,1);
elem2edge = uint32(reshape(j,NT,3));
i1(j(3*NT:-1:1)) = 3*NT:-1:1; 
i1 = i1';
k1 = ceil(i1/NT); 
k2 = ceil(i2/NT); 
t1 = i1 - NT*(k1-1);
t2 = i2 - NT*(k2-1);
ix = (i1 ~= i2); 
neighbor = uint32(accumarray([[t1(ix),k1(ix)];[t2,k2]],[t2(ix);t1],[NT 3]));
end
