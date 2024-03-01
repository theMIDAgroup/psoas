function [phi,aux_seg,c]= lsm_gmfd_3D_mod(Ib,init_mask,max_its,min_its,p,subj,dx,nu,mu,ep,c_min,mask_tot,pos)

phi = init_mask;
[yDim,xDim,zDim] = size(phi);
[xx,yy,zz] = meshgrid(1:xDim,1:yDim,1:zDim);

% figure
% isosurface(xx,yy,zz,flip(phi,3),0)
% % isosurface(xx,yy,zz,phi,0)
% % axis image
% axis([1 xDim,1 yDim,1 zDim])
% % set(gca, 'ZDir','reverse')
% 
% title('Initial front')

%-- compute velocity field
 c=calcola_c(Ib,p,dx,c_min);
c(mask_tot>0)=0;
d2x=2*dx;
dx2=dx*dx;
Dc1=(shiftL_3D(c)-shiftR_3D(c))/d2x;
Dc2=(shiftU_3D(c)-shiftD_3D(c))/d2x;
%%NOTA!! 
%sulla z probabilmente dobbiamo scegliere un nuovo step dz appropriato
Dc3=(shiftT_3D(c)-shiftB_3D(c))/d2x;
%Meglio prima o dopo il calcolo di Dc?
%c(mask_tot>0)=0;

%--main loop
f_app=zeros(size(phi));
f0=f_app;
%-- get the narrow band to compute the error
band=1.8*dx;
%     ind_f = find(phi <= band & phi >= -band);
ind_f = phi <= band & phi >= -band;
% f_app(ind_f)=phi(ind_f);
f_app(ind_f)=1;
%f_app(~ind_f)=0;
err=100;

%-- maintain the CFL condition
if ep
    dt = 0.25 * dx2;
%     tol = dx2;
else
    max_c=mu*max(max(max(c)));
    max_Dc1=nu*max(max(max(abs(Dc1))));
    max_Dc2=nu*max(max(max(abs(Dc2))));
    max_Dc3=nu*max(max(max(abs(Dc3))));
    lam_inv=3*max_c + max_Dc1+ max_Dc2 + max_Dc3;
    dt = dx/lam_inv;
%     tol=dx;
end
tol=1;
its=1;

%Video of the evolution
if pos==1
    nome_video=sprintf(['Video_sol/Video_fdgm_mod_3D_subj_',subj,'_p',num2str(p),'_nu',num2str(nu),'_s.avi']);
elseif pos==2
    nome_video=sprintf(['Video_sol/Video_fdgm_mod_3D_subj_',subj,'_p',num2str(p),'_nu',num2str(nu),'_d.avi']);
end
        
%aviobj = VideoWriter(nome_video);
%open(aviobj);

fig=figure;

% max_it=10;
% while (its <=max_it)
% while (its <=max_its && err>tol)
ind_phi=1:size(phi,1)*size(phi,2)*size(phi,3);
% while its<min_its || (its <=max_its && err>tol)
while its<min_its || (its <=max_its)
    %--
    %-- evolve the curve
    if ep
        K= get_curvature_3D(phi,ind_phi,dx);
        phi(ind_phi) = S_gm(phi,ind_phi,c,Dc1,Dc2,Dc3,nu,mu,dt,dx)+ep*dt*c(ind_phi).*K;
    else
        phi(ind_phi) = S_gm(phi,ind_phi,c,Dc1,Dc2,Dc3,nu,mu,dt,dx);
    end
    
     % stopping rule
    %-- get the narrow band to compute the error
%     ind_f = find(phi <= band & phi >= -band);
    ind_f = phi <= band & phi >= -band;
%     f0(ind_f)=phi(ind_f);
    f0(ind_f)=1; 
    f0(~ind_f)=0;
%     err=max(max(max(abs(f0-f_app))));   %max-error norm
    err=sum(sum(sum(abs(f0-f_app))));   %"l1"-error norm
    f_app=f0;
    
    if mod(its-1,10)==0
%      if mod(its-1,1)==0
       fprintf(['\nerr(',num2str(its),')=',num2str(err)]);
   %-% Print whole 3D psoas
        isosurface(xx,yy,zz,flip(phi,3),0)
%         isosurface(xx,yy,zz,phi,0)
        axis([1 xDim,1 yDim,1 zDim])
%         set(gca, 'ZDir','reverse')
        title(['AC.frame:',num2str(its)]);
        pause(0.1);
%        Fi = getframe(fig);
        %writeVideo(aviobj,Fi);
    end
    
  
    %-- Keep SDF smooth - reinitialization?
%      phi = sussman_3D(phi, .5);
    its=its+1;
end
% figure(6)
% surf(phi)
if (its >= max_its)
    fprintf('\nSuperato limite iterazioni\n');
    fprintf(['\nerr(',num2str(its),')=',num2str(err)])
    fprintf('\n');
       %-% Print whole 3D psoas
        isosurface(xx,yy,zz,flip(phi,3),0)
%         isosurface(xx,yy,zz,phi,0)
        axis([1 xDim,1 yDim,1 zDim])
%         set(gca, 'ZDir','reverse')
        title(['AC.frame:',num2str(its)]);
        pause(0.1);
        Fi = getframe(fig);
        %writeVideo(aviobj,Fi);
else
    fprintf('\nN_i=%d\n',its);
    fprintf(['\nerr(',num2str(its),')=',num2str(err)])
    fprintf('\n');
    if mod(its-2,10)~=0
         %-% Print whole 3D psoas
        isosurface(xx,yy,zz,flip(phi,3),0)
%         isosurface(xx,yy,zz,phi,0)
        axis([1 xDim,1 yDim,1 zDim])
%         set(gca, 'ZDir','reverse')
        title(['AC.frame:',num2str(its)]);
        pause(0.1);
        Fi = getframe(fig);
        %writeVideo(aviobj,Fi);
    end
end
%close(aviobj);
 aux_seg=phi<0;
%  seg(aux_seg)=phi(aux_seg);
end



%---------------------------------------------------------------------
%-- AUXILIARY FUNCTIONS ----------------------------------------------
%---------------------------------------------------------------------
    
%-- converts a mask to a SDF
% function phi = mask2phi(init_a)
%   phi = bwdist(init_a)-bwdist(1-init_a)+im2double(init_a)-.5;
% end

 %Hamiltonian and derivatives
function v = H(c,Dc1,Dc2,Dc3,p,q,r,nu,mu)
v=mu*c.*sqrt(p.*p+q.*q+r.*r)-nu*(Dc1.*p+Dc2.*q+Dc3.*r);
% v=mu*c.*sqrt(p.*p+q.*q+r.*r)+nu*(Dc1.*p+Dc2.*q+Dc3.*r);
% v=mu*c.*(sqrt(p.*p+q.*q+r.*r)-nu*(Dc1.*p+Dc2.*q+Dc3.*r));
end

function v = Hp(c,Dc1,p,q,r,nu,mu)
v=mu*c.*p./sqrt(p.*p+q.*q+r.*r+eps)-nu*Dc1;
% v=mu*c.*p./sqrt(p.*p+q.*q+r.*r+eps)+nu*Dc1;
% v=mu*c.*(p./sqrt(p.*p+q.*q+r.*r+eps)-nu*Dc1);
end

function v = Hq(c,Dc2,p,q,r,nu,mu)
v=mu*c.*q./sqrt(p.*p+q.*q+r.*r+eps)-nu*Dc2;
% v=mu*c.*q./sqrt(p.*p+q.*q+r.*r+eps)+nu*Dc2;
% v=mu*c.*(q./sqrt(p.*p+q.*q+r.*r+eps)-nu*Dc2);
end

function v = Hr(c,Dc3,p,q,r,nu,mu)
v=mu*c.*r./sqrt(p.*p+q.*q+r.*r+eps)-nu*Dc3;
% v=mu*c.*r./sqrt(p.*p+q.*q+r.*r+eps)+nu*Dc3;
% v=mu*c.*(r./sqrt(p.*p+q.*q+r.*r+eps)-nu*Dc3);
end

%Local Lax-Friedrichs scheme
function v = G_lf(c,Dc1,Dc2,Dc3,u_m,u_p,v_m,v_p,w_m,w_p,nu,mu)
H_p=abs(Hp(c,Dc1,u_m,0,0,nu,mu));
ind=abs(Hp(c,Dc1,u_p,0,0,nu,mu))>H_p;
H_p(ind)=abs(Hp(c(ind),Dc1(ind),u_p(ind),0,0,nu,mu));
H_q=abs(Hq(c,Dc2,0,v_m,0,nu,mu)); 
ind=abs(Hq(c,Dc2,0,v_p,0,nu,mu))>H_q;
H_q(ind)=abs(Hq(c(ind),Dc2(ind),0,v_p(ind),0,nu,mu));
H_r=abs(Hr(c,Dc3,0,0,w_m,nu,mu));
ind=abs(Hr(c,Dc3,0,0,w_p,nu,mu))>H_r;
H_r(ind)=abs(Hr(c(ind),Dc3(ind),0,0,w_p(ind),nu,mu));
v=H(c,Dc1,Dc2,Dc3,0.5*(u_m+u_p),0.5*(v_m+v_p),0.5*(w_m+w_p),nu,mu)-0.5*(H_p.*(u_p-u_m)+H_q.*(v_p-v_m)+H_r.*(w_p-w_m));

%  Lax-Friedrichs
% v=H(c_xy,0.5*(u_m+u_p),0.5*(v_m+v_p),0.5*(w_m+w_p))-0.5*(u_p-u_m)-0.5*(v_p-v_m)+0.5*(w_p-w_m);

%  Prova evoluzione
%  v=H(1,0.5*(u_m+u_p),0.5*(v_m+v_p),0.5*(w_m+w_p))-0.5*(u_p-u_m)-0.5*(v_p-v_m)-0.5*(w_p-w_m);
end

  
  function v = S_gm(phi,ind_phi,c,Dc1,Dc2,Dc3,nu,mu,dt,dx)
      [dimy, dimx, dimz] = size(phi);        
    [y, x, z] = ind2sub([dimy,dimx,dimz],ind_phi);                                % get subscripts

    %-- get subscripts of neighbors
    ym1 = y-1; 
    xm1 = x-1; 
    zm1 = z-1; 
    yp1 = y+1; 
    xp1 = x+1;  
    zp1 = z+1;

    %-- bounds checking  
    ym1(ym1<1) = 1; 
    xm1(xm1<1) = 1;  
    zm1(zm1<1) = 1;            
    yp1(yp1>dimy) = dimy; 
    xp1(xp1>dimx) = dimx;    
    zp1(zp1>dimz) = dimz;

%   %-- get indexes for 6 neighbors 
    izm = sub2ind(size(phi),y,x,zm1); 
    iym = sub2ind(size(phi),ym1,x,z); 
    ixm = sub2ind(size(phi),y,xm1,z);    
    i = sub2ind(size(phi),y,x,z); 
    ixp = sub2ind(size(phi),y,xp1,z);   
    iyp = sub2ind(size(phi),yp1,x,z); 
    izp = sub2ind(size(phi),y,x,zp1); 
    
%   %-- compute Hamiltonian    
    phi_xm = (phi(i)-phi(ixm))/dx;
    phi_xp = (phi(ixp)-phi(i))/dx;
    phi_ym = (phi(i)-phi(iym))/dx;
    phi_yp = (phi(iyp)-phi(i))/dx;
    phi_zm = (phi(i)-phi(izm))/dx;
    phi_zp = (phi(izp)-phi(i))/dx;
    v=phi(ind_phi) - dt.*G_lf(c(ind_phi),Dc1(ind_phi),Dc2(ind_phi),Dc3(ind_phi),phi_xm,phi_xp,phi_ym,phi_yp,phi_zm,phi_zp,nu,mu);  
%     v=G_lf(c(ind_phi),Dc1(ind_phi),Dc2(ind_phi),Dc3(ind_phi),phi_xm,phi_xp,phi_ym,phi_yp,phi_zm,phi_zp,nu,mu);  
  end 

function v=calcola_c(I,p,dx,c_min)
v=zeros(size(I));
[dimy, dimx,dimz] = size(I);
for i=1:dimy
    for j=1:dimx
        for k=1:dimz
            if((j~=dimx)&&(j~=1))
                %Ix=(I[i][j+1]-I[i][j-1])/h2x;		%DIFFERENZE CENTRATE
                Ix=abs(I(i,j+1,k)-I(i,j,k))/dx;
                app=abs(I(i,j,k)-I(i,j-1,k))/dx;
                if(app>Ix)
                    Ix=app;
                end
%             else
%                 %Ix=(I[i][1]-I[i][Nx-1])/h2x;		%PERIODICHE
%                 Ix=0;								%NEUMANN OMOGENEE
%             end
            elseif j==1
                Ix=abs(I(i,j+1,k)-I(i,j,k))/dx;    %ONE-SIDED
            else
                Ix=abs(I(i,j,k)-I(i,j-1,k))/dx;
            end
            
            if((i~=dimy)&&(i~=1))
                %Iy=(I[i+1][j]-I[i-1][j])/h2y;	
                Iy=abs(I(i+1,j,k)-I(i,j,k))/dx;
                app=abs(I(i,j,k)-I(i-1,j,k))/dx;
                if(app>Iy)
                    Iy=app;
                end
%             else
%                 %Iy=(I[1][j]-I[Ny-1][j])/h2y;		%PERIODICHE
%                 Iy=0;								%NEUMANN OMOGENEE
%             end
            elseif i==1
                Iy=abs(I(i+1,j,k)-I(i,j,k))/dx;         %ONE-SIDED
            else
                Iy=abs(I(i,j,k)-I(i-1,j,k))/dx;
            end
            if((k~=dimz)&&(k~=1))
                %Iy=(I[i+1][j]-I[i-1][j])/h2y;	
                Iz=abs(I(i,j,k+1)-I(i,j,k))/dx;
                app=abs(I(i,j,k)-I(i,j,k-1))/dx;
                if(app>Iz)
                    Iz=app;
                end
%             else
%                 %Iz=(I[1][j]-I[Ny-1][j])/h2y;		%PERIODICHE
%                 Iz=0;								%NEUMANN OMOGENEE
%             end
            elseif k==1
                Iz=abs(I(i,j,k+1)-I(i,j,k))/dx;         %ONE-SIDED
            else
                Iz=abs(I(i,j,k)-I(i,j,k-1))/dx;
            end
            mod_DI=sqrt(Ix*Ix+Iy*Iy+Iz*Iz);
            v(i,j,k)=1/(1+mod_DI^p);
%               if(mod)
              if v(i,j,k)<c_min
                  v(i,j,k)=0;
              end
        end
    end
end
end

%-- compute curvature along SDF
function curvature = get_curvature_3D(phi,idx,dx)
    [dimy, dimx, dimz] = size(phi);        
    [y, x, z] = ind2sub([dimy,dimx,dimz],idx);                                % get subscripts

    %-- get subscripts of neighbors
    ym1 = y-1; 
    xm1 = x-1; 
    zm1 = z-1; 
    yp1 = y+1; 
    xp1 = x+1;  
    zp1 = z+1;

    %-- bounds checking  
    ym1(ym1<1) = 1; 
    xm1(xm1<1) = 1;  
    zm1(zm1<1) = 1;            
    yp1(yp1>dimy) = dimy; 
    xp1(xp1>dimx) = dimx;    
    zp1(zp1>dimz) = dimz;

%     %-- get indexes for 26 neighbors
    id1 = sub2ind(size(phi),ym1,xm1,zm1);    
    id2 = sub2ind(size(phi),ym1,x,zm1); 
    id3 = sub2ind(size(phi),ym1,xp1,zm1); 
    id4 = sub2ind(size(phi),y,xm1,zm1);    
    id5 = sub2ind(size(phi),y,x,zm1); 
    id6 = sub2ind(size(phi),y,xp1,zm1); 
    id7 = sub2ind(size(phi),yp1,xm1,zm1);    
    id8 = sub2ind(size(phi),yp1,x,zm1); 
    id9 = sub2ind(size(phi),yp1,xp1,zm1); 
    
    id10 = sub2ind(size(phi),ym1,xm1,z);    
    id11 = sub2ind(size(phi),ym1,x,z); 
    id12 = sub2ind(size(phi),ym1,xp1,z); 
    id13 = sub2ind(size(phi),y,xm1,z);    
    id14 = sub2ind(size(phi),y,x,z); 
    id15 = sub2ind(size(phi),y,xp1,z); 
    id16 = sub2ind(size(phi),yp1,xm1,z);    
    id17 = sub2ind(size(phi),yp1,x,z); 
    id18 = sub2ind(size(phi),yp1,xp1,z); 
    
    id19 = sub2ind(size(phi),ym1,xm1,zp1);    
    id20 = sub2ind(size(phi),ym1,x,zp1); 
    id21 = sub2ind(size(phi),ym1,xp1,zp1); 
    id22 = sub2ind(size(phi),y,xm1,zp1);    
    id23 = sub2ind(size(phi),y,x,zp1); 
    id24 = sub2ind(size(phi),y,xp1,zp1); 
    id25 = sub2ind(size(phi),yp1,xm1,zp1);    
    id26 = sub2ind(size(phi),yp1,x,zp1); 
    id27 = sub2ind(size(phi),yp1,xp1,zp1); 
    
    dx2=dx*dx;
    d2x=2*dx;
    %-- get central derivatives of SDF at x,y,z
    phi_x  = (-phi(id13)+phi(id15))/d2x;
    phi_y  = (-phi(id11)+phi(id17))/d2x;
    phi_z = (-phi(id5)+phi(id23))/d2x;
    
    phi_xx = (phi(id13)-2*phi(idx)+phi(id15))/dx2;
    phi_yy = (phi(id11)-2*phi(idx)+phi(id17))/dx2;
    phi_zz = (phi(id5)-2*phi(idx)+phi(id23))/dx2;
    
    phi_xy = (-phi(id10)-phi(id18)...
             +phi(id12)+phi(id16))/(4*dx2);
    phi_xz = (-phi(id4)-phi(id24)...
             +phi(id6)+phi(id22))/(4*dx2);
    phi_yz = (-phi(id2)-phi(id26)...
             +phi(id8)+phi(id20))/(4*dx2);
    
%     phi_xyz = (1/(4*dx2*d2x))*(-phi(id1)-phi(id27)...
%              +phi(id3)+phi(id25)...
%              -phi(id7) - phi(id21)...
%              +phi(id9)+phi(id19));
         
    phi_x2 = phi_x.^2;
    phi_y2 = phi_y.^2;
    phi_z2 = phi_z.^2;
%     %-- compute curvature (Kappa)
%     curvature = ((phi_x2.*phi_yy + phi_y2.*phi_xx - 2*phi_x.*phi_y.*phi_xy)./...
%               (phi_x2 + phi_y2 +eps).^(3/2)).*(phi_x2 + phi_y2).^(1/2);        
    %-- compute curvature (Kappa)
%     curvature = (sqrt((phi_zz.*phi_y - phi_yy.* phi_z ).^2 + (phi_xx.* phi_z - phi_zz.*phi_x).^2 + (phi_yy.*phi_x - phi_y.*phi_xx).^2))./...
%               (phi_x2 + phi_y2 + phi_z2 + eps).^(3/2);%).*(phi_x2 + phi_y2 + phi_z2).^(1/2);
    curvature = (2*(phi_xy.*phi_x.*phi_y+phi_xz.*phi_x.*phi_z+phi_yz.*phi_z.*phi_y)-(phi_xx.*(phi_y2+phi_z2)+phi_yy.*(phi_x2+phi_z2)...
        +phi_zz.*(phi_x2+phi_y2)))./(phi_x2 + phi_y2 + phi_z2 + eps).^(3/2).*(phi_x2 + phi_y2 + phi_z2).^(1/2);
    %rispetto alla formula di Osher e Sethian il segno è diverso
%    curvature=curvature/2;
end

%-- level set re-initialization by the sussman method
function D = sussman_3D(D, dt)
  % forward/backward differences
  a = D - shiftR_3D(D); % Right on the x
  b = shiftL_3D(D) - D; % Left on the x
  c = D - shiftD_3D(D); % Down to the y
  d = shiftU_3D(D) - D; % Up to the y
  e = D - shiftB_3D(D);
  f = shiftT_3D(D) - D;
  
  a_p = a;  a_n = a; % a+ and a-
  b_p = b;  b_n = b;
  c_p = c;  c_n = c;
  d_p = d;  d_n = d;
  e_p = e;  e_n = e;
  f_p = f;  f_n = f;
  a_p(a < 0) = 0;
  a_n(a > 0) = 0;
  b_p(b < 0) = 0;
  b_n(b > 0) = 0;
  c_p(c < 0) = 0;
  c_n(c > 0) = 0;
  d_p(d < 0) = 0;
  d_n(d > 0) = 0;
  e_p(e < 0) = 0;
  e_n(e > 0) = 0;
  f_p(f < 0) = 0;
  f_n(f > 0) = 0;
  
  % vedi parismi-master file Image3D.h
  dD = zeros(size(D));
  D_neg_ind = find(D < 0);
  D_pos_ind = find(D > 0);
  dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                       + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2) + max(e_p(D_pos_ind).^2, f_n(D_pos_ind).^2)) - 1;
  dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                       + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2) + max(e_n(D_neg_ind).^2, f_p(D_neg_ind).^2)) - 1;
  
  D = D - dt .* sussman_sign(D) .* dD;
end
  
%-- whole matrix derivatives
%shift up on y
function shift = shiftU_3D(M)
  shift = [M(2:size(M,1),:,:); M(size(M,1),:,:) ];
end

% shift down on y
function shift = shiftD_3D(M)
  shift = [M(1,:,:); M(1:size(M,1)-1,:,:) ];
end

% shift left on x
function shift = shiftL_3D(M)
  shift = [ M(:,2:size(M,2),:) M(:,size(M,2),:) ];
end

%shift right on x
function shift = shiftR_3D(M)
  shift = [ M(:,1,:) M(:,1:size(M,2)-1,:) ];
end

%shift top on z
function shift = shiftT_3D(M)
  shift(:,:, 1:size(M,3)-1) =  M(:,:, 2:size(M,3));
  shift(:,:, size(M,3)) =  M(:,:, size(M,3));
end

%shift bottom on z
function shift = shiftB_3D(M)
  shift(:,:,1) =  M(:,:,1);
  shift(:,:, 2:size(M,3)) = M(:,:, 1:size(M,3)-1);
end

  
function S = sussman_sign(D)
  S = D ./ sqrt(D.^2 + 1);    
end