function ac_gmfd_2ord(folder_data)
%% ACTIVE CONTOUR 2D

%% CODE

IT_MAX = 2200;
IT_MIN = 10;

% Numerical scheme: 1 geodesic model 2 chan-vese 3 monotone (classical)
met = 1;

if ispc
    sep = '\';
else
    sep = '/';
end


% read file names in folder_data
FilesList = dir(folder_data);
while FilesList(1).name(1) == '.'
    FilesList(1) = [];
end
prompt = 'sigma = ';
sigma=input(prompt);
prompt = 'p = ';
p=input(prompt);



c_min=0;

ep = 0.05;
mu = 1;
nu = 0;

dx = 0.1;



[~,subj,X] = fileparts(folder_data);
folder_res = [folder_data, '_res'];
fd1 = [folder_res,...
    sep,'psoas_', subj, '_centers.mat']
load(fd1);


number_of_file = length(FilesList);
% inizialization of variables
dim_tot=[y2-y1+1,x2-x1+1,number_of_file];
dim_half=[y2-y1+1,0.5*(x2-x1)+1,number_of_file];

E_tot = zeros(dim_tot);
I_tot = zeros(dim_tot);
Ib_tot = zeros(dim_tot);

E_d = zeros(dim_half);
Ib_d = zeros(dim_half);

E_s = zeros(dim_half);
Ib_s = zeros(dim_half);

mask_tot_s = zeros(dim_half);
mask_tot_d = zeros(dim_half);

tic
% cycles through dicoms in the subfolder
for is = 1 : number_of_file
    %% SAVE and MODIFY image
    % save dicom and dicom info
    %         disp(is)
    I_aux  = dicomread([folder_data,sep,FilesList(is).name]);
    %         size(I_aux)
    INFO = dicominfo([folder_data,sep,FilesList(is).name]);
    
    % save position of slice on the z-axis
    z_dicom(is) = INFO.ImagePositionPatient(3);
    
    % save slope and intercept for rescaling of HU
    Resc_Slope = INFO.RescaleSlope;
    Resc_Intercept = INFO.RescaleIntercept;
    
    I_tot(:,:,is) = I_aux(y1:y2,x1:x2);
    I_aux = I_aux*Resc_Slope+Resc_Intercept;
    mask_osso_upper = I_aux>120;
    I_aux(mask_osso_upper) = 0;
    
    %Save position of non-muscle tissue
    mask_osso_upper=I_aux<5;
    
    
    mask_app=zeros(size(I_aux));
    mask_app(mask_osso_upper)=1;
    mask_tot_s(:,:,is)=mask_app(y1:y2,x1:xc);
    mask_tot_d(:,:,is)=mask_app(y1:y2,xc:x2);
    
    %           %-- Histogram equalization
    I_aux = histequaliz(I_aux);
    %         I_aux = imgaussfilt(I_aux,sigma);
    
    %-% Normalizzare il dato tra 0 e 1
    % non necessario per modello classico e GM
    I_aux=I_aux-min(min(I_aux));
    I_aux=I_aux/max(max(I_aux));
    % cut dicom image to speed up computations
    %(to modify xi,yi go to psoas center line 76:79)
    %         I = I_aux(y1:y2,x1:x2);
    Is =I_aux(y1:y2,x1:xc);
    Id =I_aux(y1:y2,xc:x2);
    
    
    if sigma
        I = imgaussfilt(I_aux,sigma);
        Is = imgaussfilt(Is,sigma);
        Id = imgaussfilt(Id,sigma);
    end
    
    Ib_s(:,:,is)=Is;
    Ib_d(:,:,is)=Id;
    Ib_tot(:,:,is)=[Is,Id(:,2:size(Id,2))];
end

Ib_s(mask_tot_s>0)=0;
Ib_d(mask_tot_d>0)=0;
fprintf("\nEnded dicom reading\n");
toc
%% Segmentation with Active Contour

% Create a circular mask with specified
% radius, center, and image size.
is=it_slice;

%% Left psoas Segmentation

% Initial mask

%Take only first input
c_sx =[x_sx(1),y_sx(1),is];
c_sx=dx*c_sx;


%--My code (more precise shape)
r=5*dx;

[dimy,dimx,dimz]=size(Ib_s);
init_mask_s=zeros(size(Ib_s));
for ii=1:dimy
    for jj=1:dimx
        for kk=1:dimz
            app_phi=(jj*dx-c_sx(1))*(jj*dx-c_sx(1))+(ii*dx-c_sx(2))*(ii*dx-c_sx(2))+(kk*dx-c_sx(3))*(kk*dx-c_sx(3));
            init_mask_s(ii,jj,kk)=sqrt(app_phi)-r;
            if init_mask_s(ii,jj,kk)>r
                init_mask_s(ii,jj,kk)=r;
            end
        end
    end
end

tic
if met == 2
    [phi_s, aux_seg,c_s]= lsm_cvfd_3D_v(double(Ib_s),init_mask_s,IT_MAX, IT_MIN,subj,dx,ep,mu,lam1,lam2,1);
elseif met == 1
    [phi_s, aux_seg,c_s]= lsm_gmfd_3D_mod(double(Ib_s),init_mask_s,IT_MAX, IT_MIN,p,subj,dx,nu,mu,ep,c_min,mask_tot_s,1);
else
    [phi_s, aux_seg,c_s]= lsm_classical_3D(double(Ib_s),init_mask_s,IT_MAX, IT_MIN,p,subj,dx,c_min,mask_tot_s,1);
end




toc

seg_s=zeros(size(phi_s));
seg_s(aux_seg)=phi_s(aux_seg);

E_s(aux_seg) = 1;


%% Right psoas segmentation
c_dx =[x_dx(1)-xc+x1,y_dx(1),is];
c_dx=dx*c_dx;

[dimy,dimx,dimz]=size(Ib_d);
init_mask_d=zeros(size(Ib_d));
for ii=1:dimy
    for jj=1:dimx
        for kk=1:dimz
            app_phi=(jj*dx-c_dx(1))*(jj*dx-c_dx(1))+(ii*dx-c_dx(2))*(ii*dx-c_dx(2))+(kk*dx-c_dx(3))*(kk*dx-c_dx(3));
            init_mask_d(ii,jj,kk)=sqrt(app_phi)-r;
            if init_mask_d(ii,jj,kk)>r
                init_mask_d(ii,jj,kk)=r;
            end
        end
    end
end


tic
if met == 2
    [phi_d, aux_seg,c_d]= lsm_cvfd_3D_v(double(Ib_d),init_mask_d,IT_MAX, IT_MIN,subj,dx,ep,mu,lam1,lam2,2);
elseif met == 1
    [phi_d, aux_seg,c_d]= lsm_gmfd_3D_mod(double(Ib_d),init_mask_d,IT_MAX, IT_MIN,p,subj,dx,nu,mu,ep,c_min,mask_tot_d,2);
else
    [phi_d, aux_seg,c_d]= lsm_classical_3D(double(Ib_d),init_mask_d,IT_MAX, IT_MIN,p,subj,dx,c_min,mask_tot_d,2);
end



toc

seg_d=zeros(size(phi_d));
seg_d(aux_seg)=phi_d(aux_seg);

E_d(aux_seg) = 1;

init_mask=cat(2,init_mask_s,init_mask_d(:,2:size(init_mask_d,2),:));
phi=cat(2,phi_s,phi_d(:,2:size(phi_d,2),:));
c=cat(2,c_s,c_d(:,2:size(c_d,2),:));
E_tot=cat(2,E_s,E_d(:,2:size(E_d,2),:));


%% 3D WHOLE PSOAS
x_dicom = INFO.ImagePositionPatient(1)+INFO.PixelSpacing(1)*x1 : ...
    INFO.PixelSpacing(1) : ...
    INFO.ImagePositionPatient(1)+INFO.PixelSpacing(1)*x2;

y_dicom = INFO.ImagePositionPatient(2)+INFO.PixelSpacing(2)*y1 : ...
    INFO.PixelSpacing(2) : ...
    INFO.ImagePositionPatient(2)+INFO.PixelSpacing(1)*y2;

[yp, xp, zp] = ind2sub(size(E_tot),find(E_tot > 0));

clear psoas psoas3d

psoas(:,1) = y_dicom(yp);
psoas(:,2) = x_dicom(xp);
psoas(:,3) = z_dicom(zp);

%% Print 3D PSOAS


figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 1]);
plot3(psoas(:,1),psoas(:,2),psoas(:,3),'.');
hold on
grid on
axis equal
[yDim,xDim,zDim] = size(phi);
[xx,yy,zz] = meshgrid(1:xDim,1:yDim,1:zDim);
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 1]);
isosurface(xx,yy,zz,flip(phi,3),0)
% isosurface(xx,yy,zz,phi,0)
% axis image
axis([1 xDim,1 yDim,1 zDim])
% set(gca, 'ZDir','reverse')
[~,subj,X] = fileparts(folder_data);
title(['Psoas-subj:',subj]);
folder_fig = [folder_data, '_fig'];
mkdir(folder_fig)
if met == 2
    saveas(gcf,[folder_fig, sep, 'psoas_' subj '_surf_met',num2str(met),'_ep',num2str(ep),'_mu',num2str(mu),'_lam',num2str(lam1),'_sig',num2str(sigma),'_ITMAX',num2str(IT_MAX), '_3D.fig'],'fig')
    saveas(gcf,[folder_fig, sep, 'psoas_' subj, '_surf_met',num2str(met),'_ep',num2str(ep),'_mu',num2str(mu),'_lam',num2str(lam1),'_sig',num2str(sigma),'_ITMAX',num2str(IT_MAX), '_3D.png'],'png')
else
    saveas(gcf,[folder_fig, sep, 'psoas_' subj, '_surf_met',num2str(met),'_ep',num2str(ep),'_mu',num2str(mu),'_nu',num2str(nu),'_p',num2str(p),'_sig',num2str(sigma),'_ITMAX',num2str(IT_MAX),'_3D.fig'],'fig')
    saveas(gcf,[folder_fig, sep, 'psoas_' subj, '_surf_met',num2str(met),'_ep',num2str(ep),'_mu',num2str(mu),'_nu',num2str(nu),'_p',num2str(p),'_sig',num2str(sigma),'_ITMAX',num2str(IT_MAX), '_3D.png'],'png')
end

%% STORE FILES
if met == 2
    addr_data = [folder_res, sep, 'psoas_' subj,'_data_3D_met',num2str(met),'_ep',num2str(ep),'_mu',num2str(mu),'_lam',num2str(lam1),'_sig',num2str(sigma),'_ITMAX',num2str(IT_MAX),'.mat'];
else
    addr_data = [folder_res, sep, 'psoas_' subj,'_data_3D_met',num2str(met),'_ep',num2str(ep),'_mu',num2str(mu),'_nu',num2str(nu),'_p',num2str(p),'_sig',num2str(sigma),'_ITMAX',num2str(IT_MAX),'.mat'];
end
if isfile(addr_data)
    save(addr_data, 'phi','c','init_mask',...
        'psoas', 'E_tot','I_tot','Ib_tot',...
        'x_dicom','y_dicom','z_dicom',...
        'Resc_Slope','Resc_Intercept','-append')
else
    save(addr_data,'phi','c','init_mask',...
        'psoas', 'E_tot','I_tot','Ib_tot',...
        'x_dicom','y_dicom','z_dicom',...
        'Resc_Slope','Resc_Intercept','-mat')
end



end

