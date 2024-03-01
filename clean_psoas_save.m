%% Clean results saved by active_contour_3D_g
function clean_psoas_save(folder_data,print_surf)
close all
if ispc
    sep = '\';
else
    sep = '/';
end

% load centers of psoas
[~,subj,X] = fileparts(folder_data);
folder_res = [folder_data, '_res'];
folder_res_clean = [folder_data, '_res_clean'];
mkdir(folder_res_clean)

fd1 = [folder_res,...
    sep,'psoas_', subj, '_centers.mat'];
load(fd1);

% load dicoms
FilesList = dir(folder_res);
while FilesList(1).name(1) == '.'
    FilesList(1) = [];
end
while FilesList(1).name(end-10:end-4) == 'centers'
    FilesList(1) = [];
end
number_of_file = length(FilesList);


for is = 1 : number_of_file
    [folder_res,sep,FilesList(is).name]
    load([folder_res,sep,FilesList(is).name]);
    E_open = E_tot;
    
    E_s=E_open(:,1:xc-x1+1,:);
    E_d=E_open(:,xc-x1+1:size(I_tot,2),:);
    c_sx =[x_sx(1),y_sx(1),it_slice];
    E_s_app = my_conncomp(E_s,c_sx);
    E_s_clean = my_conncomp2(E_s_app,c_sx);
    c_dx =[x_dx(1)-xc+x1,y_dx(1),it_slice];
    E_d_app = my_conncomp(E_d,c_dx);
    E_d_clean = my_conncomp2(E_d_app,c_dx);
    E_clean=cat(2,E_s_clean,E_d_clean(:,2:size(E_d_clean,2),:));
    E_app2=cat(2,E_s_app,E_d_app(:,2:size(E_d_app,2),:));
    E_app2=flip(E_app2,3);
    E_ap2=E_app2;
    E_cl_app=E_clean;
    E_cl_app2=E_clean;
    E_clean=flip(E_clean,3);
    [yp, xp, zp] = ind2sub(size(E_clean),find(E_clean > 0));
    clear psoas_cl        
    psoas_cl(:,1) = y_dicom(yp);
    psoas_cl(:,2) = x_dicom(xp);
    psoas_cl(:,3) = z_dicom(zp);
    %% Print 3D CLEAN PSOAS
    
    if print_surf == 1
        [yDim,xDim,zDim] = size(phi);
        [xx,yy,zz] = meshgrid(1:xDim,1:yDim,1:zDim);

        %     fprintf("\nstampa");
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 1]);
        %     scatter3(xx(:),yy(:),zz(:),[],E_clean(:));
        %     scatter3(xx(:),yy(:),zz(:),[],E_clean(:),'.');
        isosurface(xx,yy,zz,E_clean,0.1)
        axis([1 xDim,1 yDim,1 zDim])
        hold on
        %     grid on
        title(['Clean psoas 2nd step-subj: ',subj]);
    end
    
    
    %% Clean psoas from the new point if necessary
    while(input('Do you want to repeat the procedure? 0) no 1) yes\n'))
        c_sx =[x_sx(1),y_sx(1),it_slice];
        c_dx =[x_dx(1)-xc+x1,y_dx(1),it_slice];
        it_slice_new = input('Choice new initial slice\n');
        
        E_s=E_cl_app2(:,1:xc-x1+1,:);
        E_d=E_cl_app2(:,xc-x1+1:size(I_tot,2),:);
        CC = bwconncomp(E_s(:,:,it_slice_new),4);
        s = regionprops(CC,'centroid');
        centroids = cat(1,s.Centroid);
        centroids = round(centroids);
        if(size(centroids,1)>1)         %if more than one centroid we take the one from largest CC
            idc=sub2ind(size(E_s(:,:,it_slice_new)),centroids(:,2),centroids(:,1));
            numPixels = cellfun(@numel,CC.PixelIdxList);
            [biggest,idx] = max(numPixels);
            idc=idc(ismember(idc,CC.PixelIdxList{idx}));
            size(E_s(:,:,it_slice_new));
            [idc_app(2),idc_app(1)]=ind2sub(size(E_s(:,:,it_slice_new)),idc);
            centroids = idc_app;
        end
        p_sx= [centroids,it_slice_new];
        E_s_app = my_conncomp(E_s,p_sx);
        E_s_clean = my_conncomp2(E_s_app,p_sx);
        CC = bwconncomp(E_d(:,:,it_slice_new),4);
        s = regionprops(CC,'centroid');
        centroids = cat(1,s.Centroid);
        centroids = round(centroids);
        if(size(centroids,1)>1)         %if more than one centroid we take the one from largest CC
            idc=sub2ind(size(E_d(:,:,it_slice_new)),centroids(:,2),centroids(:,1));
            numPixels = cellfun(@numel,CC.PixelIdxList);
            [biggest,idx] = max(numPixels);
            idc=idc(ismember(idc,CC.PixelIdxList{idx}));
            size(E_d(:,:,it_slice_new));
            [idc_app(2),idc_app(1)]=ind2sub(size(E_d(:,:,it_slice_new)),idc);
            centroids = idc_app;
        end
        p_dx= [centroids,it_slice_new];
        E_d_app = my_conncomp(E_d,p_dx);
        E_d_clean = my_conncomp2(E_d_app,p_dx);
        E_clean=cat(2,E_s_clean,E_d_clean(:,2:size(E_d_clean,2),:));
        E_app2=cat(2,E_s_app,E_d_app(:,2:size(E_d_clean,2),:));
        E_app2=flip(E_app2,3);
        E_cl_app=E_clean;
        E_clean=flip(E_clean,3);
        %% Print 3D CLEAN PSOAS
        clear psoas_cl
        [yp, xp, zp] = ind2sub(size(E_cl_app),find(E_cl_app > 0));
            
        psoas_cl(:,1) = y_dicom(yp);
        psoas_cl(:,2) = x_dicom(xp);
        psoas_cl(:,3) = z_dicom(zp);
        if print_surf == 1

            figure
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 1]);
            plot3(psoas_cl(:,1),psoas_cl(:,2),psoas_cl(:,3),'.');
            hold on
            grid on
            axis equal
            
            %     fprintf("\nstampa");
            figure
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 1]);
            %     scatter3(xx(:),yy(:),zz(:),[],E_clean(:));
            %     scatter3(xx(:),yy(:),zz(:),[],E_clean(:),'.');
            isosurface(xx,yy,zz,E_app2,0.1)
            axis([1 xDim,1 yDim,1 zDim])
            hold on
            %         grid on
            title(['Clean psoas 1st step-subj:',subj]);
            
            figure
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 1]);
            %     scatter3(xx(:),yy(:),zz(:),[],E_clean(:));
            %     scatter3(xx(:),yy(:),zz(:),[],E_clean(:),'.');
            isosurface(xx,yy,zz,E_clean,0.1)
            axis([1 xDim,1 yDim,1 zDim])
            hold on
            %         grid on
            title(['Clean psoas 2nd step-subj :',subj]);
        end
        
        
        
        
    end
    
    
    
    E_cl_app2=flip(E_cl_app2,3);
    
    %% STORE RESULTS

    addr_data = [folder_res_clean, sep, FilesList(is).name(1:end-4),'_clean.mat'];
    
    if isfile(addr_data)
        save(addr_data,'E_tot','E_clean','E_app2','E_cl_app2','E_ap2',...
            'psoas_cl','x_dicom','y_dicom','z_dicom',...
            'Resc_Slope','Resc_Intercept','-append')
    else
        save(addr_data,'E_tot','E_clean','E_app2','E_cl_app2','E_ap2',...
            'psoas_cl','x_dicom','y_dicom','z_dicom',...
            'Resc_Slope','Resc_Intercept','-mat')
    end
    
    
end
end
function v = my_conncomp(E_app,c_x)
clear CC CC_app CC_n
E_tot=imbinarize(E_app);
toll = 9;
number_of_file=size(E_tot,3);
for is = 1:number_of_file
    E_tot(:,:,is)=imfill(E_tot(:,:,is),'holes');
end
dim_E = size(E_tot(:,:,c_x(3)));
%     v=NaN(size(E_tot));
v=zeros(size(E_tot));
idc=sub2ind(dim_E,c_x(2),c_x(1));
%     idc_full(c_x(3))=idc;
%     CC = bwconncomp(E_tot(:,:,c_x(3)));
%     CC = bwconncomp(E_tot(:,:,c_x(3)),4);
BW_app=bwareaopen(E_tot(:,:,c_x(3)),5,4);
CC = bwconncomp(BW_app,4);
index = cellfun(@(x) ismember(idc,x),CC.PixelIdxList);
idx = (index == 1);
v_app=zeros(size(E_tot(:,:,c_x(3))));
if idx == 0
   v_app = BW_app; 
else
    v_app(CC.PixelIdxList{idx})=1;   
end
v(:,:,c_x(3))=v_app;
%     v(CC.PixelIdxList{idx},c_x(3)) = 1;
if (c_x(3)<number_of_file)
    CC_app=CC;
    idx_app=idx;
    for is = c_x(3) + 1 : number_of_file
        idc_prev=idc;
        [idc_y_prev,idc_x_prev] = ind2sub(dim_E,idc_prev);
        %             CC_n = bwconncomp(E_tot(:,:,is));
        %             s = regionprops(E_tot(:,:,is),'centroid');
        BW_app=bwareaopen(E_tot(:,:,is),5,4);
        CC_n = bwconncomp(BW_app,4);
        %             CC_n = bwconncomp(E_tot(:,:,is),4);
        v_app=zeros(size(E_tot(:,:,is)));
        if (~isempty(CC_n.PixelIdxList) && ~isempty(nonzeros(idx_app)))                        %Necessario per risolvere problema slice vuota
            %             if ~isempty(CC_n.PixelIdxList)
            s = regionprops(CC_app,'centroid');
            centroids = cat(1,s.Centroid);
            centroids= round(centroids);
            idc=sub2ind(size(E_tot(:,:,is)),centroids(:,2),centroids(:,1));
            idc_app=idc(ismember(idc,CC_app.PixelIdxList{idx_app}));
            if isempty(nonzeros(idc_app))                                       %Necessario per risolvere problema buchi nelle componenti connesse
                min_app = min( abs( idc'-CC_app.PixelIdxList{idx_app}) );       %Si dovrebbe usare distanza euclidea ma qui non è così importante
                [~,icx]=min(min_app);                                           %(L'imfill dovrebbe risolverlo a prescindere)
                idc=idc(icx);
            else
                idc=idc_app;
            end
            [idc_y,idc_x] = ind2sub(dim_E,idc);
            dist=sqrt((idc_x_prev-idc_x).^2+(idc_y_prev-idc_y).^2);
            if(dist>toll)
                idc=idc_prev;
                [idc_y,idc_x] = ind2sub(dim_E,idc);
            end
            index=cellfun(@(x) ismember(idc,x),CC_n.PixelIdxList);
            if isempty(nonzeros(index))
                CC_n_tot=cat(1,CC_n.PixelIdxList{:});
                %                     [ ~, ix ] = min( abs( idc-CC_n_tot) );
                [CC_n_y,CC_n_x]=ind2sub(dim_E,CC_n_tot);
                [ ~, ix ] = min( sqrt((CC_n_x-idc_x).^2+(CC_n_y-idc_y).^2)) ;
                idc=CC_n_tot(ix);
                index=cellfun(@(x) ismember(idc,x),CC_n.PixelIdxList);
                idx_app = (index==1);
                v_app(CC_n.PixelIdxList{idx_app})=1;
            else
                idx_app = (index == 1);
                v_app(CC_n.PixelIdxList{idx_app})=1;
            end
        else
            idx_app = false;
        end
        v(:,:,is)=v_app;
        %             v(cat(1,CC_n.PixelIdxList{idx_app}),is) = 1;
        CC_app=CC_n;
    end
end
idc=sub2ind(size(E_tot(:,:,c_x(3))),c_x(2),c_x(1));
if (c_x(3)>1)
    CC_app=CC;
    idx_app=idx;
    %         idx_app=idx;
    for is = c_x(3) - 1 : -1 : 1
        idc_prev=idc;
        [idc_y_prev,idc_x_prev] = ind2sub(dim_E,idc_prev);
        %             CC_n = bwconncomp(E_tot(:,:,is));
        %             s = regionprops(E_tot(:,:,is),'centroid');
        BW_app=bwareaopen(E_tot(:,:,is),5,4);
        CC_n = bwconncomp(BW_app,4);
        %             CC_n = bwconncomp(E_tot(:,:,is),4);
        v_app=zeros(size(E_tot(:,:,is)));
        if (~isempty(CC_n.PixelIdxList) && ~isempty(nonzeros(idx_app)))                        %Necessario per risolvere problema slice vuota
            %             if ~isempty(CC_n.PixelIdxList)
            s = regionprops(CC_app,'centroid');
            centroids = cat(1,s.Centroid);
            centroids= round(centroids);
            idc=sub2ind(size(E_tot(:,:,is)),centroids(:,2),centroids(:,1));
            idc_app=idc(ismember(idc,CC_app.PixelIdxList{idx_app}));
            if isempty(nonzeros(idc_app))                                       %Necessario per risolvere problema buchi nelle componenti connesse
                min_app = min( abs( idc'-CC_app.PixelIdxList{idx_app}) );
                [~,icx]=min(min_app);
                idc=idc(icx);
            else
                idc=idc_app;
            end
            [idc_y,idc_x] = ind2sub(dim_E,idc);
            dist=sqrt((idc_x_prev-idc_x).^2+(idc_y_prev-idc_y).^2);
            if(dist>toll)
                idc=idc_prev;
                [idc_y,idc_x] = ind2sub(dim_E,idc);
            end
            index=cellfun(@(x) ismember(idc,x),CC_n.PixelIdxList);
            if isempty(nonzeros(index))
                CC_n_tot=cat(1,CC_n.PixelIdxList{:});
                %                     [ ~, ix ] = min( abs( idc-CC_n_tot) );
                [CC_n_y,CC_n_x]=ind2sub(dim_E,CC_n_tot);
                [ ~, ix ] = min( sqrt((CC_n_x-idc_x).^2+(CC_n_y-idc_y).^2)) ;
                idc=CC_n_tot(ix);
                index=cellfun(@(x) ismember(idc,x),CC_n.PixelIdxList);
                idx_app = (index==1);
                v_app(CC_n.PixelIdxList{idx_app})=1;
            else
                idx_app = (index == 1);
                v_app(CC_n.PixelIdxList{idx_app})=1;
            end
        else
            idx_app = false;
        end
        v(:,:,is)=v_app;
        %             v(cat(1,CC_n.PixelIdxList{idx_app}),is) = 1;
        CC_app=CC_n;
    end
end
end

function v = my_conncomp2(E_tot,c_x)
tol = 0.75;
number_of_file=size(E_tot,3);
v=E_tot;
dim_E = size(E_tot(:,:,c_x(3)));
idx = find(E_tot(:,:,c_x(3))==1);
[idx_y,idx_x] = ind2sub(dim_E,idx);
if (c_x(3)<number_of_file)
    size_idx = size(idx);
    pixel_tol= zeros(size_idx);
    % Assume that previous slice is "good" and compute min_tol for each
    % point of current slice
    idx_prev = find(E_tot(:,:,c_x(3)-1)==1);
    [idx_y_prev,idx_x_prev] = ind2sub(dim_E,idx_prev);
    for j = 1 : size_idx(1)
        if(size(idx_prev,1))                %Needed in case of empty previous slice
            min_tol = min(sqrt((idx_x_prev-idx_x(j)).^2+(idx_y_prev-idx_y(j)).^2));
        else
            min_tol = 0;
        end
        pixel_tol(j) = min_tol+tol;
    end
    idx_prev=idx;
    idx_x_prev=idx_x;
    idx_y_prev=idx_y;
    for is = c_x(3) + 1 : number_of_file
        idx_app = find(E_tot(:,:,is)==1);
        idx_diff = setdiff(idx_app,idx_prev);
        for j = 1 : length(idx_diff)
            [idx_y_diff,idx_x_diff] = ind2sub(dim_E,idx_diff(j));
            %                 pixel_dist = min(norm([idx_x_prev,idx_y_prev]-[idx_x_diff,idx_y_diff]));
            [pixel_dist,idx_min] = min(sqrt((idx_x_prev-idx_x_diff).^2+(idx_y_prev-idx_y_diff).^2));
            if pixel_dist > pixel_tol(idx_min)
                v(idx_y_diff,idx_x_diff,is)=0;
            end
        end
        v(:,:,is) = bwmorph(v(:,:,is),'majority');
        %             v(:,:,is)=imfill(v(:,:,is),'holes');
        idx_app = find(v(:,:,is)==1);
        [idx_y_app,idx_x_app] = ind2sub(dim_E,idx_app);
        size_idx = size(idx_app);
        pixel_tol = zeros(size_idx);
        % Compute new tolerance
        for j = 1 : size_idx(1)
            if(size(idx_prev,1))                %Needed in case of empty previous slice
                min_tol = min(sqrt((idx_x_prev-idx_x_app(j)).^2+(idx_y_prev-idx_y_app(j)).^2));
            else
                min_tol = 0;
            end
            pixel_tol(j) = min_tol+tol;
        end
        idx_prev=idx_app;
        idx_x_prev=idx_x_app;
        idx_y_prev=idx_y_app;
    end
end
if (c_x(3)>1)
    size_idx = size(idx);
    pixel_tol= zeros(size_idx);
    % Assume that previous slice is "good" and compute min_tol for each
    % point of current slice
    %         idx_prev = find(E_tot(:,:,c_x(3)+1)==1);
    idx_prev = find(v(:,:,c_x(3)+1)==1);
    [idx_y_prev,idx_x_prev] = ind2sub(dim_E,idx_prev);
    for j = 1 : size_idx(1)
        if(size(idx_prev,1))
            min_tol = min(sqrt((idx_x_prev-idx_x(j)).^2+(idx_y_prev-idx_y(j)).^2));
        else
            min_tol = 0;
        end
        pixel_tol(j) = min_tol+tol;
    end
    idx_prev=idx;
    idx_x_prev=idx_x;
    idx_y_prev=idx_y;
    for is = c_x(3) - 1 : -1 : 1
        %             is
        idx_app = find(E_tot(:,:,is)==1);
        idx_diff = setdiff(idx_app,idx_prev);
        for j = 1 : length(idx_diff)
            [idx_y_diff,idx_x_diff] = ind2sub(dim_E,idx_diff(j));
            %                 pixel_dist = min(norm([idx_x_prev,idx_y_prev]-[idx_x_diff,idx_y_diff]));
            [pixel_dist,idx_min] = min(sqrt((idx_x_prev-idx_x_diff).^2+(idx_y_prev-idx_y_diff).^2));
            if pixel_dist > pixel_tol(idx_min)
                v(idx_y_diff,idx_x_diff,is)=0;
            end
        end
        v(:,:,is) = bwmorph(v(:,:,is),'majority');
        %             v(:,:,is)=imfill(v(:,:,is),'holes');
        idx_app = find(v(:,:,is)==1);
        [idx_y_app,idx_x_app] = ind2sub(dim_E,idx_app);
        size_idx = size(idx_app);
        pixel_tol = zeros(size_idx);
        % Compute new tolerance
        for j = 1 : size_idx(1)
            if(size(idx_prev,1))
                min_tol = min(sqrt((idx_x_prev-idx_x_app(j)).^2+(idx_y_prev-idx_y_app(j)).^2));
            else
                min_tol = 0;
            end
            pixel_tol(j) = min_tol+tol;
        end
        idx_prev=idx_app;
        idx_x_prev=idx_x_app;
        idx_y_prev=idx_y_app;
    end
end
end