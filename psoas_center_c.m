function psoas_center_c(folder_data)
%% Funzione psoas_center, chiamata dalla GUI con il pushbutton3 'Center'.
% takes user input psoas center and save it
%
%
%
%% CODE

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

I = dicomread([folder_data,sep,FilesList(1).name]);

imagesc(I);                                                             % Mostra l'immagine in modo da poter interagire
colormap bone

[xc, yc]=ginput(1);

xc= floor(xc);
yc= floor(yc);



x1 = xc-100;
x2 = xc+100;
y1 = yc-75;
y2 = yc+75;
prompt='Enter number of slice for the center';
it_app=inputdlg(prompt);
it_slice=str2double(it_app);

I  = dicomread([folder_data, sep,FilesList(it_slice).name]);                % Legge il file
I_aux = I(y1:y2,x1:x2);                                                 % Seleziona gli estremi in base al valore dato in precedenza

imagesc(I_aux);                                                         % Mostra l'immagine in modo da poter interagire
colormap bone                                                           % Utilizziamo imagesc per non aver problemi con la scala di grigi
%     title(['slice ', num2str(it_slice)])


[x,y] = ginput(1);
x_sx = floor(x);
y_sx = floor(y);

[x,y] = ginput(1);
x_dx = floor(x);
y_dx = floor(y);


folder_res = [folder_data, '_res'];
mkdir(folder_res)
[~,subj,X] = fileparts(folder_data);

save([folder_res,...
    sep,'psoas_', subj, '_centers.mat'], ...
    'it_slice','x_dx','y_dx','x_sx','y_sx','x1','x2', 'xc',...
    'y1', 'y2', 'yc', '-mat')

end
