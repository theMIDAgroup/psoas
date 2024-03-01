%% La funzione viene chiamata dalla GUI con la selezione del comando Folder
%  La chiamata alla funzione richiede la selezione della cartella dove sono
%  presenti i dati dei pazienti. Successivamente crea, se non presenti di
%  default, le cartelle relative al salvataggio dei risultati, una per ogni
%  paziente. Viene inoltre creata una cartella per il salvataggio delle
%  figure generate dal codice. 

function load_data(handles)
global folder_data


folder_data=uigetdir;
set(handles.text2, 'string', folder_data)

if folder_data ==  0
else
    handles.text2 = folder_data;
    
    
end
