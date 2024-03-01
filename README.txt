How to use the code

running Gui_psoas.m the user get a GUI. 

First step:
Clicking the button "load data" the user can select 
the folder containing the DICOM images of the psoas.

Second step:
Clicking the button "centers" a prompt asks the user 
to select the image number where identify the 
psoas centers. It shold be an image containing a 
central section of the psoas. Then the user should click
in the midpoint between the two psoas and then click 
first on the center of the psoas in the 
left part of the image and then on 
the center of the psoas in the right part of the 
image one.

The GUI then create a folder, with the same name of
the one loaded by the user and suffix "_res"
where the information are stored in un .mat file.

The user than can use one of the three proposed methods,
runninig one of the following command 
ac_classical(folder_data)
ac_gmfd_1ord(folder_data)
ac_gmfd_2ord(folder_data)
where folder_data is the folder containing the DICOM images

The parameters sigma and p are asked to be prompt by the user using 
the command window.

The results are all saved in the folder with suffix "_res" 
created by the GUI.

Finally, the cleanining procedure is done runing
clean_psoas_save(folder_data,print_surf)
where folder_data is the folder containing the DICOM images and 
print_surf needd to be 0 (to not print the recognized surface) or 1 (to print
the recognized surfaces).
The user is asked to insert the number of the slice they want
to use to start the cleaning procedure that can be repeated until
the user is satisfied.

The results are all saved in the folder with suffix "_res_clean"  
 




