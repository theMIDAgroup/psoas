How to use the code

running Gui_psoas.m the user gets a GUI. 

First step:
By clicking the button "Folder" the user can select 
the folder containing the DICOM images of the psoas
(only the slices containing the psoas).

Second step:
By clicking the button "Centers" the image appears.
It should be an image containing a central section of the psoas.

Third step:
The user should click in the midpoint between the two psoas.
Then a prompt asks the user to select the image 
number where to identify the psoas centers.
Insert the slice number.

Fourth step:
Click first on the center of the psoas
in the left part of the image and then on 
the center of the psoas in the right part of the 
image. The user can now close the GUI.

The GUI then creates a folder, with the same name of
the one loaded by the user and suffix "_res",
where the information are stored in a .mat file.

Fifth step:
The user then can use one of the three proposed methods,
running one of the following commands
ac_classical(folder_data)
ac_gmfd_1ord(folder_data)
ac_gmfd_2ord(folder_data)
where folder_data is the folder containing the DICOM images,
e.g., folder_data = "/psoas/CT_001".

The parameters sigma and p are asked to be inserted by the user 
using the command window.
The user could use, for instance, sigma=1, p=2.

The results are all saved in the folder with suffix "_res" 
created by the GUI.

Sixth step:
Finally, the cleaning procedure is done running
clean_psoas_save(folder_data,print_surf)
where folder_data is the folder containing the DICOM images and 
print_surf needs to be 0 (not to print the recognized surface) or 1 (to print
the recognized surfaces).
The user is asked to insert the number of the slice they want
to use to start the cleaning procedure that can be repeated until
the user is satisfied.

The results are all saved in the folder with suffix "_res_clean".  
 




