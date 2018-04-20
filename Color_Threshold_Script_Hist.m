function [R_N_count, R_CT_count, G_N_count, G_CT_count, B_N_count, B_CT_count, total_N_count, total_CT_count, sample_name]=color_threshold(input)

input_image=imread(input,'tif'); %read in image
%split color image into RGB Channels
Red=input_image(:,:,1);
Green=input_image(:,:,2);
Blue=input_image(:,:,3);
 
%Get threshold values

nucleus_red=(98/255)*255;
nucleus_green=(82/255)*255;
nucleus_blue=(126/255)*255;

cytoplasm_red=(160/255)*255;
cytoplasm_green=(160/255)*255;
cytoplasm_blue=(160/255)*255;

%Apply threshold by RGB channel 

Red_nucleus= Red<nucleus_red;
Red_cytoplasm_tendrils= Red<cytoplasm_red;
Green_nucleus= Green<nucleus_green;
Green_cytoplasm_tendrils= Green<cytoplasm_green;
Blue_nucleus= Blue<nucleus_blue;
Blue_cytoplasm_tendrils= Blue<cytoplasm_blue;

%Make and apply threshold mask to original image 

CT_Red= uint8(Red_cytoplasm_tendrils);
CT_Green= uint8(Green_cytoplasm_tendrils);
CT_Blue= uint8(Blue_cytoplasm_tendrils);

N_Red= uint8(Red_nucleus);
N_Green= uint8(Green_nucleus);
N_Blue= uint8(Blue_nucleus);

Nucleus_Red=Red .*N_Red;
Nucleus_Green=Green .*N_Green;
Nucleus_Blue=Blue .*N_Blue;

CyT_Red=Red .*CT_Red;
CyT_Green=Green .*CT_Green;
CyT_Blue=Blue .*CT_Blue;

%generate mask and output as a .tif

input_image(:,:,1)=CyT_Red;
input_image(:,:,2)=CyT_Green;
input_image(:,:,3)=CyT_Blue;
figure
imshow(input_image);

%get total counts for nucleus values 

[R_N_counts, R_N_pixel_values]= imhist(Nucleus_Red);
[G_N_counts, G_N_pixel_values]= imhist(Nucleus_Green);
[B_N_counts, B_N_pixel_values] = imhist(Nucleus_Blue);
[R_CT_counts, R_CT_pixel_values]= imhist(CyT_Red);
[G_CT_counts, G_CT_pixel_values]= imhist(CyT_Green);
[B_CT_counts, B_CT_pixel_values] = imhist(CyT_Blue);
R_N_count=sum(R_N_counts(2:end,1));
G_N_count =sum(G_N_counts(2:end,1));
B_N_count =sum(B_N_counts(2:end,1));
R_CT_count=sum(R_CT_counts(2:end,1));
G_CT_count =sum(G_CT_counts(2:end,1));
B_CT_count =sum(B_CT_counts(2:end,1));
total_N_count=R_N_count+G_N_count+B_N_count;
total_N_pixel_counts=R_N_counts+G_N_counts+B_N_counts;
total_CT_count=R_CT_count+G_CT_count+B_CT_count;
total_CT_pixel_counts=R_CT_counts+G_CT_counts+B_CT_counts;

%make unique name for outfiles based on file name
strname=strsplit(input,'_');
sample_name=char(strcat(strname(2),{' '},strname(3)));
filename_nucleus=[sample_name '_nucleus' '.txt'];
filename_cytoplasm_tendrils= [sample_name '_cytoplasm_tendrils' '.txt'];

%open outfile and make headers
fileID=fopen(filename_nucleus,'wt');
header1='Red Pixel Counts';
header2='Green Pixel Counts';
header3='Blue Pixel Counts';
header4='Total Pixel Counts';
header5='Sum';
fprintf(fileID,[header1 '\t' header2 '\t' header3 '\t' header4 '\t' header5 '\n']);
fprintf(fileID, '%f \t %f \t %f \t %f \t %f \n',[R_N_counts, G_N_counts, B_N_counts, total_N_pixel_counts, R_N_pixel_values]');
fprintf(fileID, '%f \t %f \t %f \t %f \t %f \n', [R_N_count, G_N_count, B_N_count, total_N_count, 255]');
fclose(fileID)

fileID2=fopen(filename_cytoplasm_tendrils,'wt');
fprintf(fileID2,[header1 '\t' header2 '\t' header3 '\t' header4 '\t' header5 '\n']);
fprintf(fileID2, '%f \t %f \t %f \t %f \t %f \n',[R_CT_counts, G_CT_counts, B_CT_counts, total_CT_pixel_counts, R_CT_pixel_values]');
fprintf(fileID, '%f \t %f \t %f \t %f \t %f \n', [R_CT_count, G_CT_count, B_CT_count, total_CT_count, 255]');
fclose(fileID2)
end