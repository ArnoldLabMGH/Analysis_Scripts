function [pixel_count_total, pixel_count_background, pixel_count_auto, pixel_count_nucleus, sample_name]=segmented_histogram(input)
input_image=imread(input,'tif'); %read in image
input_grayscale=rgb2gray(input_image); %convert to grayscale
[pixel_count, pixel_values]=imhist(input_grayscale); %get histogram data of grayscale image
pixel_count_total=sum(pixel_count(2:end,1)); 
pixel_matrix_background=im2bw(input_grayscale, (203/255)); 
pixel_count_background=sum(sum(pixel_matrix_background));
pixel_matrix_auto= im2bw(input_grayscale, (160/255));

figure
imshow(pixel_matrix_auto)

pixel_count_auto= sum(sum(pixel_matrix_auto));
pixel_matrix_nucleus=im2bw(input_grayscale, (76/255));

figure
imshow(pixel_matrix_nucleus)
pixel_count_nucleus=sum(sum(pixel_matrix_nucleus));
strname=strsplit(input,'_');
sample_name=char(strcat(strname(2),{' '},strname(3)));
filename=[sample_name '.txt'];
fileID=fopen(filename,'wt');
header1='Pixel Counts';
header2='Pixel Values';
fprintf(fileID,[header1 '\t' header2 '\n']);
fprintf(fileID, '%f \t %f \n',[pixel_count, pixel_values]');
fprintf(fileID, '%f \t %f \n', [pixel_count_total, 255]');
fclose(fileID)
end