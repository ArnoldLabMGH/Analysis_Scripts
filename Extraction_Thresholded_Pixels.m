names=dir('*.tif'); %looks for files in current directory with .tif extension
name_list= {names.name}.'; %extract a transposed column file in name_list
outfile=fopen('aggregated_file.txt','wt');
header1='Total Pixel Counts';
header2='Background Pixel Counts';
header3='Auto Threshold Pixel Counts';
header4='Nucleus Pixel Counts';
header5= 'Sample_ID';
fprintf(outfile, [header1 '\t' header2 '\t' header3 '\t' header4 '\t' header5 '\n']');
for k=1:length(name_list)
    input=name_list{k,:}; %extract the first element of the structure
    [pixel_count_total, pixel_count_background, pixel_count_auto, pixel_count_nucleus, sample_name]=Threshold_Histogram(input);
    outfile=fopen('aggregated_file.txt', 'a');
    fprintf(outfile, '%f \t %f \t %f \t %f \t',[pixel_count_total, pixel_count_background, pixel_count_auto, pixel_count_nucleus]');
    fprintf(outfile, '%s \n', [sample_name]');
end
fclose(outfile)