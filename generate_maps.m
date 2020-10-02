close all; clearvars;

addpath(genpath('/scratch/FLEXL3L4/codes'));

FILE_IN = 'median_OTCI_reprojected.tif';
[A,~] = geotiffread(FILE_IN);    
FILE_OUT = strrep(FILE_IN,'tif','png');
save_png_map(A,FILE_OUT,'OTCI');
    
FILE_IN = 'mean_OTCI_reprojected.tif';
[A,~] = geotiffread(FILE_IN);    
FILE_OUT = strrep(FILE_IN,'tif','png');
save_png_map(A,FILE_OUT,'OTCI');

FILE_IN = 'GlobalJan_GEO_reprojected.tif';
[A,~] = geotiffread(FILE_IN);    
FILE_OUT = strrep(FILE_IN,'tif','png');
save_png_map(A,FILE_OUT,'OTCI');

FILE_IN = 'ValidJan_GEO_reprojected.tif';
[A,~] = geotiffread(FILE_IN);    
FILE_OUT = strrep(FILE_IN,'tif','png');
save_png_map_val(A,FILE_OUT,15);

FILE_IN = 'Incer_Jan_reprojected.tif';
[A,~] = geotiffread(FILE_IN);    
FILE_OUT = strrep(FILE_IN,'tif','png');
save_png_map_max(A,FILE_OUT,10);

FILE_IN = 'Incer_Jan2_reprojected.tif';
[A,~] = geotiffread(FILE_IN);    
FILE_OUT = strrep(FILE_IN,'tif','png');
save_png_map_max2(A,FILE_OUT,10);