%% Choosing data --DONE
close all; clearvars;
addpath(genpath(fullfile('/scratch/FLEXL3L4/codes')));

BASE_FOLDER = fullfile('/scratch/FLEXL3L4/datasets/S3_ParisTrento_05_April_01-30_2019/L2_reprojected/Semana4');

LIST_VARIABLES = {'GLOB_isnan.tif','GLOB_LQSF.tif','GLOB_OTCI.tif'};
NUM_VAR = numel(LIST_VARIABLES);
disp('Getting data...');
for i = 1:NUM_VAR
    PATH{i} = fullfile(BASE_FOLDER,LIST_VARIABLES{i});
end
[nandata_data, R] = geotiffread(PATH{1});
lqsf_data = geotiffread(PATH{2});
otci_data = geotiffread(PATH{3});
%lqsf_data_nan = lqsf_data;
otci_data_nan = otci_data;
otci_data_nan(nandata_data==1) = NaN;
%lqsf_data(nandata_data==NaN) = 0;
lqsf_data(isnan(lqsf_data))= 0;
%lqsf_data_nan(nandata_data==1) = NaN;
 %% FILTERING 
 % Initialization
 [rows,cols,NUM] = size(lqsf_data); % Get size of data
 filter_mask = uint32(zeros(rows,cols,NUM)); % Create zero matrix for the filtered values 
 bin_min = zeros(1,10); % Create zero vector for binary values
 
 disp('Filtering...');
 % Filtering loop
 for i = 1:rows
    for j = 1:cols
        for k = 1:NUM 
            data_loc = lqsf_data(i,j,k); % Local data
            bin = flip(dec2bin(data_loc)-'0'); % Data to binary compounds
            bin(length(bin_min)) = 0; % Adjust length of binary values
            if isnan(otci_data_nan(i,j,k)) % Product value is not NaN
                filter_mask(i,j,k) = 0;
            else
                if (bin(1) ~= 1 && bin(10) ~= 1) % Value is not invalid or >70º
                    if bin(2) == 1 % Value is Water
                        filter_mask(i,j,k) = 1;
                    elseif bin(3) == 1 % value is Land
                        filter_mask(i,j,k) = 1;
                    elseif bin(5)== 1 % Value is Snow_ice
                        filter_mask(i,j,k) = 1;
                    else
                        filter_mask(i,j,k) = 0;
                    end
                else
                    filter_mask(i,j,k) = 0;
                end
            end
        end
    end
    %disp(i)
end
% 
% % 
% % %imshow(mat2gray(filter_mask(:,:,1))) % To print the filter mask
% % 
 fil_data = otci_data_nan.*single(filter_mask); % Filtered Data
 fil_labels = lqsf_data.*single(filter_mask);
% % 
% % for i = 1:rows
% %     for j = 1:cols
% %         for k = 1:NUM
% %             if isnan(final_matrix(i,j,k)) % Product value is not NaN
% %                 disp("NAN")
% %             end
% %         end
% %     end
% % end
% 
 NUM_OF_VALID_SAMPLES = sum(filter_mask,3); % Addition of the valid pixels in all the 3 dim
 
 %% TEMPORAL RESAMPLING
%  %final_matrix = zeros(rows,cols); % Create final matrix of values
  final_matrix = NaN(rows,cols); % Create final matrix of values
  value = 1; % Variable to increase index
 
 disp('temporal resampling...');
 for i = 1:rows
     for j = 1:cols
         valid_index = int32(zeros(1,NUM_OF_VALID_SAMPLES(i,j))); % Initialization of vector of valid index
         set_of_values = zeros(1,NUM_OF_VALID_SAMPLES(i,j)); % Initialization of valid values
        for k = 1:NUM % Index for a valid value
            if fil_labels(i,j,k) ~= 0 % If the value is valid, add it to the vector
                valid_index(value) = k;
                set_of_values(value) = otci_data_nan(i,j,k); % Set of correct values of OGVI
                value = value + 1;
            end
        end
        value = 1;
        %num = NUM_OF_VALID_SAMPLES(i,j)
        if NUM_OF_VALID_SAMPLES(i,j) > 4
            final_matrix(i,j) = median(set_of_values); % Median in case of more than 4 samples
        else
            if NUM_OF_VALID_SAMPLES(i,j) >= 2 % Only two samples, direct tree function
                [final_matrix(i,j), index1] = STC(fil_labels(i,j,valid_index(1)), fil_labels(i,j,valid_index(2)), otci_data_nan(i,j,valid_index(1)), otci_data_nan(i,j,valid_index(2)));
                %disp(index1)
            end
            if NUM_OF_VALID_SAMPLES(i,j) >= 3
                [final_matrix(i,j), index2] = STC(fil_labels(i,j,valid_index(3)), fil_labels(i,j,valid_index(index1)), otci_data_nan(i,j,valid_index(3)), otci_data_nan(i,j,valid_index(index1)));
            end
            if NUM_OF_VALID_SAMPLES(i,j) == 4
                [final_matrix(i,j), index3] = STC(fil_labels(i,j,valid_index(4)), fil_labels(i,j,valid_index(index2)), otci_data_nan(i,j,valid_index(4)), otci_data_nan(i,j,valid_index(index2)));       
            end
            if NUM_OF_VALID_SAMPLES(i,j) == 1
                final_matrix(i,j) = otci_data_nan(i,j,valid_index);
            end
        end
    end
    %disp(i)
 end
 
 %% Incertidumbre
 disp('creating uncertainty map...');
  incer_map = incer(NUM_OF_VALID_SAMPLES, otci_data_nan, fil_labels);
  min_incer = min(incer_map(:));
  max_incer = max(incer_map(:));
  save_png_map_max2(incer_map,'Incer_Ap4.png',NUM);
  geotiffwrite('Incer_Ap4_GEO',incer_map,R); 
  disp('creating valid pixels map...');
  geotiffwrite('Valid_Ap4_GEO',NUM_OF_VALID_SAMPLES,R); 
 %% Printing
  disp('saving data and creating composite...');
  min_val = min(final_matrix(:));
  max_val = max(final_matrix(:));
  save_png_map(final_matrix,'Ap4.png','OTCI');
  filename = 'Ap4.mat';
  save(filename,'-v7.3');
  geotiffwrite('Ap4_GEO',final_matrix,R); 
    
%% Unifing the size of the products (if necessary) and rotate
function [data,min_val,max_val] = unify_size_and_rotate(PF)

    NUM_PROD = numel(PF);
    [rows,cols] = size(PF{1});
    for i=2:NUM_PROD
        [r,c] = size(PF{i});
        if r<rows; rows=r; end
        if c<cols; cols=c; end
    end
        
    for i=1:NUM_PROD
        P = PF{i};
        P = P(1:rows,1:cols); % adjunsting size
        P = flip(P,2); P = rot90(P); % flip+rotate
        data(:,:,i) = P;
    end

    min_val = min(data(:));
    max_val = max(data(:));
    
end
%% Exporting maps into PNG files
function save_png_map(iP,FILE_OUT_PNG,VARIABLE)

    if nargin<3
        error('ERROR ''write_geotif'': too few parameters! Usage: save_png_map(iP,FILE_OUT_PNG,VARIABLE)');
    end

    
    if strcmp(VARIABLE,'OTCI')

        % https://github.com/sentinel-hub/custom-scripts/blob/master/sentinel-3/otci/script.js
        ranges = [ 0, 1, 1.8, 2.5, 4, 4.5, 5];
        base_colors= [
            0, 0,   0.5; ...
            0, 0.3, 0.8; ...
            1, 0.2, 0.2; ...
            1, 0.9, 0;   ...
            0, 0.8, 0.1; ...
            0, 0.6, 0.2; ...
            %1, 1,   1        
            0, 0.25,0.25
            ];

        cmap = [];    
        NUM_COLORS = 20;
        for i=1:(size(base_colors,1)-1)
            colors = colorGradient(base_colors(i,:),base_colors(i+1,:),NUM_COLORS);
            cmap = [cmap; colors];                
        end

        clims = [ranges(1), ranges(end)];
        fig=figure;
        set(fig,'visible','off');        
        colormap(fig,cmap); imagesc(iP,'AlphaData',double(~isnan(iP)),clims); % we set transparent 'nan' values
        set(gca,'fontsize', 20);
        set(gca,'visible','off');
        set(gcf,'units','points','position',[100,100,1024,1024]);
        colorbar;    
        export_fig(FILE_OUT_PNG,'-png','-nocrop'); % '-transparent'

    end

    
end
%% Función de arbol STC-S3
function [result, index] = STC(first_point, second_point, first_product, second_product)

    bin_min = zeros(1,20);            
    first_bin = flip(dec2bin(first_point)-'0'); % Data to binary compounds
    first_bin(length(bin_min)) = 0;
    second_bin = flip(dec2bin(second_point)-'0'); % Data to binary compounds
    second_bin(length(bin_min)) = 0;
    
    if (first_bin(3) == 1 || second_bin(3) == 1) % Any pixel is LAND?
        if (first_bin(3) == 1 && second_bin(3) == 1) % Both pixels are LAND?
            if (first_bin(18) ~= 1 || first_bin(17) ~= 1 || second_bin(18) ~= 1 || second_bin(17) ~= 1) % Pixels not classified as water,cloud or ice in OGVI
                if (first_bin(18) ~= 1 && first_bin(17) ~= 1 && second_bin(18) ~= 1 && second_bin(17) ~= 1)
                    if (first_product >= second_product) %OGVI1 >= OGVI2
                        result = first_product;
                        index = 2;
                    else
                        result = second_product;
                        index = 1;
                    end
                end
            elseif (first_bin(18) ~= 1 || first_bin(17) ~= 1) % First is not water, cloud or ice in OGVI
                result = first_product;
                index = 2;
            elseif (second_bin(18) ~= 1 || second_bin(17) ~= 1) % Second is not water, cloud or ice in OGVI
                result = second_product;
                index = 1;
            elseif (first_product >= second_product) % OGVI1 >= OGVI2
                result = first_product;
                index = 2;
            else
                result = second_product;
                index = 1;
            end
        end
        if (first_bin(3) == 1 || first_bin(5) == 1) % First is water or ice
            result = second_product;
            index = 1;
        else 
            result = first_product;
            index = 2;
        end
    end
    if (first_bin(5) == 1 || second_bin(5) == 1)
        if (first_bin(5) == 1 && second_bin(5) == 1)
            if (first_product >= second_product)
                result = first_product;
                index = 2;
            else
                result = second_product;
                index = 1;
            end
        end
        if (first_bin(2) == 1)
            result = second_product;
            index = 1;
        else
            result = first_product;
            index = 2;
        end
    end
    if (first_bin(2) == 1 || second_bin(2) == 1)
        if (first_bin(18) == 1 && second_bin(18) ~= 1)
            result = first_product;
            index = 2;
        elseif (first_bin(18) ~= 1 && second_bin(18) == 1)
            result = second_product;
            index = 1;
        elseif (first_product >= second_product)
            result = first_product;
            index = 2;
        else
            result = second_product;
            index = 1;
        end
        if (first_product >= second_product)
            result = first_product;
            index = 2;
        else 
            result = second_product;
            index = 1;
        end
    end
end

%% Incertidumbre 

function [incer_map] = incer(NUM_OF_VALID_SAMPLES, otci_data_nan, fil_labels)
    [rows, cols, NUM] = size(otci_data_nan);
    t_value = [13.97 4.53 3.31 2.87 2.65 2.52 2.43 2.37 2.32 2.28 2.25 2.23 2.21 2.20 2.18 2.17 2.16 2.15 2.14 2.0];
    incer_map = NaN(rows,cols);
    dev_2 = 0;
    var2 = 0;
    value = 1;
    for i = 1:rows
        for j = 1:cols
            if NUM_OF_VALID_SAMPLES(i,j) >= 2
                valid_index = int32(zeros(1,NUM_OF_VALID_SAMPLES(i,j))); % Initialization of vector of valid index
                set_of_values = zeros(1,NUM_OF_VALID_SAMPLES(i,j)); % Initialization of valid values
                for k = 1:NUM % Index for a valid value
                    if fil_labels(i,j,k) ~= 0 % If the value is valid, add it to the vector
                        valid_index(value) = k;
                        set_of_values(value) = otci_data_nan(i,j,k); % Set of correct values of OGVI
                        value = value + 1;
                    end
                end
                value = 1;
                NUM_2 = NUM_OF_VALID_SAMPLES(i,j);
                mean_value = mean(set_of_values);
                dev_2 = 0;
                for k = 1:NUM_2
                    var = (mean_value - set_of_values(k))^2;
                    dev_2 = dev_2 + var;
                end
                dev_1 = 1/(NUM_OF_VALID_SAMPLES(i,j)-1);
                dev_end = sqrt(dev_1 * dev_2);
                dev_m = dev_end/sqrt(NUM_OF_VALID_SAMPLES(i,j));
                incer_map(i,j)=exp(-t_value(NUM_2-1)*dev_m);
                if (incer_map(i,j)>1)
                    var2 = var2+1;
                end
            end
        end
    end
    %incer_map2 = log(incer_map+1)/log(max(incer_map(:)));
    %geotiffwrite('IncerAp1',incer_map,R); 
end


           