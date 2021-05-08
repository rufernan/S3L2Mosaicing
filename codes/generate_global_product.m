%% Generating a global product with the corresponding spatial limits given a folder with latitude/longitude geolocated products in '.nc' format

function [GP_OUT, R] = generate_global_product(path_str, field_product, interpolator, field_lon, field_lat)

    if nargin<2
        error('ERROR ''generate_global_product'': too few parameters! Usage: generate_global_product(path_str,field_product[,field_lon][,field_lat])');
    end
    if nargin<4
        interpolator = 'none'; % the pixels of the products are interpolated over the global grid using the indicated algorithm
    end       
    if nargin<4
        field_lon = 'lon'; % default field name for longitude data
    end
    if nargin<5
        field_lat = 'lat'; % default field name for latitude data 
    end          

    % BORRAR %
    %clearvars;
    %path_str = '/scratch/FLEXL3L4/datasets/S3_multiOrbit/L2_reprojected/S3*';
    %field_product = 'OTCI';
    %field_lon = 'lon';
    %field_lat = 'lat';
    %%%%%%%%%%

    VALID_INTERP = {'linear','nearest','natural','cubic','v4'};  % see https://es.mathworks.com/help/matlab/ref/griddata.html#d118e532760
    
    [folder,~,~] = fileparts(path_str);    
    names_struct = dir(path_str);
    file_names = {names_struct(:).name}';
    num_products = numel(file_names);
    
    if isempty(file_names)
        error('ERROR ''find_lonlat_limits'': no detected input products!');
    end   
    
    disp('finding the geo-spatial limits of the products...');    
    % 1st product initialization
    pname = (fullfile(folder,file_names{1}));
    disp(['--processing 01/',num2str(numel(file_names),'%02d'),'...']);
    [~,~,DX,DY,lonmin,lonmax,latmin,latmax] = extract_lonlat_nc(pname, field_lon, field_lat);
    GLON_MIN = lonmin; GLON_MAX = lonmax;
    GLAT_MIN = latmin; GLAT_MAX = latmax;    
    % ... and the rest
    for i=2:num_products
        pname = (fullfile(folder,file_names{i}));
        disp(['--processing ',num2str(i,'%02d'),'/',num2str(numel(file_names),'%02d'),'...']);
        [~,~,DX,DY,lonmin,lonmax,latmin,latmax] = extract_lonlat_nc(pname, field_lon, field_lat);
        if GLON_MIN>lonmin; GLON_MIN = lonmin; end
        if GLON_MAX<lonmax; GLON_MAX = lonmax; end
        if GLAT_MIN>latmin; GLAT_MIN = latmin; end
        if GLAT_MAX<latmax; GLAT_MAX = latmax; end
    end

    % Global reference
    R = makerefmat(GLON_MIN, GLAT_MAX, DX, DY);
        
    %LON_VEC = GLON_MIN:DX:GLON_MAX+DX-eps; % we add (DX-eps) to introduce an additional column in case we reach a position close to GLON_MAX using a DX step
    %num_cols = numel(LON_VEC);
    %LAT_VEC = GLAT_MAX:DY:GLAT_MIN-DY+eps; % reverse
    %num_rows = numel(LAT_VEC);
    
    LON_VEC = GLON_MIN:DX:GLON_MAX+DX; % we add DX to guarantee that GLON_MAX is always reached
    num_cols = numel(LON_VEC);
    LAT_VEC = GLAT_MAX:DY:GLAT_MIN+DY; % idem (note that DY is negative and we want to reach GLAT_MIN in this case)
    num_rows = numel(LAT_VEC);
    
    % Output joint product
    GP_OUT = single(NaN([num_rows,num_cols,num_products]));
        
    disp('loading the data products...');
    for i=1:num_products
        
        pname = (fullfile(folder, file_names{i}));                        
        disp(['--loading ',num2str(i,'%02d'),'/',num2str(numel(file_names),'%02d'),'...']);
        P = ncread(pname, field_product);
        P = unify_size_and_rotate({P});
        [rows,cols,~] = size(P);
        [P_LON,P_LAT,~,~,~,~,~,~] = extract_lonlat_nc(pname, field_lon, field_lat);             
        
        % finding the closest initial global longitude/latitude
        dist = pdist2(LON_VEC(:),P_LON(1));
        [~,ind_lon] = min(dist(:));
        dist = pdist2(LAT_VEC(:),P_LAT(1));
        [~,ind_lat] = min(dist(:));
        

        % if required interpolating the products over the global grid 
        
        if any(contains(VALID_INTERP,interpolator))                
            
            P_LAT_vec = P_LAT(:);
            P_LON_vec = P_LON(:);
            
            [x1_grid,x2_grid] = meshgrid(P_LAT_vec,P_LON_vec);
            x1_grid = x1_grid';
            x2_grid = x2_grid';
                        
            x1_grid_vec = x1_grid(:);
            x2_grid_vec = x2_grid(:);            
            P_vec = P(:);
            
            x1_grid_vec_nn = x1_grid_vec(not(isnan(P_vec)));
            x2_grid_vec_nn = x2_grid_vec(not(isnan(P_vec)));
            P_vec_nn = P_vec(not(isnan(P_vec)));
                        
            [xq1_grid,xq2_grid] = meshgrid(LAT_VEC(ind_lat:ind_lat+rows-1),LON_VEC(ind_lon:ind_lon+cols-1));
            xq1_grid = xq1_grid';
            xq2_grid = xq2_grid';
            
            disp('----interpolating the product to the global grid...');
            vq = griddata(x1_grid_vec,x2_grid_vec,P_vec,xq1_grid,xq2_grid,interpolator);            
            disp('----done!');
            
            P = vq';
                        
        end
        
        
         GP_OUT(ind_lat:ind_lat+rows-1,ind_lon:ind_lon+cols-1,i) = P;
        
            
    end
                         
end







