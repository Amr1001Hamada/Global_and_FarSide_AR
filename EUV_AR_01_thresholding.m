

clear all
count  = 0;
count2 = 0;

load('map_304.mat');
load('TH.mat'); % [yyyymmdd tt TH]
load('seg_coord.mat'); %[yyyymmdd lat1 lat2 long1 long2]

%%
TH_smooth = smooth(TH (:,3),27*2);
TH_smooth(isnan(TH(:,3))) = NaN;

%% Start and End year
y1 = 2010; m1 = 05; d1 = 13;
y2 = 2016; m2 = 05; d2 = 14;

%% Arranging the date
D = datenum(y1,m1,d1):datenum(y2,m2,d2);
D = 10000*year(D) + 100*month(D) + day(D);
D = D';
num_days = length(D);

%% Make dir
mkdir([pwd '/TH_Stand_STEREO_JPL_07_200by500/fts/']);
mkdir([pwd '/TH_Stand_STEREO_JPL_07_200by500/png/']);
mkdir([pwd '/TH_Stand_STEREO_JPL_07_200by500/ARs/']);

%% STEREO-Date
stereo_d = datenum(y1,m1,d1):datenum(y2,m2,d2);
stereo_d = stereo_d';

%% Folder containing STEREO .fts files
folderPath = ([pwd '/TH_Stand_STEREO_JPL/fts_27/']) ;
ftsFiles = dir(fullfile(folderPath, '*.fts'));
ftsFileNames = {ftsFiles.name};

%% Loop through each STEREO .fts file
count = 0;
count_dd_tt = 0;

%% AR----AREA parameters
R_sun = 695700; % Solar radius in Km
A_sun = 4*pi * R_sun^2; % Total area of one hemesphere

A_sun  = A_sun / 2000;

A_sun = 3043.7 * 10^3; %1000 MH corresponding to 3.043,7 million square kilometers.

%  Area per pixel
num_lat_pixels = 200; % Total number of latitude pixels
num_long_pixels = 500; % Total number of longitude pixels
total_pixels = num_lat_pixels * num_long_pixels; % Total number of pixels

% Calculate the latitude values with sine transformation
latitude_values = linspace(-90, 90, num_lat_pixels); % Latitude values
latitude_values2 = linspace(-90, 90, num_lat_pixels + 1); % Latitude values

latitude_sine = sind(latitude_values); % Apply sine transformation
latitude_sine2 = sind(latitude_values2); % Apply sine transformation

% Calculate the width of latitude bands
lat_band_width = abs(diff(latitude_sine2)); % Width of each latitude band

% Calculate the area per pixel taking into account the latitude transformation
area_per_pixel = lat_band_width.* (A_sun / total_pixels) ;

Area_MOH_grid = ones(200,500).*area_per_pixel';

A = sum(Area_MOH_grid(:));

%% Create a waitbar
totalIterations = num_days*2; % Total number of iterations in your code
wait = waitbar(0, 'Processing...');

%%
for i = 1:num_days
    count = count + 1;
    file_date = num2str(D(i));

    for tt = [00 12]
        count_dd_tt = count_dd_tt+1;

        waitbar(count_dd_tt / totalIterations, wait, sprintf('Processing... %d%%', round(count_dd_tt / totalIterations * 100)));

        if tt == 0
            tt_str = '00';
        else
            tt_str = '12';
        end


        %% FIND the STEREO map match with DATA_date_time
        st_map_name = ['TH_stand_', num2str(D(i)), '_',tt_str, '*.fts']
        matchingFile = dir(fullfile(folderPath, st_map_name));

        try ftsFileName = matchingFile.name
            ftsFilePath = fullfile(folderPath, ftsFileName);

            im = fitsread(ftsFilePath);
            im = imresize(im,[200 500],'bilinear');
            im (im>0)=1;
            im = logical(im);

            im_euv = fitsread([pwd '/Stand_STEREO_JPL/fts/stand_',ftsFileName(10:end)]);
            im_euv = imresize(im_euv, [200 500]);

            im_euv_linear = fitsread([pwd '/STEREO_JPL/fts/',ftsFileName(10:end)]);
            im_euv_linear = imresize(im_euv_linear, [200 500]);
            mean_euv = nanmean(im_euv_linear(:));

            im_fqm = fitsread([pwd '/fqm/mrfqm',file_date(3:end),'t',tt_str,'00.fits']);
            im_fqm(im_fqm ~= 0) = 1;

            temp = im_fqm;
            temp = temp(any(temp, 2), :); % Extract non-zero rows
            temp = abs(temp);
            temp(temp>0)=1;
            temp (isnan(temp))=0;
            BWfqm = edge(temp,'Canny');
            [X1,Y1]= find(BWfqm == 1);

            ims = im;
            im2 = ims;

            %% Select the Best segment/TH
            AR_TH = TH_smooth (count_dd_tt,1);

            lat1 = seg_coord(count_dd_tt,3);
            lat2 = seg_coord(count_dd_tt,4);
            long1 = seg_coord(count_dd_tt,5);
            long2 = seg_coord(count_dd_tt,6);

            %% THRESHOLDING
            im_bw = im2;
            im_bw2 = im_bw;

            %% EDGING
            im_bw3 = im_bw2;
            im_bw4 = im_bw3;

            %% MASKING THE CH-FILTER ON THE ORIGINAL EIT-IMAGE
            im_edge = edge(im_bw4,'canny');
            [x_edge,y_edge]= find(im_edge == 1);

            %%
            im_x = size(im,2);
            im_y = size(im,1);

            long = linspace(0, 360,im_x);
            lat = linspace(-90, 90,im_y);
            sineLat = sind(lat);

            %% start ALLOCATING AR-groups.
            clear cc
            clear C0 C1 C2 C3_X C3_Long C4_Y C4_sinLat C4_Lat C5 C6 C7 C8 C9 C10; % Date, 'AR', 'Area(% n.of pix)', X, 'Cent.Long.', 'Y', 'Cent.SinLat.', 'Cent_Lat.', 'AR Tot_EUV_I', 'AR Mean_EUV_I', 'AR_tilt_angle','Long.Extent','Lat.Extent'};

            regionProps = regionprops(im_bw4, 'BoundingBox', 'Centroid','Orientation','PixelIdxList','BoundingBox'); % Label connected regions and get properties
            tilt_ang = cat(1,regionProps.Orientation);

            % Label connected components in the binary map
            im4labelled = bwlabel(im_bw4);


            cc = bwconncomp(im_bw4);      % returns the connected components CC found in the binary image im_bw
            AR_n = cc.NumObjects;        % Number of CH detected.
            AR_ndx = cc.PixelIdxList;    % List by CH indiceies.
            if AR_n>0
                S = regionprops(cc,'Centroid');         % Calculate centroids of the objects in the array.
                centroids = cat(1, S.Centroid);         % Concatenate structure array containing centroids into a single matrix.

                C1 = 1:AR_n; % ARs number
                C1 = C1';

                C0 = ones(size(C1,1),1) * D(count);

                C3_X = centroids(:,1);    % AR-X
                C3_Long = (C3_X - 1) * 360 / (500 - 1);    % AR-Longitude

                C4_Y = centroids(:,2);    % AR-Y
                C4_sinLat = sind( -90 + (C4_Y- 1) * (180 / (200 - 1)) );    % AR-sin(lat)
                C4_Lat = rad2deg(asin(C4_sinLat));    % AR-linear(lat)

                for F = 1:AR_n
                    A = sum(Area_MOH_grid(AR_ndx{F})); % AR-Area
                    C2(F,1) = A; % AR-%Area
                    C5(F,1) = mean(im_euv_linear(AR_ndx{F}),"omitnan")/mean_euv;          % AR-Brightness(% of the Background)
                    C6(F,1) = mean(im_euv_linear(AR_ndx{F}),"omitnan");         % AR-Mean(EUV-Intensity)
                    C7(F,1) = tilt_ang(F);         % AR-Mean(EUV-Intensity)


                    % Defining the scale of each pixel to the latitude/longitude
                    y_axis_range = 1 - (-1);
                    y_degrees_per_pixel = y_axis_range / im_y;

                    x_axis_range = 360 ;
                    x_degrees_per_pixel = x_axis_range / 500;

                    boundingBox = regionProps(F).BoundingBox;
                    LongExtent = boundingBox(3); % Width of the bounding box
                    LongExtent = LongExtent * x_degrees_per_pixel ;
                    C8(F,1) = LongExtent;

                    
                    L1 = boundingBox(2);
                    L1 = sind( -90 + (L1- 1) * (180 / (200 - 1)) );
                    L1 = rad2deg(asin(L1));

                    L2 = boundingBox(2) + boundingBox(4);
                    L2 = sind( -90 + (L2- 1) * (180 / (200 - 1)) );
                    L2 = rad2deg(asin(L2));

                    C9(F,1) = L2 - L1;

                    % Check if the AR is EarthSide or FarSide
                    centXY = round( centroids(F,:) );
                    if im_fqm (centXY(2), centXY (1)) ==1
                        C10(F,:) = 'FS';
                    else
                        C10(F,:) = 'ES';
                    end


                end


                Data = [C0 C1 round(C2) round(C3_X) round(C3_Long,1) round(C4_Y) round(C4_sinLat,2) round(C4_Lat,1) round(C5,3) round(C6,2) round(C7,1) round(C8,1) round(C9,1)];
            else
                C1=0;
                C2=0;
                C3_X=0;
                C3_Long=0;
                C4_Y=0;
                C4_sinLat=0;
                C4_Lat=0;
                C5=0;
                C6=0;
                C7=0;
                C8=0;
                C9=0;
                C10 = 'NaN';
            end

            %% Saving fits
            fitswrite(double(im_bw4), [pwd '/TH_Stand_STEREO_JPL_07_200by500/fts/TH_',ftsFileName]);

            %% PLOT
            f2 = figure('visible','off');
            % f2 = figure;
            ax1 = subplot (511);
            im_bw_show = double(im_bw);
            a = pcolor(im_bw_show); shading flat; grid on;
            colormap(ax1,gray)
            set(gca,'YDir','normal');
            set(gca,'YDir','normal');
            ylabel('Latitude[deg]', 'FontSize', 12);
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
            set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca,'TickDir','both');

            ax2 = subplot (512);
            im_bw2_show = double(im_bw2);
            im_bw2_show (isnan(im)) = NaN;
            a = pcolor(im_bw2_show); shading flat; grid on;
            colormap(ax2,gray)
            set(gca,'YDir','normal');
            set(gca,'YDir','normal');
            ylabel('Latitude[deg]', 'FontSize', 12);
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
            set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca,'TickDir','both');

            ax3 = subplot (513);
            im_bw3_show = double(im_bw4);
            im_bw3_show (isnan(im)) = NaN;
            a = pcolor(im_bw3_show); shading flat; grid on;
            colormap(ax3,gray)
            set(gca,'YDir','normal');
            set(gca,'YDir','normal');
            ylabel('Latitude[deg]', 'FontSize', 12);
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
            set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca,'TickDir','both');

            ax4 = subplot(514);
            a = pcolor(im_bw3_show); shading flat; grid on;
            colormap(ax4,gray)
            set(a,'alphaData',1-double(isnan(im)));
            set(gca,'YDir','normal');
            set(gca,'YDir','normal');
            %             title(['AR-mask [Morphological Op.] ',file_date,' - t00']);
            ylabel('Latitude[deg]', 'FontSize', 12);
            %             xlabel('Carring.Longitude[deg]', 'FontSize', 12);
            ax = gca;
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
            set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca,'TickDir','both');
            hold on;
            for n_region = 1:numel(regionProps) % Loop through each connected region
                rectangle('Position', regionProps(n_region).BoundingBox, 'EdgeColor', 'r', 'LineWidth', 2); % Draw bounding box
                centroid = regionProps(n_region).Centroid; % Label with centroid
                text(centroid(1), centroid(2), num2str(n_region), 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
            end
            hold off;

            ax5 = subplot(515);
            a = pcolor(im_euv);
            shading flat; grid on;
            set(gca,'YDir','normal');
            set(gca,'YDir','normal');
            colormap(ax5,map_304)
            caxis ([-1 2]);
            ylabel('Latitude[deg]', 'FontSize', 12);
            xlabel('Carring.Longitude[deg]', 'FontSize', 12);
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
            set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca,'TickDir','both');
            hold on
            scatter (y_edge, x_edge, 'b', '.')
            scatter(Y1, X1, 'k','.')

            screenSize = get(0, 'ScreenSize'); % Get screen size
            set(f2, 'Position', screenSize); % Set figure position to match screen size
            f2.Position(4) = 1000; % Expand the figure vertically
            f2.Position(3) = 700;  % Expand the figure horizontally
            saveas(f2,[pwd '/TH_Stand_STEREO_JPL_07_200by500/png/TH_',ftsFileName(1:end-4),'_02.png']);
            close;

            %% Extracting the AR information
            columnNames = {'Date', '#AR', 'Area(MH)', 'AR_X', 'AR.Long.', 'AR_Y', 'AR_SinLat.', 'AR_Lat', 'AR Brightness [% of Background]', 'AR Mean_EUV_I', 'AR_tilt_angle','Long.Extent','Lat.Extent','Solar.Side'};
            tableData = table(Data(:,1), Data(:,2), Data(:,3), Data(:,4), Data(:,5), Data(:,6), Data(:,7), Data(:,8), Data(:,9), Data(:,10), Data(:,11), Data(:,12), Data(:,13),C10, 'VariableNames', columnNames);
            writetable(tableData, [pwd '/TH_Stand_STEREO_JPL_07_200by500/ARs/ARs_',ftsFileName(1:end-4),'.csv']);




        catch
        end
    end
end

close(wait)