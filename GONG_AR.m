


clear all
load ('map_304.mat');

% Create directories for outputs
mkdir([pwd '/TH_Avg_fqm_10/png']);
mkdir([pwd '/TH_Avg_fqm_10/fits/']);
mkdir([pwd '/TH_Avg_fqm_10/csv/']);

%% Start and End year
y1 = 2010; m1 = 05; d1 = 13;
y2 = 2016; m2 = 05; d2 = 14;

%% Arranging the date
D = datenum(y1,m1,d1):datenum(y2,m2,d2);
D = 10000*year(D) + 100*month(D) + day(D);
D = D';
file_dates = D;   
file_dates = num2str(file_dates);

num_days = length(D);

stereo_d = datenum(y1,m1,d1):datenum(y2,m2,d2);
stereo_d = stereo_d';

%% Folder containing FarSide fqm .fts files
folderPath = ([pwd '/Avg_fqm_f6u/filtered_fits/']) ;
ftsFiles = dir(fullfile(folderPath, '*.fits'));
ftsFileNames = {ftsFiles.name};

%% Folder containing TH-EUV .fts files
folderPath_euv = ([pwd '/TH_Stand_STEREO_JPL/fts_27/']) ;
ftsFiles_euv = dir(fullfile(folderPath_euv, '*.fts'));
ftsFileNames_euv = {ftsFiles_euv.name};

%% Folder containing STAND-EUV .fts files
folderPath_euv2 = ([pwd '/Stand_STEREO_JPL/fts/']) ;
ftsFiles_euv2 = dir(fullfile(folderPath_euv2, '*.fts'));
ftsFileNames_euv2 = {ftsFiles_euv2.name};

%% Initialize variables for processing
count = 0;
count_dd_tt = 0;

%% AR----AREA 
% Solar parameters
R_sun = 695700; % Solar radius in Km
A_sun = 3043.7 * 10^3; %1000 MH corresponding to 3.043,7 million square kilometers. 

% Define pixel grid
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

Area_MOH_grid = ones(num_lat_pixels, num_long_pixels).*area_per_pixel';

A = sum(Area_MOH_grid(:));

%% Create a waitbar
totalIterations = num_days*2; % Total number of iterations in your code
wait = waitbar(0, 'Processing...');

%% Main loop for processing
for i = 1:num_days
    count = count + 1;
    file_date = file_dates(i,:);

    for tt = [00 12]
        count_dd_tt = count_dd_tt+1;

        waitbar(count_dd_tt / totalIterations, wait, sprintf('Processing... %d%%', round(count_dd_tt / totalIterations * 100)));

        if tt == 0
            tt_str = '0000';
        else
            tt_str = '1200';
        end

        %% Loading STEREO-Ar Binary map
        
        % TH EUV Map
        euv_map_name1 = ['TH_stand_', file_date, '_',tt_str(1:2), '*.fts'];
        matchingFile_euv1 = dir(fullfile(folderPath_euv, euv_map_name1));

        try ftsFileName_euv = matchingFile_euv1.name;
            ftsFilePath_euv = fullfile(folderPath_euv, ftsFileName_euv);

            st = fitsread(ftsFilePath_euv);
            if ~all(isnan(st(:)))
                st1 = imresize(st,[200,500],'bilinear');
                st1 (st1>0)=1;
                st1 = logical(st1);

                cc_st = bwconncomp(st1);      % returns the connected components CC found in the binary image im_bw
                regionProps_euv = regionprops(st1, 'BoundingBox', 'Centroid','Orientation','PixelIdxList'); % Label connected regions and get properties
                centroids_euv = cat(1, regionProps_euv.Centroid);         % Concatenate structure array containing centroids into a single matrix.

                ARnumberStereo = cc_st.NumObjects;

            end

             % EUV Map
            euv_map_name2 = ['stand_', file_date, '_',tt_str(1:2), '*.fts'];
            matchingFile_euv2 = dir(fullfile(folderPath_euv2, euv_map_name2));

            ftsFileName_euv2 = matchingFile_euv2.name;
            ftsFilePath_euv2 = fullfile(folderPath_euv2, ftsFileName_euv2);

            im_euv = fitsread(ftsFilePath_euv2);
            im_euv = imresize(im_euv,[200,500],'bilinear');

        catch
            st1 = zeros(200,500);
            im_euv = nan(200,500);

            cc_st = bwconncomp(st1);      % returns the connected components CC found in the binary image im_bw
            regionProps_euv = regionprops(st1, 'BoundingBox', 'Centroid','Orientation','PixelIdxList'); % Label connected regions and get properties
            centroids_euv = cat(1, regionProps_euv.Centroid);         % Concatenate structure array containing centroids into a single matrix.

            ARnumberStereo = cc_st.NumObjects;
        end


        %% FIND the fqm map match with DATA_date_time
        fqm_map_name = ['mrf6u', file_date(3:end), 't',tt_str, '_fs.fits'];
        matchingFile = dir(fullfile(folderPath, fqm_map_name));

        try ftsFileName = matchingFile.name
            ftsFilePath = fullfile(folderPath, ftsFileName);

            avg_sin = fitsread(ftsFilePath);
            avg_sin (avg_sin == 0) = NaN;
            if ~all(isnan(avg_sin(:)))

                % % edge the boundaries of each map
                temp1 = avg_sin;
                temp1 = temp1(any(temp1, 2), :); % Extract non-zero rows
                temp1 = abs(temp1);
                temp1(temp1>0)=1;
                temp1 (isnan(temp1))=0;
                BW1 = edge(temp1,'Canny');
                [x1,y]= find(BW1 == 1);

                im_x = size(avg_sin,2);
                im_y = size(avg_sin,1);

                fqm_th = -1* nanstd(avg_sin(:))*2;
                bw0 = avg_sin < fqm_th;

                %% Latitudinal Filttering of the areas
                cc0 = bwconncomp(bw0);      % returns the connected components CC found in the binary image im_bw
                regionProps0 = regionprops(bw0, 'BoundingBox', 'Centroid','Orientation','PixelIdxList'); % Label connected regions and get properties
                centroids0 = cat(1, regionProps0.Centroid);         % Concatenate structure array containing centroids into a single matrix.

                long = linspace(0, 360,im_x);
                lat = linspace(-90, 90,im_y);
                sineLat = sind(lat)';

                % % Setting the Latitude boundaries
                sinPlus50 = sind(50);
                sinMinus50 = sind(-50);

                % % Removing all the regions outside the Lat.Range
                [~, ndxPositive] = min(abs(sineLat - sinPlus50));
                [~, ndxNegative] = min(abs(sineLat - sinMinus50));

                areasToRemove = find(centroids0(:, 2) < ndxNegative | centroids0(:, 2) > ndxPositive);

                bw1 = bw0;
                for J = 1:numel(areasToRemove)
                    bw1(cc0.PixelIdxList{areasToRemove(J)}) = 0;
                end

                %% Area Threshold - 02 [Removing all the AR-areas < 10 MH]
                bw2 = imfill(bw1,'holes'); %  fills holes in im_bw
                
                bw3 = bw2;
                cc_fqm = bwconncomp(bw3);      % returns the connected components CC found in the binary image im_bw
                AR_n = cc_fqm.NumObjects;        % Number of CH detected.
                if AR_n > 0
                    AR_ndx = cc_fqm.PixelIdxList;    % List by CH indiceies.

                    for F = 1:AR_n
                        AR_area = sum(Area_MOH_grid(AR_ndx{F})); % AR-Area (MH)
                        if AR_area < 10
                            bw3(AR_ndx{F}) = 0;
                        end
                    end
                end

                %%
                bw4 = bw3;
                clear cc_fqm
                clear C0 C1 C2 C3 C4 C5 C6 C7 C8 C9;
                clear AR_n AR_ndx F

                cc_fqm = bwconncomp(bw4);      % returns the connected components CC found in the binary image im_bw
                regionProps_fqm = regionprops(bw4, 'BoundingBox', 'Centroid','Orientation','PixelIdxList'); % Label connected regions and get properties
                centroids_fqm = cat(1, regionProps_fqm.Centroid);         % Concatenate structure array containing centroids into a single matrix.

                AR_n = cc_fqm.NumObjects;        % Number of CH detected.
                AR_ndx = cc_fqm.PixelIdxList;    % List by CH indiceies.

                avg_sin2 = avg_sin;
                avg_sin2(avg_sin2>0)=NaN;
                if AR_n > 0
                    C1 = 1:AR_n; % ARs number
                    C1 = C1';

                    C0 = ones(size(C1,1),1) * str2num(file_date);

                    for F = 1:AR_n
                        C2(F,1) = sum(Area_MOH_grid(AR_ndx{F}));                  % AR-Area (MH)
                        C5(F,1) = sum(avg_sin2(AR_ndx{F}),"omitnan");          % Sum PhaseShift
                        C6(F,1) = mean(avg_sin2(AR_ndx{F}),"omitnan");         % Mean Phaseshift Value
                        C7(F,1) = nanmin(avg_sin2(AR_ndx{F}));             % Minimum PhaseShift value
                        C8(F,1) = nanstd(avg_sin2(AR_ndx{F}));             % Minimum PhaseShift value
                    end

                    C3_X = centroids_fqm(:,1);    % AR-X
                    C3_Long = (C3_X - 1) * 360 / (im_x - 1);    % AR-Longitude

                    C4_Y = centroids_fqm(:,2);    % AR-Y
                    C4_sinLat = sind( -90 + (C4_Y- 1) * (180 / (im_y - 1)) );    % AR-sin(lat)

                    C4_Lat = rad2deg(asin(C4_sinLat));    % AR-lat

                    Data = [C0 C1 C2 C3_X C3_Long C4_Y C4_sinLat C4_Lat C5 C6 C7 C8];

                else
                    C0 = str2num(file_date);                    
                    C1=NaN;
                    C2=NaN;
                    C3_X=NaN;
                    C3_Long=NaN;
                    C4_Y=NaN;
                    C4_sinLat=NaN;
                    C4_Lat=NaN;
                    C5=NaN;
                    C6=NaN;
                    C7=NaN;
                    C8=NaN;
                    C9=NaN;
                end


            end
        catch
            avg_sin = nan(200,500);
            bw0 = nan(200,500);
            bw1 = nan(200,500);
            bw3 = nan(200,500);
            bw4 = nan(200,500);

            C0 = str2num(file_date);
            C1=NaN;
            C2=NaN;
            C3_X=NaN;
            C3_Long=NaN;
            C4_Y=NaN;
            C4_sinLat=NaN;
            C4_Lat=NaN;
            C5=NaN;
            C6=NaN;
            C7=NaN;
            C8=NaN;
            C9=NaN;
        end

        %% PLOTTING
        yticklabels = -90:45:90;
        yticklabels_sine = -1:0.5:1;
        yticks = linspace (1,size(avg_sin,1), numel(yticklabels_sine));
        xticklabels = 0:60:360;
        xticks = linspace (1,size(avg_sin,2), numel(xticklabels));

        f = figure('visible','off');
        % f = figure
        ax1 = subplot(5,2,1);
        pcolor(avg_sin);
        shading flat;
        grid on
        colormap(ax1,'bone');
        title(['Average FS-Map - \color{black}',file_date , 't',tt_str, '\color{black}'], 'FontSize', 12);
        caxis([-0.15 0.15])
        ylabel('sin(Lat.)', 'FontSize', 14);
        ax = gca;
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels_sine)
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca(), 'XGrid','on', 'YGrid','on', 'Layer','top', 'TickDir','out')
        hold on
        line(xlim, [ndxPositive,ndxPositive],'color','red','LineWidth',1)
        line(xlim, [ndxNegative,ndxNegative],'color','red','LineWidth',1)
        text(10, 180, '(a)', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold');

        subplot(5,2,3)
        bw0 = double (bw0);
        bw0(isnan(avg_sin))= NaN;
        pcolor(bw0);
        shading flat;
        grid on
        ylabel('sin(Lat.)', 'FontSize', 12);
        ax = gca;
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels_sine)
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca(), 'XGrid','on', 'YGrid','on', 'Layer','top', 'TickDir','out')
        hold on
        line(xlim, [ndxPositive,ndxPositive],'color','red','LineWidth',1)
        line(xlim, [ndxNegative,ndxNegative],'color','red','LineWidth',1)
        text(10, 180, '(b)', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');

        subplot(5,2,5)
        bw1 = double (bw1);
        bw1(isnan(avg_sin))= NaN;
        pcolor(bw1);
        shading flat;
        grid on
        ylabel('sin(Lat.)', 'FontSize', 12);
        ax = gca;
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels_sine)
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca(), 'XGrid','on', 'YGrid','on', 'Layer','top', 'TickDir','out')
        hold on
        line(xlim, [ndxPositive,ndxPositive],'color','red','LineWidth',1)
        line(xlim, [ndxNegative,ndxNegative],'color','red','LineWidth',1)
        text(10, 180, '(c)', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');

        subplot(5,2,7)
        bw3 = double (bw3);
        bw3(isnan(avg_sin))= NaN;
        pcolor(bw3);
        shading flat
        grid on
        ylabel('sin(Lat.)', 'FontSize', 12);
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels_sine)
        set(gca(), 'XGrid','on', 'YGrid','on', 'Layer','top', 'TickDir','out')
        hold on
        line(xlim, [ndxPositive,ndxPositive],'color','red','LineWidth',1)
        line(xlim, [ndxNegative,ndxNegative],'color','red','LineWidth',1)
        text(10, 180, '(d)', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');


        subplot(5,2,9);
        bw4_double = double (bw4);
        bw4_double(isnan(avg_sin))= NaN;
        pcolor(bw4_double);
        shading interp
        colormap('bone')
        grid on
        ylabel('sin(Lat.)', 'FontSize', 12);
        xlabel('Heliographic Longitude[deg]', 'FontSize', 12);
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels_sine)
        set(gca(), 'XGrid','on', 'YGrid','on', 'Layer','top', 'TickDir','out')
        hold on
        line(xlim, [ndxPositive,ndxPositive],'color','red','LineWidth',1)
        line(xlim, [ndxNegative,ndxNegative],'color','red','LineWidth',1)


        if  ~all(isnan(bw4(:)))

            if any(st1(:))
                % % Compute bounding boxes of labeled regions
                bbox_fqm = regionprops(bw4, 'BoundingBox');
                bbox_euv = regionprops(st1, 'BoundingBox');

                % % Initialize array to store overlapped areas
                overlapped_areas = zeros(numel(bbox_fqm), 1);

                % % Iterate through each region in the fqm map
                for n_fqm = 1:numel(bbox_fqm)
                    centroid_fqm = regionProps_fqm(n_fqm).Centroid; % Label with centroid
                    distance = sqrt((centroid_fqm(1) - centroids_euv(:,1)).^2 + (centroid_fqm(2) - centroids_euv(:,2)).^2);
                    ndx = find(distance < 30, 1);
                    
                    bbox_fqm_one = bbox_fqm(n_fqm).BoundingBox;

                    % Find corresponding region in the second image based on maximum overlap
                    max_overlap_area = 0;
                    for n_euv = 1:numel(bbox_euv)
                        bbox_euv_one = bbox_euv(n_euv).BoundingBox;

                        % Compute overlap area between bounding boxes
                        overlap_area = rectint(bbox_fqm_one, bbox_euv_one);

                        % Update maximum overlap area
                        if overlap_area > max_overlap_area
                            max_overlap_area = overlap_area;
                        end
                    end
                    % Store maximum overlap area for the current region
                    overlapped_areas(n_fqm) = max_overlap_area;

                    if ~isempty(ndx) || overlapped_areas(n_fqm) > 0
                        rectangle('Position', regionProps_fqm(n_fqm).BoundingBox, 'EdgeColor', 'g', 'LineWidth', 2, 'FaceColor', 'none'); % Draw bounding box
                        text(centroid_fqm(1), centroid_fqm(2), num2str(n_fqm), 'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold');
                        C9(n_fqm,1) = 1;
                    else
                        rectangle('Position', regionProps_fqm(n_fqm).BoundingBox, 'EdgeColor', 'r', 'LineWidth', 2, 'FaceColor', 'none'); % Draw bounding box
                        text(centroid_fqm(1), centroid_fqm(2), num2str(n_fqm), 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
                        C9(n_fqm,1) = 0;
                    end
                end


            else % in case NO EUV-maps
                for n_fqm = 1:numel(regionProps_fqm) % Loop through each connected region
                    centroid_fqm = regionProps_fqm(n_fqm).Centroid; % Label with centroid
                    rectangle('Position', regionProps_fqm(n_fqm).BoundingBox, 'EdgeColor', 'y', 'LineWidth', 2, 'FaceColor', 'none'); % Draw bounding box
                    text(centroid_fqm(1), centroid_fqm(2), num2str(n_fqm), 'Color', 'y', 'FontSize', 12, 'FontWeight', 'bold');
                    C9(n_fqm,1) = NaN;
                end
            end

        else

            C9 = NaN;
        end
        hold off;
        text(10, 180, '(e)', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');


        subplot(5,2,2)
        [n_im,x_im]=hist(reshape(avg_sin,numel(avg_sin),1),256);
        plot (x_im,n_im,'k','LineWidth',2);hold on
        grid on
        axis tight
        xlabel('PhaseShift', 'FontSize', 14)
        ylabel('Num. of pix', 'FontSize', 14);
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca,'TickDir','out');
        line([fqm_th, fqm_th], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
        xlim([-0.15 0.15])
        ax8 = subplot(5,2,8);
        pcolor(im_euv);
        shading flat;
        grid on
        colormap(ax8, map_304);
        title(['EUV Synchronic Map - ',file_date ,'t', tt_str], 'FontSize', 12);
        caxis([-2 3])
        ylabel('sin(Lat.)', 'FontSize', 14);
        ax = gca;
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels_sine)
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca(), 'XGrid','on', 'YGrid','on', 'Layer','top', 'TickDir','out')
        hold on
        line(xlim, [ndxPositive,ndxPositive],'color','w','LineWidth',1)
        line(xlim, [ndxNegative,ndxNegative],'color','w','LineWidth',1)
        text(10, 180, '(a)', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold');
        
        st_bw = double(st1);
        st_bw (isnan(im_euv)) = NaN;
        ax10 = subplot(5,2,10);
        pcolor(st_bw);
        shading flat; grid on;
        colormap(ax10, gray)
        set(gca,'YDir','normal');
        set(gca,'YDir','normal');
        ylabel('sin(Latitude)', 'FontSize', 12);
        xlabel('Heliographic Longitude[deg]', 'FontSize', 12);
        ax = gca;
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels_sine)
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca,'TickDir','both');
        hold on;
        line(xlim, [ndxPositive,ndxPositive],'color','w','LineWidth',1)
        line(xlim, [ndxNegative,ndxNegative],'color','w','LineWidth',1)
 
        if  ~all(isnan(im_euv(:)))
            for n_region = 1:ARnumberStereo % Loop through each connected region
                rectangle('Position', regionProps_euv(n_region).BoundingBox, 'EdgeColor', 'g', 'LineWidth', 2); % Draw bounding box
                centroid = regionProps_euv(n_region).Centroid; % Label with centroid
                text(centroid(1), centroid(2), num2str(n_region), 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
            end
        end

        hold on
        title(['BW Synchronic Map - ',file_date ,'t', tt_str], 'FontSize', 12);
        scatter (y,x1,'w','.'); % At time t
        
        f.Position(4) = 1000; % Expand the figure vertically
        f.Position(3) = 1200;  % Expand the figure horizontally
        
        %% SAVE .png
        saveas(f,[pwd '/TH_Avg_fqm_10/png/AR_avg_fqm_',file_date,'t',tt_str,'.png']);
        close;

        %% SAVE .fits
        fitswrite(bw4_double, [pwd '/TH_Avg_fqm_10/fits/AR_avg_fqm_',file_date,'t',tt_str,'.fits']);

        %% Extracting/ SAVING the AR information
        Data = [C0 C1 C2 C3_X round(C3_Long) C4_Y round(C4_sinLat,2) round(C4_Lat) round(C5,4) round(C6,4) round(C7,4) round(C8,4) C9];

        columnNames = {'Date', '#AR', 'Area(% n.of pix)', 'AR_X', 'AR.Long.', 'AR_Y', 'AR_SinLat.', 'AR_Lat', 'SumPhaseShift', 'MeanPhaseShift', 'MinimumPhaseShift', 'stdPhaseShift', 'confusion matrix/performance evaluation [TP-FP]'};
        tableData = table(Data(:,1), Data(:,2), Data(:,3), Data(:,4), Data(:,5), Data(:,6), Data(:,7), Data(:,8), Data(:,9), Data(:,10), Data(:,11), Data(:,12), Data(:,13), 'VariableNames', columnNames);
        
        writetable(tableData, [pwd '/TH_Avg_fqm_10/csv/AR_avg_fqm_',file_date,'t',tt_str,'.csv']);
        

    end
end
close(wait)



