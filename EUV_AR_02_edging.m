%----------------------------------------------------------------------
% Title: Normalization and Thresholding of STEREO Data
% Author: Dr.AMR HAMADA [National Solar Observatory]
% Date: November 2024
% 
%----------------------------------------------------------------------

%% Initialization

clear all
count  = 0;
count2 = 0;

load('map_304.mat');

%% Start and End dares
y1 = 2010; m1 = 05; d1 = 13;
y2 = 2016; m2 = 05; d2 = 14;

%% Arranging the date
D = datenum(y1,m1,d1):datenum(y2,m2,d2);
D = 10000*year(D) + 100*month(D) + day(D);
D = D';
num_days = length(D);

%% Make dir
mkdir([pwd '/TH_Stand_STEREO_JPL/fts/']);
mkdir([pwd '/TH_Stand_STEREO_JPL/png/']);
mkdir([pwd '/TH_Stand_STEREO_JPL/ARs/']);

%% STEREO-Date
stereo_d = datenum(y1,m1,d1):datenum(y2,m2,d2);
stereo_d = stereo_d';

%% Folder containing STEREO .fts files
folderPath = ([pwd '/REPLACE_WITH_FOLDER_PATH/']) ;
ftsFiles = dir(fullfile(folderPath, '*.fts'));
ftsFileNames = {ftsFiles.name};

%% Segment Size
w_step = round(1801/180); % 1801 is the number of rows of STEREO map
W = (15:5:60)*w_step;
W_n = numel(W); % w_n is the number of windowes.

%% Loop through each STEREO .fts file
count = 0;
count_dd_tt = 0;

for i = 1:num_days
    count = count + 1;
    for tt = [00 12]
        count_dd_tt = count_dd_tt+1;
        if tt == 0
            tt_str = '00';
        else
            tt_str = '12';
        end
        
        st_map_name = ['stand_', num2str(D(i)), '_',tt_str, '*.fts'];
        
        %% FIND the STEREO map match with DATA_date_time
        matchingFile = dir(fullfile(folderPath, st_map_name));
        
        ftsFileName = matchingFile.name
        ftsFilePath = fullfile(folderPath, ftsFileName);
        
        im = fitsread(ftsFilePath);
        if ~all(isnan(im(:)))
            ims = imbilatfilt(im, 50, std(im(:),"omitmissing"));
            %         ims = imgaussfilt (im,2);     % 2D Gaussian filter
            
            %% LOAD the FS-Boundaries (x1, y1)
            FS_name = ftsFileName(9:14);
            FS_name = [FS_name, 't', tt_str,'00' ];
            try FS_boundaries = load([pwd filesep '/ES_FS_edges/', FS_name,'.mat']);
                X1 = FS_boundaries.x1;
                Y1 = FS_boundaries.y1;
            catch
                X1 = NaN;
                Y1 = NaN;
            end
            
            %% Thresholding
            count2 = 0;
            seg_win_ind = nan(W_n,5);  % TH_w, lat1, lat2, long1, long2
            
            for I = 1 : W_n  
                w = W(I);
                count2 = count2 + 1;
                
                im2=ims;
                M=ones(size(im2));
                M(~isfinite(im2))=0;
                im2(M==0)=NaN;
                
                A=ones(w,w);
                tmp1=filter2(A,im2)./filter2(A,M);
                tmp2=filter2(A,im2.^2)./filter2(A,M);
                tmp3=sqrt(tmp2-tmp1.^2);
                
                int=reshape(tmp1(ceil(w/2):end-floor(w/2),ceil(w/2):end-floor(w/2)),[],1);
                int(~isfinite(int))=NaN;
                int = real(int);
                
                cont=reshape(tmp3(ceil(w/2):end-floor(w/2),ceil(w/2):end-floor(w/2)),[],1);
                cont(~isfinite(cont))=NaN;
                cont = real(cont);
                
                Mu_C = mean(cont(:),"omitmissing");
                Sg_C = std(cont(:),"omitmissing");
                
                Mu_I = mean(int(:),"omitmissing");
                Sg_I = std(int(:),"omitmissing");
                
                C=tmp3(ceil(w/2):end-floor(w/2),ceil(w/2):end-floor(w/2));
                C(~isfinite(C))=NaN;  
                C = real(C);
                
                I_sgm=tmp1(ceil(w/2):end-floor(w/2),ceil(w/2):end-floor(w/2));
                R=(C-Mu_C)./Sg_C+2*((I_sgm-Mu_I)./Sg_I);
                
                [NN,MM]=meshgrid(1:length(R(1,:)),1:length(R(:,1)));
                
                if ~isnan(max(reshape(R,[],1),[],"omitnan"))
                    try ndx(1) = MM(R==max(reshape(R,[],1)));                         
                        ndx(2) = NN(R==max(reshape(R,[],1)));
                    catch
                    end
                else
                    ndx(1)=NaN;
                    ndx(2)=NaN;
                end
                
                if ~isnan (ndx)                    
                    lat1 = ndx(1);      
                    lat2 = lat1 +w-1;  
                    long1 = ndx(2);     
                    long2 = long1+w-1; 
                    
                    Rbest(count,count2)=R(ndx(1),ndx(2));
                end
                
                seg_win_ind (count2,2:5) = [lat1 lat2 long1 long2]; % [TH lat1 lat2 long1 long2]
                
                %%---------------- Gaussian Mixture Models--------------
                imseg = im2(lat1: lat2, long1: long2);
                imseg = (imseg);
                
                AIC = zeros(1,3);
                GMModels = cell(1,3);
                options = statset('MaxIter',500);
                for k = 2          % k is the number of component in Gaussian mixture
                    err=1;
                    while err==1
                        try
                            GMModels{k} = fitgmdist(reshape(imseg,numel(imseg),1),k,'Options',options,'Start','randSample','Replicates',10);
                            err=0;
                        catch
                            err=1;
                        end
                    end
                    AIC(k)= GMModels{k}.AIC;
                end
                [minAIC,numComponents] = min(AIC);
                numComponents = 2;
                BestModel = GMModels{numComponents};  
                amplitude=BestModel.ComponentProportion;
                
                [~,im_x]=hist((reshape(imseg,numel(imseg),1)),512);
                theta=numel(imseg)*median(diff(im_x));
                all_Sig=sqrt(squeeze(BestModel.Sigma(1,1,:)));
                all_Mu =squeeze(BestModel.mu);
                mu = sort(all_Mu);
                
                [~,imseg_x]=hist(reshape(imseg,numel(imseg),1),512);
                for ii=1:length(BestModel.mu)
                    eval(['mu' num2str(ii) '= mu(ii);']);
                    eval(['sig' num2str(ii) '= all_Sig(all_Mu == mu' num2str(ii) ');']);
                    eval(['G' num2str(ii) '= amplitude(all_Mu == mu' num2str(ii) ')*normpdf(imseg_x,mu' num2str(ii) ',sig' num2str(ii) ')*theta;']);
                end
                
                [xout,yout] = intersections(imseg_x, G1, imseg_x, G2, false);
                
                try
                    TH_w(count,count2) = xout(xout > mu1 & xout<mu2);
                catch
                    TH_w(count,count2) = mu1+3*sig1;
                end
                                
                seg_win_ind (count2,1) = TH_w(count,count2);
            end
            
            %% Select the Best segment/TH
            TH (count) = sum(TH_w(count,:).*(Rbest(count,:)==repmat(max(Rbest(count,:),[],2),1,W_n)),2)./sum((Rbest(count,:)==repmat(max(Rbest(count,:),[],2),1,W_n)),2);
            
            lat1 = seg_win_ind(seg_win_ind(:,1)==TH (count),2);
            lat2 = seg_win_ind(seg_win_ind(:,1)==TH (count),3);
            long1 = seg_win_ind(seg_win_ind(:,1)==TH (count),4);
            long2 = seg_win_ind(seg_win_ind(:,1)==TH (count),5);
            
            %% THRESHOLDING
            im_bw = im2 > TH(1,count);
            
            im_bw2 = im_bw;
            im_bw2 = imclose(im_bw2,strel('square',2)); % morphorogical operations
            im_bw2 = imopen(im_bw2,strel('square',2));  % morphorogical operations
                        
            %% EDGING
            im_bw2 = imfill(im_bw,'holes'); %  fills holes in im_bw
            im_bw2 = bwareaopen(im_bw,1000); % removes all connected components (objects) that have fewer than 10 pixels

            im_bw3 = imclose(im_bw2,strel('square',25)); % morphorogical operations
            
            %% MASKING THE CH-FILTER ON THE ORIGINAL EIT-IMAGE
            im_edge = edge(im_bw3,'canny');
            [x_edge,y_edge]= find(im_edge == 1);
            
            %% start ALLOCATING AR-groups.
            clear cc
            clear C1 C2 C3 C4 C5 C6;
            regionProps = regionprops(im_bw3, 'BoundingBox', 'Centroid'); % Label connected regions and get properties                     
            cc = bwconncomp(im_bw3);      % returns the connected components CC found in the binary image im_bw
            AR_n = cc.NumObjects;        % Number of CH detected.
            AR_ndx = cc.PixelIdxList;    % List by CH indiceies.                       
            if AR_n>0
                S = regionprops(cc,'Centroid');         % Calculate centroids of the objects in the array.
                centroids = cat(1, S.Centroid);         % Concatenate structure array containing centroids into a single matrix.                                
                
                C1 = 1:AR_n; % ARs number
                C1 = C1';

                for F = 1:AR_n
                    C2(F,1) = numel([AR_ndx{F}])/numel(im_bw3) *100; % AR-%Area
                    C5(F,1) = sum(im(AR_ndx{F}),"omitnan");                    % AR-sum(EUV-Intensity)
                    C6(F,1) = mean(im(AR_ndx{F}),"omitnan");                   % AR-Mean(EUV-Intensity)
                end

                C3 = centroids(:,1);    % AR-latitude
                C4 = centroids(:,2);    % AR-Longitude
                Data = [C1 C2 C3 C4 C5 C6];
            else
                C1=0;
                C2=0;
                C3=0;
                C4=0;
            end

            %% Saving fits
            fitswrite(double(im_bw3), [pwd '/TH_Stand_STEREO_JPL/fts/TH_',ftsFileName]);

            %% plotting.------------
            imseg = im2(lat1: lat2, long1: long2);

            f = figure('visible','off');
            xticklabels = 0 : 60 : 360;
            xticks = linspace(1, size(im, 2), numel(xticklabels));
            yticklabels = -1 : 0.5 : +1;
            yticks = linspace(1, size(im, 1), numel(yticklabels));
            
            ax_im = subplot (3,5,[1 2 3]);
            a = imagesc(im);
            shading flat; grid on;
            set(a,'alphaData',1-double(isnan(im)));
            set(gca,'YDir','normal');
            set(gca,'YDir','normal');
            colormap(ax_im,map_304)
            caxis ([-2 2]);
            h = colorbar; ylabel(h, '$Norm(I_{\mathrm{log}})$', 'Interpreter', 'latex')
            ylabel('Latitude[deg]', 'FontSize', 12);
            xlabel('Carring.Longitude[deg]', 'FontSize', 12);
            ax = gca;
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
            set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca,'TickDir','both');
            hold on
            line([long1, long1],[lat1, lat2],'LineWidth',1,'Color','k');
            line([long2, long2],[lat1, lat2],'LineWidth',1,'Color','k');
            line([long1, long2],[lat1, lat1],'LineWidth',1,'Color','k');
            line([long1, long2],[lat2, lat2],'LineWidth',1,'Color','k');
            scatter(Y1, X1, 'k','.')
            
            
            subplot (3,5,[4 5]);
            [nn,xx]=hist(reshape(im,numel(im),1),256);
            plot (xx,nn,'k');grid on;
            xlabel('$Norm(I_{\mathrm{log}})$', 'Interpreter', 'latex', 'FontSize', 12);
            ylabel('Number of pixles', 'FontSize', 12);
            axis tight
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca,'TickDir','out');
            
            
            % % SEGMENT.
            subplot(3,5,[6 7]);
            [mmm,nnn] = size (imseg);
            pcolor(imseg);
            shading flat; grid on
            colormap(map_304)
            caxis ([-1.5 2]);
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca,'TickDir','both');
            xlabel('Pix.', 'FontSize', 12);
            ylabel('Pix.', 'FontSize', 12);
            
            
            % % SEGMENT_HIST.
            subplot(3,5,[9 10]);
            AIC = zeros(1,3);
            GMModels = cell(1,3);
            options = statset('MaxIter',500);
            for k = 2          
                err=1;
                while err==1
                    try
                        GMModels{k} = fitgmdist(reshape(imseg,numel(imseg),1),k,'Options',options,'Start','randSample','Replicates',10);
                        err=0;
                    catch
                        err=1;
                    end
                end
                AIC(k)= GMModels{k}.AIC;
            end
            [minAIC,numComponents] = min(AIC);
            numComponents = 2;
            BestModel = GMModels{numComponents};  % BestModel
            amplitude=BestModel.ComponentProportion;
            
            [im_n,im_x]=hist((reshape(imseg,numel(imseg),1)),512);
            plot (im_x,im_n);
            grid on; hold on;
            theta=numel(imseg)*median(diff(im_x));
            all_Sig=sqrt(squeeze(BestModel.Sigma(1,1,:)));
            all_Mu =squeeze(BestModel.mu);
            mu = sort(all_Mu);
            [~,imseg_x]=hist(reshape(imseg,numel(imseg),1),512);
            for ii=1:length(BestModel.mu)
                eval(['mu' num2str(ii) '= mu(ii);']);
                eval(['sig' num2str(ii) '= all_Sig(all_Mu == mu' num2str(ii) ');']);
                eval(['G' num2str(ii) '= amplitude(all_Mu == mu' num2str(ii) ')*normpdf(imseg_x,mu' num2str(ii) ',sig' num2str(ii) ')*theta;']);
                plot(im_x,eval(['G' num2str(ii)]),'LineWidth',2);
            end
            plot(im_x,G1+G2,'LineWidth',2);
            lgnd = legend ('Norm(I)','G1','G2','G1+G2');
            set(lgnd,'color','none', 'Box', 'off');
            axis tight
            xlabel('$Norm(I_{\mathrm{log}})$', 'Interpreter', 'latex', 'FontSize', 12);
            ylabel('Number of pixles', 'FontSize', 12);
            line([TH(count) TH(count)], ylim,'LineStyle','-.','color','k');
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca,'TickDir','out');
            
            ax_bw = subplot (3,5,[12 13 14]);
            im_bw_show = double(im_bw);
            im_bw_show (isnan(im)) = NaN;
            a = pcolor(im_bw_show); shading flat; grid on;
            colormap(ax_bw,gray)
            set(a,'alphaData',1-double(isnan(im2)));
            set(gca,'YDir','normal');
            set(gca,'YDir','normal');
%             title(['AR-mask  ',file_date,' - t00']);
            ylabel('Latitude[deg]', 'FontSize', 12);
            xlabel('Carring.Longitude[deg]', 'FontSize', 12);
            ax = gca;
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
            set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
            set(gca,'XMinorTick','on','YMinorTick','on')
            set(gca,'TickDir','both');
                        
            screenSize = get(0, 'ScreenSize'); % Get screen size
            set(f, 'Position', screenSize); % Set figure position to match screen size            
            f.Position(4) = 1000; % Expand the figure vertically
            f.Position(3) = 1200;  % Expand the figure horizontally
            
            % Save the figure as a full-screen image
            saveas(f,[pwd '/TH_Stand_STEREO_JPL/png/TH_',ftsFileName(1:end-4),'_01.png']);
            
            %% PLOT-2  
            f2 = figure('visible','off');
            ax1 = subplot (511);
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
            im_bw3_show = double(im_bw3);
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
            ylabel('Latitude[deg]', 'FontSize', 12);
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
            a = pcolor(im);
            shading flat; grid on;
%             set(a,'alphaData',1-double(isnan(im)));
            set(gca,'YDir','normal');
            set(gca,'YDir','normal');
            colormap(ax5,map_304)
            caxis ([-2 2]);
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
            saveas(f2,[pwd '/TH_Stand_STEREO_JPL/png/TH_',ftsFileName(1:end-4),'_02.png']);
            
            %% Extracting the AR information
            columnNames = {'AR', 'Area(% n.of pix)', 'Cent.Lat.', 'Cent.Long', 'AR Tot_EUV_I', 'AR Mean_EUV_I'};
            tableData = table(Data(:,1), Data(:,2), Data(:,3), Data(:,4), Data(:,5), Data(:,6), 'VariableNames', columnNames);          
            writetable(tableData, [pwd '/TH_Stand_STEREO_JPL/ARs/ARs_',ftsFileName(1:end-4),'.csv']);


        end
    end
end

