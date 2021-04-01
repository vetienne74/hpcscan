
clear all; close all ;

%--------------------------------------------------------------------------
% start input parameters
DIR  = '.' ;

% log file name with .log extension
%FILE = 'hpcscanFD_D2Shaheen' ;
FILE = 'hpcscan.perf.FD_D2' ;
FILE = 'hpcscanFD_D2cacheBlkAnalysisNECPlconec01' ;

% indicate below the range of parameters in the log file
fdOrderInLog = [2,4,8,12,16] ;
cb1InLog     = 9999 ;
cb2InLog     = 1:32 ;
cb3InLog     = 1:32 ;

PLOT_GPOINT = 0 ;
PLOT_GFLOP  = 0 ;
PLOT_GBYTE  = 1 ;
% end input parameters
%--------------------------------------------------------------------------

% allocate tables to store perf.
nfdOrderInLog = numel(fdOrderInLog)
ncb1InLog     = numel(cb1InLog)
ncb2InLog     = numel(cb2InLog)
ncb3InLog     = numel(cb3InLog)
perfMapD2Axis1Gpoint = zeros(nfdOrderInLog, ncb1InLog, ncb2InLog, ncb3InLog) ;
perfMapD2Axis1Gflop  = zeros(nfdOrderInLog, ncb1InLog, ncb2InLog, ncb3InLog) ;
perfMapD2Axis1Gbyte  = zeros(nfdOrderInLog, ncb1InLog, ncb2InLog, ncb3InLog) ;
perfMapD2Axis2Gpoint = zeros(nfdOrderInLog, ncb1InLog, ncb2InLog, ncb3InLog) ;
perfMapD2Axis2Gflop  = zeros(nfdOrderInLog, ncb1InLog, ncb2InLog, ncb3InLog) ;
perfMapD2Axis2Gbyte  = zeros(nfdOrderInLog, ncb1InLog, ncb2InLog, ncb3InLog) ;
perfMapD2Axis3Gpoint = zeros(nfdOrderInLog, ncb1InLog, ncb2InLog, ncb3InLog) ;
perfMapD2Axis3Gflop  = zeros(nfdOrderInLog, ncb1InLog, ncb2InLog, ncb3InLog) ;
perfMapD2Axis3Gbyte  = zeros(nfdOrderInLog, ncb1InLog, ncb2InLog, ncb3InLog) ;
perfMapD2LaplaGpoint = zeros(nfdOrderInLog, ncb1InLog, ncb2InLog, ncb3InLog) ;
perfMapD2LaplaGflop  = zeros(nfdOrderInLog, ncb1InLog, ncb2InLog, ncb3InLog) ;
perfMapD2LaplaGbyte  = zeros(nfdOrderInLog, ncb1InLog, ncb2InLog, ncb3InLog) ;

% allocate table to store max val and optimal cache size
maxD2Axis1_gpoint = zeros(nfdOrderInLog,1) ;
maxD2Axis2_gpoint = zeros(nfdOrderInLog,1) ;
maxD2Axis3_gpoint = zeros(nfdOrderInLog,1) ;
maxD2Lapla_gpoint = zeros(nfdOrderInLog,1) ;
maxD2Axis1_gflop = zeros(nfdOrderInLog,1) ;
maxD2Axis2_gflop = zeros(nfdOrderInLog,1) ;
maxD2Axis3_gflop = zeros(nfdOrderInLog,1) ;
maxD2Lapla_gflop = zeros(nfdOrderInLog,1) ;
maxD2Axis1_gbyte = zeros(nfdOrderInLog,1) ;
maxD2Axis2_gbyte = zeros(nfdOrderInLog,1) ;
maxD2Axis3_gbyte = zeros(nfdOrderInLog,1) ;
maxD2Lapla_gbyte = zeros(nfdOrderInLog,1) ;
maxD2Axis1_cb = zeros(nfdOrderInLog,3) ;
maxD2Axis2_cb = zeros(nfdOrderInLog,3) ;
maxD2Axis3_cb = zeros(nfdOrderInLog,3) ;
maxD2Lapla_cb = zeros(nfdOrderInLog,3) ;

% index of parameters in log file
iFDOrder          = 9 ;
iD2Axis1Gflop     = 10 ;
iD2Axis1GpointFD  = 11 ;
iD2Axis1GpointEff = 12 ;
iD2Axis1GB        = 13 ;

iD2Axis2Gflop     = 14 ;
iD2Axis2GpointFD  = 15 ;
iD2Axis2GpointEff = 16 ;
iD2Axis2GB        = 17 ;

iD2Axis3Gflop     = 18 ;
iD2Axis3GpointFD  = 19 ;
iD2Axis3GpointEff = 20 ;
iD2Axis3GB        = 21 ;

iD2LaplaGflop     = 22 ;
iD2LaplaGpointFD  = 23 ;
iD2LaplaGpointEff = 24 ;
iD2LaplaGB        = 25 ;

iCB1              = 26 ;
iCB2              = 27 ;
iCB3              = 28 ;

% load table
pathFile = sprintf('%s/%s.log', DIR, FILE) ;
val = importdata(pathFile) ;

% build perf maps
sizeLog = size(val.data)
if (sizeLog(1) ~= numel(perfMapD2Axis1Gpoint))
    disp('*** ERROR, LOG FILE INCONSISTENT WITH PARAM DESCRIPTION ! ***')
    return
else
    disp('START BUILD PERF MAPS...')
end

for iconfig = 1:sizeLog(1)
    % get index in perf map
    idxFdOrder = find(fdOrderInLog(:) == val.data(iconfig,iFDOrder)) ;
    idxCb1     = find(cb1InLog(:) == val.data(iconfig,iCB1)) ;
    idxCb2     = find(cb2InLog(:) == val.data(iconfig,iCB2)) ;
    idxCb3     = find(cb3InLog(:) == val.data(iconfig,iCB3)) ;
    
    % store GFlop
    perfMapD2Axis1Gflop(idxFdOrder, idxCb1, idxCb2, idxCb3) = val.data(iconfig,iD2Axis1Gflop) ;
    perfMapD2Axis2Gflop(idxFdOrder, idxCb1, idxCb2, idxCb3) = val.data(iconfig,iD2Axis2Gflop) ;
    perfMapD2Axis3Gflop(idxFdOrder, idxCb1, idxCb2, idxCb3) = val.data(iconfig,iD2Axis3Gflop) ;
    perfMapD2LaplaGflop(idxFdOrder, idxCb1, idxCb2, idxCb3) = val.data(iconfig,iD2LaplaGflop) ;
    
    % store GPoint
    perfMapD2Axis1Gpoint(idxFdOrder, idxCb1, idxCb2, idxCb3) = val.data(iconfig,iD2Axis1GpointFD) ;
    perfMapD2Axis2Gpoint(idxFdOrder, idxCb1, idxCb2, idxCb3) = val.data(iconfig,iD2Axis2GpointFD) ;
    perfMapD2Axis3Gpoint(idxFdOrder, idxCb1, idxCb2, idxCb3) = val.data(iconfig,iD2Axis3GpointFD) ;
    perfMapD2LaplaGpoint(idxFdOrder, idxCb1, idxCb2, idxCb3) = val.data(iconfig,iD2LaplaGpointFD) ;
    
    % store GB
    perfMapD2Axis1Gbyte(idxFdOrder, idxCb1, idxCb2, idxCb3) = val.data(iconfig,iD2Axis1GB) ;
    perfMapD2Axis2Gbyte(idxFdOrder, idxCb1, idxCb2, idxCb3) = val.data(iconfig,iD2Axis2GB) ;
    perfMapD2Axis3Gbyte(idxFdOrder, idxCb1, idxCb2, idxCb3) = val.data(iconfig,iD2Axis3GB) ;
    perfMapD2LaplaGbyte(idxFdOrder, idxCb1, idxCb2, idxCb3) = val.data(iconfig,iD2LaplaGB) ;
end

disp('START FIND OPTIM CACHE SIZE...')
for iOrder = 1:nfdOrderInLog
    for iCb1 = 1:ncb1InLog
        for iCb2 = 1:ncb2InLog
            for iCb3 = 1:ncb3InLog
                if perfMapD2Axis1Gpoint(iOrder, iCb1, iCb2, iCb3) > maxD2Axis1_gpoint(iOrder)
                    maxD2Axis1_gpoint(iOrder) = perfMapD2Axis1Gpoint(iOrder, iCb1, iCb2, iCb3) ;
                    maxD2Axis1_gflop(iOrder) = perfMapD2Axis1Gflop(iOrder, iCb1, iCb2, iCb3) ;
                    maxD2Axis1_gbyte(iOrder) = perfMapD2Axis1Gbyte(iOrder, iCb1, iCb2, iCb3) ;
                    maxD2Axis1_cb(iOrder, 1) = iCb1 ;
                    maxD2Axis1_cb(iOrder, 2) = iCb2 ;
                    maxD2Axis1_cb(iOrder, 3) = iCb3 ;
                end
                if perfMapD2Axis2Gpoint(iOrder, iCb1, iCb2, iCb3) > maxD2Axis2_gpoint(iOrder)
                    maxD2Axis2_gpoint(iOrder) = perfMapD2Axis2Gpoint(iOrder, iCb1, iCb2, iCb3) ;
                    maxD2Axis2_gflop(iOrder) = perfMapD2Axis2Gflop(iOrder, iCb1, iCb2, iCb3) ;
                    maxD2Axis2_gbyte(iOrder) = perfMapD2Axis2Gbyte(iOrder, iCb1, iCb2, iCb3) ;
                    maxD2Axis2_cb(iOrder, 1) = iCb1 ;
                    maxD2Axis2_cb(iOrder, 2) = iCb2 ;
                    maxD2Axis2_cb(iOrder, 3) = iCb3 ;
                end
                if perfMapD2Axis3Gpoint(iOrder, iCb1, iCb2, iCb3) > maxD2Axis3_gpoint(iOrder)
                    maxD2Axis3_gpoint(iOrder) = perfMapD2Axis3Gpoint(iOrder, iCb1, iCb2, iCb3) ;
                    maxD2Axis3_gflop(iOrder) = perfMapD2Axis3Gflop(iOrder, iCb1, iCb2, iCb3) ;
                    maxD2Axis3_gbyte(iOrder) = perfMapD2Axis3Gbyte(iOrder, iCb1, iCb2, iCb3) ;
                    maxD2Axis3_cb(iOrder, 1) = iCb1 ;
                    maxD2Axis3_cb(iOrder, 2) = iCb2 ;
                    maxD2Axis3_cb(iOrder, 3) = iCb3 ;
                end
                if perfMapD2LaplaGpoint(iOrder, iCb1, iCb2, iCb3) > maxD2Lapla_gpoint(iOrder)
                    maxD2Lapla_gpoint(iOrder) = perfMapD2LaplaGpoint(iOrder, iCb1, iCb2, iCb3) ;
                    maxD2Lapla_gflop(iOrder) = perfMapD2LaplaGflop(iOrder, iCb1, iCb2, iCb3) ;
                    maxD2Lapla_gbyte(iOrder) = perfMapD2LaplaGbyte(iOrder, iCb1, iCb2, iCb3) ;
                    maxD2Lapla_cb(iOrder, 1) = iCb1 ;
                    maxD2Lapla_cb(iOrder, 2) = iCb2 ;
                    maxD2Lapla_cb(iOrder, 3) = iCb3 ;
                end
            end
        end
    end
end


for iOrder = 1:nfdOrderInLog
    %for iOrder = 5:5
    
    if (PLOT_GFLOP)
        %Gflop
        figure('Position', [100 100 1000 800'])
        
        subplot(2,2,1); hold on
        title(sprintf('FD O%d D2Axis1 Gflop \n Max %.1f / cb2=%d cb3=%d', ...
            fdOrderInLog(iOrder), maxD2Axis1_gflop(iOrder), maxD2Axis1_cb(iOrder, 2), maxD2Axis1_cb(iOrder, 3)))
        xlabel('block size axis3'); ylabel('block size axis2');
        perfMap = reshape(perfMapD2Axis1Gflop(iOrder, 1, :, :), [ncb2InLog ncb3InLog]) ;
        surf(perfMap); colorbar
        view(0,90); axis tight
        
        subplot(2,2,2); hold on
        title(sprintf('FD O%d D2Axis2 Gflop \n Max %.1f / cb2=%d cb3=%d', ...
            fdOrderInLog(iOrder), maxD2Axis2_gflop(iOrder), maxD2Axis2_cb(iOrder, 2), maxD2Axis2_cb(iOrder, 3)))
        xlabel('block size axis3'); ylabel('block size axis2');
        perfMap = reshape(perfMapD2Axis2Gflop(iOrder, 1, :, :), [ncb2InLog ncb3InLog]) ;
        surf(perfMap); colorbar
        view(0,90); axis tight
        
        subplot(2,2,3); hold on
        title(sprintf('FD O%d D2Axis3 Gflop \n Max %.1f / cb2=%d cb3=%d', ...
            fdOrderInLog(iOrder), maxD2Axis3_gflop(iOrder), maxD2Axis3_cb(iOrder, 2), maxD2Axis3_cb(iOrder, 3)))
        xlabel('block size axis3'); ylabel('block size axis2');
        perfMap = reshape(perfMapD2Axis3Gflop(iOrder, 1, :, :), [ncb2InLog ncb3InLog]) ;
        surf(perfMap); colorbar
        view(0,90); axis tight
        
        subplot(2,2,4); hold on
        title(sprintf('FD O%d Lapla Gflop \n Max %.1f / cb2=%d cb3=%d', ...
            fdOrderInLog(iOrder), maxD2Lapla_gflop(iOrder), maxD2Lapla_cb(iOrder, 2), maxD2Lapla_cb(iOrder, 3)))
        xlabel('block size axis3'); ylabel('block size axis2');
        perfMap = reshape(perfMapD2LaplaGflop(iOrder, 1, :, :), [ncb2InLog ncb3InLog]) ;
        surf(perfMap); colorbar
        view(0,90); axis tight
        
        colormap(jet)
        
        figFile = sprintf('%s_GflopFD_O%d.jpg', FILE, fdOrderInLog(iOrder)) ;
        print('-djpeg', figFile)
    end
    
    %Gpoint
    if (PLOT_GPOINT)
        figure('Position', [100 100 1000 800'])
        
        subplot(2,2,1); hold on
        title(sprintf('FD O%d D2Axis1 Gpoint \n Max %.1f / cb2=%d cb3=%d', ...
            fdOrderInLog(iOrder), maxD2Axis1_gpoint(iOrder), maxD2Axis1_cb(iOrder, 2), maxD2Axis1_cb(iOrder, 3)))
        xlabel('block size axis3'); ylabel('block size axis2');
        perfMap = reshape(perfMapD2Axis1Gpoint(iOrder, 1, :, :), [ncb2InLog ncb3InLog]) ;
        surf(perfMap); colorbar
        view(0,90); axis tight
        
        subplot(2,2,2); hold on
        title(sprintf('FD O%d D2Axis2 Gpoint \n Max %.1f / cb2=%d cb3=%d', ...
            fdOrderInLog(iOrder), maxD2Axis2_gpoint(iOrder), maxD2Axis2_cb(iOrder, 2), maxD2Axis2_cb(iOrder, 3)))
        xlabel('block size axis3'); ylabel('block size axis2');
        perfMap = reshape(perfMapD2Axis2Gpoint(iOrder, 1, :, :), [ncb2InLog ncb3InLog]) ;
        surf(perfMap); colorbar
        view(0,90); axis tight
        
        subplot(2,2,3); hold on
        title(sprintf('FD O%d D2Axis3 Gpoint \n Max %.1f / cb2=%d cb3=%d', ...
            fdOrderInLog(iOrder), maxD2Axis3_gpoint(iOrder), maxD2Axis3_cb(iOrder, 2), maxD2Axis3_cb(iOrder, 3)))
        xlabel('block size axis3'); ylabel('block size axis2');
        perfMap = reshape(perfMapD2Axis3Gpoint(iOrder, 1, :, :), [ncb2InLog ncb3InLog]) ;
        surf(perfMap); colorbar
        view(0,90); axis tight
        
        subplot(2,2,4); hold on
        title(sprintf('FD O%d Lapla Gpoint \n Max %.1f / cb2=%d cb3=%d', ...
            fdOrderInLog(iOrder), maxD2Lapla_gpoint(iOrder), maxD2Lapla_cb(iOrder, 2), maxD2Lapla_cb(iOrder, 3)))
        xlabel('block size axis3'); ylabel('block size axis2');
        perfMap = reshape(perfMapD2LaplaGpoint(iOrder, 1, :, :), [ncb2InLog ncb3InLog]) ;
        surf(perfMap); colorbar
        view(0,90); axis tight
        
        colormap(jet)
        
        figFile = sprintf('%s_GpointFD_O%d.jpg', FILE, fdOrderInLog(iOrder)) ;
        print('-djpeg', figFile)
    end
    
    %Gbyte
    if (PLOT_GBYTE)
        figure('Position', [100 100 1000 800'])
        
        subplot(2,2,1); hold on
        title(sprintf('FD O%d D2Axis1 Gbyte \n Max %.1f / cb2=%d cb3=%d', ...
            fdOrderInLog(iOrder), maxD2Axis1_gbyte(iOrder), maxD2Axis1_cb(iOrder, 2), maxD2Axis1_cb(iOrder, 3)))
        xlabel('block size axis3'); ylabel('block size axis2');
        perfMap = reshape(perfMapD2Axis1Gbyte(iOrder, 1, :, :), [ncb2InLog ncb3InLog]) ;
        surf(perfMap); colorbar
        view(0,90); axis tight
        
        subplot(2,2,2); hold on
        title(sprintf('FD O%d D2Axis2 Gbyte \n Max %.1f / cb2=%d cb3=%d', ...
            fdOrderInLog(iOrder), maxD2Axis2_gbyte(iOrder), maxD2Axis2_cb(iOrder, 2), maxD2Axis2_cb(iOrder, 3)))
        xlabel('block size axis3'); ylabel('block size axis2');
        perfMap = reshape(perfMapD2Axis2Gbyte(iOrder, 1, :, :), [ncb2InLog ncb3InLog]) ;
        surf(perfMap); colorbar
        view(0,90); axis tight
        
        subplot(2,2,3); hold on
        title(sprintf('FD O%d D2Axis3 Gbyte \n Max %.1f / cb2=%d cb3=%d', ...
            fdOrderInLog(iOrder), maxD2Axis3_gbyte(iOrder), maxD2Axis3_cb(iOrder, 2), maxD2Axis3_cb(iOrder, 3)))
        xlabel('block size axis3'); ylabel('block size axis2');
        perfMap = reshape(perfMapD2Axis3Gbyte(iOrder, 1, :, :), [ncb2InLog ncb3InLog]) ;
        surf(perfMap); colorbar
        view(0,90); axis tight
        
        subplot(2,2,4); hold on
        title(sprintf('FD O%d Lapla Gbyte \n Max %.1f / cb2=%d cb3=%d', ...
            fdOrderInLog(iOrder), maxD2Lapla_gbyte(iOrder), maxD2Lapla_cb(iOrder, 2), maxD2Lapla_cb(iOrder, 3)))
        xlabel('block size axis3'); ylabel('block size axis2');
        perfMap = reshape(perfMapD2LaplaGbyte(iOrder, 1, :, :), [ncb2InLog ncb3InLog]) ;
        surf(perfMap); colorbar
        view(0,90); axis tight
        
        colormap(jet)
        
        figFile = sprintf('%s_GbyteFD_O%d.jpg', FILE, fdOrderInLog(iOrder)) ;
        print('-djpeg', figFile)
    end
end

%end
