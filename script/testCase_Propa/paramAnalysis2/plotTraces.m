clear all ; 
close all ;

% display seismic traces
% compare traces with same scale between 2 datasets
% compute NRMS error (normalised root mean square)
% display NRMS versus time
% Vincent Etienne - AlgoDoers - July 2025

%===================
% input parameters
%===================

FIG_NAME  = 'CompareTraces' ;

icase = 3 ;

switch(icase)

    case 1
        FILE_NAME1 = './PropaEigenModeRef.proc0.grid.bin' ;
        TITLE1 = 'FP32-REF' ;
        TYPE1  = 'ieee-le' ;
        FLOAT_SIZE1 = 'float32' ;

        FILE_NAME2 = './PropaEigenModePrn.proc0.grid.bin' ;
        TITLE2 = 'FP32-FD' ;
        TYPE2  = 'ieee-le' ;
        FLOAT_SIZE2 = 'float32' ;

    case 2
        FILE_NAME1 = './PropaEigenModeRef.proc0.grid.bin.fp32' ;
        TITLE1 = 'FP32-REF' ;
        TYPE1  = 'ieee-le' ;
        FLOAT_SIZE1 = 'float32' ;

        FILE_NAME2 = './PropaEigenModePrn.proc0.grid.bin.fp32' ;
        TITLE2 = 'FP32-FD' ; 
        TYPE2  = 'ieee-le' ;
        FLOAT_SIZE2 = 'float32' ;

    case 3
        FILE_NAME1 = './PropaEigenModeRef.proc0.grid.bin.fp64' ;
        TITLE1 = 'FP64-REF' ;
        TYPE1  = 'ieee-le' ;
        FLOAT_SIZE1 = 'float64' ;

        FILE_NAME2 = './PropaEigenModePrn.proc0.grid.bin.fp64' ;
        TITLE2 = 'FP64-FD' ;
        TYPE2  = 'ieee-le' ;
        FLOAT_SIZE2 = 'float64' ;

end

NT = 10001 ;
NREC = 508 ;
DT = 0.1 ;
DREC = 2.00401 ;
MIN_REC_OFFSET = 0 ;
REC_MIN = 5 ;
REC_MAX = 504 ;
NREC2 = REC_MAX-REC_MIN+1;

PERC = 0.1 ;
FONT_SIZE  = 18 ;
MAX_NRMS   = 0.1 ;

% end of user input parameters

mkdir Fig

file1 = fopen(FILE_NAME1, 'r', TYPE1);
file2 = fopen(FILE_NAME2, 'r', TYPE2);
XIT = (1:NT) ;

%----------------------
% loop over time steps
%----------------------

nrms_time = zeros(1, NT) ;
nrms_time_cumul = zeros(1, NT) ;
cumul1 = 0 ;
cumul2 = 0 ;
trace1    = zeros(NT, NREC2) ;
trace2    = zeros(NT, NREC2) ;

for it = 1:NT

    % read data
    val1 = fread(file1, [1, NREC], FLOAT_SIZE1);
    val2 = fread(file2, [1, NREC], FLOAT_SIZE2) ;

    % copy subset of trace
    trace1(it,1:NREC2) = val1(REC_MIN:REC_MAX) ;
    trace2(it,1:NREC2) = val2(REC_MIN:REC_MAX) ;   

    % difference
    diff = val1 - val2 ;

    % compute nrms  
    sum1 = 0.0;
    sum2 = 0.0;
    for irec=1:NREC
        sum1=sum1+diff(1,irec)^2 ;
        sum2=sum2+val1(1,irec)^2 ;
    end
    if (sum2 == 0)
        nrms = 0 ;
    else
        nrms=sqrt(sum1/sum2) ;
    end
    nrms_time(it) = nrms ;
    cumul1 = cumul1 + sum1 ;
    cumul2 = cumul2 + sum2 ;
    if (cumul2 == 0)
        nrms2 = 0 ;
    else
        nrms2 = sqrt(cumul1/cumul2) ;
    end
    nrms_time_cumul(it) = nrms2 ;

end

fclose(file1);
fclose(file2);

%======================================================================
% plot error vs time (figure 2)
%======================================================================
figure
hold on
TITLE=sprintf('%s vs %s', TITLE1, TITLE2) ;
title(TITLE, 'FontSize', FONT_SIZE+1)
%plot(XIT, nrms_time, '-k', LineWidth=1.0)
plot(XIT, nrms_time_cumul, '-k', LineWidth=1.0)
plot(XIT, nrms_time_cumul, '-b', LineWidth=1.5)
xlabel('Time step', 'FontSize', FONT_SIZE)
ylabel('NRMS error', 'FontSize', FONT_SIZE)
legend('NRMS instantaneous','NRMS cumulated','Location','northwest')
grid on
%ylim([0 MAX_NRMS])

% save figure as jpeg picture
FIG_NAME2  = sprintf('Fig/NRMS-%s-%s.jpg', TITLE1, TITLE2) ;
print('-djpeg', FIG_NAME2)

nrms_time(NT)


%======================================================================
% plot traces (figure 1)
%======================================================================
figure('Position',[100 500 1600 600])

XIT  = (1:NT) ;
XREC = REC_MIN:REC_MAX ;

% 1st file
subplot(1,3,1)
hold on
xlabel('Receiver index', 'FontSize', FONT_SIZE)
ylabel('Timestep', 'FontSize', FONT_SIZE)
TITLE=sprintf('%s', TITLE1) ;
title(TITLE, 'FontSize', FONT_SIZE+1)

colormap(gray) ; box on ;
imagesc(XREC, XIT, trace1) ;
max_axis = max(max(abs(trace1))) * PERC ;
caxis([-max_axis max_axis])
axis tight; axis ij

% 2nd file
subplot(1,3,2)
hold on
xlabel('Receiver index', 'FontSize', FONT_SIZE)
ylabel('Timestep', 'FontSize', FONT_SIZE)
TITLE=sprintf('%s', TITLE2) ;
title(TITLE, 'FontSize', FONT_SIZE+1)

colormap(gray) ; box on ;
imagesc(XREC, XIT, trace2) ;
caxis([-max_axis max_axis])
axis tight; axis ij

% diff
tracediff = trace1 - trace2 ;
subplot(1,3,3)
hold on
xlabel('Receiver index', 'FontSize', FONT_SIZE)
ylabel('Timestep', 'FontSize', FONT_SIZE)
TITLE=sprintf('%s', 'Difference') ;
title(TITLE, 'FontSize', FONT_SIZE+1)

colormap(gray) ; box on ;
imagesc(XREC, XIT, tracediff) ;
caxis([-max_axis max_axis])
axis tight; axis ij

% save figure as jpeg picture
FIG_NAME2  = sprintf('Fig/TRACE-%s-%s.jpg', TITLE1, TITLE2) ;
print('-djpeg', FIG_NAME2)