clear all ; 
close all ;

% plot error vs time step
% input 1 or 2 log files

icase = 0 ;

switch(icase)

    case 0
        FILENAME1 = 'hpcscan.timestepError.Propa.log' ;
        NAME1 = 'FDO8-FP32' ;
        COLOR1 = '-k' ;
        FILENAME2 = 'none' ;
        NAME2 = 'none' ;
        COLOR2 = '-b' ;

    case 1
        FILENAME1 = 'fp64.propa.paramAnalysis2.ouessant.2026-08-13.log' ;
        NAME1 = 'FP64' ;
        COLOR1 = '-k' ;
        FILENAME2 = 'none' ;
        NAME2 = 'none' ;
        COLOR2 = '-b' ;

    case 2
        FILENAME1 = 'fp64.propa.paramAnalysis2.ouessant.2026-08-13.log' ;
        NAME1 = 'FP64' ;
        COLOR1 = '-k' ;
        FILENAME2 = 'fp32.propa.paramAnalysis2.ouessant.2026-08-13.log' ;
        NAME2 = 'FP32' ;
        COLOR2 = '-b' ;

    case 3
        FILENAME1 = 'fp32.propa.paramAnalysis2.ouessant.2026-08-13.log' ;
        NAME1 = 'FP32' ;
        COLOR1 = '-b' ;
        FILENAME2 = 'fp16.fp16.propa.paramAnalysis2.ouessant.2026-08-13.log' ;
        NAME2 = 'FP16-32' ;
        COLOR2 = '-r' ;

    case 4
        FILENAME1 = 'fp32.propa.paramAnalysis2.ouessant.2026-08-13.log' ;
        NAME1 = 'FP32' ;
        COLOR1 = '-b' ;
        FILENAME2 = 'fp16.fp16.propa.paramAnalysis2.ouessant.2026-08-13.log' ;
        NAME2 = 'FP16-16' ;
        COLOR2 = '-r' ;

    case 5
        FILENAME1 = 'fp32.propa.paramAnalysis2.ouessant.2026-08-13.log' ;
        NAME1 = 'FP32' ;
        COLOR1 = '-b' ;
        FILENAME2 = 'fp32.fp16.propa.paramAnalysis2.ouessant.2026-08-14.log' ;
        NAME2 = 'FP32-16' ;
        COLOR2 = '--r' ;

    case 6
        FILENAME1 = 'fp32.propa.paramAnalysis2.ouessant.2026-08-20.log' ;
        NAME1 = 'FP32' ;
        COLOR1 = '-k' ;
        FILENAME2 = 'none' ;
        NAME2 = 'none' ;
        COLOR2 = '-r' ;

    case 7
        FILENAME1 = 'fp32.propa.paramAnalysis2.ouessant.2026-08-20.log' ;
        NAME1 = 'FP32' ;
        COLOR1 = '-k' ;
        FILENAME2 = 'fp16.fp32.propa.paramAnalysis2.ouessant.2026-08-20.log' ;
        NAME2 = 'FP16-32' ;
        COLOR2 = '-r' ;

    case 8
        FILENAME1 = 'fp32.propa.paramAnalysis2.ouessant.2026-08-20.log' ;
        NAME1 = 'FP32' ;
        COLOR1 = '-k' ;
        FILENAME2 = 'custom_mp.propa.paramAnalysis2.ouessant.2026-08-20.log' ;
        NAME2 = 'BlockMP' ;
        COLOR2 = '-b' ;

end

MIN_ERR = 1.e-8 ;
MAX_ERR = 1.e+2 ;

TARGET_ERR = 1.e-2 ;

% end user input parameters

val1 = load(FILENAME1) ;
val1Timestep = val1(:,1) ;
val1Error    = val1(:,2) ;

if (~strcmp(NAME2, 'none'))
    val2 = load(FILENAME2) ;
    val2Timestep = val2(:,1) ;
    val2Error    = val2(:,2) ;
end

valErrX = val1Timestep ;
valErrY = val1Timestep ;
valErrY(:) = TARGET_ERR ;

FONT_SIZE = 14 ;

figure
hold on
TITLE='L1 error vs time step' ;
title(TITLE, 'FontSize', FONT_SIZE+1)
plot(valErrX, valErrY, '--k', LineWidth=1.5)
plot(val1Timestep, val1Error, COLOR1, LineWidth=1.5)
if (~strcmp(NAME2, 'none'))
    plot(val2Timestep, val2Error, COLOR2, LineWidth=1.5)
end
xlabel('Time step', 'FontSize', FONT_SIZE)
ylabel('L1 error', 'FontSize', FONT_SIZE)
legend('Err. target', NAME1,NAME2,'Location','northwest')
grid on
ylim([MIN_ERR MAX_ERR])

ax=gca ;
%ax.XScale='log'
ax.YScale='log' ;

minval1Error = min(val1Error(:))
maxval1Error = max(val1Error(:))
if (~strcmp(NAME2, 'none'))
    minval2Error = min(val2Error(:))
    maxval2Error = max(val2Error(:))
end

% save figure as jpeg picture
FIG_NAME  = sprintf('Fig/error-%s-%s.jpg', NAME1, NAME2) ;
print('-djpeg', FIG_NAME)
