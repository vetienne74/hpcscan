
clear all; close all ;

iNthread = 5 ;
iFillGB = 10 ;
iFillGpoint = 11 ;
iMaxErrGB = 12 ;
iMaxErrGpoint = 13 ;
iL1ErrGB = 14 ;
iL1ErrGpoint = 15 ;
iGetSumAbsGB = 16 ;
iGetSumAbsGpoint = 17 ;
iGetSumAbsDiffGB = 18 ;
iGetSumAbsDiffGpoint = 19 ;
iGetMaxGB = 20 ;
iGetMaxGpoint = 21 ;
iGetMinGB = 22 ;
iGetMinGpoint = 23 ;
iUpdatePressureGB = 24 ;
iUpdatePressureGpoint = 25 ;
iApplyBoundaryConditionGB = 26 ;
iApplyBoundaryConditionGpoint = 27 ;

DIR  = '.' ;
FILE = 'hpcscanGridShaheen' ;
TITLE = 'Test Case Grid' ;

pathFile = sprintf('%s/%s.log', DIR, FILE) ;
val = importdata(pathFile) ;

figure('Position',[100 100 1000 400])

subplot(1,2,1); hold on

yVal = [ val.data(1,iFillGB) val.data(2,iFillGB) ; ...
    val.data(1,iMaxErrGB) val.data(2,iMaxErrGB) ; ...
    val.data(1,iL1ErrGB) val.data(2,iL1ErrGB) ; ...
    val.data(1,iGetSumAbsGB) val.data(2,iGetSumAbsGB) ; ...
    val.data(1,iGetSumAbsDiffGB) val.data(2,iGetSumAbsDiffGB) ; ...
    val.data(1,iGetMaxGB) val.data(2,iGetMaxGB) ; ...
    val.data(1,iGetMinGB) val.data(2,iGetMinGB) ; ...
    val.data(1,iUpdatePressureGB) val.data(2,iUpdatePressureGB) ; ...
    %val.data(1,iApplyBoundaryConditionGB) val.data(2,iApplyBoundaryConditionGB) ; ...
    ]

cVal = categorical({'Fill','MaxErr','L1Err','SumAbs','SumAbsDiff','Max','Min','UptPressure'});
cVal = reordercats(cVal,{'Fill','MaxErr','L1Err','SumAbs','SumAbsDiff','Max','Min','UptPressure'});
bar(cVal,yVal)

ylabel('GByte/s'); title (TITLE)
grid on

subplot(1,2,2); hold on

yVal = [ val.data(1,iFillGpoint) val.data(2,iFillGpoint) ; ...
    val.data(1,iMaxErrGpoint) val.data(2,iMaxErrGpoint) ; ...
    val.data(1,iL1ErrGpoint) val.data(2,iL1ErrGpoint) ; ...
    val.data(1,iGetSumAbsGpoint) val.data(2,iGetSumAbsGpoint) ; ...
    val.data(1,iGetSumAbsDiffGpoint) val.data(2,iGetSumAbsDiffGpoint) ; ...
    val.data(1,iGetMaxGpoint) val.data(2,iGetMaxGpoint) ; ...
    val.data(1,iGetMinGpoint) val.data(2,iGetMinGpoint) ; ...
    val.data(1,iUpdatePressureGpoint) val.data(2,iUpdatePressureGpoint) ; ...
    %val.data(1,iApplyBoundaryConditionGpoint) val.data(2,iApplyBoundaryConditionGpoint) ; ...
    ]

bar(cVal,yVal)

ylabel('GPoint/s'); title (TITLE)
grid on

figFile = sprintf('%s.jpg', FILE) ;
print('-djpeg', figFile)

% END
