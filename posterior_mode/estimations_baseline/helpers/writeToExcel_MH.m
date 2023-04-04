% output=round(x,2);
%RUN AFTER metropolisHastings.m

  output_file='SW_Estimation_Results.csv';
 output_sheet='SW_ESTIMATION_RESULTS';
%  xlswrite(output_file,priorMean,output_sheet,'c15');
%  xlswrite(output_file,priorStd,output_sheet,'d15');
mode=round(mode,2);
posteriorMean=round(posteriorMean,2);
hpd_interval=round(hpd_interval,2);
 xlswrite(output_file,mode,output_sheet,'e3');
 xlswrite(output_file,posteriorMean,output_sheet,'f3');
 xlswrite(output_file,hpd_interval(:,1),output_sheet,'g3');
 xlswrite(output_file,hpd_interval(:,2),output_sheet,'h3');

%  xlswrite(output_file,laplace2,output_sheet,'e35');