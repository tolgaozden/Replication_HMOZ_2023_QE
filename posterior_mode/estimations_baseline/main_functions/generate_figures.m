function [] =generate_figures(file_name,location)


 fig = gcf;
 fig.PaperPositionMode = 'auto'
 fig_pos = fig.PaperPosition;
 fig.PaperSize = [fig_pos(3) fig_pos(4)];
 
 

string=(char(strcat(location,'/',file_name)));
print(fig,string,'-dpdf');  
% savefig(fig,string);  


end