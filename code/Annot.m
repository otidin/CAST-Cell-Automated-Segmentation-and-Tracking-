if doAnnot==1
hold on  
for annoti=1:length(annotvec)
yylim=ylim;
x=[annotvec(annoti),annotvec(annoti)]*expe.dt;
y=[yylim(1),yylim(2)];
plot(x,y,'--');
hold all
end
end
hold off