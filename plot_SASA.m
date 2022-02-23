####xfsdfsdfsdf

clear all


number_plot=4;
thickness_paramether=0.006; %%%0.006 or 0.005 could be fine by default
time_increment_MD_code_1=500;%% becareful!! multiply the time increment in LAMMPS by the "save increment" that used in Ovito 
time_increment_MD_code_2=100000; %% becareful!! multiply the time increment in LAMMPS by the "save increment" that used in Ovito 

% % % % red  	'#FF0000'
% % % % blue	'#0000FF'
% % % % green 	'#00FF00'
% % % % blak	'#000000'
% % % % Cyan	'#00FFFF'
% % % % magenta '#FF00FF'
% % % % orange	'#D95319'
% % % % yellow	'#FFFF00'
% % % % purple  '#7E2F8E'
% % % % brown	'#A2142F'

color_order={'#FF0000' '#0000FF' '#00FF00'  '#000000' '#00FFFF' '#FF00FF' '#7E2F8E' '#A2142F' '#FFFF00'};




filen = 'SASA_data.xlsx';

%%% dp2 1100K 0R
AA = xlsread(filen, 1)
FVVF=	AA';																																																																																																																																																																																																							
ii=1;
for ik=1:number_plot
 SASA_1=[];   
 Time_1=[];
 SASA_2=[];   
 Time_2=[];
Time_1=FVVF(ii,:);
SASA_1=FVVF(ii+1,:);


find_1=find(isnan(SASA_1));


YY_1=0;
for i=1:size(SASA_1,2)-size(find_1,2)
    if i>=2
    Time_1(1,i)=((i*time_increment_MD_code_1*1e-15)/1e-9);
    end
    
end
Time_1(1,1)=Time_1(1,2);
SASA_1(1,1)=SASA_1(1,2);
semilogx(Time_1,SASA_1/(SASA_1(1,1)),'Linewidth',2,'Color' , color_order{1,ik} )
set(gcf,'color','w');
xlabel('Time,{\it t}, ns');
ylabel('Normalized surface area, {\it a/a_0}');
ylim([ .79 1.008]);
hold on

Time_2=FVVF(ii+2,:);																																																																																																																																																																																																						
SASA_2=FVVF(ii+3,:);





YY_1=0;
for i=1:size(SASA_2,2)
   
    if i>=2
    Time_2(1,i)=Time_1(1,(size(SASA_1,2)-size(find_1,2)))+((i*time_increment_MD_code_2*1e-15)/1e-9);
    end
     Time_2(1,1)=Time_1(1,(size(SASA_1,2)-size(find_1,2)));
end

YY_1=0;

for j=1:size(SASA_2,2)-1

     
     
    y_variation=(SASA_2(1,j)/(SASA_1(1,1)))-(SASA_2(1,j+1)/(SASA_1(1,1))); 

for i=1:20
   
     r = -1+2*rand(1,1);
 
SASA_2_extnended(1,i,j)=((SASA_2(1,j)/(SASA_1(1,1)))*r*thickness_paramether)+(SASA_2(1,j)/(SASA_1(1,1)));
end


 
end

SASA_2_normalized_extended=[];
SASA_2_normalized=SASA_2/(SASA_1(1,1));
k=1;
for j=1:size(SASA_2,2)-1
   A=SASA_2_extnended(1,:,j);

 SASA_2_normalized_extended=[SASA_2_normalized_extended SASA_2_normalized(1,j) A ];
 k=k+20+1;
end

Time_2_extended=[];



Time_2_extended=Time_2(1,1):(Time_2(1,end)-Time_2(1,1))/(size(SASA_2_normalized_extended,2)-1):Time_2(1,end);

loglog(Time_2_extended,SASA_2_normalized_extended,'Linewidth',2,'Color' , color_order{1,ik})

hold on



i=i+4;
ii=ii+5;
end
1