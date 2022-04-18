clear all
close all
clc



atom_number=2444;
%-------------------------------     atom_number=1760;
frame_number=1;
Au_lattice_constant=4.078;
Ni_lattice_constant=3.52;
Au_mass = 196.97; 
Ni_mass=58.693;
lattice_constant=Ni_lattice_constant; % 2A is more than the true value in MD file
mass=Ni_mass;
CN = 12; 
L = 6; % value of spherical harmonic degree for fcc structure

textline_per_frame=9;
atom_radius=0.25*sqrt(2*(lattice_constant^2));  %just for FCC
probe_radius=2.25;
% % % % % box_length=8.5826300000000003;
box_length=13;
total_mass = mass*atom_number;

length_bx=132;






box_max=1.3e2;
box_min=-2e0;

length_box= box_max-box_min;
%     !!!!!!!!!!!!!!!     Format Should be like this length_box=1.3e2-(-2e0);
ratio_atoms_alighn=0.02;
atom_number_alighn=round(atom_number*ratio_atoms_alighn);




%Disorder_parameter_total = zeros(atom_number,1);
%Disorder_parameter = zeros(natoms,(2*L+1),i);



nbins = 1000;
dr = 0.1*length_bx/nbins;



rpeak1=lattice_constant*sqrt(2)/2;
rpeak2=lattice_constant;




fid=fopen('BB');
%%%%%%fid=fopen('MD_Ni_NPs_ROOT2_1200K_R15_gap2_rotation60_09042021_1316pm');
A=textscan(fid, '%s', 'delimiter', '\n');


%%%   Please NOTE that the MD file should be created sumetrically and the particle should be located at the origin of coordination system 




n=0;
for i=1:frame_number
    for j=1:atom_number
   
DD=A{1}{n+textline_per_frame+j};
text_per_line=textscan(DD, '%s', 'delimiter', ' ');

atom_coordinate(j,1,i)=str2double(text_per_line{1}{1});
atom_coordinate(j,2,i)=str2double(text_per_line{1}{2});
atom_coordinate(j,3,i)=str2double(text_per_line{1}{3});
atom_coordinate(j,4,i)=str2double(text_per_line{1}{4});
atom_coordinate(j,5,i)=str2double(text_per_line{1}{5});



    end
    n=n+atom_number+textline_per_frame;
end



for i=1:frame_number
      for j=1:atom_number

find_id_1=find(atom_coordinate(:,1,i)==atom_coordinate(j,1,1));      
atom_coordinate_trimmed(j,1,i)=atom_coordinate(find_id_1,1,i);
atom_coordinate_trimmed(j,2,i)=atom_coordinate(find_id_1,2,i);
atom_coordinate_trimmed(j,3,i)=atom_coordinate(find_id_1,3,i);
atom_coordinate_trimmed(j,4,i)=atom_coordinate(find_id_1,4,i);
atom_coordinate_trimmed(j,5,i)=atom_coordinate(find_id_1,5,i);


 end
      end
atom_coordinate=[];
atom_coordinate=atom_coordinate_trimmed;




hist = zeros(nbins,1);

D_frame_saved=[];
D_frame_saved_new=[];

fid = fopen('MD_3_trimmed_steinhardt_test_aug.txt', 'w'); 
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

for i=1:frame_number
 %atom_coordinate_new=NPs_rotation_PBC_adjustment(atom_coordinate,i,atom_number,box_min,box_max,atom_number_alighn);
% % atom_coordinate_new=NPs_rotation_PBC_adjustment_newversion(atom_coordinate,i,atom_number,box_min,box_max,atom_number_alighn);
% % atom_coordinate=atom_coordinate_new;   
    
    
           q_lm = zeros(atom_number,2*L+1);
           Nb = zeros(atom_number,1);
           Ylm_Nb = zeros(atom_number,2*L+1);
           phi = zeros(atom_number-1,1);
           theta = zeros(atom_number-1,1);   
           r = zeros(nbins,1);
           NN_index = zeros(atom_number-1,atom_number-1);
           atoms_distance=zeros(atom_number-1,atom_number-1);
           
  
        
        for k=1:nbins
            r(k,1) = (k-0.1)*dr;
        

        end  
  
        
        
    
        [c1 indpeak1] = min(abs(r-rpeak1));
        [c2 indpeak2] = min(abs(r-rpeak2));    
        indmin=(indpeak1+indpeak2)/2;    
        
        r_cut=1*r(round(indmin));

%         if i<=15
%             r_cut=1*r(round(indmin));
%         end
%         if i>15
%             r_cut=1*r(round(indmin));
%         end
%         
for j = 1:atom_number
    
    
for jk = 1:atom_number
    
                dx = atom_coordinate(j,3,i)-atom_coordinate(jk,3,i);
                dy = atom_coordinate(j,4,i)-atom_coordinate(jk,4,i);
                dz = atom_coordinate(j,5,i)-atom_coordinate(jk,5,i);

                
                
                
                
                
               atoms_distance(j,jk)= sqrt(dx.^2+dy.^2+dz.^2);
               [phi_cur,theta_cur,r_coor] = cart2sph(dx,dy,dz);
               
               
               
               if (jk~=j)
                    if (atoms_distance(j,jk)<r_cut)
                        Nb(j,1) = Nb(j,1) + 1;
                        NN_index(j,jk) = jk;
                        theta(jk,1) = theta_cur + pi/2;
                        phi(jk,1) = phi_cur;
                    else
                    
                        NN_index(j,jk) = 0;
                        theta(jk,1) = 1e130;
                        phi(jk,1) = 1e130;
                    end
                else
               
                    theta(jk,1) = 1e130;
                    phi(jk,1) = 1e130;               
               
end

end


            if (Nb(j,1)==0)
                disp(j)
            end
            
            theta(theta>1e100) = [];
            phi(phi>1e100) = [];
            
            
            
            if (Nb(j,1)~=0)
                % Calculate the spherical harmonics  for each nearest neighbor of
                % each reference particle i and determine the complex quantity
                for jnn = 1:Nb(j,1)
                    THETA = theta(jnn,1);
                    PHI = phi(jnn,1);
                    for k = 1:(2*L+1)
                        M = k-(L+1);
                        [Ylm]=sph_harm(L,M,THETA,PHI);
                        Ylm_Nb(j,jnn,k) = Ylm;
                    end
                end
            end            
            
            
            
end


        for kl = 1:(2*L+1)
            for imn = 1:atom_number
                for jnn = 1:Nb(imn,1)
                    q_lm(imn,kl) = q_lm(imn,kl) + Ylm_Nb(imn,jnn,kl);
                end
            end
        end
        
        for imn = 1:atom_number
            q_lm(imn,:) = q_lm(imn,:)/sum(Nb(imn,1));
        end
        
        
count_NaN=0;
        D = zeros(atom_number,1);
        Dk = zeros(atom_number,(2*L+1));
        %     for k=1:(2*L+1)
        for ipo=1:atom_number
            jnn_in = find(NN_index(ipo,:)~=0);
            for jnn = 1:(Nb(ipo,1))
                jnn_in2 = jnn_in(jnn);
                Dk(ipo,:) = Dk(ipo,:) + abs(q_lm(ipo,:)-q_lm(jnn_in2,:)).^2;
            end
            D(ipo,1) = sum(Dk(ipo,:));
            D(ipo,1) = D(ipo,1)/Nb(ipo,1);
            if any(D(ipo,1))==0
               D(ipo,1)=0; 
               count_NaN=count_NaN+1;
            end
        end        
        
     
      
   D_ave = sum(D(:,1)./((atom_number-count_NaN)-Nb(:,1)))  
   D_ave_new = sum(D(:,1)./(atom_number-count_NaN))      
  D_ave_base = sum(D(:,1)./((atom_number-Nb(:,1)))) 
        
        D_frame_saved=[D_frame_saved,D_ave];
         D_frame_saved_new=[D_frame_saved_new,D_ave_new];       
        
        
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
        fprintf(fid,"ITEM:"+" "+"TIMESTEP"+" "+"\n"+frame_number+"\n");
fprintf(fid,"ITEM:"+" "+"NUMBER OF ATOMS"+" "+"\n"+ atom_number+"\n");
fprintf(fid,"ITEM: BOX BOUNDS pp pp pp"+"\n");
fprintf(fid,-box_length+"e+01"+" "+box_length+"e+01"+"\n");
fprintf(fid,-box_length+"e+01"+" "+box_length+"e+01"+"\n");
fprintf(fid,-box_length+"e+01"+" "+box_length+"e+01"+"\n");
fprintf(fid,"ITEM:"+" "+"ATOMS"+" "+"id"+" "+"type"+" "+"x"+" "+"y"+" "+"z"+" "+"c_3[1]"+"\n");
       
for atom_id=1:atom_number
         fprintf(fid,atom_id +" "+"1"+" "+atom_coordinate(atom_id,3,i)+" "+atom_coordinate(atom_id,4,i)+" "+atom_coordinate(atom_id,5,i)+" "+D(atom_id,1)+"\n");


end
        
 %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    i  
 %#$%$%TTDSETHBH     D_avarage=sum(sum(D_frame_saved)/atom_number)/size(D_frame_saved,2)   
end


 %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
mm=[];
D_avarage=[];
D_avarage=sum(D_frame_saved)/atom_number;
time_increment=1e-15;
frame_interval=100;
time_ns=0;
for i=1:3500
    
    time_ns(1,i)=((i*frame_interval*time_increment)/1e-9);


end




semilogx(time_ns,D_frame_saved,'b', 'Linewidth', 4)
xlabel('Time (ns)');
ylabel('Average Disorder Variable, D');
title('Ni/T=1200k/dp=3nm')

hold on;

filen = 'AA.xlsx';
AA = xlsread(filen)
FVVF=	AA';																																																																																																																																																																																																							
semilogx(FVVF(1,:),FVVF(2,:),'r', 'Linewidth', 4)


hold on;

filen = 'AA.xlsx';
AA = xlsread(filen)
FVVF=	AA';																																																																																																																																																																																																							
hold on;
semilogx(FVVF(5,:),FVVF(6,:),'r', 'Linewidth', 4)



for  i=1: frame_number
  mm=[mm, D_frame(1,i)];  
  for j=1:size(mm,2)  
    D_avarage=D_avarage+mm(1,j);

end

end

 %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
hold on
ds=[];
DDDSSA=0;
D_avarage_average=[];
for  i=1: size(D_frame_saved,2)
  ds=[D_frame_saved(1,1:i)];  
  for j=1:size(ds,2)  
     DDDSSA= DDDSSA+ds(1,j);
    

  end
DDDSSA= DDDSSA/size(ds,2);
D_avarage_average=[D_avarage_average,DDDSSA];
DDDSSA=0;
end

semilogx(time_ns,D_avarage_average,'g', 'Linewidth', 4)








correlation_mean=[];
aa=0;
for i=2: size(correlation_norm,1)
mm=[mm, correlation_norm(i,1)];
 for j=1:size(mm,2)
   aa=aa+mm(1,j);

 end
 correlation_mean=[correlation_mean;aa/size(mm,2)];
aa=0;
end
 










hold on

% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % for i=1:frame_number
% % % % % % % % % % % %            q_lm = zeros(atom_number,2*L+1);
% % % % % % % % % % % %            Nb = zeros(atom_number,1,1);
% % % % % % % % % % % %            Ylm_Nb = zeros(atom_number,2*L+1);
% % % % % % % % % % % %            phi = zeros(atom_number-1,1);
% % % % % % % % % % % %            theta = zeros(atom_number-1,1);   
% % % % % % % % % % % %            r = zeros(nbins,1);
% % % % % % % % % % % %            NN_index = zeros(atom_number-1,atom_number-1);
% % % % % % % % % % % %            
% % % % % % % % % % % %            
% % % % % % % % % % % % % Calculation of the position vectors    
% % % % % % % % % % % %   for im = 1:(atom_number-1)
% % % % % % % % % % % %             for ip = (im+1):atom_number
% % % % % % % % % % % %                 distx = atom_coordinate(im,3,i)-atom_coordinate(ip,3,i);
% % % % % % % % % % % %                 disty = atom_coordinate(im,4,i)-atom_coordinate(ip,4,i);
% % % % % % % % % % % %                 distz = atom_coordinate(im,5,i)-atom_coordinate(ip,5,i);
% % % % % % % % % % % % 
% % % % % % % % % % % %                 
% % % % % % % % % % % %                 dist= sqrt(distx.^2+disty.^2+distz.^2);
% % % % % % % % % % % %                 
% % % % % % % % % % % % 
% % % % % % % % % % % %                 if (im~=ip)
% % % % % % % % % % % %                     BIN = round(dist/dr)+1;
% % % % % % % % % % % %                     if (BIN <= nbins)
% % % % % % % % % % % %                         hist(BIN,1) = hist(BIN,1) + 1;
% % % % % % % % % % % %                     end
% % % % % % % % % % % %                 end
% % % % % % % % % % % %             end
% % % % % % % % % % % %   end
% % % % % % % % % % % %   
% % % % % % % % % % % %   
% % % % % % % % % % % %   
% % % % % % % % % % % %         
% % % % % % % % % % % %         for k=1:nbins
% % % % % % % % % % % %             r(k,1) = (k-0.1)*dr;
% % % % % % % % % % % %         
% % % % % % % % % % % % 
% % % % % % % % % % % %         end  
% % % % % % % % % % % %   
% % % % % % % % % % % %         
% % % % % % % % % % % %         
% % % % % % % % % % % %     
% % % % % % % % % % % %         [c1 indpeak1] = min(abs(r-rpeak1));
% % % % % % % % % % % %         [c2 indpeak2] = min(abs(r-rpeak2));    
% % % % % % % % % % % %         indmin=(indpeak1+indpeak2)/2;    
% % % % % % % % % % % %         r_cut=r(round(indmin));
% % % % % % % % % % % %         
% % % % % % % % % % % %         
% % % % % % % % % % % % for j = 1:atom_number
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % % for jk = 1:atom_number
% % % % % % % % % % % %     
% % % % % % % % % % % %                 dx = atom_coordinate(j,3,i)-atom_coordinate(jk,3,i);
% % % % % % % % % % % %                 dy = atom_coordinate(j,4,i)-atom_coordinate(jk,4,i);
% % % % % % % % % % % %                 dz = atom_coordinate(j,5,i)-atom_coordinate(jk,5,i);
% % % % % % % % % % % % 
% % % % % % % % % % % %                atoms_distance(j,jk)= sqrt(dx.^2+dy.^2+dz.^2);
% % % % % % % % % % % %                [phi_cur,theta_cur,r_coor] = cart2sph(dx,dy,dz);
% % % % % % % % % % % %                
% % % % % % % % % % % %                
% % % % % % % % % % % %                
% % % % % % % % % % % %                if (jk~=j)
% % % % % % % % % % % %                     if (atoms_distance(j,jk)<r_cut)
% % % % % % % % % % % %                         Nb(j,1) = Nb(j,1) + 1;
% % % % % % % % % % % %                         NN_index(j,jk) = jk;
% % % % % % % % % % % %                         theta(jk,1) = theta_cur + pi/2;
% % % % % % % % % % % %                         phi(jk,1) = phi_cur;
% % % % % % % % % % % %                     else
% % % % % % % % % % % %                     
% % % % % % % % % % % %                         NN_index(j,jk) = 0;
% % % % % % % % % % % %                         theta(jk,1) = 1e130;
% % % % % % % % % % % %                         phi(jk,1) = 1e130;
% % % % % % % % % % % %                     end
% % % % % % % % % % % %                 else
% % % % % % % % % % % %                
% % % % % % % % % % % %                     theta(jk,1) = 1e130;
% % % % % % % % % % % %                     phi(jk,1) = 1e130;               
% % % % % % % % % % % %                
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % %             if (Nb(j,1)==0)
% % % % % % % % % % % %                 disp(j)
% % % % % % % % % % % %             end
% % % % % % % % % % % %             
% % % % % % % % % % % %             theta(theta>1e100) = [];
% % % % % % % % % % % %             phi(phi>1e100) = [];
% % % % % % % % % % % %             
% % % % % % % % % % % %             
% % % % % % % % % % % %             
% % % % % % % % % % % %             if (Nb(j,1)~=0)
% % % % % % % % % % % %                 % Calculate the spherical harmonics  for each nearest neighbor of
% % % % % % % % % % % %                 % each reference particle i and determine the complex quantity
% % % % % % % % % % % %                 for jnn = 1:Nb(j,1)
% % % % % % % % % % % %                     THETA = theta(jnn,1);
% % % % % % % % % % % %                     PHI = phi(jnn,1);
% % % % % % % % % % % %                     for k = 1:(2*L+1)
% % % % % % % % % % % %                         M = k-(L+1);
% % % % % % % % % % % %                         [Ylm]=sph_harm(L,M,THETA,PHI);
% % % % % % % % % % % %                         Ylm_Nb(i,jnn,k) = Ylm;
% % % % % % % % % % % %                     end
% % % % % % % % % % % %                 end
% % % % % % % % % % % %             end            
% % % % % % % % % % % %             
% % % % % % % % % % % %             
% % % % % % % % % % % %             
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % %         for kl = 1:(2*L+1)
% % % % % % % % % % % %             for imn = 1:atom_number
% % % % % % % % % % % %                 for jnn = 1:Nb(imn,1)
% % % % % % % % % % % %                     q_lm(imn,kl) = q_lm(imn,kl) + Ylm_Nb(imn,jnn,kl);
% % % % % % % % % % % %                 end
% % % % % % % % % % % %             end
% % % % % % % % % % % %         end
% % % % % % % % % % % %         
% % % % % % % % % % % %         for imn = 1:atom_number
% % % % % % % % % % % %             q_lm(imn,:) = q_lm(imn,:)/sum(Nb(imn,1));
% % % % % % % % % % % %         end
% % % % % % % % % % % %         
% % % % % % % % % % % %         
% % % % % % % % % % % % 
% % % % % % % % % % % %         D = zeros(atom_number,1);
% % % % % % % % % % % %         Dk = zeros(atom_number,(2*L+1));
% % % % % % % % % % % %         %     for k=1:(2*L+1)
% % % % % % % % % % % %         for ipo=1:atom_number
% % % % % % % % % % % %             jnn_in = find(NN_index(ipo,:)~=0);
% % % % % % % % % % % %             for jnn = 1:(Nb(ipo,1))
% % % % % % % % % % % %                 jnn_in2 = jnn_in(jnn);
% % % % % % % % % % % %                 Dk(ipo,:) = Dk(ipo,:) + abs(q_lm(ipo,:)-q_lm(jnn_in2,:)).^2;
% % % % % % % % % % % %             end
% % % % % % % % % % % %             D(ipo,1) = sum(Dk(ipo,:));
% % % % % % % % % % % %             D(ipo,1) = D(ipo,1)/Nb(ipo,1);
% % % % % % % % % % % %         end        
% % % % % % % % % % % %         
% % % % % % % % % % % %         
% % % % % % % % % % % %         
% % % % % % % % % % % %         
% % % % % % % % % % % %         
% % % % % % % % % % % %         
% % % % % % % % % % % %          D_frame=[D_frame;D];
% % % % % % % % % % % %         
% % % % % % % % % % % %         
% % % % % % % % % % % %         
% % % % % % % % % % % %         
% % % % % % % % % % % %        
% % % % % % % % % % % %         
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % xx=x(:,1);
% % % % % % % % % % % % yy=y(:,1);
% % % % % % % % % % % % zz=z(:,1);
% % % % % % % % % % % % figure
% % % % % % % % % % % % S = repmat([80],numel(x),1);
% % % % % % % % % % % % s = S(:);
% % % % % % % % % % % % h=scatter3(xx,yy,zz,s,D,'fill');
% % % % % % % % % % % % set(h,'MarkerEdgeColor','k')
% % % % % % % % % % % % colorbar
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % mean_aggl=zeros(frame_number,3);
% % % % % % % % % % % % 
% % % % % % % % % % % % for i=1:frame_number
% % % % % % % % % % % %       for j=1:atom_number
% % % % % % % % % % % %           
% % % % % % % % % % % % mean_aggl(i,1)=mean_aggl(i,1)+atom_coordinate(j,3,i);
% % % % % % % % % % % % mean_aggl(i,2)=mean_aggl(i,2)+atom_coordinate(j,4,i);    
% % % % % % % % % % % % mean_aggl(i,3)=mean_aggl(i,3)+atom_coordinate(j,5,i);
% % % % % % % % % % % % 
% % % % % % % % % % % %       end
% % % % % % % % % % % %  
% % % % % % % % % % % %   
% % % % % % % % % % % % end
% % % % % % % % % % % %      mean_agglomerate=mean_aggl/atom_number;
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % for i=1:frame_number
% % % % % % % % % % % %       for j=1:atom_number
% % % % % % % % % % % % 
% % % % % % % % % % % %         find_id_1=find(atom_coordinate(:,1,i)==atom_coordinate(j,1,1))  ;      
% % % % % % % % % % % % atom_coordinate_trimmed(j,1,i)=atom_coordinate(find_id_1,1,i);
% % % % % % % % % % % % atom_coordinate_trimmed(j,2,i)=atom_coordinate(find_id_1,2,i);
% % % % % % % % % % % % atom_coordinate_trimmed(j,3,i)=atom_coordinate(find_id_1,3,i);
% % % % % % % % % % % % atom_coordinate_trimmed(j,4,i)=atom_coordinate(find_id_1,4,i);
% % % % % % % % % % % % atom_coordinate_trimmed(j,5,i)=atom_coordinate(find_id_1,5,i);
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % %  end
% % % % % % % % % % % %       end
% % % % % % % % % % % %  
% % % % % % % % % % % % atom_coordinate=[];  
% % % % % % % % % % % % atom_coordinate=atom_coordinate_trimmed;
% % % % % % % % % % % %      
% % % % % % % % % % % %      
% % % % % % % % % % % %      
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % atom_id_saved=[];
% % % % % % % % % % % % fid = fopen('MD_1_trimmed5.txt', 'w'); 
% % % % % % % % % % % % for i=1:frame_number
% % % % % % % % % % % %     
% % % % % % % % % % % % Atom=[];    
% % % % % % % % % % % % shp=[];    
% % % % % % % % % % % % polar_coord=[];
% % % % % % % % % % % % ijk=1;
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % coord_min= min(atom_coordinate(:,:,i));
% % % % % % % % % % % % coord_max=max(atom_coordinate(:,:,i));
% % % % % % % % % % % % 
% % % % % % % % % % % % min_x= coord_min(1,3);
% % % % % % % % % % % % min_y= coord_min(1,4);
% % % % % % % % % % % % min_z= coord_min(1,5);
% % % % % % % % % % % % 
% % % % % % % % % % % % max_x= coord_max(1,3);
% % % % % % % % % % % % max_y= coord_max(1,4);
% % % % % % % % % % % % max_z= coord_max(1,5);  
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % centroid_x=0;
% % % % % % % % % % % % centroid_y=0;
% % % % % % % % % % % % centroid_z=0;
% % % % % % % % % % % % for mno=1:atom_number
% % % % % % % % % % % %   
% % % % % % % % % % % %     centroid_x=centroid_x+atom_coordinate(mno,3,i);
% % % % % % % % % % % %     centroid_y=centroid_y+atom_coordinate(mno,4,i);
% % % % % % % % % % % %     centroid_z=centroid_z+atom_coordinate(mno,5,i);
% % % % % % % % % % % % end
% % % % % % % % % % % %     centroid_x=centroid_x/atom_number;
% % % % % % % % % % % %     centroid_y=centroid_y/atom_number;
% % % % % % % % % % % %     centroid_z=centroid_z/atom_number;
% % % % % % % % % % % % 
% % % % % % % % % % % % for mno=1:atom_number
% % % % % % % % % % % %     
% % % % % % % % % % % %     atom_coordinate(mno,3,i)=atom_coordinate(mno,3,i)-centroid_x;
% % % % % % % % % % % %     atom_coordinate(mno,4,i)=atom_coordinate(mno,4,i)-centroid_y;
% % % % % % % % % % % %     atom_coordinate(mno,5,i)=atom_coordinate(mno,5,i)-centroid_z;
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % Atom=[atom_coordinate(:,3,i),atom_coordinate(:,4,i),atom_coordinate(:,5,i)];
% % % % % % % % % % % % 
% % % % % % % % % % % % AA(:,1,i)=Atom(:,1);
% % % % % % % % % % % % AA(:,2,i)=Atom(:,2);
% % % % % % % % % % % % AA(:,3,i)=Atom(:,3);
% % % % % % % % % % % % 
% % % % % % % % % % % % shp=[];
% % % % % % % % % % % % shp = alphaShape(Atom(:,1),Atom(:,2),Atom(:,3),3);
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % %    for nn=1:atom_number
% % % % % % % % % % % % 
% % % % % % % % % % % %        
% % % % % % % % % % % %        
% % % % % % % % % % % %        
% % % % % % % % % % % %        [polar_coord(1,nn),polar_coord(2,nn),polar_coord(3,nn)] = cart2sph(atom_coordinate(nn,3,i),atom_coordinate(nn,4,i),atom_coordinate(nn,5,i));
% % % % % % % % % % % % %         polar_coord(1,nn)=mod(cart2sph(relative_coord_y(1,nn),relative_coord_z(1,nn)),2*pi);
% % % % % % % % % % % % %         polar_coord(2,nn)=mod(cart2sph(relative_coord_y(1,nn),relative_coord_z(1,nn)),2*pi);        
% % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % % % % % %         polar_coord(1,nn) =rad2deg(polar_coord(1,nn));
% % % % % % % % % % % % % % % % % % % % % %         polar_coord(2,nn) =rad2deg(polar_coord(2,nn));
% % % % % % % % % % % %         
% % % % % % % % % % % %         
% % % % % % % % % % % % % [azimuth,elevation,r] = cart2sph
% % % % % % % % % % % % 
% % % % % % % % % % % % trial_movment_r=polar_coord(3,nn)+atom_radius;
% % % % % % % % % % % % 
% % % % % % % % % % % % [trial_movment_x,trial_movment_y,trial_movment_z]=sph2cart(polar_coord(1,nn),polar_coord(2,nn),trial_movment_r);
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % if inShape(shp,trial_movment_x,trial_movment_y,trial_movment_z)==0
% % % % % % % % % % % %     
% % % % % % % % % % % %    surface_atom_id(ijk,i)=atom_coordinate(nn,1,i); 
% % % % % % % % % % % %     ijk=ijk+1;
% % % % % % % % % % % %     
% % % % % % % % % % % % end
% % % % % % % % % % % %    end
% % % % % % % % % % % % 
% % % % % % % % % % % %   inner_atom_id=[];
% % % % % % % % % % % %    for imn=1:atom_number
% % % % % % % % % % % %     count_not_inner=0;
% % % % % % % % % % % %        
% % % % % % % % % % % %        for rst=1:size(surface_atom_id,1)
% % % % % % % % % % % %            
% % % % % % % % % % % %            if  atom_coordinate(imn,1,i)~=surface_atom_id(rst,i)
% % % % % % % % % % % %             count_not_inner=count_not_inner+1;   
% % % % % % % % % % % % 
% % % % % % % % % % % %        
% % % % % % % % % % % %            end
% % % % % % % % % % % %    
% % % % % % % % % % % %        
% % % % % % % % % % % %        end
% % % % % % % % % % % %          if count_not_inner==size(surface_atom_id(:,i),1)
% % % % % % % % % % % %               inner_atom_id=[inner_atom_id;atom_coordinate(imn,1,i)];
% % % % % % % % % % % %        
% % % % % % % % % % % %          end
% % % % % % % % % % % %    
% % % % % % % % % % % %    end
% % % % % % % % % % % %    
% % % % % % % % % % % %    
% % % % % % % % % % % %    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    for ij=1:size(surface_atom_id,1)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         ABC=find(atom_coordinate==surface_atom_id(ij,i));   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %        plot3(atom_coordinate(surface_atom_id(ij,i),3,i),atom_coordinate(surface_atom_id(ij,i),4,i),atom_coordinate(surface_atom_id(ij,i),5,i),'.','MarkerSize', 20,'color','b')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %        hold on
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    end
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % fprintf(fid,"ITEM:"+" "+"TIMESTEP"+" "+"\n"+frame_number+"\n");
% % % % % % % % % % % % fprintf(fid,"ITEM:"+" "+"NUMBER OF ATOMS"+" "+"\n"+ atom_number+"\n");
% % % % % % % % % % % % fprintf(fid,"ITEM: BOX BOUNDS pp pp pp"+"\n");
% % % % % % % % % % % % fprintf(fid,-box_length+"e+01"+" "+box_length+"e+01"+"\n");
% % % % % % % % % % % % fprintf(fid,-box_length+"e+01"+" "+box_length+"e+01"+"\n");
% % % % % % % % % % % % fprintf(fid,-box_length+"e+01"+" "+box_length+"e+01"+"\n");
% % % % % % % % % % % % fprintf(fid,"ITEM:"+" "+"ATOMS"+" "+"id"+" "+"type"+" "+"x"+" "+"y"+" "+"z"+"\n");
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % for atom_id=1:atom_number
% % % % % % % % % % % % 
% % % % % % % % % % % %     number_not=0;
% % % % % % % % % % % %     for in=1: size(surface_atom_id,1)
% % % % % % % % % % % %         
% % % % % % % % % % % %     find_id=surface_atom_id(in,i)==atom_coordinate(atom_id,1,i)  ;    
% % % % % % % % % % % %     
% % % % % % % % % % % %      
% % % % % % % % % % % %      if find_id==1
% % % % % % % % % % % %              fprintf(fid,atom_id +" "+"1"+" "+atom_coordinate(atom_id,3,i)+" "+atom_coordinate(atom_id,4,i)+" "+atom_coordinate(atom_id,5,i)+"\n");
% % % % % % % % % % % %     
% % % % % % % % % % % %            
% % % % % % % % % % % %      end
% % % % % % % % % % % %      
% % % % % % % % % % % %  
% % % % % % % % % % % % 
% % % % % % % % % % % %     end
% % % % % % % % % % % %     
% % % % % % % % % % % %      for in=1: size(inner_atom_id,1)
% % % % % % % % % % % %         
% % % % % % % % % % % %     find_id=inner_atom_id(in,1)==atom_coordinate(atom_id,1,i)  ;    
% % % % % % % % % % % %     
% % % % % % % % % % % %      
% % % % % % % % % % % %      if find_id==1
% % % % % % % % % % % %              fprintf(fid,atom_id +" "+"2"+" "+atom_coordinate(atom_id,3,i)+" "+atom_coordinate(atom_id,4,i)+" "+atom_coordinate(atom_id,5,i)+"\n");
% % % % % % % % % % % %     
% % % % % % % % % % % %            
% % % % % % % % % % % %      end
% % % % % % % % % % % %      
% % % % % % % % % % % %  
% % % % % % % % % % % % 
% % % % % % % % % % % %     end
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % % % % % % %     for in=1: size(surface_atom_id,1)
% % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % %     find_id=surface_atom_id(in,i)==atom_id  ;    
% % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % %      
% % % % % % % % % % % % % % % % %      if find_id==0
% % % % % % % % % % % % % % % % %          number_not=number_not+1;
% % % % % % % % % % % % % % % % %              
% % % % % % % % % % % % % % % % %   
% % % % % % % % % % % % % % % % %       
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %              
% % % % % % % % % % % % % % % % %      end
% % % % % % % % % % % % % % % % %      
% % % % % % % % % % % % % % % % %    if number_not== size(surface_atom_id,1) 
% % % % % % % % % % % % % % % % %      fprintf(fid,atom_id +" "+"2"+" "+atom_coordinate(atom_id,3,i)+" "+atom_coordinate(atom_id,4,i)+" "+atom_coordinate(atom_id,5,i)+"\n");
% % % % % % % % % % % % % % % % %    end
% % % % % % % % % % % % % % % % %      
% % % % % % % % % % % % % % % % %    
% % % % % % % % % % % % % % % % %  
% % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % %     end
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % % 
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % end
% % % % % % % % % % % % fclose(fid);
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % [X,Y,Z] = meshgrid(min_x:1:max_x, min_y:1:max_y,min_z:1:max_z);
% % % % % % % % % % % % in = inShape(shp,X,Y,Z);
% % % % % % % % % % % % plot(shp)
% % % % % % % % % % % % 
% % % % % % % % % % % % hold on 
% % % % % % % % % % % % plot3(X(in),Y(in),Z(in),'r.')
% % % % % % % % % % % % plot3(X(~in),Y(~in),Z(~in),'b.')
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % clear all
% % % % % % % % % % % % % % % % % % close all
% % % % % % % % % % % % % % % % % % clc
% % % % % % % % % % % % % % % % % % atom_number=2552;
% % % % % % % % % % % % % % % % % % frame_number=1;
% % % % % % % % % % % % % % % % % % textline_per_frame=9;
% % % % % % % % % % % % % % % % % % x_increment=.1;
% % % % % % % % % % % % % % % % % % tolerance_theta=40;
% % % % % % % % % % % % % % % % % % tolerance_x=2;
% % % % % % % % % % % % % % % % % % lattice_constant=3.52; % 2A is more than the true value in MD file
% % % % % % % % % % % % % % % % % % atom_radius=0.25*sqrt(2*(lattice_constant^2));  %just for FCC
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % fid=fopen('MD_1');
% % % % % % % % % % % % % % % % % % A=textscan(fid, '%s', 'delimiter', '\n');
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % n=0;
% % % % % % % % % % % % % % % % % % for i=1:frame_number
% % % % % % % % % % % % % % % % % %     for j=1:atom_number
% % % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % % D=A{1}{j+textline_per_frame+n};
% % % % % % % % % % % % % % % % % % text_per_line=textscan(D, '%s', 'delimiter', ' ');
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % atom_coordinate(j,1,i)=str2double(text_per_line{1}{1});
% % % % % % % % % % % % % % % % % % atom_coordinate(j,2,i)=str2double(text_per_line{1}{2});
% % % % % % % % % % % % % % % % % % atom_coordinate(j,3,i)=str2double(text_per_line{1}{3});
% % % % % % % % % % % % % % % % % % atom_coordinate(j,4,i)=str2double(text_per_line{1}{4});
% % % % % % % % % % % % % % % % % % atom_coordinate(j,5,i)=str2double(text_per_line{1}{5});
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % % % %     n=n+atom_number;
% % % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % mean_aggl=zeros(frame_number,3);
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % for i=1:frame_number
% % % % % % % % % % % % % % % % % %       for j=1:atom_number
% % % % % % % % % % % % % % % % % %           
% % % % % % % % % % % % % % % % % % mean_aggl(i,1)=mean_aggl(i,1)+atom_coordinate(j,3,i);
% % % % % % % % % % % % % % % % % % mean_aggl(i,2)=mean_aggl(i,2)+atom_coordinate(j,4,i);    
% % % % % % % % % % % % % % % % % % mean_aggl(i,3)=mean_aggl(i,3)+atom_coordinate(j,5,i);
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % %       end
% % % % % % % % % % % % % % % % % %  
% % % % % % % % % % % % % % % % % %   
% % % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % % %      mean_agglomerate=mean_aggl/atom_number;
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % atom_id_saved=[];
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % for ii=1:frame_number
% % % % % % % % % % % % % % % % % % coord_min= min(atom_coordinate(:,:,ii));
% % % % % % % % % % % % % % % % % % coord_max=max(atom_coordinate(:,:,ii));
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % min_x= coord_min(1,3);
% % % % % % % % % % % % % % % % % % min_y= coord_min(1,4);
% % % % % % % % % % % % % % % % % % min_z= coord_min(1,5);
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % max_x= coord_max(1,3);
% % % % % % % % % % % % % % % % % % max_y= coord_max(1,4);
% % % % % % % % % % % % % % % % % % max_z= coord_max(1,5);  
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % moving_x=max_x;
% % % % % % % % % % % % % % % % % % number_increment=round((max_x-min_x)/x_increment);
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % surface_atom_id_all_slice=[];
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % for x_track=1:number_increment
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % %     for j=1:atom_number
% % % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % %                  if moving_x-tolerance_x<=atom_coordinate(j,3,ii) && atom_coordinate(j,3,ii)<=moving_x+tolerance_x
% % % % % % % % % % % % % % % % % %                    
% % % % % % % % % % % % % % % % % %                    atom_id_saved=[atom_id_saved,atom_coordinate(j,1,ii)];  
% % % % % % % % % % % % % % % % % %                      
% % % % % % % % % % % % % % % % % %                  end
% % % % % % % % % % % % % % % % % %                  
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % %    atom_id_saved_condensed=unique(atom_id_saved); 
% % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % %    mn_y=0;
% % % % % % % % % % % % % % % % % %    mn_z=0;
% % % % % % % % % % % % % % % % % %    
% % % % % % % % % % % % % % % % % %    
% % % % % % % % % % % % % % % % % % if  size(atom_id_saved_condensed,1) >=1   
% % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % %     for k=1:size(atom_id_saved_condensed,2)
% % % % % % % % % % % % % % % % % %        ABC=find(atom_coordinate==atom_id_saved_condensed(1,k));   
% % % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % %         mn_y=mn_y+atom_coordinate(ABC,4,ii);
% % % % % % % % % % % % % % % % % %         mn_z=mn_z+atom_coordinate(ABC,5,ii);
% % % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % % % %     centroid_slice_y=mn_y/size(atom_id_saved_condensed,2);
% % % % % % % % % % % % % % % % % %     centroid_slice_z=mn_z/size(atom_id_saved_condensed,2);
% % % % % % % % % % % % % % % % % %  
% % % % % % % % % % % % % % % % % %     for k=1:size(atom_id_saved_condensed,2)
% % % % % % % % % % % % % % % % % %      DEF=find(atom_coordinate==atom_id_saved_condensed(1,k));    
% % % % % % % % % % % % % % % % % %      atom_dist(1,k)=sqrt(((atom_coordinate(DEF,4,ii)-centroid_slice_y)^2)+((atom_coordinate(DEF,5,ii)-centroid_slice_z)^2)); 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % %     end  
% % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % surface_atom_id=surface_atoms(atom_coordinate,atom_id_saved_condensed,atom_dist,atom_radius,atom_number,i,tolerance_theta,centroid_slice_y,centroid_slice_z)
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % surface_atom_id_all_slice=[surface_atom_id_all_slice,surface_atom_id];
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % remember to delete the zeroes in the matrix
% % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % %     atom_dist=[];
% % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % %     atom_id_saved=[];
% % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % mean_agglomerate=0;    
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % %     moving_x=moving_x-x_increment
% % % % % % % % % % % % % % % % % %   
% % % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % % % surface_atom_id_all_slice=unique(surface_atom_id_all_slice);
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % surface_atom_id_all_frame=[surface_atom_id_all_frame;surface_atom_id_all_slice];
% % % % % % % % % % % % % % % % % % surface_atom_id_all_slice=[];
% % % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % for j=1:atom_number
% % % % % % % % % % % % % % % % % % plot3(atom_coordinate(j,3,1),atom_coordinate(j,4,1),atom_coordinate(j,5,1),'.')
% % % % % % % % % % % % % % % % % % hold on
% % % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % for j=1:size(surface_atom_id_all_slice,2)
% % % % % % % % % % % % % % % % % %     ABC=find(atom_coordinate==surface_atom_id_all_slice(1,j));    
% % % % % % % % % % % % % % % % % % plot3(atom_coordinate(ABC,3,1),atom_coordinate(ABC,4,1),atom_coordinate(ABC,5,1),'.','MarkerSize', 30)
% % % % % % % % % % % % % % % % % % hold on
% % % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % 
