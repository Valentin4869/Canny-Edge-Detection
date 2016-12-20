function [fres]=canny2(img,sigma,lowT,highT)

c=3;
Gauss= @(v) exp(-(-c*v:c*v).^2./(2*v.^2))/(v.*sqrt(2*pi)) ;
dGauss=@(v) (-(-c*v:c*v)./(v.^2)).*Gauss(v);



def_T=false;

 if nargin==2
     def_T=true;   
 end



if (ndims(img)==3)
   
  img =double(rgb2gray(img));
end

img = imgaussfilt(img,sigma/2);
[m,n]=size(img);



fx = conv2(img, dGauss(sigma)');
fy = conv2(img, dGauss(sigma));




step=floor(sigma.^2);
%step=size(fy,2)-n;
fy=fy(:,(step+1):(n+step));
fx=fx((step+1):(m+step),:);

Theta = atan2( fy, fx) * (180.0/pi);
Magnitude = sqrt( fx.^2 + fy.^2 );



%eliminate negative angles and discretize angles to 4 ranges [0 45 90 135]

Theta = Theta + 180*(Theta<0);

Theta_d=0;
for i=0:3
Theta_d=Theta_d+(Theta>=45*i & Theta<(45*(i+1)))*(45*i);
end

I_N=Magnitude;

for i=1:m
    
    for j=1:n
        dir_vec=round([cos(Theta_d(i,j)) sin(Theta_d(i,j))]);
        pixel=[i,j];
        pixel_p= pixel+dir_vec;
        pixel_n= pixel-dir_vec;
        
        if sum(pixel_p==0+pixel_n==0)==0
            if Magnitude(i,j)<Magnitude(pixel_p(1),pixel_p(2)) || Magnitude(i,j)<Magnitude(pixel_n(1),pixel_n(2))
                I_N(i,j)=0;
            else
                I_N(i,j)=Magnitude(i,j);
            end
            
        else
            
            try
                m1=Magnitude(pixel_p(1),pixel_p(2));
            catch
                m1=0;
            end
            
            try
                m2=Magnitude(pixel_n(1),pixel_n(2));
            catch
                m2=0;
            end
            
             if Magnitude(i,j)<m1 || Magnitude(i,j)<m2
                I_N(i,j)=0;
            else
            
            end
            
        end
    end
end


if def_T
 lowT  = 0.4*mean2(Magnitude);
 highT = 2*lowT;
end




nl=zeros(size(I_N));
nh=zeros(size(I_N));

x=find(I_N>=lowT);

nl(x)=1;

x=find(I_N>=highT);

nh(x)=1;


fres=zeros(size(nl));

nc=[-1 0 1 -1 1 -1 0 1];
nr=[-1 -1 -1 0 0 1 1 1];
for i=1:m
    
    for j=1:n
       
      
        p=nh(i,j);
        for k=1:8
           
            if  nh(i,j)==1
            
            fres(i,j)=1;
            try
            pn=nl(i+nr(k),j+nc(k));
            catch
            continue;
            end

                if nl(i+nr(k),j+nc(k))==1
                    fres(i+nr(k),j+nc(k))=1;
                    
              
                end     
      
            end
        end
    end
   
   
end



end