clear;
load phi_3d_500up_30000;
load PM_normal600 ;
data_org=PM_normal;
addpath(genpath('..'));
denoiser1='BM3D'; 
% pno=1400;
iters=80;
ksize_row=16;
ksize_col=16;
ksize_time=300; 
timenum=5;
lblocksize=ksize_row;mblocksize=ksize_col;nbtimesize=50;
a=5;b=5;c=50;
%挑出所需人数2倍的经过候选区域的人。
samplerate=0.37;
pno1=round(samplerate*ksize_row*ksize_col*ksize_time);
  %分配空间
dataphi(2*pno1,700,3)=0;
real=1;
ind=randperm(30000);
  for x=1:30000
    i=ind(x);
    count1=0;
     for j=1:700
       if phi(i,j,1)>=a && phi(i,j,1)<(a+ksize_row) && phi(i,j,2)>=b && phi(i,j,2)<(b+ksize_col)&&phi(i,j,3)>=c && phi(i,j,3)<(c+ksize_time)
            count1=count1+1;
       end  %第一块
     end
     if count1>0
        dataphi(real,:,:)=phi(i,:,:);real=real+1;
     end
     if real==2*pno1+1
         break; 
     end
  %  e1=mean(e);
  end
%预分配好多内存

e(timenum,17)=0;ejue(timenum,17)=0;
for ro=11:12
    pno=round(samplerate*ksize_row*ksize_col*nbtimesize*ro);%采样人数
    data_reshape=zeros(lblocksize*mblocksize*nbtimesize,ro);
    num=1;
    for n=1:ro
      data=data_org(1:lblocksize,1:mblocksize,1+(n-1)*nbtimesize:n*nbtimesize);
      data=reshape(data,lblocksize*mblocksize*nbtimesize,1);
      data_reshape(:,num)=data;
      num=num+1;
    end
    data_reshape1=reshape(data_reshape,ksize_row*ksize_col*nbtimesize*ro,1);
    T_gui=data_reshape/415*255;%%PM25数据需归一化

 
 for times=1:timenum
    clear M MM MM_pinv y data_hat ;
    ind=randperm(real-1);
    M=zeros(ro,8000,lblocksize*mblocksize*nbtimesize);
    countm=1; k=1;
   for x=1:real
          i=ind(x);
          one1=zeros(ksize_row,ksize_col,nbtimesize*ro);
          for j=1:700
              if dataphi(i,j,1)>=a && dataphi(i,j,1)<(a+ksize_row) && dataphi(i,j,2)>=b && dataphi(i,j,2)<(b+ksize_col)&&dataphi(i,j,3)>=c && dataphi(i,j,3)<(c+nbtimesize*ro)
                    one1(dataphi(i,j,1)-a+1,dataphi(i,j,2)-b+1,dataphi(i,j,3)-c+1)=1;
              end  
          end
          one=zeros(lblocksize,mblocksize,nbtimesize);
          num=0;
       for n=1:ro
        one=one1(1:lblocksize,1:mblocksize,1+(n-1)*nbtimesize:n*nbtimesize);
        num=num+1;
        if length(find(one==1))>1
         one=reshape(one,1,lblocksize*mblocksize*nbtimesize);
         M(num,k,:)=one;countm=countm+1;
        end
       end
        k=k+1;
      if countm>pno
          break;
      end
   end 
     MM=zeros(countm-1,ksize_row*ksize_col*nbtimesize*ro);
     MM_pinv=zeros(ksize_row*ksize_col*nbtimesize*ro,countm-1);
     M_row=0;
     for mm=1:ro
         m1=squeeze(M(mm,:,:));
         m1(all(m1==0,2),:)=[];
         m1=bsxfun(@rdivide,m1,sqrt(sum(m1,2))); 
         m1_pinv=pinv(m1);
%          y=m1*data_reshape(:,mm);
          y=m1*T_gui(:,mm);%%pm数据
         data_hatblock = DAMP(y,iters,lblocksize*mblocksize,nbtimesize,denoiser1,m1,m1_pinv);
         data_hat(:,mm)=data_hatblock;
     end
        
        data_hat=reshape(data_hat,ksize_row*ksize_col*nbtimesize*ro,1);
        data_hat=data_hat*415/255;%pm数据
        e(times,ro)=norm(data_reshape1-data_hat,2)/norm(data_reshape1,2);
        ejue(times,ro)=norm(data_reshape1-data_hat,1)/(ksize_row*ksize_col*nbtimesize*ro);
        fprintf('blockround037—%d=%.4f\n',ro,e(times,ro)); 
 end
 filename=strcat('.\pmround_notfull\',num2str(samplerate),'_',num2str(ro),'_e.mat');
 save(filename,'e');
end

%%sum(countm1)/timenum