
%分块重构 分块求伪逆 固定人数和固定采样率用50平均

clear all;
load ('.\trajectory\phi_3d_300up_30000.mat');%加载轨迹矩阵
load ('.\database\PM_normal.mat');%加载数据集
data_org=PM_normal;
addpath(genpath('..'));
denoiser1='BM3D'; 
% pno=1400;
iters=80;
ksize_row=16;
ksize_col=16;
ksize_time=300; 
timenum=3;%重复实验次数
lblocksize=16;mblocksize=16;nbtimesize=50;
blocknum=ksize_row*ksize_col*ksize_time/(lblocksize*mblocksize*nbtimesize);
a=5;b=5;c=50;
%预分配内存
data_reshape(lblocksize*mblocksize*nbtimesize,blocknum)=0;
num=1;
for l=1:ksize_row/lblocksize   
     for m=1:ksize_col/mblocksize
        for n=1:ksize_time/nbtimesize
        data=data_org(1+(l-1)*lblocksize:l*lblocksize,1+(m-1)*mblocksize:m*mblocksize,1+(n-1)*nbtimesize:n*nbtimesize);
        data=reshape(data,lblocksize*mblocksize*nbtimesize,1);
        data_reshape(:,num)=data;
        num=num+1;
        end
    end
end
data_reshape1=reshape(data_reshape,ksize_row*ksize_col*ksize_time,1);
T_gui=data_reshape/415*255;%pm数据事前需要进行归一化，而temp不需要。
sampleratearray=[0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4];
e(timenum,19)=0;ejue(timenum,19)=0;
for rate=1:19
    samplerate=sampleratearray(rate);
    pno=round(samplerate*ksize_row*ksize_col*ksize_time);
  %分配空间
   dataphi(2*pno,700,3)=0;
   M(blocknum,8000,lblocksize*mblocksize*nbtimesize)=0;
    %挑出所需人数2倍的经过候选区域的人。
   real=1;
   ind=randperm(30000);
  for x=1:30000
    i=ind(x);
    count1=0;
     for j=1:700
       if phi(i,j,1)>=a && phi(i,j,1)<(a+ksize_row) && phi(i,j,2)>=b && phi(i,j,2)<(b+ksize_col)&&phi(i,j,3)>=c && phi(i,j,3)<(c+ksize_time)
            count1=count1+1;
       end 
     end
     if count1>0
        dataphi(real,:,:)=phi(i,:,:);real=real+1;
     end
     if real==2*pno+1
         break; 
     end
   end
 for times=1:timenum
    clear M MM MM_pinv y data_hat ;
    ind=randperm(real-1);
    M(blocknum,8000,lblocksize*mblocksize*nbtimesize)=0;
    data_hat(lblocksize*mblocksize*nbtimesize,blocknum)=0;
    countm=1; k=1;
   for x=1:pno
          i=ind(x);
          one1=zeros(ksize_row,ksize_col,ksize_time);
          for j=1:700
              if dataphi(i,j,1)>=a && dataphi(i,j,1)<(a+ksize_row) && dataphi(i,j,2)>=b && dataphi(i,j,2)<(b+ksize_col)&&dataphi(i,j,3)>=c && dataphi(i,j,3)<(c+ksize_time)
                    one1(dataphi(i,j,1)-a+1,dataphi(i,j,2)-b+1,dataphi(i,j,3)-c+1)=1;
              end  
          end
          one=zeros(lblocksize,mblocksize,nbtimesize);
          num=0;
          for l=1:ksize_row/lblocksize   
             for m=1:ksize_col/mblocksize
               for n=1:ksize_time/nbtimesize
                one=one1(1+(l-1)*lblocksize:l*lblocksize,1+(m-1)*mblocksize:m*mblocksize,1+(n-1)*nbtimesize:n*nbtimesize);
                num=num+1;
                if length(find(one==1))>1
                 one=reshape(one,1,lblocksize*mblocksize*nbtimesize);
                 M(num,k,:)=one;countm=countm+1;
                end
               end
            end
          end
        k=k+1;
      if countm>pno
          break;
      end
   end 
     MM=zeros(countm-1,ksize_row*ksize_col*ksize_time);
     MM_pinv=zeros(ksize_row*ksize_col*ksize_time,countm-1);
     M_row=0;
     for mm=1:blocknum
         m1=squeeze(M(mm,:,:));
         m1(all(m1==0,2),:)=[];
         m1=bsxfun(@rdivide,m1,sqrt(sum(m1,2))); 
         m1_pinv=pinv(m1);
         y=m1*T_gui(:,mm);
         data_hatblock = DAMP(y,iters,lblocksize*mblocksize,nbtimesize,denoiser1,m1,m1_pinv);
         data_hat(:,mm)=data_hatblock;
     end
     data_hat=reshape(data_hat,ksize_row*ksize_col*ksize_time,1);
     x_hat1=data_hat*415/255;%pm数据需要反归一化，而temp数据不需要
     e(times,rate)=norm(data_reshape1-x_hat1,2)/norm(data_reshape1,2);
     ejue(times,rate)=norm(data_reshape1-x_hat1,1)/(ksize_row*ksize_col*ksize_time);
     fprintf('block%.2f_rate=%.4f\n',samplerate,e(times,rate)); 
      
  %将重构后的图像reshape成三维图形并以热度图显示
%  num=1;
%  dataorgback(ksize_row,ksize_col,ksize_time)=0;
%  data_hatorg=reshape(data_hat,[lblocksize*mblocksize*nbtimesize,ksize_row*ksize_col*ksize_time/(lblocksize*mblocksize*nbtimesize)]);
% for l=1:ksize_row/lblocksize   
%      for m=1:ksize_col/mblocksize
%         for n=1:ksize_time/nbtimesize
%             databack1=reshape(data_hatorg(:,num),lblocksize*mblocksize,nbtimesize);
%             dataorgback(1+(l-1)*lblocksize:l*lblocksize,1+(m-1)*mblocksize:m*mblocksize,1+(n-1)*nbtimesize:n*nbtimesize)=reshape(databack1,lblocksize,mblocksize,nbtimesize);
%             num=num+1;
%         end
%     end
% end
% imagesc(squeeze(dataorgback(:,:,80)));
 
  end
end