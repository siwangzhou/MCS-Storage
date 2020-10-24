clear;
load phi_3d_500up_30000;
load T_normal;
data_org=T_normal;
addpath(genpath('..'));
denoiser1='BM3D'; 
% pno=1400;
iters=80;
ksize_row=16;
ksize_col=16;
ksize_time=300; 
timenum=50;
lblocksize=ksize_row;mblocksize=ksize_col;nbtimesize=50;
blocknum=ksize_row*ksize_col*ksize_time/(lblocksize*mblocksize*nbtimesize);
a=5;b=5;c=50;

%预分配好多内存
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
data_reshape=reshape(data_reshape,ksize_row*ksize_col*ksize_time,1);
sampleratearray=[0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14];
e(timenum,13)=0;ejue(timenum,13)=0;
for rate=1:13
    samplerate=sampleratearray(rate);
%      samplerate=0.205;
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
       end  %第一块
     end
     if count1>0
        dataphi(real,:,:)=phi(i,:,:);real=real+1;
     end
     if real==2*pno+1
         break; 
     end
  %  e1=mean(e);
   end
 for times=1:timenum
    clear M MM MM_pinv y data_hat ;
    ind=randperm(real-1);
    M(blocknum,8000,lblocksize*mblocksize*nbtimesize)=0;
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
         [mcol,mrow]=size(m1);
         m1=bsxfun(@rdivide,m1,sqrt(sum(m1,2))); 
         m1_pinv=pinv(m1);
         MM((M_row+1):M_row+mcol,((mm-1)*mrow+1):mm*mrow)=m1;
         MM_pinv(((mm-1)*mrow+1):mm*mrow,(M_row+1):M_row+mcol)=m1_pinv;
         M_row=M_row+mcol;
     end
      y=MM*data_reshape;
      data_hat = DAMPblock(y,iters,lblocksize,mblocksize,nbtimesize,ksize_row,ksize_col,ksize_time,denoiser1,MM,MM_pinv);
      e(times,rate)=norm(data_reshape-data_hat,2)/norm(data_reshape,2);
      ejue(times,rate)=norm(data_reshape-data_hat,1)/(ksize_row*ksize_col*ksize_time);
      fprintf('500_%.4f_rate=%.4f\n',samplerate,e(times,rate)); 
  
  %  e1=mean(e);
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
