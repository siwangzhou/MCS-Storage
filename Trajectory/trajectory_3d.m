%构造三维轨迹矩阵30*30*700
clear;
siterow=30;sitecol=30;
height=700;
sample=30000;
sitematrix=zeros(siterow,sitecol);
phi=zeros(sample,height,3);
a=1;
for i=1:siterow
    for j=1:sitecol
        sitematrix(i,j)=a;
        a=a+1;
    end
end
for k=1:sample
    m=0;n=0;m1=zeros(1,1);n1=zeros(1,1);
    M=zeros(siterow*siterow,height);
    startsite=randperm(siterow*sitecol,1);
    step=400+randperm(height-400-1,1);%给随机步数制定下限
    starttime=randperm(height-step,1);
%     step=700;
%     starttime=1;
    ID=startsite;
    for  i=starttime:(starttime+step-1)     
        M(ID,i)=1;
        [idi,idj]=find(sitematrix==ID);
         x=rand();
        if idi==1&&idj>=2&&idj<=sitecol-1 %% 第一行特殊情况
            ID1(x<=1/4)=sitematrix(idi,idj-1);
            ID1(x>1/4&&x<=2/4 )=sitematrix(idi,idj+1);
            ID1(x>2/4&&x<=3/4)=sitematrix(idi+1,idj);
            ID1(x>3/4)=sitematrix(idi,idj);
        else if idi==siterow && idj>=2&&idj<=sitecol-1 %% 最后一行特殊情况
            ID1(x<=1/4)=sitematrix(idi,idj-1);
            ID1(x>1/4&&x<=2/4 )=sitematrix(idi,idj+1);
            ID1(x>2/4&&x<=3/4)=sitematrix(idi-1,idj);
            ID1(x>3/4)=sitematrix(idi,idj);
        else if  idj==1&& idi>=2&&idi<=siterow-1 %%第一列特殊情况
            ID1(x<=1/4)=sitematrix(idi+1,idj);
            ID1(x>1/4&&x<=2/4 )=sitematrix(idi-1,idj);
            ID1(x>2/4&&x<=3/4)=sitematrix(idi,idj+1); 
            ID1(x>3/4)=sitematrix(idi,idj);
        else if  idj==sitecol&& idi>=2&&idi<=siterow-1 %%最后一列特殊情况
            ID1(x<=1/4)=sitematrix(idi+1,idj);
            ID1(x>1/4&&x<=2/4 )=sitematrix(idi-1,idj);
            ID1(x>2/4&&x<=3/4)=sitematrix(idi,idj-1); 
            ID1(x>3/4)=sitematrix(idi,idj);
        else if  idi==1&&idj==1
              ID1(x<=1/3)=sitematrix(idi+1,idj);
              ID1(x>1/3&&x<=2/3 )=sitematrix(idi,idj+1); 
              ID1(x>2/3)=sitematrix(idi,idj);
        else if  idi==1&&idj==sitecol
              ID1(x<=1/3)=sitematrix(idi+1,idj);
              ID1(x>1/3&&x<=2/3 )=sitematrix(idi,idj-1);  
              ID1(x>2/3)=sitematrix(idi,idj);
        else if  idi==siterow&&idj==1
              ID1(x<=1/3)=sitematrix(idi-1,idj);
              ID1(x>1/3&&x<=2/3 )=sitematrix(idi,idj+1);
              ID1(x>2/3)=sitematrix(idi,idj);
        else if  idi==siterow&&idj==sitecol
              ID1(x<=1/3)=sitematrix(idi-1,idj);
              ID1(x>1/3&&x<=2/3 )=sitematrix(idi,idj-1);
              ID1(x>2/3)=sitematrix(idi,idj);
            else
            ID1(x<=1/5)=sitematrix(idi,idj-1);
            ID1(x>1/5&&x<=2/5 )=sitematrix(idi,idj+1);
            ID1(x>2/5&&x<=3/5)=sitematrix(idi+1,idj);
            ID1(x>3/5&&x<=4/5)=sitematrix(idi-1,idj);
            ID1(x>4/5)=sitematrix(idi,idj);
            end
            end
            end
            end
            end
            end
            end
        end
       ID=ID1; 
    end
 [r,c]=find(M==1);
 for w=1:length(r)
     [m,n]=find(sitematrix==r(w));
     m1(w)=m;n1(w)=n;
 end
    m1=m1';n1=n1';
    m1=[zeros(c(1)-1,1);m1;zeros(height-length(r)-c(1)+1,1)];
    n1=[zeros(c(1)-1,1);n1;zeros(height-length(r)-c(1)+1,1)];
    c=[zeros(c(1)-1,1);c;zeros(height-length(c)-c(1)+1,1)];
    m1=m1';c=c'; n1=n1';
    
    Q=cat(3,m1,n1);Q=cat(3,Q,c);
    
    phi(k,:,:)=Q;
end
save('phi_3d_400up_30000.mat','phi');