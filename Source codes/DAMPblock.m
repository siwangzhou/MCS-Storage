function x_hat = DAMPblock(y,iters,lblocksize,mblocksize,nbtimesize,ksize_row,ksize_col,ksize_time,denoiser,M_func,M_pinv)


%M_pinv=M_func';
denoi=@(noisy,sigma_hat) denoiseblock(noisy,sigma_hat,lblocksize,mblocksize,nbtimesize,ksize_row,ksize_col,ksize_time,denoiser);

n=ksize_row*ksize_col*ksize_time;
m=length(y);

z_t=y;
x_t=zeros(n,1);
%M_pinv=pinv(M_func);

for i=1:iters
    pseudo_data=M_pinv*z_t+x_t;
    sigma_hat=sqrt(1/m*sum(abs(z_t).^2));
    x_t=denoi(pseudo_data,sigma_hat);
    eta=randn(1,n);
    epsilon=max(pseudo_data)/1000+eps;
    div=eta*((denoi(pseudo_data+epsilon*eta',sigma_hat)-x_t)/epsilon);
    z_t=y-M_func*x_t+1/m.*z_t.*div;
   % if mod(i,10)==0
  %      fprintf('iters=%d\n',i);
 %   end
end
%x_hat=reshape(x_t,[width height]);
x_hat=x_t;
end