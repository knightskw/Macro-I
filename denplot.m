% Shai Knight-Winnig 2016
% DENPLOT.m: procedure to compute and plot nonparametric kernel density 
%            estimate of data vector */

function denplot(vname,tstring,argmin,argmax);


dat=evalin('base',vname);
%argmin=floor(min(dat));
%argmax=ceil(max(dat));
%fprintf(1,'%s',name1);

if (~exist('pos'));
pos=.75;
end;

x=argmin:(argmax-argmin)/1000:argmax; 
n=size(dat,1);
h=1.06/(n^(.2));
s=std(dat);
mu=mean(dat);
npden=zeros(size(x,2),1);
%To use Gaussian Kernel, uncomment below
% for i=1:size(x,2)
%       npden(i,1)=sum(exp((-((x(1,i)-dat).*(x(1,i)-dat)))/(2*h*h*s*s)))/(s*sqrt(2*pi)*h*n);
% end

% %To use Epanechnikov Kernel, uncomment below
for i=1:size(x,2)
   npden(i,1)=sum((abs((x(1,i)-dat)/(h))<=sqrt(5)).*0.75.*(1-(x(1,i)-dat).*(x(1,i)-dat)/(h*h*5))/(sqrt(5)*n*h));
end


plot(x,npden);
axis([argmin argmax 0 max(npden)]);
title(['Distribution of ' tstring]);
ylabel('Density');
xlabel('x');
text(argmin+pos*(argmax-argmin),max(npden)*.95,['Mean    ' num2str(mu)]);
text(argmin+pos*(argmax-argmin),max(npden)*.9,['Median  ' num2str(median(dat))]);
text(argmin+pos*(argmax-argmin),max(npden)*.85,['Minimum ' num2str(min(dat))]);
text(argmin+pos*(argmax-argmin),max(npden)*.8,['Maximum ' num2str(max(dat))]);
text(argmin+pos*(argmax-argmin),max(npden)*.75,['Std dev ' num2str(s)]);
text(argmin+pos*(argmax-argmin),max(npden)*.7,['N       ' num2str(n)]);

