x=[0.0200 	0.0400 	0.0600 	0.0800 	0.1000 	0.1200 	0.1400 	0.1600 	0.1800 	0.2000 	0.2200 	0.2400 	0.2600 	0.2800 	0.3000 	0.3200 	0.3400 	0.3600 	0.3800 	0.4000];
y1=[0.5017 	0.4084 	0.3548 	0.3045 	0.2733 	0.2435 	0.2197 	0.1995 	0.1799 	0.1658 	0.1503 	0.1370 	0.1255 	0.1153 	0.1059 	0.0981 	0.0907 	0.0850 	0.0789 	0.0743];
y2=[0.5050 	0.4362 	0.3810 	0.3321 	0.2979 	0.2731 	0.2498 	0.2299 	0.2126 	0.1930 	0.1765 	0.1624 	0.1490 	0.1370 	0.1256 	0.1162 	0.1059 	0.0988 	0.0913 	0.0861];
y3=[1.2195 	0.7691 	0.5317 	0.4074 	0.3493 	0.3069 	0.2976 	0.2898 	0.2777 	0.2731 	0.2598 	0.2602 	0.2463 	0.2420 	0.2324 	0.2244 	0.2081 	0.2112 	0.1985 	0.1866];
y4=[0.4746 	0.3368 	0.2990 	0.2837 	0.2646 	0.2488 	0.2289 	0.2097 	0.1914 	0.1729 	0.1561 	0.1383 	0.1241 	0.1123 	0.1011 	0.0919 	0.0844 	0.0784 	0.0729 	0.0692];
y5=[1.0436 	0.9308 	0.7336 	0.6312 	0.6112 	0.6043 	0.6046 	0.6020 	0.6011 	0.5990 	0.5985 	0.5988 	0.5975 	0.5975 	0.5970 	0.5972 	0.5957 	0.5958 	0.5963 	0.5950];

% plot(x,y1,'-p',x,y2,'o-b',x,y3,'*-r',x,y4,'*-k',x,y5,'^-g');
% legend('mcs','mcs-block','ST-CNC','cdp-st','global');

plot(x,y5,'-ob',x,y3,'-dk',x,y1,'-xm',x,y2,'-pr');
axis([0.03 0.37 0 0.8]);
legend('CNS-WSNs','ST-CNC-WSNs','The proposed (without PDR)','The proposed');
xlabel('Decoding ratio');
ylabel('The relative square error');