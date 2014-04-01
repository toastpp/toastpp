function  readcmfidat

fid = fopen('threelayer20iso.int','r')
[d,n] = fscanf(fid,'%5d %e');
for i = 1:n/2
    t(i) = d(i*2-1);
    y0(i) = d(i*2);
end

fid = fopen('threelayer20iso.cin','r')
[d,n] = fscanf(fid,'%5d %e %e %e %e %e %e %e %e %e %e');
for i = 1:n/11
    t(i) = d(i*11-10);
    x1(i,1) = d(i*11-9);
    x1(i,2) = d(i*11-8);
    x2(i,1) = d(i*11-7);
    x2(i,2) = d(i*11-6);
    x3(i,1) = d(i*11-5);
    x3(i,2) = d(i*11-4);
    x4(i,1) = d(i*11-3);
    x4(i,2) = d(i*11-2);
    x5(i,1) = d(i*11-1);
    x5(i,2) = d(i*11);
end
y1=complex(x1(:,1),x1(:,2));
y2=complex(x2(:,1),x2(:,2));
y3=complex(x3(:,1),x3(:,2));
y4=complex(x4(:,1),x4(:,2));
y5=complex(x5(:,1),x5(:,2));


figure;
plot(t,log(abs(y0)),'g-');
hold on
plot(t,log(abs(y1)));
plot(t,log(abs(y2)));
plot(t,log(abs(y3)));
plot(t,log(abs(y4)));
plot(t,log(abs(y5)));
title('complex intensity')

figure
plot(t,angle(y1),'r-');
hold on
plot(t,angle(y2),'r-');
plot(t,angle(y3),'r-');
plot(t,angle(y4),'r-');
plot(t,angle(y5),'r-');
title('phase angle rad')



