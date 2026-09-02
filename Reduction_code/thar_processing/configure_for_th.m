function [ufix,vfix] = configure_for_th(data,uc,vc,blue_data_chop_value)

image(data,'CDataMapping','scaled');colormap(gray);caxis([0 500])
hold on
plot(uc,vc-blue_data_chop_value,'r+')
disp('you will be asked to mark the difference between the thoriums and ')
disp('their expected position for 9 different parts of the spectrum ...')
disp(' ')
sz=size(data);
cnt=0;
for s1=1:3
    for s2=1:3
cnt=cnt+1;
axis([s1*sz(2)/4-250 s1*sz(2)/4+250 s2*sz(1)/4-250 s2*sz(1)/4+250])
disp('click first on + then on ThAr line')
[x1,y1]=ginput(2);
x(cnt)=x1(1);
y(cnt)=y1(1);
zx(cnt)=diff(x1);
zy(cnt)=diff(y1);
    end
end
p1=polynomial_fit_2d(x,y,zx,2,2);
p2=polynomial_fit_2d(x,y,zy,2,2);
ufix=polynomial_val_2d(p1,uc,vc);
vfix=polynomial_val_2d(p2,uc,vc);
ucfixed=uc+ufix;
vcfixed=vc+vfix;
plot(ucfixed,vcfixed-blue_data_chop_value,'ro')
