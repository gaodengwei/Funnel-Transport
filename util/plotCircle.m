function h=plotCircle(c,r)
theta = 0:pi/100:2*pi;
x=c(1)+r*cos(theta);
y=c(2)+r*sin(theta);
h=plot(x,y,'LineWidth',2,'LineStyle','--');

end