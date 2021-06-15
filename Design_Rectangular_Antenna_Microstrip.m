%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Design Rectangular Microstrip Antenna by Transmission Line Model %%%
%%%                  Luiz Guilherme Martinho, 15-2021                   %%%
%%%                   https://github.com/guiimartinho                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
format long

disp('Design Rectangular Microstrip Antenna by Transmission Line Model Using MATLAB');
disp('=============================================================================');
er=input('Enter the dielectric constant:');
h=input('Enter the substrate thickness (in mm):');
f=input('Enter the frequency (GHz):');
z=input('Enter the input impedance (ohm):');
disp('Calculating. Please wait.');

f=f*1e9;

%calculate the width
wid=(3e8/(sqrt((er+1)/2)*2*f))*1000;        %in mm

%calculate the effective dielectric constant
e_eff=((er+1)/2)+ (((er-1)/2)* (1+((12*h)/wid))^-0.5);

%calculate the extension of length L
del_l=(((e_eff+0.3)*((wid/h)+0.264))/((e_eff-0.258)*((wid/h)+0.8)))*(0.412*h);  %in mm

%calculate the effective length
l_eff=(3e8/(2*f*sqrt(e_eff)))*1000;

%calculate the actual length
L=l_eff-(2*del_l);


la=(3e8/f)*1000;
k=(2*pi)/la;
x=k*(wid);
i1=-2+cos(x)+(x*sinint(x))+(sin(x)/x);
g1=i1/(120*pi*pi);          %Conductance

a=@(th)(((sin((x./2).*cos(th))./cos(th)).^2).*(besselj(0,(k.*L.*sin(th)))).*(sin(th)).^3);
a1=integral(a,0,pi);
g12=a1/(120*pi*pi);     %in siemens
r_in=1/(2*(g1+g12));    %in ohms

inset=(L/pi)*(acos(sqrt(z/r_in)));        %in mm

Lg_min=6*h+L;   %based on Punit S. Nakar (2004) book
Wg_min=6*h+wid;

B=60*pi*pi/(z*sqrt(er)); 
m1=2*B-1;
m=log(m1);
n1=B-1;
n=log(n1);

W=(2*h/pi)*(B-1-m+(((er-1)/(2*er))*(n+(0.39*0.61)/er)));   %based on Liang J. (2004) book 

g = (3e8*4.65e-9)/(sqrt(2*e_eff)*f*10^-9); %based on Journak by M A Matin (2010)

disp(' ');
disp('RECTANGULAR PATCH:')
    disp(['The width of the patch (Wp) is: ',num2str(wid),' mm'])
    disp(['The length of the patch (Lp) is: ',num2str(L),' mm'])
    disp(['The inset feed point (Fi) is: ',num2str(inset),' mm'])
    disp(['The width of the feed line (Wf) is: ',num2str(W),' mm'])
    disp(['The gap of the feed line (Gpf) is: ',num2str(g),' mm'])
    disp(['The minimum length of ground plane is: ',num2str(Lg_min),' mm'])
    disp(['The minimum width of ground plane is: ',num2str(Wg_min),' mm'])
