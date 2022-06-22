%{
 This MATLAB code presents Method of characteristics for two-dimensional,supersonic nozzle design.
 Minimum length of the supersonic nozzle has been calculated for the optimum Mach
 number at the nozzle exit with uniform flow at the diverging section of the nozzle
 Numerical solution is established for the two dimensional, steady, in viscid, irrotational
 and supersonic flow.  The design considerations are concentrated at
 the diverging section.
 Code developed by Khalid Rafiq.
 Kashmir University.
%}
format short
clear all
clc
%%Conversions.
RTOD = 180/pi;
DTOR = pi/180;
%% Gas properties. 
gma=input('Enter value of Gamma : ');
%% Definition of Prandtl Mayers function, v.
Prandtl_Meyer = @(Me,gma) sqrt((gma+1)/(gma-1))*atan(sqrt((gma-1)/(gma+1)*(Me^2-1)))-atan(sqrt(Me^2-1)) ;
v = Prandtl_Meyer ;
%% Numerical inverse Prandtl Meyer number.
M_n = @(v) (1+0.7860*v^0.666+0.0321*v^1.333-0.0988*v^2)/(1-0.3883*v^0.666-0.1049*v^1.333) ;
%% Definition of Mach angle.
Mu_n = @(M_n) asind(1/M_n);
%% Design parameters.
Me = input('Enter Exit Mach no. : ');
%% Max Nozzle wall angle.
theta_max= 0.5 * Prandtl_Meyer(Me,gma)*RTOD ;
%% Input characteristic lines %%
n_c = input('Enter number of characteristics lines originating from throat : ');
%% Pre- allocation of parameters.
nth= 0.5*n_c*(n_c + 3); % total number of nodes in the mesh.
M=zeros(1,nth);        % pre-allocation of Mach numbers.
mu=zeros(1,nth);        % pre-allocation of Mach angles
v=zeros(1,nth);        % pre-allocation of PMF.
LRC=zeros(1,nth);      % pre-allocation of Left running characteristics.
RRC=zeros(1,nth);      % pre-allocation of Right running characteristics.
theta=zeros(1,nth);    % pre-allocation of flow angles.
%% Choosing the initial step angle.Starts Numeric computation.
if n_c > 10
     theta_0=  theta_max/n_c ;
else
    theta_0 = 0.3733 ; % Any value close to 0 would suffice.
end
%% Generating the flow parametres at each grid point.
nth= 0.5*n_c*(n_c + 3); % total number of nodes in the mesh.
dtheta=(theta_max-theta_0)/(n_c-1); % Ref. Anderson 11.1.
for i=1:n_c  % First LRC has n_c grid points in addition to wall point.
    theta(i)=theta_0+(i-1)*dtheta;
    v(i)=theta(i);
    LRC(i)=theta(i)-v(i); % theta-v=constant.
    RRC(i)=theta(i)+v(i); % theta+v=constant.
end
i=n_c+1; % The wall point.
theta(i)=theta(i-1); % node at wall has same flow parameter as one prior to it.
v(i)=v(i-1);         % node at wall has same flow parameter as one prior to it.
LRC(i)=LRC(i-1);     % node at wall has same flow parameter as one prior to it.
RRC(i)=RRC(i-1);     % node at wall has same flow parameter as one prior to it.
a=2;   
b=n_c+2;
for k=1:n_c-1    % Jumps to the second LRC,2nd point on centrline.
    j=a;
    h=b;
    RRC(h)=RRC(j);
    theta(h)=0;
    v(h)=RRC(j)-theta(h);
    LRC(h)=theta(h)-v(h);
    j=j+1;
    for i=h+1:n_c-a+b  % second LRC,points in bet C.L and wall.
        RRC(i)=RRC(j);
        LRC(i)=LRC(i-1);
        theta(i)=0.5*(LRC(i)+RRC(i)); % theta = 1/2 *((S+) + S(-))
        v(i)=0.5*(RRC(i)-LRC(i));     % theta = 1/2 *((S+) + S(-))
        j=j+1;
    end
    if i==n_c-a+b
        h=i+1;
    else
        h=h+1;
    end
    theta(h)=theta(h-1);   % node at wall has same flow parameter as one prior to it.
    v(h)=v(h-1);           % node at wall has same flow parameter as one prior to it.
    LRC(h)=LRC(h-1);       % node at wall has same flow parameter as one prior to it.
    RRC(h)=RRC(h-1);       % node at wall has same flow parameter as one prior to it.
    a=a+1;
    b=h+1;
end
%% Mach number and Mach angle at each node
v=v*DTOR;
Mu = zeros(1,nth);     % pre-allocation of variable.
for i=1:nth
    M(i)=M_n(v(i));    % Respective Mach number at each node.
    Mu(i)=Mu_n(M(i));  % Respective Mach angle at each node
end
           % Here we complete out our flow paramaters at all the grid
           % points taking each at a time. However nothing is yet known
           %about the shape of the nozzle.
%% Creating the table of the flow parameters at the Nodes.
pts = 1:nth ;
disp('Flow Parameters at the Nodes : '),disp('    Point      RRC(K_)   LRC(K+)    Theta     New(v)     M      Mu   ');
disp([pts',RRC',LRC',theta',v',M',Mu'])
%% Grid generator
figure(1)
D = input('Enter height of throat: ') ; 
i=1;
x=zeros(1,nth);
y=zeros(1,nth);
wall=theta_max;
while (i<=n_c+1)
    if i==1
        x(i)=-D/(tand(theta(i)-Mu(i)));
        y(i)=0;
        plot([0 x(i)],[D 0]);
    

        hold on
    else if i==n_c+1
            x(i)=(y(i-1)-D-x(i-1)*tand((theta(i-1)+theta(i)+Mu(i-1)+Mu(i))*0.5))/...
                (tand(0.5*(wall+theta(i)))-tand((theta(i-1)+theta(i)+Mu(i-1)+Mu(i))*0.5));
            y(i)=D+x(i)*tand(0.5*(wall+theta(i)));
            plot([x(i-1) x(i)],[y(i-1) y(i)]);
            hold on
            plot([0 x(i)],[D y(i)],'r','Linewidth',1.5);
           
            hold on
        else
            x(i)=(D-y(i-1)+x(i-1)*tand(0.5*(Mu(i-1)+theta(i-1)+Mu(i)+theta(i))))/...
                (tand(0.5*(Mu(i-1)+theta(i-1)+Mu(i)+theta(i)))-tand(theta(i)-Mu(i)));
            y(i)=tand(theta(i)-Mu(i))*x(i)+D;
            plot([x(i-1) x(i)],[y(i-1) y(i)]);
            hold on
            plot([0 x(i)],[D y(i)]);
            hold  on
        end
    end
    i=i+1;
    hold on
end
h=i;
k=0;
i=h;
for j=1:n_c-1
    while (i<=h+n_c-k-1)
        if (i==h)
            x(i)=x(i-n_c+k)-y(i-n_c+k)/(tand(0.5*(theta(i-n_c+k)+theta(i)-Mu(i-n_c+k)-Mu(i))));
            y(i)=0;
            plot([x(i-n_c+k) x(i)],[y(i-n_c+k) y(i)]);
            hold on
        else if (i==h+n_c-k-1)
                x(i)=(x(i-n_c+k)*tand(0.5*(theta(i-n_c+k)+theta(i)))-y(i-n_c+k)+y(i-1)-x(i-1)...
                    *tand((theta(i-1)+theta(i)+Mu(i-1)+Mu(i))*0.5))/(tand(0.5*(theta(i-n_c+k)...
                    +theta(i)))-tand((theta(i-1)+theta(i)+Mu(i-1)+Mu(i))*0.5));
                y(i)=y(i-n_c+k)+(x(i)-x(i-n_c+k))*tand(0.5*(theta(i-n_c+k)+theta(i)));
                plot([x(i-1) x(i)],[y(i-1) y(i)]);
                hold on
                plot([x(i-n_c+k) x(i)],[y(i-n_c+k) y(i)],'r','Linewidth',1.5);
                hold on
            else
                s1= tand(0.5*(theta(i)+theta(i-1)+Mu(i)+Mu(i-1)));
                s2= tand(0.5*(theta(i)+theta(i-n_c+k)-Mu(i)-Mu(i-n_c+k)));
                x(i)=(y(i-n_c+k)-y(i-1)+s1*x(i-1)-s2*x(i-n_c+k))/(s1-s2);
                y(i)=y(i-1)+(x(i)-x(i-1))*s1;
                plot([x(i-1) x(i)],[y(i-1) y(i)]);
                hold on
               h5 =  plot([x(i-n_c+k) x(i)],[y(i-n_c+k) y(i)]);
                hold on
            end
        end
        i=i+1;
    end
    k=k+1;
    h=i;
    i=h;
    hold on
end
h1 = get(h5,'Parent');
set(h1,'Fontsize',14,'Linewidth',2);
title(sprintf('Characteristic lines=%d for Mach=%d and Cp/Cv=%d',n_c,Me,gma))
xlabel('Length of Expansion duct :');
ylabel('Height of duct :','Rotation',90);
axis equal
xlim([0 x(nth)+0.5])
ylim([0 y(nth)+0.5])
arearatio = y(nth)/D
hold on
%% Variation of mach. no along each char line.
x = x' ;  y = y' ;  m = M' ;
nx = 500 ; ny = 500 ;
[X,Y] = meshgrid(linspace(min(x),max(x),nx),linspace(min(y),max(y),ny)) ;
M =griddata(x,y,m,X,Y) ;
figure(2)
contourf(X,Y,M,200,'LineStyle','none') ;

xlim([0 x(nth)+0.5])
ylim([0 y(nth)+0.5])
title('Variation of mach. no along each meshline: ');
xlabel('Length of Expansion duct :');
ylabel('Height of duct :','Rotation',90);
axis equal
colorbar
hold on 
%% Variation of pressure(p/p.stag) along each mach line.

a = (gma+1)/ 2 ;
b = (gma-1)/ 2 ;
pr = (1 + b*m.^2).^(-gma/2*b) ;
[A,B] = meshgrid(linspace(min(x),max(x),nx),linspace(min(y),max(y),ny)) ;
Pr =griddata(x,y,pr,A,B) ;

figure(3)
contourf(A,B,Pr,200,'LineStyle','none') ;
xlim([0 x(nth)+0.5])
ylim([0 y(nth)+0.5])
title('Variation of pressure ratio along each meshline: ');
xlabel('Length of Expansion duct :');
ylabel('Height of duct :','Rotation',90);
axis equal
colorbar
hold on 

%% %% Variation of pressure(T/T.stag) along each mach line.

tr = (1 + b*m.^2).^(-1) ;
[C,D] = meshgrid(linspace(min(x),max(x),nx),linspace(min(y),max(y),ny)) ;
Tr =griddata(x,y,tr,C,D) ;
figure(4)
contourf(A,B,Tr,200,'LineStyle','none') ;
xlim([0 x(nth)+0.5])
ylim([0 y(nth)+0.5])
title('Variation of temperature ratio along each meshline: ');
xlabel('Length of Expansion duct :');
ylabel('Height of duct :','Rotation',90);
axis equal
colorbar
hold on 

%% %% Variation of pressure(T/T.stag) along each mach line.

rho = (pr).^(gma) ;
[E,F] = meshgrid(linspace(min(x),max(x),nx),linspace(min(y),max(y),ny)) ;
Rho =griddata(x,y,rho,E,F) ;
figure(5)
contourf(A,B,Rho,200,'LineStyle','none') ;
xlim([0 x(nth)+0.5])
ylim([0 y(nth)+0.5])
title('Variation of density ratio along each meshline: ');
xlabel('Length of Expansion duct :');
ylabel('Height of duct :','Rotation',90);
axis equal
colorbar