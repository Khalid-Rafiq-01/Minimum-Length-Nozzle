%{
    This code is basically an extension of the nozzle geometry code but it
    adds control parameter of the flow regiem by controlling limits of the 
    back pressure at the end of the divergent duct of the CD nozzle or the 
    beginning of the test section. As usual the flow is considered compressible,
    invicid, two-dimensional and irrotational.
    Khalid Rafiq 
    Kashmir University
 %}
    

format short
clear all
clc
%%Conversions.
RTOD = 180/pi;
DTOR = pi/180;
%% Gas properties. 
gma= input('Enter value of Gamma : ');
%% Definition of Prandtl Mayers function, v.
Prandtl_Meyer = @(Me,gma) sqrt((gma+1)/(gma-1))*atan(sqrt((gma-1)/(gma+1)*(Me^2-1)))-atan(sqrt(Me^2-1)) ;
v = Prandtl_Meyer ;
%% Numerical inverse Prandtl Meyer number.
M_n = @(v) (1+0.7860*v^0.666+0.0321*v^1.333-0.0988*v^2)/(1-0.3883*v^0.666-0.1049*v^1.333) ;
%% Definition of Mach angle.
Mu_n = @(M_n) asind(1/M_n);
%% Design parameters.
Me =  input('Enter Exit Mach no. : ');
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
xlim([0 x(nth)+0.5]);
ylim([0 y(nth)+0.5]);
arearatio = y(nth)/D;


%% Using bisection method to find subsonic and supersonic mach numbers.
% primary purpose is to check various back pressure conditions at which
% shocks may occur.

a1 = (gma+1)/ 2 ;
b1 = (gma-1)/ 2 ;
Msub1 = zeros(1,length(arearatio));
 
f = @(Ms) (1./Ms).*((1/a1 *( 1+b1*Ms.^2)).^(a1/(2*b1))) - arearatio ;
a = 0.0000001;
b = 1;
e = 10^-5;
n = ceil((log(b-a)-log(e))/log(2));
Ms = zeros(length(n)) ;
 for i = 1:n 
 Ms(i) = (a+b)/2 ;
 if   f(a)*f(Ms(i))< 0
     b = Ms(i);
 else if   f(a)*f(Ms(i))> 0
         a = Ms(i);
     end
 
 Ms;
 end
 Msub = Ms(1,end);
 end

 Msup = M(1,end) ;

 %% Variation of pressure ratio supersonic (p/p.stag) along each mach line.
  prsup = (1 + b1*Msup.^2).^(-gma/(2*b1)) ;
  
  %% Variation of pressure ratio subsonic (p/p.stag) along each mach line.
  prsub = ((1 + b1*Msub.^2).^(-gma/(2*b1)))' ;
  
  %% If the normal shock is formed at the exit of the supersonic nozzle, it must follow 
   % the nozrmal shock isentropic relations.
   p2to1 = 1 + (2*gma/(gma+1))*(Msup.^2 - 1) ;
   
   % now calculating ratio of Pexit(downstream normal shock) to the stagnation pressure
   prns = p2to1.*prsup ;
   disp([prsub',prns',prsup'])
  %% Control of the supersonic flow and conditions of shock formation.
  fprintf('Flow will be fully subsonic if ratio of backpressure to stagnation pressure is above %d\n',prsub )
  fprintf('Normal shock developed in divergent duct if backpressure to stagnation pressure ratio is between %d and %d\n',prsub,prns )
  fprintf('Normal shock developed in divergent duct at exit plane if backpressure to stagnation pressure ratio is %d\n',prns )
  fprintf('An over expanded nozzle is obtained if backpressure to stagnation pressure is between %d and %d\n',prns,prsup)
  fprintf('An under expanded nozzle is obtained if backpressure to the stagnation pressure ratio is below %d\n',prsup)
 


 
 