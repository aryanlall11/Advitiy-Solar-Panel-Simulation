%Constants
iter=40000;
N_Cell = 2;
Vt = 25.7e-3; % Thermal voltage at temperature 25 degree C.

% Values at STC 
V_OC_STC = 4.66; % in Volts
I_SC_STC = 0.517; % in Amperes
I_MP_STC = 0.5; % in Amperes
V_MP_STC = 4.1; % in Volts
m_p=[];

%% Initial Guesses 

n_initial = (2*V_MP_STC - V_OC_STC)/(N_Cell*Vt*(log((I_SC_STC-I_MP_STC)/I_SC_STC) + (I_MP_STC/(I_SC_STC-I_MP_STC))));
Rs_initial = (V_MP_STC/I_MP_STC) - ((2*V_MP_STC - V_OC_STC)/(I_SC_STC-I_MP_STC))/(log((I_SC_STC-I_MP_STC)/I_SC_STC) + (I_MP_STC/(I_SC_STC-I_MP_STC))); 
Rsh_initial = sqrt(Rs_initial/((I_SC_STC/(n_initial*N_Cell*Vt))*(exp((Rs_initial*I_SC_STC - V_OC_STC)/(n_initial*N_Cell*Vt)))));


%% Solving 3 non linear equations 

initial_guess = zeros(1,3);
initial_guess(1) = Rs_initial;
initial_guess(2) = Rsh_initial;
initial_guess(3) = n_initial;

SDMeqn = @SDM_3eqn;

options=optimset('TolFun',1e-6,'MaxIter',100000,'MaxFunEvals',100000);

x = fsolve(SDMeqn,initial_guess,options);

%Solutions of the above equations
Rs_STC = x(1);
Rsh_STC = x(2);
n_STC = x(3);

%Calculating Isat and Iph from the above results
Isat_STC = (I_SC_STC - (V_OC_STC - I_SC_STC*Rs_STC)/Rsh_STC)*(exp(-V_OC_STC/(n_STC*N_Cell*Vt)));

Iph_STC = Isat_STC*exp(V_OC_STC/(n_STC*N_Cell*Vt)) + V_OC_STC/Rsh_STC;

%% Plotting the IV curve at STC 

figure();
f = @(V,I) I - Iph_STC + Isat_STC*(exp((V+I*Rs_STC)/(n_STC*N_Cell*Vt))-1) + ((V+I*Rs_STC)/Rsh_STC);
ezplot(f, [0 5 0 1]);  %%x axis 0-10 y-axis 0-50
title('IV Curve of the Module at STC');
ylabel('Current (in A)');
xlabel('Voltage (in V)');

%% Calculating the five parameters at different values of irradiance and temperatures

%Constants
Temp_STC = 25+273;  % in Kelvin
Temp_NOTC = 44+273; % in Kelvin
G_STC = 1000; % W per m2
G_NOTC = 800; % W per m2
I_SC_NOTC = 0.45; % in Amperes
alpha = 5e-4;
q = 1.60217662e-19;
E_G = 1.1; % for silicon 
K_B = 1.38064852e-23;

K = log((I_SC_NOTC - alpha*(Temp_NOTC-Temp_STC))/I_SC_STC)/(log(G_NOTC/G_STC));

%Initializing variables
I_SC = zeros(iter,4);
Rs = zeros(iter,4);
Rsh = zeros(iter,4);
Iph = zeros(iter,4);
Isat = zeros(iter,4);
Voc = zeros(iter,4);
Vmp = zeros(iter,4);
Imp = zeros(iter,4);
Pmp = zeros(iter,4);
FF = zeros(iter,4);
G=zeros(iter,4);
figure();
I=[-1,-1,-1,-1];
power=[];
voltage=[];
current=[];
n=0;
temp = 273 + 50;
%Calculating the parameters at different values of irradiance and temperatures
for i = 1:iter
    G(i,1:4)=side(i,1:4);
    for j = 1:4
        
      if(G(i,j)~=I(j))
        I_SC(i,j) = ((G(i,j)/G_STC)^K)*(I_SC_STC + alpha*(temp-Temp_STC));
        Rs(i,j) = Rs_STC;
        if(G(i,j)~=0)
         Rsh(i,j) = (G_STC/G(i,j))*Rsh_STC;
        else
         Rsh(i,j) = (G_STC)*Rsh_STC; 
        end
        
        Iph(i,j) = I_SC(i,j)*(1 + (Rs(i,j)/Rsh(i,j)));
        Isat(i,j) = Isat_STC*(temp/Temp_STC)^3*exp(((q*E_G)/(n_STC*K_B))*(1/Temp_STC - 1/temp));
        Vt = K_B*temp/q; % Considering the variation of thermal voltage with temperature
        
        fun1 = @(V) Iph(i,j) - Isat(i,j)*(exp(V/(n_STC*N_Cell*Vt))-1) - V/Rsh(i,j);
        Voc(i,j) = fzero(fun1,4);
        
        n=1;
      else
        I_SC(i,j) = I_SC(i-1,j);
        Rs(i,j) = Rs(i-1,j);
        Rsh(i,j) = Rsh(i-1,j);
        Iph(i,j) = Iph(i-1,j);
        Isat(i,j) = Isat(i-1,j);
        Voc(i,j) = Voc(i-1,j);
        Vt = K_B*temp/q;
      end
    end
    if(n==1)
     [x,y,z] =Test(Iph(i,1:4),n_STC,N_Cell,Vt,Isat(1,1),Rs(1,1),G(i,1:4),I_SC(i,1:4));
     %[x,y,z] =allp(Iph(i,1:4),n_STC,N_Cell,Vt,Isat(1,1),Rs(1,1),G(i,1:4),I_SC(i,1:4));
     %[x,y,z] =s3p1(Iph(i,1:4),n_STC,N_Cell,Vt,Isat(1,1),Rs(1,1),G(i,1:4),I_SC(i,1:4),Voc(i,1:4));
     %[x,y,z] =p2sp2(Iph(i,1:4),n_STC,N_Cell,Vt,Isat(1,1),Rs(1,1),G(i,1:4),I_SC(i,1:4),Voc(i,1:4));
      m_p=[m_p;x];
     voltage=[voltage;y];
     current=[current;z];
    else
      m_p=[m_p;m_p(end)];
      voltage=[voltage;voltage(end)];
      current=[current;current(end)];
    end
    I=G(i,1:4);
    n=0;
end

%% IV Curves for different irradiances and temperatures
% temperature = 273 + 50;
%  for temp = 1:4
%     for i = 1:1500
%         Iph_1 = Iph(i,temp);
%         Isat_1 = Isat(i,temp);
%         Rs_1 = Rs(i,temp);
%         Rsh_1 = Rsh(i,temp);
%         Vt = K_B*temperature/q; % Considering the variation of thermal voltage with temperature
%         f = @(V,I) I - Iph_1 + Isat_1*(exp((V+I*Rs_1)/(n_STC*N_Cell*Vt))-1) + ((V+I*Rs_1)/Rsh_1);
%         ezplot(f, [0 5 0 1]);
%         hold on;
%     end
%     title(['IV Curve for different irradiance at ' num2str(temperature-273) ' Degree Celsius']);
%     ylabel('Current (in A)');
%     xlabel('Voltage (in V)');
%     %legend({'400','500','600','700','800','900','1000','1100','1200','1300'});
%     figure();
%  end
 
 
x=zeros(1,iter); 
y=m_p(1:iter);
k=0.1;
for i =1:iter
    x(1,i)=k;
    k=k+0.1;
end
xlabel ('Time(s)'); ylabel ( 'Power(W)'); %power plot
plot(x,y); 

y=voltage(1:iter);
xlabel ('Time(s)'); ylabel ( 'Volatge(V)'); %voltage plot
plot(x,y); 

y=current(1:iter);
xlabel ('Time(s)'); ylabel ( 'Current(A)'); %current plot
plot(x,y);

current=double(current);
m_p=double(m_p);
% 
% y=Pmp(1:20,2);
% xlabel ('Time(s)'); ylabel ( 'Power(W)');
% plot(x,y); 
% 
% y=Pmp(1:20,3);
% xlabel ('Time(s)'); ylabel ( 'Power(W)');
% plot(x,y); 
% 
% y=Pmp(1:20,4);
% xlabel ('Time(s)'); ylabel ( 'Power(W)');
% plot(x,y); 
 
 
%% Plotting contour plots for different irradiances and temperatures
% temper = 20:5:80;
% irr = 400:100:1300;
% 
% figure();
% surf(temper,irr,I_SC);
% view(2);
% title('Contour Plot of Short Circuit Current');
% ylabel('Irradiance in W per m2');
% xlabel('Temperature in degree Celsius');
% 
% figure();
% surf(temper,irr,Voc);
% view(2);
% title('Contour Plot of Open Circuit Voltage');
% ylabel('Irradiance in W per m2' );
% xlabel('Temperature in degree Celsius');
% 
% figure();
% surf(temper,irr,Pmp);
% view(2);
% title('Contour Plot of Maximum Power');
% ylabel('Irradiance in W per m2' );
% xlabel('Temperature in degree Celsius');
% 
% figure();
% surf(temper,irr,FF);
% view(2);
% title('Contour Plot of Fill Factor');
% ylabel('Irradiance in W per m2' );
% xlabel('Temperature in degree Celsius');
% 
% figure();
% surf(temper,irr,Rsh);
% view(2);
% title('Contour Plot of Shunt Resistance');
% ylabel('Irradiance in W per m2');
% xlabel('Temperature in degree Celsius');
% 
% figure();
% surf(temper,irr,Iph);
% view(2);
% title('Contour Plot of Iph');
% ylabel('Irradiance in W per m2');
% xlabel('Temperature in degree Celsius');
% 
% figure();
% surf(temper,irr,Isat);
% view(2);
% title('Contour Plot of Isat');
% ylabel('Irradiance in W per m2');
% xlabel('Temperature in degree Celsius');
