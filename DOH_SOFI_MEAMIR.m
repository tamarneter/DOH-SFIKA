clear all;
close all;
clc;

%כיול מתמר הלחץ
rho = 997.1; %kg/m^3
rho_err=0.36;
g = 9.81; %m/s^2
height = [1.09, 1.30, 1.54, 1.68, 1.90, 2.16]; %m
height_err=0.02;
delt_p = rho*g*height; %pa
delt_p_err=sqrt((g*height*rho_err).^2+(rho*g*height_err)^2);
volt = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1];
volt_err=ones(1,6)*0.05;

f = @(x) (x+0.8812)/(9.476*10^-5);
figure(1);
hold on;
fplot(f);
grid on;
errorbar(volt,delt_p,delt_p_err,delt_p_err,volt_err,volt_err,'o');
ylim([10000,22000]);
xlim([0,1.2]);
%title('pressure to voltage');
ylabel('pressure [Pa]');
xlabel('voltage [Volt]');
legend('estimated curve', 'mesured curve');
hold off;

diameter = 0.19;%m
A=pi*diameter^2/4;

%חריר
H_mesured_harir = [0.2,0.4,0.2,0.4,0.4];%m
H_mesured_harir_err=0.03;
time_harir = [50.67,53.12,22.55,27.67,22.75]; %sec
time_harir_err=1;
for i=[1:1:5]
   Q_harir(i) = H_mesured_harir(i)*A/time_harir(i);
   Q_harir_err(i)=sqrt((A/time_harir(i)*H_mesured_harir_err)^2+(A*H_mesured_harir(i)/(time_harir(i)^2)*time_harir_err)^2);
end

%נחיר
H_mesured_nozzle = [0.1,0.2,0.3,0.3,0.3]; %m
H_mesured_nozzle_err=0.03;
time_nozzle = [29.15,26.17,26.57,21.07,16.66]; %sec
time_nozzle_err=1;
for i=[1:1:5]
   Q_nozzle(i) = H_mesured_nozzle(i)*A/time_nozzle(i);
   Q_nozzle_err(i)=sqrt((A/time_nozzle(i)*H_mesured_nozzle_err)^2+(A*H_mesured_nozzle(i)/(time_nozzle(i)^2)*time_nozzle_err)^2);
end

%ונטורי
H_mesured_venturi = [0.1,0.2,0.3,0.3,0.4]; %m
H_mesured_venturi_err=0.03;
time_venturi = [25.44,27.11,27.00,20.61,22.13]; %sec
time_venturi_err=1;
for i=[1:1:5]
   Q_venturi(i) = H_mesured_venturi(i)*A/time_venturi(i);
   Q_venturi_err(i)=sqrt((A/time_venturi(i)*H_mesured_venturi_err)^2+(A*H_mesured_venturi(i)/(time_venturi(i)^2)*time_venturi_err)^2);
end


R = @(x) 0.9696*x+6.181*10^-6;
rotmeter_mesurment = [1,2,3,4,5]*10^-4; %m^3/s
rotmeter_mesurment_err=(1/3)*10^(-4);
figure(2);
hold on;
grid on;
fplot(R);
% scatter(rotmeter_mesurment,Q_harir,'filled');
errorbar(rotmeter_mesurment,Q_harir,Q_harir_err,Q_harir_err,ones(1,5)*rotmeter_mesurment_err,ones(1,5)*rotmeter_mesurment_err,'o');
xlim([0*10^-4,6*10^-4]);
ylim([0*10^-4,6*10^-4]);
%title('real flow rate to displayed flow rate');
xlabel('displayed flow rate [m^3/sec]');
ylabel('real flow rate [m^3/sec]');
legend('estimated curve', 'mesured curve');
hold off;

bob_height = [0.045,0.072,0.078,0.107,0.123];
bob_height_err=0.005;
R2 = @(x) 0.00508*x-0.0001329;
figure(3);
hold on;
grid on;
fplot(R2);
%scatter(bob_height,Q_harir,'filled');
errorbar(bob_height,Q_harir,Q_harir_err,Q_harir_err,ones(1,5)*bob_height_err,ones(1,5)*bob_height_err,'o');
xlim([0.03,0.13]);
ylim([0.5*10^-4,6*10^-4]);
%title('bob height to displayed flow rate');
xlabel('bob height [m]');
ylabel('real flow rate [m^3/sec]');
legend('estimated curve', 'mesured curve');
hold off;

%חיישן מערבולות
MTM = @(x) (1.253*10^-5)*x - 6.305*10^-6;
w = [9.6,17.6,20.31,33,40.5]; %Hz
w_err=0.5;
figure(4);
hold on;
grid on;
fplot(MTM);
%scatter(w,Q_harir,'filled');
errorbar(w,Q_harir,Q_harir_err,Q_harir_err,w_err*ones(1,5),w_err*ones(1,5),'o');
xlim([5,45]);
ylim([0.5*10^-4,6*10^-4]);
%title('real flow rate to frequency');
xlabel('frequency [Hz]');
ylabel('real flow rate [m^3/sec]');
legend('estimated curve', 'mesured curve');
hold off;

harir_p1=([0.7,1.0,1.1,1.9,2.5]+0.8812)/(9.476*10^-5);
harir_p1_err=1.055*(10^4)*0.05*ones(1,5);
harir_p2=([0.7,0.9,1.0,1.6,2.0]+0.8812)/(9.476*10^-5);
nozzle_p1=([0.7,0.8,1.1,1.3,1.7]+0.8812)/(9.476*10^-5);
nozzle_p1_err=1.055*(10^4)*0.05*ones(1,5);
nozzle_p2=([0.7,0.8,0.9,1.1,1.4]+0.8812)/(9.476*10^-5);
venturi_p1=([0.8,1.0,1.2,1.6,2.0]+0.8812)/(9.476*10^-5);
venturi_p1_err=1.055*(10^4)*0.05*ones(1,5);
venturi_p2=([0.7,0.7,0.7,0.8,0.8]+0.8812)/(9.476*10^-5);

%6.2.3
D1 = 24*10^-3; %m
D2 = 11.5*10^-3; %m
A1 = pi*((D1)^2)/4; %m^2
A2 = pi*((D2)^2)/4; %m^2
Q_theo_harir1 = (A2/sqrt(1-(A2/A1)^2))*sqrt(2*harir_p1/rho);
Q_theo_harir1_err=A2./(sqrt((1-(A2/A1))^2*2*rho*harir_p1)).*harir_p1_err;
Q_theo_harir2 = (A2/sqrt(1-(A2/A1)^2))*sqrt(2*harir_p2/rho);
Q_theo_nozzle1 = (A2/sqrt(1-(A2/A1)^2))*sqrt(2*nozzle_p1/rho);
Q_theo_nozzle1_err=A2./(sqrt((1-(A2/A1))^2*2*rho*harir_p1)).*harir_p1_err;
Q_theo_nozzle2 = (A2/sqrt(1-(A2/A1)^2))*sqrt(2*nozzle_p2/rho);
Q_theo_venturi1 = (A2/sqrt(1-(A2/A1)^2))*sqrt(2*venturi_p1/rho);
Q_theo_venturi1_err=A2./(sqrt((1-(A2/A1))^2*2*rho*harir_p1)).*harir_p1_err;
Q_theo_venturi2 = (A2/sqrt(1-(A2/A1)^2))*sqrt(2*venturi_p2/rho);

miu = 0.9*10^-3; %pa*s
%מספרי ריינולדסזזזזז
RE_harir1 = (rho*(Q_harir/A1)*D1)/miu;
RE_nozzle1 = (rho*(Q_nozzle/A1)*D1)/miu;
RE_venturi2 = (rho*(Q_venturi/A2)*D2)/miu;

%מקדם פריקה חריר
beta_harir = 11.5/24;
MC_harir = 0.625;
C_harir = MC_harir*sqrt(1-(A2/A1)^2);
%מקדם פריקה ונטורי
C_venturi = [0.96,0.97,0.97,0.97,0.98];
%מקדם פריקה נחיר
C_nozzle = [0.94,0.95,0.96,0.96,0.96];

%תיקון ספיקה לפי C
Q_harir_fixed_c = Q_theo_harir1*C_harir;
Q_harir_fixed_c_err=Q_theo_harir1_err*C_harir;
for i=[1:1:5]
    Q_nozzle_fixed_c(i) = Q_theo_nozzle1(i)*C_nozzle(i);
    Q_nozzle_fixed_c_err(i)=Q_theo_nozzle1_err(i).*C_nozzle(i);
end
for i=[1:1:5]
    Q_venturi_fixed_c(i) = Q_theo_venturi1(i)*C_venturi(i);
    Q_venturi_fixed_c_err(i)=Q_theo_venturi1_err(i).*C_nozzle(i);
end

%תיקון ספיקה לפי שטות אחרת
L_venturi = 58*10^-3 + 263.5*10^-3; %m
L_nozzle = 47*10^-3 + 228.5*10^-3; %m 
L_harir = 36*10^-3 + 236.5*10^-3; %m
l_venturi = 58*10^-3; %m
l_nozzle = 47*10^-3; %m
l_harir = 36*10^-3; %m
%הפרשי לחץ שטות אחרת
P_friction_venturi = (l_venturi/L_venturi)*venturi_p2; %pa
P_friction_nozzle = (l_nozzle/L_nozzle)*nozzle_p2; %pa
P_friction_harir = (l_harir/L_harir)*harir_p2; %pa

Q_venturi_fixed_fric_1 = (A2/sqrt(1-(A2/A1)^2))*sqrt(2*(venturi_p1-P_friction_venturi)/rho);
Q_venturi_fixed_fric_1_err=Q_theo_venturi1_err.*sqrt(venturi_p1./(venturi_p1-P_friction_venturi));
Q_nozzle_fixed_fric_1 = (A2/sqrt(1-(A2/A1)^2))*sqrt(2*(nozzle_p1-P_friction_nozzle)/rho);
Q_nozzle_fixed_fric_1_err=Q_theo_nozzle1_err.*sqrt(nozzle_p1./(nozzle_p1-P_friction_nozzle));
Q_harir_fixed_fric_1 = (A2/sqrt(1-(A2/A1)^2))*sqrt(2*(harir_p1-P_friction_harir)/rho);
Q_harir_fixed_fric_1_err=Q_theo_harir1_err.*sqrt(harir_p1./(harir_p1-P_friction_harir));

%alla graphim 1
%     %venturi
figure(5);
hold on;
grid on;
f_venturi_1 = fit(Q_venturi',Q_theo_venturi1','poly1');
f_venturi_2 = fit(Q_venturi',Q_venturi_fixed_c','poly1');
f_venturi_3 = fit(Q_venturi',Q_venturi_fixed_fric_1','poly1');
%scatter(Q_venturi,Q_theo_venturi1,'filled');
errorbar(Q_venturi,Q_theo_venturi1,Q_theo_venturi1_err,Q_theo_venturi1_err,Q_venturi_err,Q_venturi_err,'o');
plot(f_venturi_1,'g');
%scatter(Q_venturi,Q_venturi_fixed_c,'filled');
errorbar(Q_venturi,Q_venturi_fixed_c,Q_venturi_fixed_c_err,Q_venturi_fixed_c_err,Q_venturi_err,Q_venturi_err,'o');
plot(f_venturi_2,'b');
%scatter(Q_venturi,Q_venturi_fixed_fric_1,'filled');
errorbar(Q_venturi,Q_venturi_fixed_fric_1,Q_venturi_fixed_fric_1_err,Q_venturi_fixed_fric_1_err,Q_venturi_err,Q_venturi_err,'o');
plot(f_venturi_3,'c');
legend('ideal flow','idael flow linear fit','discharge coefficient compensation','discharge coefficient compensation linear fit','pressure compensation','pressure compensation linear fit');
ylabel('calculated flow [m^3/s]');
xlabel('real flow [m^3/s]');
title('venturi- calculated flow to real flow');
hold off;

    %harir
figure(6);
hold on;
grid on;
f_harir_1 = fit(Q_harir',Q_theo_harir1','poly1');
f_harir_2 = fit(Q_harir',Q_harir_fixed_c','poly1');
f_harir_3 = fit(Q_harir',Q_harir_fixed_fric_1','poly1');
% scatter(Q_harir,Q_theo_harir1,'filled');
errorbar(Q_harir,Q_theo_harir1,Q_theo_harir1_err,Q_theo_harir1_err,Q_harir_err,Q_harir_err,'o');
plot(f_harir_1,'g');
% scatter(Q_harir,Q_harir_fixed_c,'filled');
errorbar(Q_harir,Q_harir_fixed_c,Q_harir_fixed_c_err,Q_harir_fixed_c_err,Q_harir_err,Q_harir_err,'o');
plot(f_harir_2,'b');
% scatter(Q_harir,Q_harir_fixed_fric_1,'filled');
errorbar(Q_harir,Q_harir_fixed_fric_1,Q_harir_fixed_fric_1_err,Q_harir_fixed_fric_1_err,Q_harir_err,Q_harir_err,'o');
plot(f_harir_3,'c');
title('harir- calculated flow to real flow');
ylabel('calculated flow [m^3/s]');
xlabel('real flow [m^3/s]');
legend('ideal flow','idael flow linear fit','discharge coefficient compensation','discharge coefficient compensation linear fit','pressure compensation','pressure compensation linear fit');
hold off;

    %nozzle
figure(7);
hold on;
grid on;
f_nozzle_1 = fit(Q_nozzle',Q_theo_nozzle1','poly1');
f_nozzle_2 = fit(Q_nozzle',Q_nozzle_fixed_c','poly1');
f_nozzle_3 = fit(Q_nozzle',Q_nozzle_fixed_fric_1','poly1');
% scatter(Q_nozzle,Q_theo_nozzle1,'filled');
errorbar(Q_nozzle,Q_theo_nozzle1,Q_theo_nozzle1_err,Q_theo_nozzle1_err,Q_nozzle_err,Q_nozzle_err,'o');
plot(f_nozzle_1,'g');
% scatter(Q_nozzle,Q_nozzle_fixed_c,'filled');
errorbar(Q_nozzle,Q_nozzle_fixed_c,Q_nozzle_fixed_c_err,Q_nozzle_fixed_c_err,Q_nozzle_err,Q_nozzle_err,'o');
plot(f_nozzle_2,'b');
% scatter(Q_nozzle,Q_nozzle_fixed_fric_1,'filled');
errorbar(Q_nozzle,Q_nozzle_fixed_fric_1,Q_nozzle_fixed_fric_1_err,Q_nozzle_fixed_fric_1_err,Q_nozzle_err,Q_nozzle_err,'o');
plot(f_nozzle_3,'c');
title('nozzle- calculated flow to real flow');
ylabel('calculated flow [m^3/s]');
xlabel('real flow [m^3/s]');
legend('ideal flow','idael flow linear fit','discharge coefficient compensation','discharge coefficient compensation linear fit','pressure compensation','pressure compensation linear fit');
hold off;

%graphei kiool
p1_venturi_sqr = sqrt(venturi_p1);
p1_venturi_sqr_err=venturi_p1_err./(2*p1_venturi_sqr);
p1_nozzle_sqr = sqrt(nozzle_p1);
p1_nozzle_sqr_err=nozzle_p1_err./(2*p1_nozzle_sqr);
p1_harir_sqr = sqrt(harir_p1);
p1_harir_sqr_err=harir_p1_err./(2*p1_harir_sqr);
    %venturi
figure(8);
hold on;
grid on;
%scatter(p1_venturi_sqr,Q_venturi,'filled');
errorbar(p1_venturi_sqr,Q_venturi,Q_venturi_err,Q_venturi_err,p1_venturi_sqr_err,p1_venturi_sqr_err,'o');
f_venturi_kiool = fit(p1_venturi_sqr',Q_venturi','poly1');
plot(f_venturi_kiool);
legend('measured results','measured results linear fit');
title('venturi- real flow to sqare delta(P)');
xlabel('(p1-p2)^(0.5) [Pa^(0.5)]');
ylabel('real flow [m^3/s]');
xlim([130,180]);
ylim([0.5,6]*10^-4);
hold off;

    %nozzle
figure(9);
hold on;
grid on;
%scatter(p1_nozzle_sqr,Q_nozzle,'filled');
errorbar(p1_nozzle_sqr,Q_nozzle,Q_nozzle_err,Q_nozzle_err,p1_nozzle_sqr_err,p1_nozzle_sqr_err,'o');
f_nozzle_kiool = fit(p1_nozzle_sqr',Q_nozzle','poly1');
plot(f_nozzle_kiool);
legend('measured results','measured results linear fit');
title('nozzle- real flow to sqare delta(P)');
xlabel('(p1-p2)^(0.5) [Pa^(0.5)]');
ylabel('real flow [m^3/s]');
xlim([125,170]);
ylim([5*10^-5,6*10^-4]);
hold off;

    %harir
figure(10);
hold on;
grid on;
%scatter(p1_harir_sqr,Q_harir,'filled');
errorbar(p1_harir_sqr,Q_harir,Q_harir_err,Q_harir_err,p1_harir_sqr_err,p1_harir_sqr_err,'o');
f_harir_kiool = fit(p1_harir_sqr',Q_harir','poly1');
plot(f_harir_kiool);
legend('measured results','measured results linear fit');
title('harir- real flow to sqare delta(P)');
xlabel('(p1-p2)^(0.5) [Pa^(0.5)]');
ylabel('real flow [m^3/s]');
xlim([120,200]);
ylim([0.5,6]*10^-4);
hold off;

    %HAKOL
figure(11);
hold on;
grid on;
xlim([120,190]);
ylim([8*10^-5,6*10^-4]);
plot(f_venturi_kiool,'g');
plot(f_nozzle_kiool,'b');
plot(f_harir_kiool,'c');
title('real flow to sqare delta(P)');
xlabel('(p1-p2)^(0.5) [Pa^(0.5)]');
ylabel('real flow [m^3/s]');
legend('venturi linear calibration line','nozzle linear calibration line','harir linear calibration line');
hold off;

    %ANOTHA ONE
    
    for i=[1:1:5]
      Q_venturi_sqare(i) = (Q_venturi(i))^2;
      Q_venturi_sqare_err(i)=2*Q_venturi(i)*Q_venturi_err(i);
    end
    
     for i=[1:1:5]
      Q_nozzle_sqare(i) = (Q_nozzle(i))^2;
      Q_nozzle_sqare_err(i)=2*Q_nozzle(i)*Q_nozzle_err(i);
     end
    
      for i=[1:1:5]
      Q_harir_sqare(i) = (Q_harir(i))^2;
      Q_harir_sqare_err(i)=2*Q_harir(i)*Q_harir_err(i);
      end
    
figure(12);
hold on;
grid on;
xlim([0,3.5]*10^(-7));
g_venturi = fit(Q_venturi_sqare',venturi_p1','poly1');
g_nozzle = fit(Q_nozzle_sqare',nozzle_p1','poly1');
g_harir = fit(Q_harir_sqare',harir_p1','poly1');
% scatter(Q_venturi_sqare,venturi_p1,'filled');
errorbar(Q_venturi_sqare,venturi_p1,venturi_p1_err,venturi_p1_err,Q_venturi_sqare_err,Q_venturi_sqare_err,'o');
plot(g_venturi,'g');
% scatter(Q_nozzle_sqare,nozzle_p1,'filled');
errorbar(Q_nozzle_sqare,nozzle_p1,nozzle_p1_err,nozzle_p1_err,Q_nozzle_sqare_err,Q_nozzle_sqare_err,'o');
plot(g_nozzle,'b');
% scatter(Q_harir_sqare,harir_p1,'filled');
errorbar(Q_harir_sqare,harir_p1,harir_p1_err,harir_p1_err,Q_harir_sqare_err,Q_harir_sqare_err,'o');
plot(g_harir,'c');
title('pressure drop to squared real flow');
xlabel('squared real flow [m^6/s^2]');
ylabel('p1-p2 [Pa]');
legend('venturi','venturi linear fit','nozzle','nozzle linear fit','harir','harir linear fit');
hold off;

%ANALIZA BILTI NISBELET
d_pump = 0.18; %m
%N = 1.093 * V
    %160V
P_VOLT_err=2.5 ;
N_err=1.093*P_VOLT_err;
P_PSI_err=0.25;
P_FINAL_PA_err=P_PSI_err*6894.75729;
P_INIT_M_err=0.02;
Q_err=1.253*10^(-5)*w_err;

P_FINAL_PSI_160 = [6.5,6.5,6,4.5,4.5]; %PSI
P_FINAL_PA_160 = P_FINAL_PSI_160*6894.75729; %PA
P_INIT_M_160 = [0.233,0.241,0.233,0.236,0.236]+0.02; %m
P_INIT_PA_160 = rho*g*P_INIT_M_160; %Pa
DELTA_P_PA_160 = P_FINAL_PA_160 - P_INIT_PA_160; %Pa
TEDER_HZ_160 = [9.7,12.4,15.8,22.3,25.5]; %Hz הערך הנמוך

P_INIT_PA_160_err=sqrt((g*P_INIT_M_160*rho_err).^2+(rho*g*P_INIT_M_err).^2);
DELTA_P_PA_160_err=sqrt(P_FINAL_PA_err.^2 + P_INIT_PA_160_err.^2);

Q_160 = MTM(TEDER_HZ_160);
N_160 = 1.093*160;
MISPAR_P_160 = DELTA_P_PA_160/(rho*(N_160*d_pump)^2);
MISPAR_Q_160 = Q_160/(N_160*d_pump^3);
MISPAR_P_160_err=sqrt((DELTA_P_PA_160_err/(rho*(N_160*d_pump)^2)).^2+(2*DELTA_P_PA_160*N_err/(rho*(N_160^3)*(d_pump^2)).^2));
MISPAR_Q_160_err=sqrt((Q_err/(N_160*(d_pump^3))).^2+(Q_160*N_err/((N_160^2)*(d_pump^3))).^2);
    %180
P_FINAL_PSI_180 = [9,8.5,7.5,7,6.5,6.5]; %PSI
P_FINAL_PA_180 = P_FINAL_PSI_180*6894.75729; %PA
P_INIT_M_180 = [0.336,0.336,0.336,0.336,0.336,0.336]-0.1; %m
P_INIT_PA_180 = rho*g*P_INIT_M_180; %Pa
DELTA_P_PA_180 = P_FINAL_PA_180 - P_INIT_PA_180; %Pa
TEDER_HZ_180 = [12,17.1,22.6,26.1,28.7,29.9]; %Hz הערך הנמוך
Q_180 = MTM(TEDER_HZ_180);
N_180 = 1.093*180;

P_INIT_PA_180_err=sqrt((g*P_INIT_M_180*rho_err).^2+(rho*g*P_INIT_M_err).^2);
DELTA_P_PA_180_err=sqrt(P_FINAL_PA_err.^2 + P_INIT_PA_180_err.^2);

MISPAR_P_180 = DELTA_P_PA_180/(rho*(N_180*d_pump)^2);
MISPAR_P_180_err=sqrt((DELTA_P_PA_180_err/(rho*(N_180*d_pump)^2)).^2+(2*DELTA_P_PA_180*N_err/(rho*(N_180^3)*(d_pump^2)).^2));
MISPAR_Q_180 = Q_180/(N_180*d_pump^3);
MISPAR_Q_180_err=sqrt((Q_err/(N_180*(d_pump^3))).^2+(Q_180*N_err/((N_180^2)*(d_pump^3))).^2);

    %200
P_FINAL_PSI_200 = [11,10.5,10.5,9.5,9.5,8.5]; %PSI
P_FINAL_PA_200 = P_FINAL_PSI_200*6894.75729; %PA
P_INIT_M_200 = [0.336,0.337,0.337,0.336,0.333,0.335]-0.1; %m
P_INIT_PA_200 = rho*g*P_INIT_M_200; %Pa
DELTA_P_PA_200 = P_FINAL_PA_200 - P_INIT_PA_200; %Pa
TEDER_HZ_200 = [14.2,17.5,19.4,24.9,28.0,34.8]; %Hz הערך הנמוך
Q_200 = MTM(TEDER_HZ_200);
N_200 = 1.093*200;

P_INIT_PA_200_err=sqrt((g*P_INIT_M_200*rho_err).^2+(rho*g*P_INIT_M_err).^2);
DELTA_P_PA_200_err=sqrt(P_FINAL_PA_err.^2 + P_INIT_PA_200_err.^2);

MISPAR_P_200 = DELTA_P_PA_200/(rho*(N_200*d_pump)^2);
MISPAR_Q_200 = Q_200/(N_200*d_pump^3);
MISPAR_P_200_err=sqrt((DELTA_P_PA_200_err/(rho*(N_200*d_pump)^2)).^2+(2*DELTA_P_PA_200*N_err/(rho*(N_200^3)*(d_pump^2)).^2));
MISPAR_Q_200_err=sqrt((Q_err/(N_200*(d_pump^3))).^2+(Q_200*N_err/((N_200^2)*(d_pump^3))).^2);
    %220
P_FINAL_PSI_220 = [12.5,12.5,12,11.5,11,10]; %PSI
P_FINAL_PA_220 = P_FINAL_PSI_220*6894.75729; %PA
P_INIT_M_220 = [0.237,0.237,0.236,0.235,0.234,0.232]; %m
P_INIT_PA_220 = rho*g*P_INIT_M_220; %Pa
DELTA_P_PA_220 = P_FINAL_PA_220 - P_INIT_PA_220; %Pa
TEDER_HZ_220 = [11.6,15.1,22.4,29.0,35.6,45.0]; %Hz הערך הנמוך
Q_220 = MTM(TEDER_HZ_220);
N_220 = 1.093*220;

P_INIT_PA_220_err=sqrt((g*P_INIT_M_220*rho_err).^2+(rho*g*P_INIT_M_err).^2);
DELTA_P_PA_220_err=sqrt(P_FINAL_PA_err.^2 + P_INIT_PA_220_err.^2);


MISPAR_P_220 = DELTA_P_PA_220/(rho*(N_220*d_pump)^2);
MISPAR_Q_220 = Q_220/(N_220*d_pump^3);
MISPAR_P_220_err=sqrt((DELTA_P_PA_220_err/(rho*(N_220*d_pump)^2)).^2+(2*DELTA_P_PA_220*N_err/(rho*(N_220^3)*(d_pump^2)).^2));
MISPAR_Q_220_err=sqrt((Q_err/(N_220*(d_pump^3))).^2+(Q_220*N_err/((N_220^2)*(d_pump^3))).^2);

figure(13);
hold on;
grid on;
% scatter(Q_160,DELTA_P_PA_160,'filled');
errorbar(Q_160,DELTA_P_PA_160,DELTA_P_PA_160_err,DELTA_P_PA_160_err,Q_err*ones(1,5),Q_err*ones(1,5),'o');
% scatter(Q_180,DELTA_P_PA_180,'filled');
errorbar(Q_180,DELTA_P_PA_180,DELTA_P_PA_180_err,DELTA_P_PA_180_err,Q_err*ones(1,6),Q_err*ones(1,6),'o');
% scatter(Q_200,DELTA_P_PA_200,'filled');
errorbar(Q_200,DELTA_P_PA_200,DELTA_P_PA_200_err,DELTA_P_PA_200_err,Q_err*ones(1,6),Q_err*ones(1,6),'o');
% scatter(Q_220,DELTA_P_PA_220,'filled');
errorbar(Q_220,DELTA_P_PA_220,DELTA_P_PA_220_err,DELTA_P_PA_220_err,Q_err*ones(1,6),Q_err*ones(1,6),'o');

title('delta pressure to flow rate');
xlabel('flow rate [m^3/s]');
ylabel('delta pressure [Pa]');
legend('160V','180V','200V','220V');
hold off;

figure(14);
hold on;
grid on;
% scatter(MISPAR_Q_160,MISPAR_P_160,'filled');
errorbar(MISPAR_Q_160,MISPAR_P_160,MISPAR_P_160_err,MISPAR_P_160_err,MISPAR_Q_160_err,MISPAR_Q_160_err,'o');
% scatter(MISPAR_Q_180,MISPAR_P_180,'filled');
errorbar(MISPAR_Q_180,MISPAR_P_180,MISPAR_P_180_err,MISPAR_P_180_err,MISPAR_Q_180_err,MISPAR_Q_180_err,'o');
% scatter(MISPAR_Q_200,MISPAR_P_200,'filled');
errorbar(MISPAR_Q_200,MISPAR_P_200,MISPAR_P_200_err,MISPAR_P_200_err,MISPAR_Q_200_err,MISPAR_Q_200_err,'o');
% scatter(MISPAR_Q_220,MISPAR_P_220,'filled');
errorbar(MISPAR_Q_220,MISPAR_P_220,MISPAR_P_220_err,MISPAR_P_220_err,MISPAR_Q_220_err,MISPAR_Q_220_err,'o');

title('MISPAR_P TO MISPAR_Q');
xlabel('MISPAR_Q');
ylabel('MISPAR_P');
legend('160V','180V','200V','220V');
hold off;
table_160=[P_INIT_PA_160;P_FINAL_PA_160;DELTA_P_PA_160;TEDER_HZ_160;Q_160;MISPAR_P_160;MISPAR_Q_160];
table_180=[P_INIT_PA_180;P_FINAL_PA_180;DELTA_P_PA_180;TEDER_HZ_180;Q_180;MISPAR_P_180;MISPAR_Q_180];
table_200=[P_INIT_PA_200;P_FINAL_PA_200;DELTA_P_PA_200;TEDER_HZ_200;Q_200;MISPAR_P_200;MISPAR_Q_200];
table_220=[P_INIT_PA_220;P_FINAL_PA_220;DELTA_P_PA_220;TEDER_HZ_220;Q_220;MISPAR_P_220;MISPAR_Q_220];
