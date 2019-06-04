%We start by declaring the constants of the exercise
Vtn=0.4; %V
Vs=0;   %V
Eo=8.85E-12; %F/M
h=0.02; %Ve-1
ut=26e-3; %V
n=1.2;
Eox=3.9*Eo; %F/M
W=180e-9;
L=90e-9;
tox=3.5e-9;
u=180e-4; %cm^2/(V·s)
Vgs=[0,0.3,0.6,0.9,1.2,1.5,1.8];

%Now we calculate the Voltage overdrive and other parameter
Vod=Vgs-Vtn
Cox=Eox/tox;
Beta=u*Cox*W/L;
K=Beta/2;

%Now we proceed with analyzing the different equations of Id
Vds=linspace(0,2,2000);

%----------PROCEDURE TO SOLVE THE FIRST EXERCISE---------------%
%( As we see from the results of the overdrive voltage..
%.. we only need to apply the weak inversion model on the first..
%.. and second value of Vgs() )
Id0=2*n*Beta*(ut^2)*exp(-Vtn/(n*ut));
Idcut1=Id0*exp(Vgs(1)/(n*ut))*(exp(Vs/ut)-exp(-Vds/ut));
Idcut2=Id0*exp(Vgs(2)/(n*ut))*(exp(Vs/ut)-exp(-Vds/ut));
Idcut3=0;
Idcut4=0;
Idcut5=0;
Idcut6=0;
Idcut7=0;

Idsat=K*(Vod.*Vod);

%Before start ploting the results, we need to apply..
%the channel length modulation for the saturation region
%Id1oh=Beta*((Vgs(1)-Vtn).*Vds-((Vds.*Vds)/2));
%Id1sat=Idsat(1)*(1+h.*Vds);

%Id2oh=Beta*((Vgs(2)-Vtn).*Vds-((Vds.*Vds)/2));
%Id2sat=Idsat(2)*(1+h.*Vds);

Id3oh=Beta*((Vgs(3)-Vtn).*Vds-((Vds.*Vds)/2));
Id3sat=Idsat(3)*(1+h.*Vds);

Id4oh=Beta*((Vgs(4)-Vtn).*Vds-((Vds.*Vds)/2));
Id4sat=Idsat(4)*(1+h.*Vds);

Id5oh=Beta*((Vgs(5)-Vtn).*Vds-((Vds.*Vds)/2));
Id5sat=Idsat(5)*(1+h.*Vds);

Id6oh=Beta*((Vgs(6)-Vtn).*Vds-((Vds.*Vds)/2));
Id6sat=Idsat(6)*(1+h.*Vds);

Id7oh=Beta*((Vgs(7)-Vtn).*Vds-((Vds.*Vds)/2));
Id7sat=Idsat(7)*(1+h.*Vds);

%Now we plot the obtained curves
%For doiing it well, we first discriminate the regions of the current
%curve1=(Vs>Vod(1)&(Vds>Vod(1))).*Idcut1+((Vds<Vod(1))&(Vs<Vod(1))).*Id1oh+(Vds>Vod(1)&Vs<Vod(1)).*Id1sat;
%curve2=(Vs>Vod(2)&(Vds>Vod(2))).*Idcut2+((Vds<Vod(2))&(Vs<Vod(2))).*Id2oh+(Vds>Vod(2)&Vs<Vod(2)).*Id2sat;
curve1=Idcut1;
curve2=Idcut2;
curve3=(Vs>Vod(3)&(Vds>Vod(3))).*Idcut3+((Vds<Vod(3))&(Vs<Vod(3))).*Id3oh+(Vds>Vod(3)&Vs<Vod(3)).*Id3sat;
curve4=(Vs>Vod(4)&(Vds>Vod(4))).*Idcut4+((Vds<Vod(4))&(Vs<Vod(4))).*Id4oh+(Vds>Vod(4)&Vs<Vod(4)).*Id4sat;
curve5=(Vs>Vod(5)&(Vds>Vod(5))).*Idcut5+((Vds<Vod(5))&(Vs<Vod(5))).*Id5oh+(Vds>Vod(5)&Vs<Vod(5)).*Id5sat;
curve6=(Vs>Vod(6)&(Vds>Vod(6))).*Idcut6+((Vds<Vod(6))&(Vs<Vod(6))).*Id6oh+(Vds>Vod(6)&Vs<Vod(6)).*Id6sat;
curve7=(Vs>Vod(7)&(Vds>Vod(7))).*Idcut7+((Vds<Vod(7))&(Vs<Vod(7))).*Id7oh+(Vds>Vod(7)&Vs<Vod(7)).*Id7sat;
%Vds=Vod
curve10=K*(Vds.*Vds);


figure(1)
subplot(1,2,1)
plot(Vds,curve1,Vds,curve2,Vds,curve3,Vds,curve4,Vds,curve5,Vds,curve6,Vds,curve7,Vds,curve10,'k--');
title('I_D vs. V_{DS} for different values of V_{GS}');
xlabel('V_{DS} (V)');
ylabel('I_D (A)');
ylim([0 4e-4]);
grid on;
legend('V_{GS}=0 V','V_{GS}=0.3 V','V_{GS}=0.6 V','V_{GS}=0.9 V','V_{GS}=1.2 V','V_{GS}=1.5 V','V_{GS}=1.8 V','V_{DS}=V_{GS}-V_{TH}');

subplot(1,2,2)
plot(Vds,curve1,Vds,curve2,Vds,curve3,Vds,curve4,Vds,curve5,Vds,curve6,Vds,curve7,Vds,curve10,'k--');
title('log(I_D) vs. V_{DS} for different values of V_{GS}');
xlabel('V_{DS} (V)');
ylabel('I_D (A)');
ylim([0 4e-4]);
grid on;
set(gca, 'YScale', 'log')
legend('V_{GS}=0 V','V_{GS}=0.3 V','V_{GS}=0.6 V','V_{GS}=0.9 V','V_{GS}=1.2 V','V_{GS}=1.5 V','V_{GS}=1.8 V','V_{DS}=V_{GS}-V_{TH}');

%--------------------PROCEDURE TO SOLVE THE SECOND EXERCISE----------%
%We begin by determining a high enough value of Vds...
%...and we create a Vgs linspace variable for the curve
Vds_ex2=2.5;
Vgs_ex2=linspace(0,1.8,1000);

%With the previous values, we calculate the different...
%...current drains depending on the weak or strong inversion

%For the weak inversion model
Id0_ex2=2*n*Beta*(ut^2)*exp(-Vtn/(n*ut));
Id_ex2weak=Id0_ex2*exp(Vgs_ex2./(n*ut))*(exp(Vs/ut)-exp(-Vds_ex2/ut));

%For the strong inversion model
Vod_ex2=(Vgs_ex2)-Vtn;

Id=K*(Vod_ex2.*Vod_ex2);
Id_ex2sat=Id*(1+h.*Vds_ex2);

Id_ex2ohm=Beta*((Vgs_ex2-Vtn).*Vds_ex2-((Vds_ex2.*Vds_ex2)/2));

%Now we plot the obtained curves
curve8=(Vod_ex2<-0.2).*Id_ex2weak+((Vod_ex2>0.2)&(Vds_ex2>Vod_ex2)&(Vs<Vod_ex2)).*Id_ex2sat+((Vod_ex2>0.2)&(Vds_ex2<Vod_ex2)&(Vs<Vod_ex2)).*Id_ex2ohm;
curve9=log(curve8);

figure(2)
subplot(1,2,1)
plot(Vgs_ex2,curve8);
title('I_D vs. V_{GS} for V_{DS}=2.5 V');
xlabel('V_{GS} (V)');
ylabel('I_D (A)');
xlim([0 1.8]);
grid on;

subplot(1,2,2)
plot(Vgs_ex2,curve8);
title('log(I_D) vs. V_{GS} for V_{DS}=2.5 V');
xlabel('V_{GS} (V)');
ylabel('I_D (A)');
xlim([0 1.8]);
set(gca, 'YScale', 'log')
grid on;
