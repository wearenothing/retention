clear all; clc;
tic
shot=95122;
% shot=100000;
filepath='D:\li\zhipu\';
zhipu=textread([filepath num2str(shot),'.txt'],'','headerlines',1);
time_b=abs(zhipu(:,1)-59);%Ѱ�������зŵ翪ʼʱ�䣬ͨ������ֵ��60�룬�õ���ӽ���ʼ�ŵ�ʱ�̵ļ�¼ֵ
time_e=abs(zhipu(:,1)-300);%Ѱ�������зŵ���ֹʱ�䣬ͨ������ֵ��70�룬�õ���ӽ������ŵ�ʱ�̵ļ�¼ֵ
[b1,b2]=find(time_b==min(time_b));%ȷ���ڵ�������зŵ翪ʼʱ�̵�ָ��λ��
[b3,b4]=find(time_e==min(time_e));%ȷ���ڵ�������зŵ����ʱ�̵�ָ��λ��
DH=(zhipu(b1:1:b3,4).*2+zhipu(b1:1:b3,3))./((zhipu(b1:1:b3,2)+zhipu(b1:1:b3,3)+zhipu(b1:1:b3,4)).*2);
timeDH=zhipu(b1:1:b3,1)-60;  
P=polyfit(timeDH,DH,5); % P=polyfit(timeDH,DH,7)

% V_vessel = 110;
Pam2P = 2.41e20;

if shot>88344
    V_SMBI=2.0431e-4; 
else
    V_SMBI=3984.31.*(1e-6); %SMBI3�����������m-3
end
%2021����ѹ�޽������飬ȡ��KHG1��JHG6,JHG3�����ڳ�Ar
%���ڼ��ϵĹ���
V_JHG1 = 3.118e-4;%PUFFϵͳ���������m-3
V_DHG1=2.919E-4;%D����E���Ӳ����߳��������m-3 
V_JHG2 = 2.91e-4;%J����ˮƽ������λ��ƫ�������������m-3
V_JHG4 = 2.99e-4;% N���ڱ�ICRF���߼��ϳ��������m-3
V_JHG5 = 2.997e-4; % N������ICRF���߼��ϳ��������m-3
%V_KHG1=2.919E-4;%K����N���Ӳ����߳��������m-3
%V_JHG6 = 2.949e-4;% J������ICRF���߼��ϳ��������m-3

%����ƫ���������������
V_OUG1 = 1.4687e-3;% J����ƫ�������������m-3
V_ODG1 = 1.4601e-3;% J����ƫ�������������m-3
V_CDG1 = 6.792e-3;% J����ƫ�������������m-3
V_HDG1 = 6.75e-3;% J����ƫ�������������m-3


mdsconnect('mds.ipp.ac.cn');
%���ϵ��
KH=2.9;
mdsopen('east', shot);%��ʱ��east_1
G107 = 10.^(1.667*mdsvalue('\G107')-9.333).*KH;
G105 = 10.^(1.667*mdsvalue('\G105')-9.333).*KH;
G109 = 10.^(1.667*mdsvalue('\G109')-9.333).*KH;
G106 = 10.^(1.667*mdsvalue('\G106')-9.333).*KH;
VJHpev3=mdsvalue('\VJHpev3');
time_VJHpev3=mdsvalue('dim_of(\VJHpev3)');
index=find(VJHpev3>0.3);
%-------------------------
jhg1=mdsvalue('\jhg1');
jhg2=mdsvalue('\jhg2');
dhg1=mdsvalue('\Dhg1');
jhg4=mdsvalue('\jhg4');
jhg5=mdsvalue('\jhg5');

JHG1 = (mdsvalue('\jhg1')-1)*2.5e4; % JHG1Ϊ�����źţ���λPa
JHG2 = (mdsvalue('\jhg2')-1)*2.5e4; % 
DHG1 = (mdsvalue('\Dhg1')-1)*2.5e4; %D���ڵ��Ӳ����߳���
JHG4 = (mdsvalue('\jhg4')-1)*2.5e4; % B���ڱ�ICRF���߼���
JHG5 = (mdsvalue('\jhg5')-1)*2.5e4; % J������ICRF���߼���

%�������ݳ���
JHG1=JHG1(1:10:end,1);
JHG2=JHG2(1:10:end,1);
DHG1=DHG1(1:10:end,1);
JHG4=JHG4(1:10:end,1);
JHG5=JHG5(1:10:end,1);
%-------------------------
%ƫ��������
oug1=mdsvalue('\OUG1');
odg1=mdsvalue('\ODG1');
cdg1=mdsvalue('\CDG1');
hdg1=mdsvalue('\HDG1');
OUG1 = (mdsvalue('\OUG1')-1)*2.5e4; % J����ƫ����
ODG1 = (mdsvalue('\ODG1')-1)*2.5e4; % J����ƫ����
CDG1 = (mdsvalue('\CDG1')-1)*2.5e4; % J����ƫ����
HDG1 = (mdsvalue('\HDG1')-1)*2.5e4; % J����ƫ����
%�������ݳ���
OUG1=OUG1(1:10:end,1);
ODG1=ODG1(1:10:end,1);
CDG1=CDG1(1:10:end,1);
HDG1=HDG1(1:10:end,1);
%�������ݳ���

NBI1LHI=0.8.*(mdsvalue('\NBI1LHI'));%NBI ע������ź�
NBI1RHI=0.8.*(mdsvalue('\NBI1RHI'));
NBI2LHI=0.8.*(mdsvalue('\NBI2LHI'));
NBI2RHI=0.8.*(mdsvalue('\NBI2RHI'));
time_NBI=mdsvalue('dim_of(\NBI2RHI)');

%�������ݳ���
G107=G107(1:10:end,1);
G105=G105(1:10:end,1);
G109=G109(1:10:end,1);
G106=G106(1:10:end,1);
%�������ݳ���
%-----------------------------------
ipm = mdsvalue('\ipm');
time = mdsvalue('dim_of(\G107)');
time=time(1:10:end,1);%�ɱ仯 
time1 = mdsvalue('dim_of(\ipm)');
time1=time1(1:10:end,1);%�ɱ仯 
length_time1=length(time);
length_time2=length(time1);
time_puffing=mdsvalue('dim_of(\jhg1)');
length_timepuffing=length(time_puffing);

background = mean(G107(1:100));
G107 = G107-background;
%-------------------------------------

%��������

%LHW ����
background_DHG1 = mean(DHG1(1:5000));
if max(dhg1)-min(dhg1)<0.3
    DHG1=zeros(length_time1,1);
    puff_DHG=V_DHG1.*DHG1;
else
    for j=1:length_time1
        puff_DHG(j,1)= V_DHG1*(background_DHG1-DHG1(j));
    end
end
puff_DHG;
puff_LHW=puff_DHG;
puff_LHW=smooth(puff_LHW,1000,'sgolay');
[M_LHW,N_LHW]=find(puff_LHW==max(puff_LHW));
M_LHW=min(M_LHW);
if length(M_LHW)>10
    puff_LHWlong=puff_LHW;
else
puff_LHWlong(1:1:M_LHW,1)=puff_LHW(1:1:M_LHW,1);
puff_LHWlong((M_LHW+1):1:length_time1,1)=ones(length_time1-M_LHW,1)*puff_LHW(M_LHW,1);
end
puff_LHWlong=smooth(puff_LHWlong,100,'sgolay');

%GAS PUFF��������JHG1+JHG2
background_JHG1 = mean(JHG1(1:500));
JHG1=smooth(JHG1,100,'sgolay');
for j=1:length_time1
        puff_JHG1(j,1)= V_JHG1*(background_JHG1-JHG1(j));  
end
puff_JHG1;
background_JHG2 = mean(JHG2(1:500));
JHG2=smooth(JHG2,100,'sgolay');
for j=1:length_time1
        puff_JHG2(j,1)= V_JHG2*(background_JHG2-JHG2(j));  
end
puff_GAS=puff_JHG1+puff_JHG2;
puff_GAS=smooth(puff_GAS,1000,'sgolay');%��ǰƽ����������ŵ��������ֶ������������

% [M_GAS,N_GAS]=find(puff_GAS==max(puff_GAS(1:length(time1),1)));
% if length(M_GAS)>1
%     M_GAS=max(M_GAS);
% end
M_GAS=floor(max(index)./10);
puff_GASlong(1:1:M_GAS,1)=puff_GAS(1:1:M_GAS,1);
puff_GASlong((M_GAS+1):1:length_time1,1)=ones(length_time1-M_GAS,1)*puff_GASlong(M_GAS,1);

   
%ICRF��������
background_JHG4 = mean(JHG4(1:5000));
background_JHG5 = mean(JHG5(1:5000));

if max(jhg4)-min(jhg4)<0.3
    JHG4=zeros(length_time1,1);
    puff_JHG4=V_JHG4.*JHG4;
else
    for j=1:length_time1
        puff_JHG4(j,1)= V_JHG4*(background_JHG4-JHG4(j));
    end
end
puff_JHG4;

if max(jhg5)-min(jhg5)<0.3
    JHG5=zeros(length_time1,1);
    puff_JHG5=V_JHG5.*JHG5;
else
    for j=1:length_time1
        puff_JHG5(j,1)= V_JHG5*(background_JHG5-JHG5(j));
    end
end
puff_JHG5;
puff_ICRF=puff_JHG4+puff_JHG5;
puff_ICRF=smooth(puff_ICRF,100,'sgolay');
[M_ICRF,N_ICRF]=find(puff_ICRF==max(puff_ICRF));
if length(M_ICRF)>10
    puff_ICRFlong=puff_ICRF;
else
    puff_ICRFlong(1:1:M_ICRF,1)=puff_ICRF(1:1:M_ICRF,1);
    puff_ICRFlong((M_ICRF+1):1:length_time1,1)=ones((length_time1-M_ICRF),1)*puff_ICRF(M_ICRF,1);
end
puff_ICRFlong=smooth(puff_ICRFlong,100,'sgolay');

%ƫ������������
background_OUG1 = mean(OUG1(1:5000));
background_ODG1 = mean(ODG1(1:5000));
background_CDG1= mean(CDG1(1:5000));
background_HDG1= mean(HDG1(1:5000));

if max(oug1)-min(oug1)<0.3
    OUG1=zeros(length_time1,1);
    puff_OUG1=V_OUG1.*OUG1;
else
    for j=1:length_time1
        puff_OUG1(j,1)= V_OUG1*(background_OUG1-OUG1(j));
    end
end
puff_OUG1;
if max(cdg1)-min(cdg1)<0.3
    CDG1=zeros(length_time1,1);
    puff_CDG1=V_CDG1.*CDG1;
else
    for j=1:length_time1
        puff_CDG1(j,1)= V_CDG1*(background_CDG1-CDG1(j));
    end
end
puff_CDG1;
if max(odg1)-min(odg1)<0.3
    ODG1=zeros(length_time1,1);
    puff_ODG1=V_ODG1.*ODG1;
else
    for j=1:length_time1
        puff_ODG1(j,1)= V_ODG1*(background_ODG1-ODG1(j));
    end
end
puff_ODG1;
if max(hdg1)-min(hdg1)<0.3
    HDG1=zeros(length_time1,1);
    puff_HDG1=V_HDG1.*HDG1;
else
    for j=1:length_time1
        puff_HDG1(j,1)= V_HDG1*(background_HDG1-HDG1(j));
    end
end
puff_HDG1;
puff_divertor=puff_OUG1+puff_CDG1+puff_ODG1+puff_HDG1;
puff_divertor=smooth(puff_divertor,100,'sgolay');
[M_divertor,N_divertor]=find(puff_divertor==max(puff_divertor));
if length(M_divertor)>10
    puff_divertorlong=puff_divertor;
else
    puff_divertorlong(1:1:M_divertor,1)=puff_divertor(1:1:M_divertor,1);
    puff_divertorlong((M_divertor+1):1:length_time1,1)=ones(length_time1-M_divertor,1)*puff_divertor(M_divertor,1);
end
puff_divertorlong=smooth(puff_divertorlong,100,'sgolay');


mdsopen('east', shot);
%SMBI���� SMBI2�ź�ΪPAS,SMBI3�ź�ΪPJS
%SMBI2����
SMBI2=mdsvalue('\SMBI2');
time_SMBI001=mdsvalue('dim_of(\SMBI2)');
%SMBI_zs=mean(SMBI2(1:70000));
MAXSMBI=max(SMBI2);
PJS203=(mdsvalue('\PJS203').*(4*10^5));%data link and rectify
time_PJS203=mdsvalue('dim_of(\PJS203)');
background2 = mean(PJS203(6970:7000));%SMBI���ӳ�ʼѹǿ
SMBI2g= V_SMBI*(background2-(PJS203)); %����������Ħ�����仯
%SMBI2g=SMBI2g(1:10:length_time1,1);%����̫�࣬ȡԭ����ʮ��֮һ������

if V_SMBI==3984.31*(1e-6) %SMBI3�����������m-3
    SMBI2g=smooth(SMBI2g,5000,'sgolay');
else
    SMBI2g=smooth(SMBI2g,100,'sgolay');
end

for i=1:10*length_time1
   if SMBI2g(i)<0
        SMBI2g(i)=0;
   end
end
[M1,N1]=find(SMBI2g==max(SMBI2g));
if length(M1)>1
    M1=max(M1);
end
SMBIlong2g(1:1:M1,1)=SMBI2g(1:1:M1,1);
SMBIlong2g((M1+1):1:10*length_time1,1)=ones(10*length_time1-M1,1)*SMBIlong2g(M1,1);

if MAXSMBI<1.4
    SMBIlong2g=zeros(10*length_time1,1);
end
SMBI_zs2=mean(SMBI2(1:70000));
SMBI002=SMBI2-SMBI_zs2;
SMBIlong2=cumtrapz(time_SMBI001, SMBI002);
SMBIlong2=SMBIlong2(1:10:end);
%SMBIlong2=SMBIlong2(1:length_time1);

for i=1:1:length(time)
    rate(i,1)=2*Pam2P*SMBIlong2g(i,1)./SMBIlong2(i,1);
    if SMBIlong2(i,1)==0
        rate(i,1)=0;
    end
end
coeff=abs(mean(rate(11000:1:end)));
SMBIlong22=coeff*SMBIlong2;
SMBIlong22=SMBIlong2g;


%SMBI3����
SMBI3=mdsvalue('\SMBI3');
time_SMBI002=mdsvalue('dim_of(\SMBI3)');
MAXSMBI3=max(SMBI3);
PAS103=(mdsvalue('\PAS103').*(4*10^5));%data link and rectify
time_PAS103=mdsvalue('dim_of(\PAS103)');
background3 = mean(PAS103(69300:70000));%SMBI���ӳ�ʼѹǿ
SMBI3g= V_SMBI*(background3-(PAS103)); %����������Ħ�����仯
%SMBI3g=SMBI3g(1:10:(10*length_time1),1);%����̫�࣬ȡԭ����ʮ��֮һ

if V_SMBI==3984.31*(1e-6) %SMBI3�����������m-3
    SMBI3g=smooth(SMBI3g,5000,'sgolay');
else
    SMBI3g=smooth(SMBI3g,100,'sgolay');
end
%SMBI=SMBI2+SMBI3;
for i=1:10*length_time1
   if SMBI3g(i)<0
        SMBI3g(i)=0;
   end
end
%SMBI=smooth(SMBI,5000,'sgolay');
[M2,N2]=find(SMBI3g==max(SMBI3g));
SMBIlong3g(1:1:M2,1)=SMBI3g(1:1:M2,1);
SMBIlong3g((M2+1):1:10*length_time1,1)=ones(10*length_time1-M2,1)*SMBIlong3g(M2,1);

if MAXSMBI3<1.4
    SMBIlong3g=zeros(10*length_time1,1);
end

SMBI_zs3=mean(SMBI3(1:70000));
SMBI003=SMBI3-SMBI_zs3;
SMBIlong3=cumtrapz(time_SMBI002,SMBI003);

SMBIlong3=SMBIlong3(1:10:end);
%SMBIlong3=SMBIlong3(1:length_time1);

for i=1:1:length(time)
    rate(i,1)=2*Pam2P*SMBIlong3g(i,1)./SMBIlong3(i,1);
    if SMBIlong3(i,1)==0
        rate(i,1)=0;
    end
end
coeff=abs(mean(rate(11000:1:end)));
SMBIlong33=coeff*SMBIlong3;
SMBIlong33=SMBIlong3g;
SMBIlong=2*Pam2P*(SMBIlong33+SMBIlong22);
SMBIlong=SMBIlong(1:10:end);

%NBI ע�벿��
background_NBI1LHI=mean(NBI1LHI(1:1500));
NBI1LHI=NBI1LHI-background_NBI1LHI;
background_NBI1RHI=mean(NBI1RHI(1:1500));
NBI1RHI=NBI1RHI-background_NBI1RHI;
background_NBI2RHI=mean(NBI2RHI(1:1500));
NBI2RHI=NBI2RHI-background_NBI2RHI;
background_NBI2LHI=mean(NBI2LHI(1:1500));
NBI2LHI=NBI2LHI-background_NBI2LHI;
NBI_add= NBI1LHI+NBI1RHI+NBI2LHI+NBI2RHI;
NBI_add= cumtrapz(time_NBI, NBI_add);
%NBI_add2=NBI_add(1:5:end);
NBI_add2=interp1(time_NBI,NBI_add,time1,'linear');
NBI_addlong(1:length(time1),1)=NBI_add2;
NBI_addlong(length(time1)+1:1:length(time),1)=ones((length(time)-length(time1)),1).*NBI_add2(length(time1),1);
%%
%�ܳ�����
%���Ӳ���̨���±ó���
for i=1:length_time1
    if G107(i,1)>1.3*10^(-3)
        pump_LHWlow(i,1)=(1.32*(G107(i,1)).^(-0.2757)-0.3481); %У׼��Ĺ�ϵʽ��ax^b+c,abc 1.32 -0.2757 -0.3481)
    else%vif G107(i,1)<1.3*10^(-3)&&G107(i,1)==1.3*10^(-3)
        pump_LHWlow(i,1)=(1.32*(1.3*10^(-3)).^(-0.2757)-0.3481); 
    end
end
%���з��ӱã�����4̨+���Ӳ���̨+���ӻ���+�ڴ���1̨��
for i=1:length_time1
    if G107(i,1)>7*10^(-4)
        pump_midallcry(i,1)=(11.56*(G107(i,1)).^(-0.1736)-9.375); %У׼��Ĺ�ϵʽ��ax^b+c,abc 11.56 -0.1736 9.375)
    else%if G107(i,1)<7*10^(-4)&&G107(i,1)==7*10^(-4)
        pump_midallcry(i,1)=(11.56*(7*10^(-4)).^(-0.1736)-9.375); 
    end
end
%��ƽ��5̨���±ó���
for i=1:length_time1
    if G107(i,1)>7*10^(-4)
        pump_middw(i,1)=(5.667*(G107(i,1)).^(-0.2915)); %У׼��Ĺ�ϵʽ��ax^b,ab 5.667 -0.2915)
    else %G107(i,1)<7*10^(-4)&&G107(i,1)==7*10^(-4)
        pump_middw(i,1)=(5.667*(7*10^(-4)).^(-0.2915)); 
    end
end
pump_mid=(pump_LHWlow+pump_midallcry+pump_middw).*G107;

%��ƫ��������
%��ƫ��������ó���
for i=1:length_time1
     if G109(i,1)>6.2*10^(-3)
        pump_upin(i,1)=(0.829*(G109(i,1)).^(-0.7144)+19.71); %У׼��Ĺ�ϵʽ��ax^b+c)
     else %if G109(i,1)<6.2^(-3)&&G109(i,1)==6.2^(-3)
        pump_upin(i,1)=(0.829*(6.2*10^(-3)).^(-0.7144)+19.71);
         
     end
end
%������ó���
for i=1:length_time1
     if G109(i,1)>5*10^(-3)
        pump_upout(i,1)=(0.8307*(G109(i,1)).^(-0.4968)+3.69); %У׼��Ĺ�ϵʽ��ax^b+c,a=0.8307,b=-0.4968,c=3.69)
     else %if G109(i,1)<5*10^(-3)&&G109(i,1)==5*10^(-3)
        pump_upout(i,1)=(0.8307*(5*10^(-3)).^(-0.4968)+3.69);
        
     end
end
pump_up=(pump_upin+pump_upout).*G109;
%��ƫ��������
%��ƫ��������ó���
for i=1:length_time1
    if G106(i,1)>6.9*10^(-3)
        pump_lowin(i,1)=(2.062*(G106(i,1)).^(-0.49)+21.26);%(abc,2.062,-0.49,21.26)
    else%if G106(i,1)<6.9*10^(-3)&&G106(i,1)==6.9*10^(-3)
        pump_lowin(i,1)=(2.062*(6.9*10^(-3)).^(-0.49)+21.26);
    end
end
%�������3̨����
for i=1:length_time1
    if G106(i,1)>5.1*10^(-3)
        pump_lowout(i,1)=(0.2514*(G106(i,1)).^(-0.4809)+2.462);%(abc,2.062,-0.49,21.26)
    else%if G106(i,1)<5.1*10^(-3)&&G106(i,1)==5.1*10^(-3)
        pump_lowout(i,1)=(0.2514*(5.1*10^(-3)).^(-0.4809)+2.462);
    end
end
pump_low=(pump_lowout+pump_lowin).*G106;  
pump = pump_mid+pump_up+pump_low;   
pump = cumtrapz(time, pump);
%%


mdsopen('efit_east', shot);
c = mdsvalue('\volume');
b = mdsvalue('dim_of(\volume)');

for i = 1:length_time1
    if time(i)<max(b) && time(i)>min(b)
        V_plasma(i) = spline(b, c, time(i));
    else
        V_plasma(i) = 0;
    end
    if V_plasma(i)<0
       V_plasma(i)=0;
    end
end
V_plasma = V_plasma';

mdsopen('pcs_east', shot);
B = mdsvalue('\dfsdev');
A = mdsvalue('dim_of(\dfsdev)');
% B = y';
% A=x';
%mdsclose;
for i = 1:length_time1
    if time(i)<max(A) && time(i)>min(A)
        Ne(i) = spline(A, B, time(i));
    else
        Ne(i) = 0;
    end
   % if Ne(i) < 0
   %     Ne(i) = 0;
    %end
end
Ne = Ne';
P_plasma = 1e19.*Ne.*V_plasma./1e6;

%����������
PH=polyval(P,time);%D��������
% PH = 0.5;
V=60;%����������
R=8.31;
QVV=2.*PH.*G107.*V./(300*R) *6.02*10^23;%������ڵ�����������;
%QVV=QVV(1:length_time1,:);
puff_all=puff_LHWlong+puff_GASlong+puff_ICRFlong+puff_divertorlong;
%puff_all=smooth(puff_all,100,'sgolay');

add_all=2*Pam2P*(puff_all)+SMBIlong+NBI_addlong.*6.25*10^(18);
%add_all=2*Pam2P*(puff_all);
add_all=smooth(add_all,500,'sgolay');

retention =add_all-2*Pam2P.*PH.*pump-QVV;
retention_wall = retention-P_plasma;
retention_wall=smooth(retention_wall,100,'sgolay');
[m,n]=find(retention_wall==max(retention_wall));
for i=1:1:length(retention_wall)
    if i<m+1
        recover_wall(i,1)=0;
    else
        recover_wall(i,1)=max(retention_wall)-retention_wall(i,1);
    end
end
recover_wall;
retention_ratio=retention_wall./(add_all);

writematrix([time,retention_ratio],'retention_ratio_ph.csv')

PUMP=2*Pam2P.*PH.*pump;
PP=max(PUMP);
mdsclose;
figure(1);
plot(time,retention_wall/1e22,time,P_plasma/1e22,time,2*Pam2P*puff_all/1e22,time,SMBIlong/1e22,time,2.*Pam2P.*PH.*pump/1e22,time,NBI_addlong.*6.25*10^(18)/1e22,time,add_all/1e22,'linewidth',3)
xlabel(('Time(s)'));
xlim([0 110]);
ylabel({'Particles','(10^{22})'});
legend('retention','plasma','puff','SMBI','pump','NBI','add all');
set(get(gca,'xlabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(get(gca,'ylabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(gca,'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(gca,'linewidth',3);

figure(2);
subplot(2,1,1)
plot(time,retention_wall/1e22,time,add_all/1e22,'linewidth',3)
xlabel(('Time(s)'));
ylabel({'Particles','(10^{22})'});
legend('Retention','Puff');
set(get(gca,'xlabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(get(gca,'ylabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(gca,'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(gca,'linewidth',3); 
subplot(2,1,2)
plot(time,retention_ratio,'linewidth',3)
xlabel(('Time(s)'));
ylabel({'Retention ratio','(a.u.)'});
legend(num2str(shot));
legend('boxoff');
set(get(gca,'xlabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(get(gca,'ylabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(gca,'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(gca,'linewidth',3);

figure(3)
retention_rate=diff(retention_wall).*1000;
retention_rate=smooth(retention_rate,1000,'sgolay');
puff_rate=diff(add_all).*1000;
puff_rate=smooth(puff_rate,1000,'sgolay');
puffing_rate=diff(puff_all).*1000;
puffing_rate=2*Pam2P*smooth(puffing_rate,1000,'sgolay');
SMBI_rate=diff(SMBIlong).*1000;
SMBI_rate=smooth(SMBI_rate,1000,'sgolay');
pump_rate=diff(PUMP).*1000;
plot(time(1:end-1),retention_rate/1e21,time(1:end-1),puff_rate/1e21,time(1:end-1),puffing_rate/1e21,time(1:end-1),SMBI_rate/1e21,time(1:end-1),pump_rate/1e21,'linewidth',3)
xlabel(('Time(s)'));
ylabel({'Rate','(10^{21})'});
legend('Retention','Puff_rate','puffing_rate','SMBI_rate','pump');
set(get(gca,'xlabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(get(gca,'ylabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(gca,'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(gca,'linewidth',3);


figure(4)
subplot(3,1,1)
set(gca,'position',[0.2 0.65 0.5 0.25]);
plot(time,retention_wall/1e22,time,add_all/1e22,time,recover_wall/1e22,time,PUMP/1e22,'linewidth',3)
xlim([0 30]);
xlabel(('Time(s)'));
ylabel({'Particles','(10^{22})'});
legend('Retention','Puffing','recovery','pumping');
set(get(gca,'xlabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(get(gca,'ylabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(gca,'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(gca,'XTickLabel','','linewidth',3);
set(gca,'linewidth',3);
subplot(3,1,2)
set(gca,'position',[0.2 0.4 0.5 0.25]);
plot(time,retention_ratio,'linewidth',3)
xlim([0 30]);% ʱ�䶼Ҫ��0ʱ�̿�ʼ��
ylabel({'Retention ratio','(a.u.)'});
legend(num2str(shot));
legend('boxoff');
set(get(gca,'xlabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(get(gca,'ylabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(gca,'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(gca,'XTickLabel','','linewidth',3);
set(gca,'linewidth',3);
subplot(3,1,3)
set(gca,'position',[0.2 0.15 0.5 0.25]);
plot(time(1:end-1),retention_rate/1e21,time(1:end-1),puff_rate/1e21,'linewidth',3)
xlim([0 30]);
xlabel(('Time(s)'));
ylabel({'Retention rate','(10^{21}D/s)'});
legend('Retention rate');
set(get(gca,'xlabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(get(gca,'ylabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(gca,'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
set(gca,'linewidth',3);
% figure(5)
% a=2*Pam2P.*PH;
% plot(time,a.*pump_lowin,time,a.*pump_upin,time,a.*pump_lowout,time,a.*pump_upout,'linewidth',3);
% xlabel(('Time(s)'));
% ylabel({'Pump Rate','(s^{-1})'});
% legend('pump_{lowin}','pump_{upin}','pump_{lowout}','pump_{upout}');
% set(get(gca,'xlabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
% set(get(gca,'ylabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
% set(gca,'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
% set(gca,'linewidth',2);
% figure(6)
% plot(time,a.*pump_lowout,time,a.*pump_upout,'linewidth',3);
% xlabel(('Time(s)'));
% ylabel({'Pump Rate','(m^{3}s^{-1})'});
% legend('pump_{lowout}','pump_{upout}');
% set(get(gca,'xlabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
% set(get(gca,'ylabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
% set(gca,'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
% set(gca,'linewidth',2);


% figure(5)
% pump_rate=diff(PH.*pump).*1000;
% plot(time(1:end-1),pump_rate,'linewidth',3)
% xlabel(('Time(s)'));
% ylabel({'Pump Rate','(Pa m^{-3}s^{-1})'});
% set(get(gca,'xlabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
% set(get(gca,'ylabel'),'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
% set(gca,'Fontname','Times New Roman','FontSize',24,'Fontweight','bold');
% set(gca,'linewidth',3);


toc
