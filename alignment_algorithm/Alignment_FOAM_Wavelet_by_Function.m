% function [attitue vec] = opti_Align(pos, vel, imu, vec)
clear all
clc
%glvs
close all
global glv

% Defining the parameters
%=============================================================================
% Time related parameters
fs_gps=5;
fs_imu= 100;         % The gyroscope and accelerometer sampling rate is 100 Hz.


t=1/fs_imu;          % [t=0.0100] , update time of gyro output
T=1/fs_gps;          % [T=1.000] ,  update time of rotation of the geographical frame in second
k=T/t;         % [k=10] , update ratio of T and t



% Earth parameters
wie=7.292115e-5;       % rad/sec  omega of the earth
Re=6371000;            % radius of the earth
e = 0.0818191908425;   % eccentricity of WGS 84 model
R0 = 6378137.0;
Rp = 6356752.3142;     % polar radius
f = 1 / 298.257223563; % the flattening of the ellipsoid, f.
gn=[0 0 -1];           % NED navigation frame
meu=3.986004418*10^14;

% A degree of longitude is widest at the equator with a distance of 111.321 kilometers.
% Each degree of latitude is approximately 111 kilometers apart.
one_degLatLong = (2*pi*Re)/360;         % one degree latitute is approsimately equal
one_m_dis    = 1 /one_degLatLong;       %  degree lat or long in one meter
MonEnd=1;

% uploading the data and defining it by different variable
% uploading the data and defining it by different variable
% data_imu=importdata('dataIMU_unfilter.txt');

% data_imu=importdata('E:\Research Work\IMU data and Code\Working\data\dataIMU.txt');
% data_GPS=importdata('E:\Research Work\IMU data and Code\Working\data\dataGPS.txt');
data_imu=importdata('dataIMU.txt');
data_GPS=importdata('dataGPS.txt');

StrTime=49;
endTime=349;       % data time in second
num=endTime/T;

strtTime_imu=find(data_imu(:,1)==StrTime);
strtTime_gps=find(data_GPS(:,end)==StrTime);
endTime_imu=find(data_imu(:,1)==endTime);     % Number of sample data
endTime_gps=find(data_GPS(:,end)==endTime);      % Number of sample data


size(data_GPS)
Timu=data_imu(strtTime_imu:endTime_imu,1);                % Data time array in sec
Wib_b_nois=data_imu(strtTime_imu:endTime_imu,2:4);
Fib_b_nois=data_imu(strtTime_imu:endTime_imu,5:7);
DataAngle_Ref=data_GPS(strtTime_gps:endTime_gps,8:10);

% Wib_b=Wib_Filter(Wib_b_nois);
% Fib_b=Fib_Filter(Fib_b_nois);

Wib_b=Wib_b_nois;
Fib_b=Fib_b_nois;

% data_GPS=importdata('dataGPS_unfilter.txt');
Tgps=data_GPS(strtTime_gps:endTime_gps,7);                % Data time array in sec
Vn=data_GPS(strtTime_gps:endTime_gps+1,1:3);
% Vn0=data_GPS(strtTime_gps-1,1:3);
if(StrTime>=1),Vn0=data_GPS(strtTime_gps-1,1:3); end
% Vn0=data_imu(strtTime_imu-1+99,11:13) ;

if(StrTime<1), Vn0=[0 0 0]; end
lat=data_GPS(strtTime_gps:endTime_gps,4);
long=data_GPS(strtTime_gps:endTime_gps,5);
h=data_GPS(strtTime_gps:endTime_gps,6);

[Tval_gps,Tpos_gps]=intersect(Tgps,Timu);


% Initialization
qntm_1_no=[1 0 0 0];
qbtm_1_bo=[1 0 0 0];
delta_theta_1=[0 0 0];
alpha1=[0 0 0]';
delta_Vel_1=[0 0 0];
Beta_bar1=[0 0 0]';
B=zeros(3);
Z=[0 0 0]';
w=1;
% T=2*t;
ni=0;
ki=0;
K = zeros(4);
Xmat = @(xx)[0 -xx(3) xx(2);
    xx(3) 0 -xx(1);
    -xx(2) xx(1) 0];
Ib=[0 0 0*0.85]'; % x=0.000, y=0.000, z=0.850 m (x-right, y-fwd, z-up)
nn=2;
len=length(Wib_b);
% for m=1:1:endTime_gps-strtTime_gps




for k=1:nn:len-nn+1
    k1 = k+nn-1;
    Tsyn = Timu(k);
    
    
    
    %% Computing the time varing Body Frame Rotation qbt_bo
    %==================================================================
    
    % computing rotation vector
    delta_theta_1= Wib_b(k,1:3);
    delta_theta_2= Wib_b(k1,1:3);
    
    
    
    Wib_b2=[delta_theta_1 ; delta_theta_2];
    qbtm_bo = AttitudeUpdate(Wib_b2, qbtm_1_bo);
    Cbtm_1_bo = q2mat(qbtm_1_bo);
    qbtm_1_bo=qbtm_bo;
    
    
    %% Compute alpha using the gyroscope/accelerometer outputs by
    % ==========================================================================
    % computing angular increasement
    delta_Vel_1= Fib_b(k,1:3);
    delta_Vel_2=Fib_b(k1,1:3);
    
    Wib_b_X= Xmat(Wib_b(k,1:3));%+Wib_b((m-1)*k+n,1:3));
    Wib_b_X_0= Xmat(Wib_b(1,1:3));
    
    LeverArmTerm=(Cbtm_1_bo*Wib_b_X-Wib_b_X_0)*Ib;
    
    alphaMth=Cbtm_1_bo*(delta_Vel_1+delta_Vel_2+1/2*cross((delta_theta_1+delta_theta_2),(delta_Vel_1+delta_Vel_2))...
        +2/3*(cross(delta_theta_1,delta_Vel_2)+cross(delta_Vel_1,delta_theta_2)))'+LeverArmTerm;
    alpha2= alpha1+alphaMth;
    alpha1=alpha2;
    
    
    
    %%
    if mod(Tsyn,0.2)==0
        
        ki=ki+1;
        ni=Tpos_gps(ki);
        angle_qbt_bo(ki,:)= q2att(qbtm_bo');
        
        alphaMthArray(1:3,ki)=alphaMth;
        alpha(1:3,ki)=alpha2;%./norm(alpha2);
        alphaM=alpha2;%./norm(alpha2);
        alphaUnNorm(1:3,ki)=alpha2;
        
        %% Computing the time varing navigation Frame Rotation
        %=========================================================================
        % computing the earth parameter
        Wie_n=[0 wie*cos(lat(ki)) wie*sin(lat(ki))];         % ENU frame
        Rn_L=R0*(1-e^2)/(1-e^2*(sin(lat(ki)))^2)^(3/2);     % finding the earth radius north side
        Re_L=R0/(1-e^2*(sin(lat(ki)))^2)^(1/2);             % finding the earth radius east side
        
        Wen_n=[-Vn(ki,2)/(Rn_L+h(ki)) Vn(ki,1)/(Re_L+h(ki))  Vn(ki,1)*tan(lat(ki))/(Re_L+h(ki)) ];
        Win_n=Wie_n+Wen_n;
        
        go_L = 9.7803253359*(1 + 0.001931853*(sin(lat(ki)))^2)/sqrt(1-e^2*sin(lat(ki))^2)*gn';
        gbn_Lh=go_L*(1-2/R0*(1+f+wie^2*R0^2*Rp/meu)*h(ki)+3/R0^2*h(ki)^2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % My computed Win_n and gn is wrong
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        pos0=[lat(ni) long(ni) h(ni)];
        eth = earth(pos0,Vn(ni,:));
        Wen_n2=eth.wnen;
        Win_n2=eth.wnin';
        
        Win_nArray2(ki,1:3)=Win_n2;
        Wen_nArray2(ki,1:3)=Wen_n2;
        gbn_Lh=eth.gn;
        GnArray(ki,1:3)=gbn_Lh';
        
        
        Win_n_X= Xmat(Win_n2);
        Wie_n_X= Xmat(Wie_n);
        
        %%  computing the rotation vector of navigation frame
        qntm_no=AttitudeUpdate(Win_n2,qntm_1_no,T);
        Cntm_1_no = q2mat(qntm_1_no);
        Cntm_no = q2mat(qntm_no);
        qntm_1_no=qntm_no;
        
        qntm_noArray(ki,:)=qntm_no;
        angle_qno_nt(ki,:)= q2att(qntm_no');
        
        %% Compute beta using the aided velocity
        %=========================================================================
        
        Beta_barMth= Cntm_1_no*((T/2*eye(3)+T^2/6*Win_n_X)*Wie_n_X*Vn(ni,1:3)'...
            +(T/2*eye(3)+T^2/3*Win_n_X)*Wie_n_X*Vn(ni+1,1:3)'-(T*eye(3)+T^2/2*Win_n_X)*gbn_Lh);
        Beta_bar2= Beta_bar1+Beta_barMth;
        
        Beta_barMthArray(1:3,ki)=Cntm_no*Vn(ni,1:3)'-Vn0'+Beta_barMth;
        Beta2=Cntm_no*Vn(ni,1:3)'-Vn0'+ Beta_bar2;
        Beta_bar1=Beta_bar2;
        
        Beta(1:3,ki)=Beta2;%./norm(Beta2);
        BetaM=Beta2;%./norm(Beta2);
        BetaUnNorm(1:3,ki)=Beta2;
        
        
        %%   Computing the constant body to navigation frame Cbn0
        %=========================================================================
        
        
        qbn0= QMethod(alpha,Beta);
        
        qbn0Array(ki,:)=qbn0;
        angle_qbn0(ki,:) = q2att(qbn0);
        Cbn0=quat2dcm(qbn0);
        ErrCbn0(1:3,ki)=Beta(1:3,ki)-Cbn0*alpha(1:3,ki);
        ErrCbn0UnNorm(1:3,ki)=BetaUnNorm(1:3,ki)-Cbn0*alphaUnNorm(1:3,ki);
        Cbn0Alpha(1:3,ki)=Cbn0*alpha(1:3,ki);
        
        
        %% Computing Product chain Rule
        %=========================================================================
        Qbn=quatmultiply(quatconj(qntm_no),quatmultiply(qbn0,qbtm_bo));   % quaternion multiplication
        angle_Qbn(ki,:) = q2att(Qbn');
        
        
        % computing the error
        angle_qbn_Ref_rad(ki,:)=DataAngle_Ref(ki,:);
        qbn_Ref=a2qua(DataAngle_Ref(ki,:));
        DataAngle_RefCheck(ki,:)=q2att(qbn_Ref);
        
        ErrDownCode(ki,:)=qq2phi(Qbn',qbn_Ref)*180/pi;
        
        dQ = quatmultiply(qbn_Ref', quatconj(Qbn));
        dQang = q2att(dQ);
        Error_quat(ki,:)=dQang'*180/pi;
        
    end
    
end


% WahbaEq=1/2*norm(sum(ErrCbn0'.^2));
angle_qbn0Deg=angle_qbn0*180/pi;
angle_qbn_Ref=DataAngle_RefCheck*180/pi;
Err=(angle_qbn_Ref-angle_Qbn*180/pi);
angle_QbnDeg=angle_Qbn*180/pi;

TimeAxisIMU=(1:1:endTime_imu-strtTime_imu)*t;
TimeAxisGPS=(1:1:endTime_gps-strtTime_gps)*T;

angle_qbt_bo_deg=angle_qbt_bo*180/pi;
%% ploting the graphs
% =========================================================================

figure
plot(Error_quat)

size(TimeAxisGPS)
size(Error_quat)
size(TimeAxisIMU)


ResultData=[angle_Qbn*180/pi angle_qbn_Ref Error_quat ];
% save ResultData.txt ResultData -ascii -double
save ResultDataFilter.txt ResultData -ascii -double

figure, plot(TimeAxisGPS,Error_quat,'LineWidth',2), title('Error in degree'); grid on, xlabel('Time /s'); ylabel('Error using quat in degree'),legend('Roll','Pitch','yaw'), %axis([0,t*nn_gps,-0.05,0.05])

figure,
subplot(3,1,1), plot(TimeAxisGPS,angle_qbt_bo*180/pi,'LineWidth',2), grid on, title({'All three decomposed angles computeded by Quaternion' ; 'qbo-bt'}); xlabel('Time /s'); ylabel('Angle /degree'),legend('Roll','Pitch','yaw')
subplot(3,1,2), plot(TimeAxisGPS,angle_qbn0*180/pi,'LineWidth',2),   grid on, title('angle-qbn0'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Roll','Pitch','yaw')
subplot(3,1,3), plot(TimeAxisGPS,angle_qno_nt*180/pi,'LineWidth',2), grid on, title('qno-nt'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Roll','Pitch','yaw')
iii
figure, plot(TimeAxisGPS, angle_Qbn*180/pi,'LineWidth',2), grid on, title('Qbn');
figure, plot(TimeAxisGPS,angle_qbn_Ref,'LineWidth',2), grid on,title('angle_qbn_Ref');



%%  Ploting the graph of all three attitude
figure, plot(TimeAxisGPS,angle_qbt_bo*180/pi,'LineWidth',2), grid on, title('Body frame attitude'); xlabel('Time /s'); ylabel('Angle /degree'), legend('Roll','Pitch','yaw')
figure, plot(TimeAxisGPS,angle_qbn0*180/pi,'LineWidth',2),   grid on, title('angle-qbn0'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Roll','Pitch','yaw'),% axis([0,t*nn_gps,-0.3,0.3])
figure, plot(TimeAxisGPS,angle_qno_nt*180/pi,'LineWidth',2), grid on, title('qno-nt'); xlabel('Time /s'); ylabel('Angle /degree'), legend('Roll','Pitch','yaw')


figure
subplot(311), plot(TimeAxisGPS,angle_qbt_bo*180/pi,'LineWidth',2), grid on, title({'All three decomposed angles computeded by Quaternion' ; 'qbo-bt'}); xlabel('Time /s'); ylabel('Angle /degree'),legend('Roll','Pitch','yaw')
subplot(312), plot(TimeAxisGPS,angle_qbn0*180/pi,'LineWidth',2),   grid on,title('angle-qbn0'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Roll','Pitch','yaw')
subplot(313), plot(TimeAxisGPS,angle_qno_nt*180/pi,'LineWidth',2), grid on, title('qno-nt');xlabel('Time /s'); ylabel('Angle /degree'),legend('Roll','Pitch','yaw')


%%   Computing the degree error
figure, plot(TimeAxisGPS,Error_quat,'LineWidth',2),  grid on, title('Error in degree');      xlabel('Time /s'); ylabel('Error using quat in degree'), legend('Roll','Pitch','yaw'),% axis([0,t*nn_gps,-0.05,0.15])
figure, plot(TimeAxisGPS,Err,'LineWidth',2),         grid on, title('Error using Eular in degree'); xlabel('Time /s'); ylabel('Error in degree'),legend('Roll','Pitch','yaw')
figure, plot(TimeAxisGPS,ErrDownCode,'LineWidth',2), grid on, title('Error using DownCode in degree'); xlabel('Time /s'); ylabel('Error in degree'),legend('Roll','Pitch','yaw')
% iii
figure,
subplot(3,1,1), plot(TimeAxisGPS,Error_quat(:,1),'LineWidth',2), grid on, title({'Error in degree' ; 'Roll angle '}); xlabel('Time /s'); ylabel('Error in degree'), %axis([0,t*nn_gps,-0.05,0.5])
subplot(3,1,2), plot(TimeAxisGPS,Error_quat(:,2),'LineWidth',2), grid on, title('Pitch angle'); xlabel('Time /s'); ylabel('Error in degree'),%axis([0,t*nn_gps,-0.05,0.2])
subplot(3,1,3), plot(TimeAxisGPS,Error_quat(:,3),'LineWidth',2), grid on, title('Yaw angle '); xlabel('Time /s'); ylabel('Error in degree'), %axis([0,t*nn_gps,-0.05,0.1])


figure
subplot(3,1,1), plot(TimeAxisGPS,Err(:,1),'LineWidth',2), grid on, title({'Error in degree ' ; 'Roll angle '}); xlabel('Time /s'); ylabel('Error in degree '),%axis([0,t*nn_gps,-0.1,0.1])
subplot(3,1,2), plot(TimeAxisGPS,Err(:,2),'LineWidth',2), grid on, title('Pitch angle'); xlabel('Time /s'); ylabel('Error in degree ')%,axis([0,t*nn_gps,-0.15,0.15])
subplot(3,1,3), plot(TimeAxisGPS,Err(:,3),'LineWidth',2), grid on, title('Yaw angle '); xlabel('Time /s'); ylabel('Error in degree ')%, axis([0,t*nn_gps,-0.15,0.15])


%%  comparasion b/w true attitude and computed attitude Qbn
figure,
subplot(3,1,1), plot(TimeAxisGPS,angle_Qbn(:,1)*180/pi,'r','LineWidth',2), hold on, plot(TimeAxisGPS,angle_qbn_Ref(:,1),'b--','LineWidth',2)
grid on, title({'comparasion b/w true attitude and computed attitude by quat Qbn' ; 'Roll angle Qbn'}); xlabel('Time /s'); ylabel('Angle /degree'), legend('Estimated Attitude quat ','True Attitude')

subplot(3,1,2), plot(TimeAxisGPS,angle_Qbn(:,2)*180/pi,'r','LineWidth',2), hold on, plot(TimeAxisGPS,angle_qbn_Ref(:,2),'b--','LineWidth',2)
grid on,title('Pitch angle bn'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Estimated Attitude quat ','True Attitude')

subplot(3,1,3), plot(TimeAxisGPS,angle_Qbn(:,3)*180/pi,'r','LineWidth',2),hold on,plot(TimeAxisGPS,angle_qbn_Ref(:,3),'b--','LineWidth',2)
grid on,title('Yaw angle Qbn'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Estimated Attitude quat ','True Attitude')


figure, plot(TimeAxisGPS,angle_Qbn(:,1)*180/pi,'r','LineWidth',2),hold on, plot(TimeAxisGPS,angle_qbn_Ref(:,1),'b--','LineWidth',2)
grid on, title({'comparasion b/w true attitude and computed attitude by quat Qbn' ; 'Roll angle Qbn'}); xlabel('Time /s'); ylabel('Angle /degree'),legend('Estimated Attitude quat ','True Attitude')

figure, plot(TimeAxisGPS,angle_Qbn(:,2)*180/pi,'r','LineWidth',2), hold on,plot(TimeAxisGPS,angle_qbn_Ref(:,2),'b--','LineWidth',2)
grid on,title('Pitch angle bn'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Estimated Attitude quat ','True Attitude')

figure, plot(TimeAxisGPS,angle_Qbn(:,3)*180/pi,'r','LineWidth',2), hold on,plot(TimeAxisGPS,angle_qbn_Ref(:,3),'b--','LineWidth',2)
hold on,grid on, title('Yaw angle Qbn'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Estimated Attitude quat ','True Attitude')


figure
plot(Wib_b)
figure
plot(Fib_b)
%% Graph of  simulated data
for i=1:1
    %% Sensor Data
    
    figure
    subplot(311), plot(TimeAxisIMU,Wib_b(:,1),'LineWidth',2), grid on,     title('Angular velocity '); xlabel('Time /s');      ylabel('rad /s')
    subplot(312), plot(TimeAxisIMU,Wib_b(:,2),'LineWidth',2), grid on,     title('Angular velocity '); xlabel('Time /s');      ylabel('rad /s')
    subplot(313), plot(TimeAxisIMU,Wib_b(:,3),'LineWidth',2), grid on,     title('Angular velocity '); xlabel('Time /s');      ylabel('rad /s')
    
    figure
    subplot(311), plot(TimeAxisIMU,Fib_b(:,1),'LineWidth',2), grid on, title('Specific Force'); xlabel('Time /s');  ylabel('m/s^2')
    subplot(312), plot(TimeAxisIMU,Fib_b(:,2),'LineWidth',2), grid on, title('Specific Force'); xlabel('Time /s');  ylabel('m/s^2')
    subplot(313), plot(TimeAxisIMU,Fib_b(:,3),'LineWidth',2), grid on, title('Specific Force'); xlabel('Time /s');  ylabel('m/s^2')
    
    % ploting the graph gps output
    figure
    subplot(3,1,1), plot(TimeAxisGPS,long*180/pi,'LineWidth',2); grid on;  title('Longtitude-time'); xlabel('time/s'); ylabel('Longitude/degree');
    subplot(3,1,2), plot(TimeAxisGPS,lat*180/pi,'LineWidth',2);  grid on;  title('latitude-time');   xlabel('time/s'); ylabel('latitude/degree');
    subplot(3,1,3), plot(long*180/pi,lat*180/pi,'LineWidth',2);  grid on;  title('Longtitude-latitude'); xlabel('Long/deg');  ylabel('lati/deg');
    
    figure
    subplot(1,2,1), plot(long(1)*180/pi,lat(1)*180/pi, 'rp'); hold on,  plot(long*180/pi,lat*180/pi);           grid on;    title('Lat-Long in degree'); xlabel('lat / deg');    ylabel('long / deg');
    subplot(1,2,2), plot(0,0, 'rp'); hold on, plot((long(:)-long(1))*Re*cos(lat(1)), (lat(:)-lat(1))*Re);grid on;   title('Distance East-North in meter'); xlabel('East / meter');    ylabel('North / meter');
    
    figure,  plot(TimeAxisGPS,h);                  grid on;    xlabel('time/s'), ylabel('height/m'), title('height-time');
    
    
    figure,  plot(TimeAxisGPS,Vn(1:length(h),1:3),'LineWidth',2),grid on;       xlabel('Time / s');    ylabel('m / s'),   legend('East direction', 'North direction','Up direction')
    figure,  plot(TimeAxisGPS,angle_qbn_Ref), grid on, title('True attitude');  xlabel('Time /s');   ylabel('Angle /degree'),    legend('Roll','Pitch','yaw')
    
    %
end

