% This algorithm is introduced by YUANXIN WU and XIANFEI PAN in Research
% paper named "Velocity/Position Integration

% This algorithm is implemented by Mahmood ul hassan in 2-2019.
% Devenport Q Method
clc,clear all;
close all;
% format long

glvs
global glv

% Defining the parameters
%=============================================================================
% Time related parameters
fs_imu= 100;         % The gyroscope and accelerometer sampling rate is 100 Hz.
t=1/fs_imu;          % [t=0.0100] , update time of gyro output
fs_gps=5;
T=1/fs_gps;          % [T=1.000] ,  update time of rotation of the geographical frame in second
k=T/t;         % [k=10] , update ratio of T and t
time=300+20;       % data time in second
nn_imu=fs_imu*time;      % Number of sample data
nn_gps=fs_gps*time;      % Number of sample data
num=time/T;

StrTime=50;
strtTime_imu=(fs_imu*StrTime+1);
strtTime_gps=(fs_gps*StrTime+1);

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
% data_imu=importdata('dataIMU_unfilter.txt');
data_imu=importdata('dataIMU.txt');
Timu=data_imu(strtTime_imu:nn_imu,1);                % Data time array in sec
Wib_b_nois=data_imu(strtTime_imu:nn_imu,2:4);
Fib_b_nois=data_imu(strtTime_imu:nn_imu,5:7);
DataAngle_Ref=data_imu(strtTime_imu:nn_imu,8:10);

%Wib_b=Wib_Filter(Wib_b_nois);
%Fib_b=Fib_Filter(Fib_b_nois);

 Wib_b=Wib_b_nois;
 Fib_b=Fib_b_nois;

% data_GPS=importdata('dataGPS_unfilter.txt');
data_GPS=importdata('dataGPS.txt');
Tgps=data_GPS(strtTime_gps:nn_gps,7);                % Data time array in sec
Vn=data_GPS(strtTime_gps:nn_gps+1,1:3);
% Vn0=data_GPS(strtTime_gps-1,1:3);
if(StrTime>=1),Vn0=data_GPS(strtTime_gps-1,1:3); end
% Vn0=data_imu(strtTime_imu-1+99,11:13) ;

if(StrTime<1), Vn0=[0 0 0]; end
lat=data_GPS(strtTime_gps:nn_gps,4);
long=data_GPS(strtTime_gps:nn_gps,5);
h=data_GPS(strtTime_gps:nn_gps,6);


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
K = zeros(4);
Xmat = @(xx)[0 -xx(3) xx(2);
    xx(3) 0 -xx(1);
    -xx(2) xx(1) 0];
kkk=0;
Ib=[0 0 0*0.85]'; %x=0.000, y=0.000, z=0.850 m (x-right, y-fwd, z-up)


for m=1:1:num-strtTime_gps+1
        ni=ni+1;
 if(m<=0)
     ni=1;
 end

 
for n=2:2:k
    
    
    %% Computing the time varing Body Frame Rotation qbt_bo
    %==================================================================
    
    % computing rotation vector
    delta_theta_1= Wib_b((m-1)*k+n-1,1:3);
    delta_theta_2= Wib_b((m-1)*k+n,1:3);
    delta_theta_2Array(ni,:)=delta_theta_2;
    delta_theta_1Array(ni,:)=delta_theta_1;
    
    zeta_b=delta_theta_1+delta_theta_2+cross(2/3*delta_theta_1,delta_theta_2);   %Rotation Vector
    zeta_bArray(ni,1:3)=zeta_b;
    
    % Quat update
    qbtm_btm_1=([cos(norm(0.5*zeta_b)) sin(norm(0.5*zeta_b))/norm(0.5*zeta_b)*0.5*zeta_b]);
    qbtm_bo=quatmultiply(qbtm_1_bo,qbtm_btm_1);
    qbtm_bo = qbtm_bo/sqrt(qbtm_bo*qbtm_bo');
    
    Cbtm_1_bo = q2mat(qbtm_1_bo);
    qbtm_1_bo=qbtm_bo;
    
    
    %% Compute alpha using the gyroscope/accelerometer outputs by
    %==========================================================================
    % computing angular increasement
    delta_Vel_1= Fib_b((m-1)*k+n-1,1:3);
    delta_Vel_2=Fib_b((m-1)*k+n,1:3);
    
    Wib_b_X= Xmat(Wib_b((m-1)*k+n-1,1:3));%+Wib_b((m-1)*k+n,1:3));
    Wib_b_X_0= Xmat(Wib_b(1,1:3));
    
    LeverArmTerm=(Cbtm_1_bo*Wib_b_X-Wib_b_X_0)*Ib;
    
    alphaMth=Cbtm_1_bo*(delta_Vel_1+delta_Vel_2+1/2*cross((delta_theta_1+delta_theta_2),(delta_Vel_1+delta_Vel_2))...
        +2/3*(cross(delta_theta_1,delta_Vel_2)+cross(delta_Vel_1,delta_theta_2)))'+LeverArmTerm;
    alpha2= alpha1+alphaMth;
    alpha1=alpha2;
    
    alphaMthArray(1:3,m)=alphaMth;
    alpha(1:3,m)=alpha2;%./norm(alpha2);
    alphaM=alpha2;%./norm(alpha2);
    alphaUnNorm(1:3,m)=alpha2;
end
        angle_qbt_bo(ni,:)= q2att(qbtm_bo');
% Timu((m-1)*k+n)
% Tgps(m)
    %% Computing the time varing navigation Frame Rotation
    %=========================================================================
    % computing the earth parameter
    Wie_n=[0 wie*cos(lat(m)) wie*sin(lat(m))];         % ENU frame
    Rn_L=R0*(1-e^2)/(1-e^2*(sin(lat(m)))^2)^(3/2);     % finding the earth radius north side
    Re_L=R0/(1-e^2*(sin(lat(m)))^2)^(1/2);             % finding the earth radius east side
    
    Wen_n=[-Vn(m,2)/(Rn_L+h(m)) Vn(m,1)/(Re_L+h(m))  Vn(m,1)*tan(lat(m))/(Re_L+h(m)) ];
    Win_n=Wie_n+Wen_n;
    
    go_L = 9.7803253359*(1 + 0.001931853*(sin(lat(m)))^2)/sqrt(1-e^2*sin(lat(m))^2)*gn';
    gbn_Lh=go_L*(1-2/R0*(1+f+wie^2*R0^2*Rp/meu)*h(m)+3/R0^2*h(m)^2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % My computed Win_n and gn is wrong
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pos0=[lat(m) long(m) h(m)];
    eth = earth(pos0,Vn(m,:));
    Wen_n2=eth.wnen;
    Win_n2=eth.wnin';
    Win_nArray2(ni,1:3)=Win_n2;
    Wen_nArray2(ni,1:3)=Wen_n2;
    gbn_Lh=eth.gn;
    GnArray(ni,1:3)=gbn_Lh';
    
    
    Win_n_X= Xmat(Win_n2);
    Wie_n_X= Xmat(Wie_n);
    
    %%  computing the rotation vector of navigation frame
    
    zeta_n=Win_n2*T;
    zeta_n_X=Xmat(zeta_n);
    
    qntm_ntm_1=([cos(norm(0.5*zeta_n)) sin(norm(0.5*zeta_n))/norm(0.5*zeta_n)*0.5*zeta_n]);
    qntm_no=quatmultiply(qntm_1_no, qntm_ntm_1);   % quaternion updating
    qntm_no = qntm_no/sqrt(qntm_no*qntm_no');
    angle_qno_nt(ni,:)= q2att(qntm_no');
    
    Cntm_1_no = q2mat(qntm_1_no);
    Cntm_no = q2mat(qntm_no);
    
            qntm_noArray(ni,:)=qntm_no;

    qntm_1_no=qntm_no;
    
    %% Compute beta using the aided velocity
    %=========================================================================
    
    Beta_barMth= Cntm_1_no*((T/2*eye(3)+T^2/6*Win_n_X)*Wie_n_X*Vn(m,1:3)'...
        +(T/2*eye(3)+T^2/3*Win_n_X)*Wie_n_X*Vn(m+1,1:3)'-(T*eye(3)+T^2/2*Win_n_X)*gbn_Lh);
    Beta_bar2= Beta_bar1+Beta_barMth;
    
    Beta_barMthArray(1:3,m)=Cntm_no*Vn(m,1:3)'-Vn0'+Beta_barMth;
    Beta2=Cntm_no*Vn(m,1:3)'-Vn0'+ Beta_bar2;
    Beta_bar1=Beta_bar2;
    
    Beta(1:3,m)=Beta2;%./norm(Beta2);
    BetaM=Beta2;%./norm(Beta2);
    BetaUnNorm(1:3,m)=Beta2;
    
    
    %%   Computing the constant body to navigation frame Cbn0
    %=========================================================================
    if(m>=0)
        %     computing 4X4 K matrix
             
%             B=(w*alpha(1:3,1:m))*Beta(1:3,1:m)';
%             Z=sum(cross(w*alpha(1:3,1:m)',Beta(1:3,1:m)',2),1) ;
%             S=B+B';
%             K=[trace(B) Z
%                 Z' S-trace(B)*eye(3)    ];
%         
%         %     computing the biggest eign value
%             [V,D] = eig(K,'nobalance');
%             [d, id] = max(diag(D));
%             qbn0  = V(:,id)';
%         qbn0 = qbn0/sqrt(qbn0*qbn0');
%         
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         dM = rq2m([0;alphaM])-lq2m([0;BetaM]);
%                 K = 0.99991*K + dM'*dM*T;
%                 [V, D] = eig(K);
%                  qbn0 = V(:,1)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      
%         
%         dM= (lq2m([0;BetaM])-rq2m([0;alphaM]))'*(lq2m([0;BetaM])-rq2m([0;alphaM]))*T;
%         K=K+dM;
%         [V, D] = eig(K);
%         [d, id] = min(diag(D));
%         qbn0  = V(:,id)';
%         qbn0 = qbn0/sqrt(qbn0*qbn0');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        qp=Beta(1:3,1:m);
        pp=alpha(1:3,1:m);
        B=qp*pp';
        
        detB=det(B);
        froBsq=norm(B, 'fro')^2;
        adjB=det(B)*inv(B);
        froadjBsq=norm(adjB, 'fro')^2;
        
        lam=0.5*(sum(sum(qp.^2)) + sum(sum(pp.^2))); % 0.5*(trace(qp*qp') + trace(pp*pp'));
        lamprev=0.0;
        
        % compute lam_max with Newton-Raphson
        while abs((lam-lamprev)/lam) >= 1E-12
            lamprev=lam;
            tmp=lam^2-froBsq;
            lam=lam - (tmp^2 - 8*lam*detB - 4*froadjBsq)/(4*tmp*lam - 8*detB);
        end
        
        % note that adj(B')=adj(B)'
        R=((lam^2 + froBsq)*B + 2*lam*adjB' - 2*B*B'*B)/(lam*(lam^2-froBsq) - 2*detB);
        qbn0=m2qua(R)';


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        qbn0Array(ni,:)=qbn0;
        angle_qbn0(ni,:) = q2att(qbn0);
        Cbn0=quat2dcm(qbn0);
        ErrCbn0(1:3,ni)=Beta(1:3,ni)-Cbn0*alpha(1:3,ni);
        ErrCbn0UnNorm(1:3,ni)=BetaUnNorm(1:3,ni)-Cbn0*alphaUnNorm(1:3,ni);
        Cbn0Alpha(1:3,ni)=Cbn0*alpha(1:3,ni);
        
        
        %% Computing Product chain Rule
        %=========================================================================
        Qbn=quatmultiply(quatconj(qntm_no),quatmultiply(qbn0,qbtm_bo));   % quaternion multiplication
        angle_Qbn(ni,:) = q2att(Qbn');
        
        
        km=k*(m);
        
        % computing the error
        angle_qbn_Ref_rad(ni,:)=DataAngle_Ref(km,:);
        qbn_Ref=a2qua(DataAngle_Ref(km,:));
        DataAngle_RefCheck(ni,:)=q2att(qbn_Ref);
        
        ErrDownCode(ni,:)=qq2phi(Qbn',qbn_Ref)*180/pi;
        
        dQ = quatmultiply(qbn_Ref', quatconj(Qbn));
        dQang = q2att(dQ);
        Error_quat(ni,:)=dQang'*180/pi;
        
    end
    
    
    
    
    
end
WahbaEq=1/2*norm(sum(ErrCbn0'.^2));
angle_qbn0Deg=angle_qbn0*180/pi;
angle_qbn_Ref=DataAngle_RefCheck*180/pi;
Err=(angle_qbn_Ref-angle_Qbn*180/pi);
angle_QbnDeg=angle_Qbn*180/pi;

TimeAxisIMU=(0:1:nn_imu-strtTime_imu)*t;
TimeAxisGPS=(0:1:nn_gps-strtTime_gps)*T;


%% ploting the graphs
% =========================================================================

size(TimeAxisGPS)
size(Error_quat)

ResultData=[angle_Qbn*180/pi angle_qbn_Ref Error_quat];
% save ResultData.txt ResultData -ascii -double
save ResultDataFilter.txt ResultData -ascii -double

figure, plot(TimeAxisGPS,Error_quat,'LineWidth',2), title('Error in degree'); grid on, xlabel('Time /s'); ylabel('Error using quat in degree'),legend('Roll','Pitch','yaw'), %axis([0,t*nn_gps,-0.05,0.05])
iii
figure,
subplot(3,1,1), plot(TimeAxisGPS,angle_qbt_bo*180/pi), grid on, title({'All three decomposed angles computeded by Quaternion' ; 'qbo-bt'}); xlabel('Time /s'); ylabel('Angle /degree'),legend('Roll','Pitch','yaw')
subplot(3,1,2), plot(TimeAxisGPS,angle_qbn0*180/pi),   grid on, title('angle-qbn0'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Roll','Pitch','yaw')
subplot(3,1,3), plot(TimeAxisGPS,angle_qno_nt*180/pi), grid on, title('qno-nt'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Roll','Pitch','yaw')

figure, plot(TimeAxisGPS, angle_Qbn*180/pi), grid on, title('Qbn'); 
figure, plot(TimeAxisGPS,angle_qbn_Ref), grid on,title('angle_qbn_Ref');



%%  Ploting the graph of all three attitude
figure, plot(TimeAxisGPS,angle_qbt_bo*180/pi), grid on, title('Body frame attitude'); xlabel('Time /s'); ylabel('Angle /degree'), legend('Roll','Pitch','yaw')
figure, plot(TimeAxisGPS,angle_qbn0*180/pi),   grid on, title('angle-qbn0'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Roll','Pitch','yaw'),% axis([0,t*nn_gps,-0.3,0.3])
figure, plot(TimeAxisGPS,angle_qno_nt*180/pi), grid on, title('qno-nt'); xlabel('Time /s'); ylabel('Angle /degree'), legend('Roll','Pitch','yaw')


figure
subplot(311), plot(TimeAxisGPS,angle_qbt_bo*180/pi), grid on, title({'All three decomposed angles computeded by Quaternion' ; 'qbo-bt'}); xlabel('Time /s'); ylabel('Angle /degree'),legend('Roll','Pitch','yaw')
subplot(312), plot(TimeAxisGPS,angle_qbn0*180/pi),   grid on,title('angle-qbn0'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Roll','Pitch','yaw')
subplot(313), plot(TimeAxisGPS,angle_qno_nt*180/pi), grid on, title('qno-nt');xlabel('Time /s'); ylabel('Angle /degree'),legend('Roll','Pitch','yaw')


%%   Computing the degree error
figure, plot(TimeAxisGPS,Error_quat),  grid on, title('Error in degree');      xlabel('Time /s'); ylabel('Error using quat in degree'), legend('Roll','Pitch','yaw'),% axis([0,t*nn_gps,-0.05,0.15])
figure, plot(TimeAxisGPS,Err),         grid on, title('Error using Eular in degree'); xlabel('Time /s'); ylabel('Error in degree'),legend('Roll','Pitch','yaw')
figure, plot(TimeAxisGPS,ErrDownCode), grid on, title('Error using DownCode in degree'); xlabel('Time /s'); ylabel('Error in degree'),legend('Roll','Pitch','yaw')
iii
figure, 
subplot(3,1,1), plot(TimeAxisGPS,Error_quat(:,1)), grid on, title({'Error in degree' ; 'Roll angle '}); xlabel('Time /s'); ylabel('Error in degree'), %axis([0,t*nn_gps,-0.05,0.5])
subplot(3,1,2), plot(TimeAxisGPS,Error_quat(:,2)), grid on, title('Pitch angle'); xlabel('Time /s'); ylabel('Error in degree'),axis([0,t*nn_gps,-0.05,0.2])
subplot(3,1,3), plot(TimeAxisGPS,Error_quat(:,3)), grid on, title('Yaw angle '); xlabel('Time /s'); ylabel('Error in degree'), axis([0,t*nn_gps,-0.05,0.1])


figure
subplot(3,1,1), plot(TimeAxisGPS,Err(:,1)), grid on, title({'Error in degree ' ; 'Roll angle '}); xlabel('Time /s'); ylabel('Error in degree '),%axis([0,t*nn_gps,-0.1,0.1])
subplot(3,1,2), plot(TimeAxisGPS,Err(:,2)), grid on, title('Pitch angle'); xlabel('Time /s'); ylabel('Error in degree '),axis([0,t*nn_gps,-0.15,0.15])
subplot(3,1,3), plot(TimeAxisGPS,Err(:,3)), grid on, title('Yaw angle '); xlabel('Time /s'); ylabel('Error in degree '), axis([0,t*nn_gps,-0.15,0.15])


%%  comparasion b/w true attitude and computed attitude Qbn
figure, 
subplot(3,1,1), plot(TimeAxisGPS,angle_Qbn(:,1)*180/pi,'r'), hold on, plot(TimeAxisGPS,angle_qbn_Ref(:,1),'b--')
grid on, title({'comparasion b/w true attitude and computed attitude by quat Qbn' ; 'Roll angle Qbn'}); xlabel('Time /s'); ylabel('Angle /degree'), legend('Estimated Attitude quat ','True Attitude')

subplot(3,1,2), plot(TimeAxisGPS,angle_Qbn(:,2)*180/pi,'r'), hold on, plot(TimeAxisGPS,angle_qbn_Ref(:,2),'b--')
grid on,title('Pitch angle bn'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Estimated Attitude quat ','True Attitude')

subplot(3,1,3), plot(TimeAxisGPS,angle_Qbn(:,3)*180/pi,'r'),hold on,plot(TimeAxisGPS,angle_qbn_Ref(:,3),'b--')
grid on,title('Yaw angle Qbn'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Estimated Attitude quat ','True Attitude')


figure, plot(TimeAxisGPS,angle_Qbn(:,1)*180/pi,'r'),hold on, plot(TimeAxisGPS,angle_qbn_Ref(:,1),'b--')
grid on, title({'comparasion b/w true attitude and computed attitude by quat Qbn' ; 'Roll angle Qbn'}); xlabel('Time /s'); ylabel('Angle /degree'),legend('Estimated Attitude quat ','True Attitude')

figure, plot(TimeAxisGPS,angle_Qbn(:,2)*180/pi,'r'), hold on,plot(TimeAxisGPS,angle_qbn_Ref(:,2),'b--')
grid on,title('Pitch angle bn'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Estimated Attitude quat ','True Attitude')

figure, plot(TimeAxisGPS,angle_Qbn(:,3)*180/pi,'r'), hold on,plot(TimeAxisGPS,angle_qbn_Ref(:,3),'b--')
hold on,grid on, title('Yaw angle Qbn'); xlabel('Time /s'); ylabel('Angle /degree'),legend('Estimated Attitude quat ','True Attitude')



%% Graph of  simulated data
for i=1:1
    %% Sensor Data
    
    figure
    subplot(311), plot(TimeAxisIMU,Wib_b(:,1)), grid on,     title('Angular velocity '); xlabel('Time /s');      ylabel('rad /s')
    subplot(312), plot(TimeAxisIMU,Wib_b(:,2)), grid on,     title('Angular velocity '); xlabel('Time /s');      ylabel('rad /s')
    subplot(313), plot(TimeAxisIMU,Wib_b(:,3)), grid on,     title('Angular velocity '); xlabel('Time /s');      ylabel('rad /s')
       
    figure
    subplot(311), plot(TimeAxisIMU,Fib_b(:,1)), grid on, title('Specific Force'); xlabel('Time /s');  ylabel('m/s^2')
    subplot(312), plot(TimeAxisIMU,Fib_b(:,2)), grid on, title('Specific Force'); xlabel('Time /s');  ylabel('m/s^2')
    subplot(313), plot(TimeAxisIMU,Fib_b(:,3)), grid on, title('Specific Force'); xlabel('Time /s');  ylabel('m/s^2')
   
    % ploting the graph gps output
    figure
    subplot(3,1,1), plot(TimeAxisGPS,long*180/pi); grid on;  title('Longtitude-time'); xlabel('time/s'); ylabel('Longitude/degree');
    subplot(3,1,2), plot(TimeAxisGPS,lat*180/pi);  grid on;  title('latitude-time');   xlabel('time/s'); ylabel('latitude/degree');
    subplot(3,1,3), plot(long*180/pi,lat*180/pi);  grid on;  title('Longtitude-latitude'); xlabel('Long/deg');  ylabel('lati/deg');
    
      figure
subplot(1,2,1), plot(long(1)*180/pi,lat(1)*180/pi, 'rp'); hold on,  plot(long*180/pi,lat*180/pi);           grid on;    title('Lat-Long in degree'); xlabel('lat / deg');    ylabel('long / deg');
subplot(1,2,2), plot(0,0, 'rp'); hold on, plot((long(:)-long(1))*Re*cos(lat(1)), (lat(:)-lat(1))*Re);grid on;   title('Distance East-North in meter'); xlabel('East / meter');    ylabel('North / meter');

    figure,  plot(TimeAxisGPS,h);                  grid on;    xlabel('time/s'), ylabel('height/m'), title('height-time');
    

    figure,  plot(TimeAxisGPS,Vn(1:length(h),1:3),'LineWidth',2),grid on;       xlabel('Time / s');    ylabel('m / s'),   legend('East direction', 'North direction','Up direction')
    figure,  plot(TimeAxisGPS,angle_qbn_Ref), grid on, title('True attitude');  xlabel('Time /s');   ylabel('Angle /degree'),    legend('Roll','Pitch','yaw')

    %
end

