
clear all
clc
close all
Re=6371000;     % radius of the earth

glvs
ts = 0.01;       % sampling interval
% avp0 = avpset2([0;0;0], [ 0 0 0], glv.pos0); % init avp
avp0 = avpset2([0;0;0], [ 0 0 0], [131 21 0]); % init avp
% trajectory segment setting
xxx = [];
seg = trjsegment2(xxx, 'init',         0);

 seg = trjsegment2(seg, 'accelerate',   20, xxx, 0.50000);
seg = trjsegment2(seg, 'uniform',      5);
% % seg = trjsegment2(seg, 'coturnleft',   350, 1, xxx, 4);
% 
seg = trjsegment2(seg, 'turnleft',   350, 1, xxx, 2);



% seg = trjsegment2(seg, 'accelerate',   20, xxx, 0.500000);
% seg = trjsegment2(seg, 'uniform',      200);
% 
% seg = trjsegment2(seg, 'coturnleft',   20, 4.5, xxx, 4);
% seg = trjsegment2(seg, 'deaccelerate',   20, xxx,0.20000);
% seg = trjsegment2(seg, 'uniform',      80);
% seg = trjsegment2(seg, 'accelerate',   30, xxx, 0.10000);
% seg = trjsegment2(seg, 'uniform',      100);
% 
% seg = trjsegment2(seg, 'coturnleft',   20, 4.5, xxx, 4);
% seg = trjsegment2(seg, 'accelerate',   20, xxx, 0.500000);
% seg = trjsegment2(seg, 'uniform',      100);
% 
% seg = trjsegment2(seg, 'coturnleft',   20, 4.5, xxx, 4);
% seg = trjsegment2(seg, 'uniform',      100);




% seg = trjsegment2(seg, 'deaccelerate',   300, xxx, .10000);
% seg = trjsegment2(seg, 'uniform',      330);
% seg = trjsegment2(seg, 'turnright',   1, 90, 1, 1);



% generate, save & plot
trj = trjsimu2(avp0, seg.wat, ts, 1);
trjfile('trj10ms.mat', trj);

% gpsVnPos = gpssimu(avp, dvn, dpos, tau, lever, imu_delay, isplot)
% gpsVnPos = gpssimu2(trj.avp, 0, 0);
% imuerr = imuerrset(0.0, 0, 0, 00, 0,0,0,0, 0, 0);

% According Wu Youin 
% gpsVnPos = gpssimu2(trj.avp, 0.1, 2 );
% imuerr = imuerrset(0.01, 50, 0.1, 500, 0,1,0,0, 0, 0);
% imuerr = imuerrset(0.01, 50, 0.1, 500, 0,1,0,0, 0, 0);

f=1.25;
% f=1;


% With filter
% gpsVnPos = gpssimu2(trj.avp, 0.1/f, 0.5/f );
% imuerr = imuerrset(0.05/f, 50/f, 0.1/f, 500/f, 0,1,0,0, 0, 0);

% With out filter
gpsVnPos = gpssimu2(trj.avp, 0.1, 0 );
imuerr = imuerrset(0.01/f, 50/f, 0.01/f, 50/f, 0,1,0,0, 0, 0);

% imuerr = imuerrset(eb, db, web, wdb, sqrtR0G, TauG, sqrtR0A, TauA, dKGii, dKAii, dKGij, dKAij, KA2)
%     eb - gyro constant bias (deg/h)
%     db - acc constant bias (ug)
%     web - angular random walk (deg/sqrt(h))
%     wdb - velocity random walk (ug/sqrt(Hz))



imu = imuadderr(trj.imu, imuerr);


data=[trj.avp(:,10) trj.avp(:,1:3) trj.avp(:,4:6) trj.avp(:,7:9) trj.imu(:,1:3) trj.imu(:,4:6)];
save data.txt data -ascii -double

dataIMU=[imu(:,7) imu(:,1:3) imu(:,4:6) trj.avp(:,1:3) trj.avp(:,4:6) trj.avp(:,7:9)];
save dataIMU.txt dataIMU -ascii -double
dataGPS=gpsVnPos;
save dataGPS.txt dataGPS  -ascii -double

% trjfile('trj10ms.mat', trj);
insplot(trj.avp);
imuplot(trj.imu);
% pos2gpx('trj_SINS_gps', trj.avp(1:round(1/trj.ts):end,7:9)); % to Google Earth
% % profile viewer



%     angle_qbn_Ref=trj.avp(:,1:3);
%     Vn=trj.avp(:,4:6);
%     lat=trj.avp(:,7);
%     long=trj.avp(:,8);
%     h=trj.avp(:,9);
%     Wib_b=trj.imu(:,1:3)/ts;
%     Fib_b=trj.imu(:,4:6)/ts;
%     tt=trj.imu(:,7);


Tgps=dataGPS(:,7);                % Data time array in sec
Vn=dataGPS(:,1:3);
lat=dataGPS(:,4);
long=dataGPS(:,5);
h=dataGPS(:,6);

Timu=dataIMU(:,1);                % Data time array in sec
Wib_b=dataIMU(:,2:4)/ts;
Fib_b=dataIMU(:,5:7)/ts;
angle_qbn_Ref=dataIMU(:,8:10);



fs_imu= 100;         % The gyroscope and accelerometer sampling rate is 100 Hz.
t=1/fs_imu;          % [t=0.0100] , update time of gyro output
fs_gps=5;
T=1/fs_gps;          % [T=1.000] ,  update time of rotation of the geographical frame in second
k=T/t;         % [k=10] , update ratio of T and t
time=370;       % data time in second
nn_imu=fs_imu*time;      % Number of sample data
nn_gps=fs_gps*time;      % Number of sample data
num=time/T;


TimeAxisIMU=(0:1:nn_imu-1)*t;
TimeAxisGPS=(0:1:nn_gps-1)*T;

size(TimeAxisIMU)
size(Wib_b)

figure
plot(TimeAxisIMU,angle_qbn_Ref*180/pi),     grid on,    title('Attitude');     xlabel('Time /s');     ylabel('Angle /degree') ,    legend('Roll','Pitch','heading (Yaw)')

figure
subplot(311),  plot(TimeAxisIMU,angle_qbn_Ref(:,1)*180/pi),     grid on,    title('Attitude');     xlabel('Time /s');     ylabel('Angle /degree') ,    legend('Roll','Pitch','heading (Yaw)')


figure
subplot(311),    plot(TimeAxisIMU,Wib_b(:,1)),      grid on, title('Angular velocity Wx in rad'),    xlabel('Time /s'),    ylabel('rad /s')
subplot(312),    plot(TimeAxisIMU,Wib_b(:,2)),      grid on, title('Angular velocity Wy in rad'),    xlabel('Time /s'),    ylabel('rad /s')
subplot(313),    plot(TimeAxisIMU,Wib_b(:,3)),      grid on, title('Angular velocity Wz in rad'),    xlabel('Time /s'),    ylabel('rad /s')

figure
subplot(311),    plot(TimeAxisIMU,Wib_b(:,1)*180/pi),      grid on, title('Angular velocity Wx in deg'),    xlabel('Time /s'),    ylabel('deg /s')
subplot(312),    plot(TimeAxisIMU,Wib_b(:,2)*180/pi),      grid on, title('Angular velocity Wy in deg'),    xlabel('Time /s'),    ylabel('deg /s')
subplot(313),    plot(TimeAxisIMU,Wib_b(:,3)*180/pi  ),      grid on, title('Angular velocity Wz in deg'),    xlabel('Time /s'),    ylabel('deg /s')

figure
subplot(311),    plot(TimeAxisIMU,Fib_b(:,1)),      grid on,  title('Specific Force Fx'),    xlabel('Time / s'),    ylabel('m/s^2')
subplot(312),    plot(TimeAxisIMU,Fib_b(:,2)),      grid on,  title('Specific Force Fy'),    xlabel('Time / s'),    ylabel('m/s^2')
subplot(313),    plot(TimeAxisIMU,Fib_b(:,3)),      grid on,  title('Specific Force Fz'),    xlabel('Time / s'),    ylabel('m/s^2')

figure
plot(TimeAxisGPS,Vn,'LineWidth',2),     grid on,    title('Velocity in Navigation Frame ENU');     xlabel('Time / s');     ylabel('m / s'),     legend('East direction', 'North direction','Up direction')

figure
subplot(3,1,1), plot(TimeAxisGPS,long*180/pi);     grid on,    title('Long-time');    xlabel('time/s');     ylabel('Long / deg');
subplot(3,1,2), plot(TimeAxisGPS,lat*180/pi);      grid on;    title('lat-time');      xlabel('time/s');     ylabel('lat / deg ');
subplot(3,1,3), plot(TimeAxisGPS,h);               grid on;    title('height-time');        xlabel('time/s');     ylabel('height / m')
avp0(7:8)

  figure
subplot(1,2,1), plot(long(1)*180/pi,lat(1)*180/pi, 'rp'); hold on,  plot(long*180/pi,lat*180/pi);           grid on;    title('Lat-Long in degree'); xlabel('lat / deg');    ylabel('long / deg');
subplot(1,2,2), plot(0,0, 'rp'); hold on, plot((long(:)-long(1))*Re*cos(lat(1)), (lat(:)-lat(1))*Re);grid on;   title('Distance East-North in meter'); xlabel('East / meter');    ylabel('North / meter');


figure
plot(long(1)*180/pi,lat(1)*180/pi, 'rp'); hold on,  plot(long*180/pi,lat*180/pi ,'LineWidth',2);           grid on;    title('Trajectory of the vehicle'); xlabel('Lat / deg');    ylabel('Long / deg');
