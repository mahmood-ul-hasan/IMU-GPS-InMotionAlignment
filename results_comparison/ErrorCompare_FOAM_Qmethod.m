ResultData=importdata('ResultData.txt');
ResultDataFilter=importdata('ResultDataFilter.txt');

% ResultData=[angle_Qbn*180/pi angle_qbn_Ref];

angle_Qbn=ResultData(:,1:3);
angle_qbn=ResultData(:,4:6);
Error=ResultData(:,7:9);
timeAttitude=ResultData(:,10);

angle_Qbn_Filter=ResultDataFilter(:,1:3);
angle_qbn_Ref_Filter=ResultDataFilter(:,4:6);
Error_Filter=ResultDataFilter(:,7:9);
timeAttitude_Filter=ResultDataFilter(:,10);

figure, plot(TimeAxisGPS,Error_quat,'LineWidth',2), title('Error in degree'); grid on, xlabel('Time /s'); ylabel('Error using quat in degree'),legend('Roll','Pitch','yaw'), %axis([0,t*nn_gps,-0.05,0.05])

size(Error)
size(Error_Filter)
size(TimeAxisGPS)

figure, 
subplot(3,1,1), plot(TimeAxisGPS,Error(:,1),'LineWidth',2), hold on, plot(TimeAxisGPS,Error_Filter(:,1),'r','LineWidth',2), grid on,  ylabel('Roll Error (deg)'), legend( 'Traditional OBA' , 'Proposed Method')%axis([0,t*nn_gps,-0.05,0.5])
subplot(3,1,2), plot(TimeAxisGPS,Error(:,2),'LineWidth',2), hold on, plot(TimeAxisGPS,Error_Filter(:,2),'r','LineWidth',2), grid on,  ylabel('Pitch Error (deg)')
subplot(3,1,3), plot(TimeAxisGPS,Error(:,3),'LineWidth',2), hold on, plot(TimeAxisGPS,Error_Filter(:,3),'r','LineWidth',2), grid on,  xlabel('Time /s'); ylabel('Yaw Error (deg)')


figure, plot(TimeAxisGPS,angle_qbn_Ref(:,3),'g','LineWidth',2), hold on, plot(TimeAxisGPS,angle_Qbn_Filter(:,3),'r','LineWidth',2), hold on, plot(TimeAxisGPS,angle_Qbn(:,3),'b','LineWidth',2), 
hold on,grid on, xlabel('Time /s'); ylabel('Yaw angle (deg)'), legend( 'True Attitude','Estimated by Proposed Method', 'Estimated by Traditional OBA' )

figure, plot(TimeAxisGPS*5,timeAttitude,'LineWidth',2), hold on, plot(TimeAxisGPS*5,timeAttitude_Filter,'r','LineWidth',2);  grid on, xlabel('Number of Observation vectors'); ylabel('Execution Time(s)'),legend('Q Method','FOAM'), %axis([0,t*nn_gps,-0.05,0.05])

