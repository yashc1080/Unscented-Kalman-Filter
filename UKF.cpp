#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <cmath>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

int n_aug = 7;
double glambda;
double weight1;
double weightn;
double dt;
VectorXd x = VectorXd(5);//state vector
VectorXd x_aug(7);//augmented vector
MatrixXd P_aug(7,7);//augmented covariance matrix
MatrixXd P = MatrixXd(5, 5);//covariance matrix
MatrixXd sigma_aug(7,15);//augmented sigma point matrix
MatrixXd predict_sigma_aug(5,15);//predicted augmented sigma point matrix
MatrixXd measured_predict_sigma_aug(3,15);// predicted sigma points in measurement space
VectorXd stateMean(5);
MatrixXd stateCovariance(5,5);
VectorXd measurementMean(3);
MatrixXd measurementCovariance(3,3);
MatrixXd crossCorrelationMatrix(5,3);
MatrixXd kalmanGain;
VectorXd z(3);
int ccount = 0;
float m_rho,m_pi,m_rdot;
 //Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a = 0.2;

  //Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd = 0.2;
 
MatrixXd  R_radar_(3,3);

 // Radar measurement noise standard deviation radius in m
 double std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
 double std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
 double  std_radrd_ = 0.3;

    
void calculateError(VectorXd state,VectorXd groundTruth){
    double error;
    for(int i=0;i<3;i++){
        error += abs(state(i)) - abs(groundTruth(i));
    }
   //print result
  cout<<"-------------------------------------------------------------"<<endl;
  cout << "Updated state x: " << endl << x << endl;
  cout << "Updated state covariance P: " << endl << P << endl;
  cout << "Error : " << endl << error << endl;
  cout << "Ground Truth: " << endl << groundTruth << endl;
   cout<<"-------------------------------------------------------------"<<endl;
}

void updateState(){
  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
   double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }


  MatrixXd Xsig_pred = predict_sigma_aug;
  MatrixXd Zsig = measured_predict_sigma_aug;
  VectorXd z_pred = measurementMean;
  MatrixXd S = MatrixXd(n_z,n_z);
 S.fill(0.0);

  for (int i = 0; i < 2*n_aug +1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
   while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R_radar_;
  //create example vector for incoming radar measurement

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x = stateMean + K * z_diff;
  P = stateCovariance - K*S*K.transpose();

}


VectorXd radar2State(VectorXd m_radar){
   double r_m_rho = m_radar(0);
   double r_m_pi = m_radar(1);
   double r_m_rdot = m_radar(2);
   double r_m_px = r_m_rho*cos(r_m_pi);
        if ( r_m_px < 0.0001 ) {
        r_m_px = 0.0001;
      }
        double r_m_py = r_m_rho*sin(r_m_pi);
        if(r_m_py< 0.0001){
          r_m_py = 0.0001;
        }  
        double r_m_vx = r_m_rdot*cos(r_m_pi);
        double r_m_vy = r_m_rdot*sin(r_m_pi);
        double r_m_v  = sqrt(r_m_vx*r_m_vx + r_m_vy*r_m_vy);
  VectorXd r_x(5);
  r_x << r_m_px,r_m_py,r_m_v,0,0;
  return r_x;
}

void NormalizeAngleOnComponent(VectorXd vector, int index) {
  while (vector(index)> M_PI) vector(index)-=2.*M_PI;
  while (vector(index)<-M_PI) vector(index)+=2.*M_PI;
}


void convertIntoMeasurementSpace(){
    double rho,theta,r_rho;
    for(int i=0;i<predict_sigma_aug.cols();i++){
        rho = sqrt((predict_sigma_aug(0,i)*predict_sigma_aug(0,i))+(predict_sigma_aug(1,i)*predict_sigma_aug(1,i)));
        theta = atan2(predict_sigma_aug(1,i),predict_sigma_aug(0,i));
         while (theta> M_PI) theta-=2.*M_PI;
        while (theta<-M_PI) theta+=2.*M_PI;
        r_rho = (predict_sigma_aug(0,i)*cos(predict_sigma_aug(3,i))*predict_sigma_aug(2,i) + predict_sigma_aug(1,i)*sin(predict_sigma_aug(3,i))*predict_sigma_aug(2,i))/rho;
        measured_predict_sigma_aug.col(i) << rho,theta,r_rho;
    }
   
}

void predictMeanAndCovariance(MatrixXd x,int distinction){
    VectorXd extractedMean(x.rows());
    extractedMean.fill(0.0);
    MatrixXd extractedCovariance(x.rows(),x.rows());
    extractedCovariance.fill(0.0);
    weight1 = glambda/(glambda+n_aug);
    weightn = 0.5/(glambda+n_aug);
    VectorXd x_diff(x.rows());
    for(int i = 0;i<x.cols();i++){
        if(i==0){
            extractedMean += weight1*x.col(i);
        }
        else{
             extractedMean += weightn*x.col(i);
        }
    }
    

    for(int i=0;i<x.cols();i++){
        int select;
        if(x.rows()>3){
            select = 3;
        }
        else{
            select = 1;
        }
        if(i==0){
             while (x_diff(select)> M_PI) x_diff(select)-=2.*M_PI;
             while (x_diff(select)<-M_PI) x_diff(select)+=2.*M_PI;
            x_diff = x.col(i)-extractedMean;
            extractedCovariance += weight1*x_diff*(x_diff.transpose());
        }
        else{
             while (x_diff(select)> M_PI) x_diff(select)-=2.*M_PI;
             while (x_diff(select)<-M_PI) x_diff(select)+=2.*M_PI;
             x_diff = x.col(i)-extractedMean;
             extractedCovariance += weightn*x_diff*(x_diff.transpose());
        }
    }
    if(distinction==0){
        stateMean = extractedMean;
        stateCovariance = extractedCovariance;
    }
    else{
        measurementMean = extractedMean;
        measurementCovariance = extractedCovariance;
    }

}

void predictSigmaPoints(double dt){
    for(int i=0;i<sigma_aug.cols();i++){
        VectorXd state(5); 
        state = sigma_aug.col(i).head(5);
        VectorXd diffrential(5);
        VectorXd noise(5);
        noise(0) = 0.5*dt*dt*cos(state(3))*sigma_aug(5,i);
        noise(1) = 0.5*dt*dt*sin(state(3))*sigma_aug(5,i);
        noise(2) = dt*sigma_aug(5,i);
        noise(3) = 0.5*dt*dt*sigma_aug(6,i);
        noise(4) = dt*sigma_aug(6,i);
        if(state(4)!=0 || state(4)!=0.0){
            diffrential(0) = (state(2)/state(4))*(sin(state(3) + (state(4)*dt)) - sin(state(3)));
            diffrential(1) = (state(2)/state(4))*(cos(state(3)) - cos(state(3) + (state(4)*dt)) );
            diffrential(2) = 0;
            diffrential(3) = dt*state(4);
            diffrential(4) = 0;
            state  = state + diffrential + noise;
           
        }
        else{
            diffrential(0) = state(2)*cos(state(3))*dt;
            diffrential(1) = state(2)*sin(state(3))*dt;
            diffrential(2) = 0;
            diffrential(3) = dt*state(4);
            diffrential(4) = 0;
            state  = state + diffrential + noise;
        }
        predict_sigma_aug.col(i) = state;
    }
}

void generateAugmentedMatrix(VectorXd x,MatrixXd P){
  x_aug.head(5) = x;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P;
  P_aug(5,5) = std_a*std_a;
  P_aug(6,6) = std_yawdd*std_yawdd;

}

MatrixXd genarateSigmaPoints(VectorXd state,MatrixXd P){
    int nsigma = 2*state.rows() + 1;
    int lambda = 3 - n_aug;
    glambda = lambda;
    MatrixXd sigmaPoints(state.rows(),nsigma);
    double coeff = lambda + state.rows();
    MatrixXd stage1 = P*coeff;
    MatrixXd stage2 = stage1.llt().matrixL();
    MatrixXd stage3_1(stage2.rows(),stage2.cols());
    MatrixXd stage3_2(stage2.rows(),stage2.cols());
    for(int i=0;i<stage2.cols();i++){
        stage3_1.col(i) = state + stage2.col(i);
        stage3_2.col(i) = state - stage2.col(i); 
    }
    MatrixXd stage3(stage3_1.rows(),state.cols()+stage3_1.cols()+stage3_2.cols());
    stage3 << state,stage3_1,stage3_2;
    return stage3;
    
}

int main(){
predict_sigma_aug.fill(0.0);
measured_predict_sigma_aug.fill(0.0);
stateMean.fill(0.0);
stateCovariance.fill(0.0);
measurementMean.fill(0.0);
measurementCovariance.fill(0.0);
kalmanGain.fill(0.0);

 R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;

ifstream reader("C:\\DEV\\website\\data\\obj_pose-laser-radar-synthetic-input.txt");
   string line;
     long long c_timestamp;
  long long p_timestamp;
  char type;
  int count = 0;//for initializing
  double gt_x;//ground truth
  double gt_y;//ground truth
  double gt_vy;//ground truth
  double gt_vx;//ground truth 
   double gt_v;

   while(getline(reader,line)){
    istringstream my_stream(line);
    my_stream>>type;
    if(type=='R'){
         my_stream>>m_rho;
      my_stream>>m_pi;
      my_stream>>m_rdot;
      z<< m_rho,m_pi,m_rdot;
       my_stream>>c_timestamp;
        if(count==0){
        x  = radar2State(z);
         P <<  1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

        p_timestamp = c_timestamp;
        count++;
        continue;
      }
       dt = (c_timestamp - p_timestamp)/1000000.0;
      p_timestamp = c_timestamp;
       my_stream>>gt_x;
      my_stream>>gt_y;
      my_stream>>gt_vx;
      my_stream>>gt_vy;
      gt_v = sqrt(gt_vx*gt_vx + gt_vy*gt_vy);
       VectorXd gt(3);
       gt << gt_x,gt_y,gt_v;
        generateAugmentedMatrix(x,P);//converts vector of size 5 to 7 and matrix of 5x5 to 7x7
    sigma_aug = genarateSigmaPoints(x_aug,P_aug);//outputs 7x15 matrix
    predictSigmaPoints(0.1);//outputs 5x15 matrix
    predictMeanAndCovariance(predict_sigma_aug,0);//outputs a vector of size 5 and matrix of size 5x5
    convertIntoMeasurementSpace();//outputs 3x15 matrix from 5x15 matrix
    predictMeanAndCovariance(measured_predict_sigma_aug,1);//mean and covariance in measurement space
    updateState();
    calculateError(x,gt);
    }
   }
    return 0;
}