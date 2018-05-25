#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    ekf_.H_ = MatrixXd(3, 4);
    ekf_.R_ = MatrixXd(3, 3);

  // Measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
                0, 0.0225;

  // Measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;

  // H Laser
      H_laser_ << 1, 0, 0, 0,
                  0, 1, 0, 0;

  // State Transition Matrix
      ekf_.F_ = MatrixXd(4, 4);
      ekf_.F_ <<  1, 0, 1, 0,
                  0, 1, 0, 1,
                  0, 0, 1, 0,
                  0, 0, 0, 1;

  // State covariance matrix P
      ekf_.P_ = MatrixXd(4, 4);
      ekf_.P_ <<  1, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0, 1000, 0,
                  0, 0, 0, 1000;

  // Acceleration noise components
      noise_ax = noise_ay = 25;

      /** :

      ----------------------
          Fusion
      ----------------------
          Noise = 100:: <0.0885, 0.0952, 0.4853, 0.4599>
          Noise = 50 :: <0.0882, 0.0889, 0.4564, 0.4128>
          Noise = 40 :: <0.0885, 0.0872, 0.4503, 0.4053>
          
          Noise = 30 :: <0.0891, 0.0853, 0.4446, 0.4008>
          
          Noise = 25 :: <0.0897, 0.0843, 0.4424, 0.4009> -- sum = 1.01719
          Noise = 25 :: <0.0736, 0.0867, 0.4135, 0.4545> -- Data Set 2

          Noise = 25 :: <0.0880, 0.0859, 0.3524, 0.3962> -- Radar first measurement
          
          Noise = 23 :: <0.0900, 0.0839, 0.4417, 0.4016> -- sum = 1.01724
          Noise = 19 :: <0.0891, 0.0853, 0.4446, 0.4008> -- sum = 1.01980
          
          Noise = 18 :: <0.0916, 0.0831, 0.4411, 0.4081>
          Noise = 17 :: <0.0920, 0.0831, 0.4414, 0.4101>
          Noise = 16 :: <0.0924, 0.0831, 0.4419, 0.4125>
          Noise = 15 :: <0.0924, 0.0831, 0.4419, 0.4125>
          Noise = 14 :: <0.0929, 0.0832, 0.4425, 0.4153>
          Noise = 13 :: <0.0935, 0.0833, 0.4434, 0.4186>
          Noise = 12 :: <0.0942, 0.0836, 0.4447, 0.4226>
          Noise = 11 :: <0.0951, 0.0840, 0.4463, 0.4274>
          Noise =  9 :: <0.0973, 0.0855, 0.4513, 0.4399>


      -----------------------
          Radar Only
      -----------------------
          Noise = 50 :: <0.1976, 0.2714, 0.5205, 0.6596>
          Noise = 45 :: <0.1975, 0.2713, 0.5189, 0.6610>
          Noise = 40 :: <0.1976, 0.2715, 0.5182, 0.6638>
          Noise = 25 :: <0.2003, 0.2785, 0.5239, 0.6884>
          Noise =  9 :: <0.2302, 0.3464, 0.5835, 0.8040>

      -----------------------
          Laser Only
      -----------------------
          Noise = 45 :: <0.1310, 0.1092, 0.6121, 0.4824>
          Noise = 25 :: <0.1322, 0.1067, 0.6012, 0.4764>
          Noise =  9 :: <0.1473, 0.1153, 0.6383, 0.5346>

      **/

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
      if (!is_initialized_) 
      {
        /**
        TODO:
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
          // first measurement
            ekf_.x_ = VectorXd(4);
            ekf_.x_ << 1, 1, 1, 1;

        //---------------------
        //    Laser
        //---------------------
          if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
           // Convert radar from polar to cartesian coordinates and initialize state.
              double rho      = measurement_pack.raw_measurements_[0];
              double phi      = measurement_pack.raw_measurements_[1];
              double rho_dot  = measurement_pack.raw_measurements_[2];

              double x = rho * cos(phi);
              double y = rho * sin(phi);

              double vx = rho_dot * cos(phi);
              double vy = rho_dot * sin(phi);

          // Assign to state matrix
              ekf_.x_(0) = x;
              ekf_.x_(1) = y;
              ekf_.x_(2) = vx;
              ekf_.x_(3) = vy;

              is_initialized_ = true;
          }

        //---------------------
        //    Radar
        //---------------------
          else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) 
          {
            return;
            // // Initialize State
            //     double x = measurement_pack.raw_measurements_[0];
            //     double y = measurement_pack.raw_measurements_[1];

            //     double vx = 5.19;
            //     double vy = 0;

            //   // Assign to state matrix
            //     ekf_.x_(0) = x;
            //     ekf_.x_(1) = y;
            //     ekf_.x_(2) = vx;
            //     ekf_.x_(3) = vy;
          }

           previous_timestamp_ = measurement_pack.timestamp_;
           
        return;
      }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

    // Time
        float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; // dt - expressed in seconds
        previous_timestamp_ = measurement_pack.timestamp_;

        float dt_2 = dt * dt;
        float dt_3 = dt * dt_2;
        float dt_4 = dt * dt_3;

      // Modify the F matrix so that the time is integrated
        ekf_.F_(0, 2) = dt;
        ekf_.F_(1, 3) = dt;

      // Set the process covariance matrix Q
        ekf_.Q_ = MatrixXd(4,4);
        ekf_.Q_ << (dt_4 / 4 * noise_ax),    0 , (dt_3 / 2 * noise_ax) , 0,
                    0, (dt_4 / 4 * noise_ay), 0, ( dt_3 / 2*noise_ay) ,
                    dt_3/2*noise_ax, 0, dt_2 * noise_ax, 0,
                    0, dt_3 / 2 * noise_ay, 0, dt_2 * noise_ax;

      // Predict
        ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/
      
      //---------------------
      //    RADAR
      //---------------------
          if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
          {
              ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
              ekf_.R_ = R_radar_;

              ekf_.UpdateEKF(measurement_pack.raw_measurements_);
          } 

      //---------------------
      //    LASER
      //---------------------
          else 
          {
              ekf_.H_ = H_laser_;
              ekf_.R_ = R_laser_;

              ekf_.Update(measurement_pack.raw_measurements_);
          }

  // Output
    // cout << "\n\nx_ = " << ekf_.x_ << endl;
    // cout << "\nP_ = " << ekf_.P_ << endl;
}
