# Sparse Identification of Nonlinear Dynamics for Model Predictive Control in the Low-Data Limit

Sparse identification of nonlinear dynamics with control (SINDYc) is combined with model predictive control (MPC). This framework learns nonlinear dynamical models affected by an exogenous control variable from few measurements. The resulting SINDYc models have the ability to enhance the performance of model predictive control (MPC), based on limited, noisy data. SINDYc models are parsimonious, identifying the fewest terms in the model needed to explain the data, making them interpretable and generalizable. We show that the resulting SINDY-MPC framework has higher performance, requires significantly less data, and is more computationally efficient and robust to noise than neural network models, making it viable for online training and execution in response to rapid system changes. SINDY-MPC also shows improved performance over linear data-driven models, although linear models may provide a stopgap until enough data is available for SINDY.
SINDY-MPC is demonstrated on a variety of dynamical systems with different challenges including the chaotic Lorenz system, a simple model for flight control of a F8 aircraft, and an HIV model incorporating drug treatment.

The publication "Sparse identification of nonlinear dynamics for model predictive control in the low-data limit"
by E. Kaiser, J. N. Kutz and S. L. Brunton. is available on [arXiv](https://arxiv.org/abs/1711.05501).

## Installation

1. Clone this repository to your desktop.
2. Add path to `SINDY-MPC/utils` folder to Matlab search path using `addpath('<path to mds>/SINDY-MPC/utils')`.

## Dependencies
No additional dependencies.

## Getting Started

The code for each example
See folder `/EX_XXXXX` 

demonstrating the approach on a stochastically forced linear system with known low-rank dynamics and an artificially inflated state dimension (Example 4.1 in [1]). Just execute this file in MatLab and it will generate the model, figures, and compute the reconstruction error.

## Organization

The algorithms are in the `SINDY-MPC/` directory.
The folder `SINDY-MPC/utils/` contains helper functions. Example specific functions are in the respective folder `SINDY-MPC/EX_YYYY/` for example `YYYY`.
Each example folder has a set of analogous functions:

	ConstraintFCN_models : Evaluates MPC constraint functions for different models
    ObjectiveFCN_models  : Computes control objective for different models
    ConstraintFCN		: Evaluates MPC constraint functions using true system model
    ObjectiveFCN		 : Computes control objective using true system model
    EX_YYYY_SI_DelayDM   : (Delay)DMDc System Identification (SI) for example YYYY
    EX_YYYY_SI_NARX  	: Neural Network System Identification (SI) for example YYYY
    EX_YYYY_SI_SINDYc    : SINDYc System Identification (SI) for example YYYY
    getTrainingData	  : Runs simulation to collect training data
    getValidationData	: Runs simulation to collect validation data
    getMPCparams 		: Defines parameters for MPC
    MPC_F8 			  : Runs MPC on the true system model
    MPC_F8_ModelComparison  : Comparison of all considered models regarding prediction on training and validation data and MPC
    getTrainingData_Ensemble : Runs the simulation for an ensemble of initial conditions
    
    
    
    EX_F8_SI_NARX_Ensemble : Neural Network System Identification (SI) for example YYYY 
    
    VIZ_3D_MODELvsTRUTH 
    VIZ_SI_Validation_MPC_ensemble
    VIZ_SI_Validation_MPC
    VIZ_SI_Validation

    
## License (MIT license)

See the [LICENSE file](LICENSE) for details.
