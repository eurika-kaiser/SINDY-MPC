# Sparse Identification of Nonlinear Dynamics for Model Predictive Control in the Low-Data Limit

Sparse identification of nonlinear dynamics with control (SINDYc) is combined with model predictive control (MPC). This framework learns nonlinear dynamical models affected by an exogenous control variable from few measurements. The resulting SINDYc models have the ability to enhance the performance of model predictive control (MPC), based on limited, noisy data. SINDYc models are parsimonious, identifying the fewest terms in the model needed to explain the data, making them interpretable and generalizable. We show that the resulting SINDY-MPC framework has higher performance, requires significantly less data, and is more computationally efficient and robust to noise than neural network models, making it viable for online training and execution in response to rapid system changes. SINDY-MPC also shows improved performance over linear data-driven models, although linear models may provide a stopgap until enough data is available for SINDY.
SINDY-MPC is demonstrated on a variety of dynamical systems with different challenges including the Lotka-Volterra system, the chaotic Lorenz system, a simple model for flight control of a F8 aircraft, and an HIV model incorporating drug treatment.

The publication "Sparse identification of nonlinear dynamics for model predictive control in the low-data limit"
by E. Kaiser, J. N. Kutz and S. L. Brunton. is available on [arXiv](https://arxiv.org/abs/1711.05501).

## Installation

1. Clone this repository to your desktop.
2. Add path to `SINDY-MPC/utils` folder to Matlab search path using `addpath('<path to SINDY-MPC>/SINDY-MPC')`.

## Dependencies
No additional dependencies.

## Getting Started

The code for each example `YYYY` is in the corresponding example folder `/EX_YYYY`. Code used for all examples can be found in `SINDY-MPC/utils`, example-specific code, e.g. for plotting, will be in the corresponding example folder.

1. Go into an example folder `SINDY-MPC/EX_YYYY`.
2. Run scripts for SINDYc system identification, e.g. `EX_YYYY_SI_SINDYc.m`. To train other models replace `SINDYc` with `NARX` for a neural network or `DMDc` for a linear system. The trained models are saved in the folder `SINDY-MPC/DATA/`.
3. Run model predictive control by choosing one of these options:
	a. Execute `MPC_YYYY_SINDYc.m` to use SINDYc in MPC.
    b. Run MPC for all models (e.g. SINDYc, neural network, linear system via DMDc) by executing `MPC_LOTKA_ModelComparison.m`.
4. Saved figures may be found in `SINDY-MPC/FIGURES/YYYY/`.

It is easier to get started with the examples `FLIGHT_CONTROL_F8` and `HIV_THERAPY`.
The models for the HIV system are included in the folder `SINDY-MPC/DATA/HIV/`. So in order to obtain the control results, one can immediately start by executing `MPC_HIV_ModelComparison.m` without prior computation of the models.

## Organization

The algorithms are in the `SINDY-MPC/` directory. The folder `SINDY-MPC/utils/` contains helper functions. Example specific functions are in the respective folder `SINDY-MPC/EX_YYYY/` for example `YYYY`. Each example folder contains a similar set of functions. The most important ones are:

    EX_YYYY_SI_SINDYc    : SINDYc System Identification (SI) for example YYYY
	EX_YYYY_SI_DelayDM   : (Delay)DMDc System Identification (SI) for example YYYY
    EX_YYYY_SI_NARX  	: Neural Network System Identification (SI) for example YYYY
    getTrainingData	  : Runs simulation to collect training data
    getValidationData	: Runs simulation to collect validation data
    getMPCparams 		: Defines parameters for MPC
	ConstraintFCN_models : Evaluates MPC constraint functions for different models
    ObjectiveFCN_models  : Computes control objective for different models
    ConstraintFCN		: Evaluates MPC constraint functions using true system model
    ObjectiveFCN		 : Computes control objective using true system model    
    MPC_YYYY_ModelComparison  : Comparison of all considered models regarding prediction on training and validation data and MPC
    
## License ([CiteMe OSS](https://github.com/cite-me/oss))

See the [LICENSE file](LICENSE) for details.
