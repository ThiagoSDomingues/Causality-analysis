#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

class HydroConditions {
private:
    static constexpr double DEFAULT_EPS = 1e-13;

public:
    // Pre-conditions
    static bool pre_condition_A1(double tau_PI, double tau_pi, double eta, double zeta, 
                                  double tau_pipi, double delta_PIPI, double lambda_PIpi, 
                                  double delta_pipi, double lambda_piPI, double cs2, 
                                  double eps = DEFAULT_EPS) {
        return (tau_pi > eps) && (tau_pi > eps) && (eta >= eps) && (zeta >= eps) && 
               (tau_pipi >= eps) && (delta_PIPI >= eps) && (lambda_PIpi >= eps) && 
               (delta_pipi >= eps) && (lambda_piPI >= eps) && (cs2 >= eps);
    }

    static bool pre_condition_A2(double e, double p, double bulk, double eps = DEFAULT_EPS) {
        return (e > eps) && (p >= eps) && (e + p + bulk > eps);
    }

    static bool pre_condition_A3(double e, double p, double bulk, double lambda_a, 
                                  double eps = DEFAULT_EPS) {
        return e + p + bulk + lambda_a > eps;
    }

    // Freeze-out temperature (requires interpolation function T_interp)
    static double temperature_freeze(double e_freeze) {
        // Note: T_interp function needs to be implemented
        // double fz_temperature = T_interp(e_freeze);
        // return fz_temperature;
        return 0.0; // Placeholder
    }

    // Necessary conditions
    static bool necessary_condition_4a(double eta, double lambda_piPI, double bulk, 
                                        double tau_pipi, double lambda_1) {
        return (2*eta + lambda_piPI*bulk) - 0.5*tau_pipi*std::abs(lambda_1) >= 0;
    }

    static bool necessary_condition_4b(double e, double p, double bulk, double tau_pi, 
                                        double eta, double lambda_piPI, double tau_pipi, 
                                        double lambda_3) {
        return e + p + bulk - (1.0/(2*tau_pi))*(2*eta + lambda_piPI*bulk) - 
               tau_pipi*lambda_3 >= 0;
    }

    static bool necessary_condition_4c(double tau_pi, double eta, double lambda_piPI, 
                                        double bulk, double tau_pipi, double lambda_a, 
                                        double lambda_d) {
        return (1.0/(2*tau_pi))*(2*eta + lambda_piPI*bulk) + 
               (tau_pipi/(4*tau_pi))*(lambda_a + lambda_d) >= 0;
    }

    static bool necessary_condition_4d(double e, double p, double bulk, double lambda_a, 
                                        double tau_pi, double eta, double lambda_piPI, 
                                        double tau_pipi, double lambda_d) {
        return e + p + bulk + lambda_a - (1.0/(2*tau_pi))*(2*eta + lambda_piPI*bulk) - 
               (tau_pipi/(4*tau_pi))*(lambda_d + lambda_a) >= 0;
    }

    static bool necessary_condition_4e(double tau_pi, double eta, double lambda_piPI, 
                                        double bulk, double tau_pipi, double lambda_d, 
                                        double delta_pipi, double delta_PIPI, 
                                        double lambda_PIpi, double tau_PI, double p, 
                                        double cs2, double zeta, double e) {
        return (1.0/(2*tau_pi))*(2*eta + lambda_piPI*bulk) + 
               (tau_pipi/(2*tau_pi))*lambda_d + 
               (1.0/(6*tau_pi))*(2*eta + lambda_piPI*bulk + (6*delta_pipi - tau_pipi)*lambda_d) + 
               (zeta + delta_PIPI*bulk + lambda_PIpi*lambda_d)/tau_PI + 
               (e + p + bulk + lambda_d)*cs2 >= 0;
    }

    static bool necessary_condition_4f(double e, double p, double bulk, double lambda_d, 
                                        double tau_pi, double eta, double lambda_piPI, 
                                        double tau_pipi, double delta_pipi, double delta_PIPI, 
                                        double zeta, double cs2, double tau_PI, 
                                        double lambda_PIpi) {
        return e + p + bulk + lambda_d - (1.0/(2*tau_pi))*(2*eta + lambda_piPI*bulk) - 
               (tau_pipi*lambda_d)/(2*tau_pi) - 
               (1.0/(6*tau_pi))*(2*eta + lambda_piPI*bulk + (6*delta_pipi - tau_pipi)*lambda_d) - 
               (zeta + delta_PIPI*bulk + lambda_PIpi*lambda_d)/tau_PI - 
               (e + p + bulk + lambda_d)*cs2 >= 0;
    }

    // Sufficient conditions
    static bool sufficient_condition_5a(double e, double p_QCD, double qcd_bulk, 
                                         double lambda_1, double eta, double lambda_piPI, 
                                         double tau_pi, double tau_pipi, double lambda_3) {
        return (e + p_QCD + qcd_bulk - std::abs(lambda_1)) - 
               (2*eta + lambda_piPI*qcd_bulk)/(2*tau_pi) - 
               tau_pipi*lambda_3/(2*tau_pi) >= 0;
    }

    static bool sufficient_condition_5b(double eta, double lambda_piPI, double bulk, 
                                         double tau_pipi, double lambda_1) {
        return 2*eta + lambda_piPI*bulk - tau_pipi*std::abs(lambda_1) > 0;
    }

    static bool sufficient_condition_5c(double tau_pipi, double delta_pipi) {
        return tau_pipi <= 6 * delta_pipi;
    }

    static bool sufficient_condition_5d(double lambda_PIpi, double tau_PI, double cs2, 
                                         double tau_pipi, double tau_pi) {
        return (lambda_PIpi/tau_PI) + cs2 - (tau_pipi/(12*tau_pi)) >= 0;
    }

    static bool sufficient_condition_5e(double e, double p, double bulk, double lambda_1, 
                                         double lambda_3, double eta, double lambda_piPI, 
                                         double tau_pi, double delta_pipi, double tau_pipi, 
                                         double zeta, double cs2, double delta_PIPI, 
                                         double lambda_PIpi, double tau_PI) {
        double denominator = e + p + bulk - std::abs(lambda_1) - 
                            (2*eta + lambda_piPI*bulk)/(2*tau_pi) - 
                            tau_pipi*lambda_3/(2*tau_pi);
        
        double lhs = (1.0/(3*tau_pi))*(4*eta + 2*lambda_piPI*bulk + (3*delta_pipi + tau_pipi)*lambda_3) + 
                     (zeta + delta_PIPI*bulk + lambda_PIpi*lambda_3)/tau_PI + 
                     std::abs(lambda_1) + lambda_3*cs2 + 
                     (((12*delta_pipi - tau_pipi)/(12*tau_pi))*(lambda_PIpi/tau_PI + cs2 - tau_pipi/(12*tau_pi))*
                     std::pow(lambda_3 + std::abs(lambda_1), 2))/denominator;
        
        return lhs <= (e + p + bulk)*(1 - cs2);
    }

    static bool sufficient_condition_5f(double tau_pi, double tau_pipi, double delta_pipi, 
                                         double bulk, double e, double p, double eta, 
                                         double lambda_piPI, double lambda_PIpi, double lambda_1, 
                                         double zeta, double cs2, double delta_PIPI, double tau_PI) {
        return (1.0/(6*tau_pi))*(2*eta + lambda_piPI*bulk + (tau_pipi - 6*delta_pipi)*std::abs(lambda_1)) + 
               (zeta + delta_PIPI*bulk - lambda_PIpi*std::abs(lambda_1))/tau_PI + 
               (e + p + bulk - std::abs(lambda_1))*cs2 >= 0;
    }

    static bool sufficient_condition_5g(double tau_pipi, double tau_PI, double delta_pipi, 
                                         double lambda_PIpi, double tau_pi, double cs2, 
                                         double lambda_1, double lambda_3, double eta, 
                                         double lambda_piPI, double bulk) {
        double numerator = ((12*delta_pipi - tau_pipi)/(12*tau_pi))*
                          (lambda_PIpi/tau_PI + cs2 - tau_pipi/(12*tau_pi))*
                          std::pow(lambda_3 + std::abs(lambda_1), 2);
        
        double denominator = std::pow((1.0/(2*tau_pi))*(2*eta + lambda_piPI*bulk) - 
                                      tau_pipi*std::abs(lambda_1)/(2*tau_pi), 2);
        
        return numerator/denominator <= 1;
    }

    static bool sufficient_condition_5h(double tau_pipi, double delta_pipi, double tau_pi, 
                                         double eta, double lambda_piPI, double zeta, 
                                         double delta_PIPI, double lambda_PIpi, double e, 
                                         double p, double bulk, double lambda_1, double cs2, 
                                         double lambda_2, double lambda_3, double tau_PI) {
        double lhs = (1.0/(3*tau_pi))*(4*eta + 2*lambda_piPI*bulk - (3*delta_pipi + tau_pipi)*std::abs(lambda_1)) + 
                     (zeta + delta_PIPI*bulk - lambda_PIpi*std::abs(lambda_1))/tau_PI + 
                     (e + p + bulk - std::abs(lambda_1))*cs2;
        
        double rhs = (((e + p + bulk + lambda_2)*(e + p + lambda_3))/
                     (3*(e + p + bulk - std::abs(lambda_1))))*
                     (1 + 2*((1.0/(2*tau_pi))*(2*eta + lambda_piPI*bulk) + 
                     (tau_pipi/(2*tau_pi))*lambda_3)/(e + p + bulk - std::abs(lambda_1)));
        
        return lhs >= rhs;
    }

    // Linear causality inequality
    static double linear_causality_inequality(double relax_shear_factor, double relax_bulk_factor) {
        double a = 0.14976611173885748;
        double x = relax_shear_factor;
        double y = relax_bulk_factor;
        
        return a + (4.0/3)*(1.0/x) + (1.0/y)*std::pow(1.0/3 - a, 2);
    }

    // Classification function for grid points
    static std::vector<std::vector<std::string>> classify_grid_points(
        const std::vector<std::vector<bool>>& pre_conditions,
        const std::vector<std::vector<bool>>& necessary_conditions,
        const std::vector<std::vector<bool>>& sufficient_conditions,
        size_t rows, size_t cols) {
        
        std::vector<std::vector<std::string>> classification(rows, std::vector<std::string>(cols));
        
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                // Check all pre-conditions
                bool pre_check = true;
                for (const auto& cond : pre_conditions) {
                    if (!cond[i * cols + j]) {
                        pre_check = false;
                        break;
                    }
                }
                
                // Check necessary conditions
                bool necessary_check = pre_check;
                if (pre_check) {
                    for (const auto& cond : necessary_conditions) {
                        if (!cond[i * cols + j]) {
                            necessary_check = false;
                            break;
                        }
                    }
                }
                
                // Check sufficient conditions
                bool sufficient_check = necessary_check;
                if (necessary_check) {
                    for (const auto& cond : sufficient_conditions) {
                        if (!cond[i * cols + j]) {
                            sufficient_check = false;
                            break;
                        }
                    }
                }
                
                // Classify
                if (!pre_check) {
                    classification[i][j] = "black";
                } else if (sufficient_check) {
                    classification[i][j] = "blue";
                } else if (necessary_check) {
                    classification[i][j] = "purple";
                } else {
                    classification[i][j] = "red";
                }
            }
        }
        
        return classification;
    }
};
