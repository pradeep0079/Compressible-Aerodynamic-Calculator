#include<iostream>
#include<cmath>
#include <limits>
using namespace std;

class IsentropicFlowRelations{
public:
    float gamma;
    float Mach;

    IsentropicFlowRelations(float gamma, float Mach){
        this->gamma = gamma;
        this->Mach = Mach;
    }

    float staticByTotalPressure(){
        return pow(1 +  (gamma - 1)/2 * pow(Mach,2), -((gamma)/(gamma-1)));
    }

    float staticByTotalTemperature(){
        return pow(1 +  (gamma - 1)/2 * pow(Mach, 2), -1);
    }

    float staticByTotalDensity(){
        return pow(1 +  (gamma - 1)/2 * pow(Mach,2), (-1/(gamma-1)));
    }

    float staticByChokeArea(){
        return  pow(1 + pow(Mach,2) * (gamma - 1)/2 , (gamma + 1)/(2*(gamma-1))) * pow((gamma + 1)/2, -((gamma + 1)/(2*(gamma -1))))/Mach;
    }

    float machAngle(){
        cout<<"Mach Angle is in Degree."<<endl;
        return asin(1/Mach)* (180/M_PI);
    }

    float prandtlMeyerAngle(){
        float a = sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1) / (gamma + 1) * (pow(Mach,2)- 1)));
        float b = atan(sqrt(pow(Mach,2) - 1));
        return (a - b) * (180 / M_PI);
    }

};

class NormalShockRelations : public IsentropicFlowRelations{
public:

    NormalShockRelations(float gamma, float Mach) : IsentropicFlowRelations(gamma, Mach){}
    float secondMach() {
        float a = ((gamma - 1) * pow(Mach,2)) + 2;
        float b = 2 * gamma * pow(Mach,2) - (gamma - 1);
        return sqrt(a / b);
    }

    float AfterDensity() {
        return ((gamma + 1) * pow(Mach,2)) / ((gamma - 1) *pow(Mach,2) + 2);
    }

    float AfterTemperature() {
        float a = (2 * pow(Mach,2) * gamma - (gamma - 1)) * ((gamma - 1) * pow(Mach,2) + 2);
        float b = (gamma + 1) * (gamma + 1) * pow(Mach,2);
        return a / b;
    }

    float AfterPressure() {
        float a = 2 * gamma * pow(Mach,2) - (gamma - 1);
        float b = gamma + 1;
        return a / b;
    }

    float AfterTotalPressure() {
        float a = pow(AfterDensity(), gamma / (gamma - 1));
        float b = pow(1 / AfterPressure(), 1 / (gamma - 1));
        return a * b;
    }
};

class ObliqueShockRelations : public IsentropicFlowRelations{
public:
    float turnAngle;

    ObliqueShockRelations(float gamma, float Mach, float turnAngle) : IsentropicFlowRelations(gamma, Mach) {
        this->turnAngle = turnAngle;
    }

    float theta_beta_m_relation(float beta) {
        // Theta-Beta-Mach relation equation.
        float beta_rad = beta * M_PI / 180.0;
        float term1 = 2 * (1 / tan(beta_rad)) * (Mach * Mach * pow(sin(beta_rad), 2) - 1);
        float term2 = Mach * Mach * (gamma + cos(2 * beta_rad)) + 2;
        return tan(turnAngle * M_PI / 180.0) - (term1 / term2);
    }

    float theta_beta_m_relation_derivative(float beta) {
        // Derivative of Theta-Beta-Mach relation equation with respect to beta.
        float beta_rad = beta * M_PI / 180.0;
        float term1 = 2 * (1 / tan(beta_rad)) * (Mach * Mach * pow(sin(beta_rad), 2) - 1);
        float term2 = Mach * Mach * (gamma + cos(2 * beta_rad)) + 2;
        float derivative_term1 = -2 * (1 / (sin(beta_rad) * sin(beta_rad))) * (Mach * Mach * pow(sin(beta_rad), 2) - 1);
        float derivative_term2 = 4 * Mach * Mach * sin(beta_rad) * cos(beta_rad);
        return -(derivative_term1 + term2 * derivative_term2) / (term2 * term2);
    }

    float find_shock_angle() {
        // Initial guess in degree
        float beta_guess = 45.0;
        float beta_solution = beta_guess;
        float tolerance = 1e-2;
        int max_iterations = 10;
        int iteration = 0;

        while (fabs(theta_beta_m_relation(beta_solution)) > tolerance && iteration < max_iterations) {
            beta_solution -= theta_beta_m_relation(beta_solution) / theta_beta_m_relation_derivative(beta_solution);
            iteration++;
        }

        if (iteration == max_iterations) {
            cout << "Solution did not converge" << endl;
            return std::numeric_limits<float>::quiet_NaN();
        }

        return beta_solution;
    }

    float deflection_angle() {
        float beta = find_shock_angle();
        return turnAngle;
    }

    float changeInDensity() {
        float beta_rad = find_shock_angle() * M_PI / 180.0;
        float a = (gamma + 1) * Mach * Mach * pow(sin(beta_rad), 2);
        float b = (gamma - 1) * Mach * Mach * pow(sin(beta_rad), 2) + 2;
        return a / b;
    }

    float changeInPressure() {
        float beta_rad = find_shock_angle() * M_PI / 180.0;
        float a = 2 * gamma * Mach * Mach * pow(sin(beta_rad), 2) - (gamma - 1);
        return a / (gamma + 1);
    }

    float changeInTemperature() {
        return changeInPressure() * (1 / changeInDensity());
    }

    float totalPressureChange() {
        float a = pow(changeInDensity(), gamma / (gamma - 1));
        float b = pow(1 / changeInPressure(), 1 / (gamma - 1));
        return a * b;
    }

    float AfterShockMachNumber() {
        float beta_rad = find_shock_angle() * M_PI / 180.0;
        float a = 2 + (gamma - 1) * Mach * Mach * pow(sin(beta_rad), 2);
        float b = 2 * gamma * Mach * Mach * pow(sin(beta_rad), 2) - (gamma - 1);
        float c = sqrt(a / b);
        return (1 / (sin(beta_rad - (turnAngle * M_PI / 180.0)))) * (c);
    }
};

int main() {
    IsentropicFlowRelations p1(1.4, 2);
    cout << "P/Po: " << p1.staticByTotalPressure() << endl;
    cout << "T/To: " << p1.staticByTotalTemperature() << endl;
    cout << "rho/rho_o: " << p1.staticByTotalDensity() << endl;
    cout << "A/A*: " << p1.staticByChokeArea() << endl;
    cout << "Mach Angle: " << p1.machAngle() << endl;
    cout << "P-M Angle: " << p1.prandtlMeyerAngle() << endl;

    NormalShockRelations p2(1.4, 2);
    cout << "After Normal Shock M: " << p2.secondMach() << endl;
    cout << "Density change: " << p2.AfterDensity() << endl;
    cout << "Temperature change: " << p2.AfterTemperature() << endl;
    cout << "Pressure change: " << p2.AfterPressure() << endl;
    cout << "Total Pressure change: " << p2.AfterTotalPressure() << endl;

    ObliqueShockRelations p3(1.4, 2, 20);
    cout << "Shock Angle: " << p3.find_shock_angle() << endl;
    cout << "Change in Density: " << p3.changeInDensity() << endl;
    cout << "Change in Pressure: " << p3.changeInPressure() << endl;
    cout << "Change in Temperature: " << p3.changeInTemperature() << endl;
    cout << "Total Pressure Change: " << p3.totalPressureChange() << endl;
    cout << "After Shock Mach Number: " << p3.AfterShockMachNumber() << endl;

    return 0;
}

