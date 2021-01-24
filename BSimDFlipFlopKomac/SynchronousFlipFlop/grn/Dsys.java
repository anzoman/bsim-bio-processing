package BSimDFlipFlopKomac.SynchronousFlipFlop.grn;

import bsim.ode.BSimOdeSystem;

/*
     * Representation of the ODE system
     */
public class Dsys implements BSimOdeSystem {
    int numEq = 4;				// System of 4 equations

    private double D = 0;				// External D chemical level
    private double CLK = 0;				// External CLK chemical level
    // Parameters from the paper:
    private double a1 = 0.8508;//s^-1
    private double a2 = 1.5299;//s^-1
    private double a3 = 0.3431;//s^-1
    private double a4 = 1.5299;//s^-1

    private double Kd1 = 99.0481;//nM
    private double Kd2 = 12.4672*100;//nM rescaled for chen oscillator lower period
    private double Kd3 = 34.9188;//nM
    private double Kd4 = 99.0481;//nM
    private double Kd5 = 14.6698*100;//nM rescaled for chen oscillator lower period
    private double Kd6 = 11.7473;//nM
    private double Kd7 = 99.8943;//nM

    private double dt1 = 0.0036;//s^-1
    private double dt2 = 0.0036;//s^-1

    private double unitStep(double inp){
        return (inp < 0)? 0 : 1;
    }

    public double[] derivativeSystem(double t, double[] y) {
        double[] dy = new double[numEq];


        //a
        dy[0] =  a1 * unitStep(D - Kd1) * unitStep(Kd2 - CLK) + a2 * unitStep(Kd3 - y[1]) - dt1 * y[0];
        //ac
        dy[1] =  a1 * unitStep(Kd1 - D) * unitStep(Kd2 - CLK) + a2 * unitStep(Kd3 - y[0]) - dt1 * y[1];

        //q
        dy[2] =  a3 * unitStep(y[0] - Kd4) * unitStep(CLK - Kd5) * unitStep(Kd7 - y[2]) + a4 * unitStep(Kd6 - y[3]) * unitStep(Kd7 - y[2]) - dt2 * y[2];
        //qc
        dy[3] =  a3 * unitStep(y[1] - Kd4) * unitStep(CLK - Kd5) * unitStep(Kd7 - y[3]) + a4 * unitStep(Kd6 - y[2]) * unitStep(Kd7 - y[3]) - dt2 * y[3];

		dy[0] *= t;
		dy[1] *= t;
		dy[2] *= t;
		dy[3] *= t;

        return dy;
    }

    // Set up external chemical level
    public void setExternalLevel(double _d, double _clk){
        D = _d;
        CLK = _clk;
    }

    // Initial conditions for the ODE system
    public double[] getICs() {
        double[] ics = new double[numEq];

        ics[0] = 0.0;
        ics[1] = 0.0;
        ics[2] = 0.0;
        ics[3] = 0.0;
        return ics;
    }

    public int getNumEq() {
        return numEq;
    }


}
