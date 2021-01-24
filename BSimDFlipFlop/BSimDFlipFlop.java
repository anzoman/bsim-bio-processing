package BSimDFlipFlop;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.draw.BSimDrawer;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;
import bsim.export.BSimMovExporter;
import bsim.export.BSimPngExporter;
import bsim.ode.BSimOdeSolver;
import bsim.ode.BSimOdeSystem;
import bsim.particle.BSimBacterium;
import processing.core.PGraphics3D;

import javax.vecmath.Vector3d;
import java.awt.*;
import java.util.Calendar;
import java.util.Random;
import java.util.Vector;

/**
 * Simulation of bacteria with D flip-flop GRNs coupled by a chemical field.</br>
 * <p>
 * ODE system for the D flip-flop taken from:
 * "Modeling a synthetic multicellular clock: Repressilators coupled by quorum sensing"
 * PNAS 2004 101 (30) 10955-10960; doi:10.1073/pnas.0307095101
 */
public class BSimDFlipFlop {

    // Initial conditions of the GRNs - Used for convenience
    public static int ICS_RANDOM = 1;
    public static int ICS_UNIFORM = 2;

    /*********************************************************
     * Simulation Definition
     */
    public static void main(String[] args) {


        /*********************************************************
         * Create a new directory for the simulation results
         */
        String filePath = BSimUtils.generateDirectoryPath("results/" + BSimUtils.timeStamp() + "/");

        /*********************************************************
         * Create a new simulation object
         */
        BSim sim = new BSim();            // New Sim object
        sim.setDt(0.01);                // Simulation Timestep
        sim.setSimulationTime(100);        // Simulation Length (6000 = 100 minutes)
        sim.setTimeFormat("0.00");        // Time Format for display
        sim.setBound(100, 100, 100);        // Simulation Boundaries


        /*********************************************************
         * Set up some constants
         */
        // diffusivity of AI in (microns^2/sec)? (BSim 1.0: diffusion coeff. of AHL = 0.23 per second)
        final double diffusivity = 100;    // Diffusivity of AHL
        final double decayRate = 0.01 / 60;    // Decay Rate (0.1666 = 10 per minute, 0.01/60 = 1e-2 per min)

        final double cellWallDiffusivity = 2.0;        // Cell wall diffusivity (Should be the same as eta in the GRN?)
        final int theInitialConditions = ICS_RANDOM;    // What initial conditions do we want to use?

        // Set up the chemical field for AHL:
        final BSimChemicalField field = new BSimChemicalField(sim, new int[]{25, 25, 25}, diffusivity, decayRate);

        BSimChemicalField hField = new BSimChemicalField(sim, new int[]{(int) sim.getBound().x, (int) sim.getBound().y, 1}, diffusivity, decayRate);
        BSimChemicalField iField = new BSimChemicalField(sim, new int[]{(int) sim.getBound().x, (int) sim.getBound().y, 1}, diffusivity, decayRate);

        BSimChemicalField dField = new BSimChemicalField(sim, new int[]{(int) sim.getBound().x, (int) sim.getBound().y, 1}, diffusivity, decayRate);
        BSimChemicalField qField = new BSimChemicalField(sim, new int[]{(int) sim.getBound().x, (int) sim.getBound().y, 1}, diffusivity, decayRate);
        BSimChemicalField qcField = new BSimChemicalField(sim, new int[]{(int) sim.getBound().x, (int) sim.getBound().y, 1}, diffusivity, decayRate);

        // ICs as in original DDEs
        hField.setConc(10.0);
        iField.setConc(10.0);


        /*********************************************************
         * BSimBacterium with a D flip-flop inside
         */
        class BSimDFlipFlopBacterium extends BSimBacterium {
            protected QuorumDFlipFlop repGRN;          // Instance of internal class
            protected double[] y, yNew;                // Local values of ODE variables

            final double cellWallDiffusivity = 2.0;  // Cell wall diffusivity, taken from other implementations using BSimCapsuleBacterium
            BSimChemicalField hField;
            BSimChemicalField iField;
            BSimChemicalField dField;
            BSimChemicalField qField;
            BSimChemicalField qcField;

            public BSimDFlipFlopBacterium(BSim sim, Vector3d position, BSimChemicalField h, BSimChemicalField i, BSimChemicalField d, BSimChemicalField q, BSimChemicalField qc) {
                super(sim, position);

                // Create the parameters and initial conditions for the ODE system
                repGRN = new QuorumDFlipFlop();
                repGRN.generateBeta();
                y = repGRN.getICs();

                hField = h;
                iField = i;
                dField = d;
                qField = q;
                qcField = qc;
            }

            /*
             * Action each time step: iterate ODE system for one time-step.
             */
            @Override
            public void action() {

                // Movement
                super.action();

                // Chemical fields concentrations
                double externalChemD;    // External data chem. field
                double externalChemCLK;    // External clock chem. field
                double externalChemQ;    // External q chem. field
                double externalChemQc;    // External qc chem. field
                double deltaChemQ;        // Change in q chemical quantity
                double deltaChemQc;        // Change in qc chemical quantity

                // external chemical level at position of the bacterium
                externalChemD = dField.getConc(position);
                externalChemCLK = iField.getConc(position);
                externalChemQ = qField.getConc(position);
                externalChemQc = qcField.getConc(position);

                // Get the external chemical field level for the GRN ode system later on:
                /* Qc for reverse! */
                repGRN.setExternalLevel(externalChemQc, externalChemCLK);

                // re-scaled time units
                yNew = BSimOdeSolver.rungeKutta45(repGRN, sim.getTime() / 60, y, sim.getDt() / 60);
                y = yNew;

                // Adjust the external chemical field
                deltaChemQ = externalChemQ - y[2];
                deltaChemQc = externalChemQc - y[3];

                // Changing external concentration
                qField.addQuantity(position, cellWallDiffusivity * (-deltaChemQ));
                qcField.addQuantity(position, cellWallDiffusivity * (-deltaChemQc));
            }

            /*
             * Representation of the D flip-flop ODE system with quorum coupling
             */
            class QuorumDFlipFlop implements BSimOdeSystem {
                int numEq = 4;                       // System of 4 equations

                private double D = 0;                // External D chemical level
                private double CLK = 0;              // External CLK chemical level
                private Random r = new Random();     // Random number generator
                private double beta;                 // beta parameter

                // Parameters are from the paper: https://www.sciencedirect.com/science/article/abs/pii/S1877750316303866
                private double a1 = 0.8508;
                private double a2 = 1.5299;
                private double a3 = 0.3431;
                private double a4 = 1.5299;

                private double Kd1 = 99.0481;
                private double Kd2 = 12.4672 * 100;
                private double Kd3 = 34.9188;
                private double Kd4 = 99.0481;
                private double Kd5 = 14.6698 * 100;
                private double Kd6 = 11.7473;
                private double Kd7 = 99.8943;

                private double dt1 = 0.0036;
                private double dt2 = 0.0036;

                private double unitStep(double inp) {
                    return (inp < 0) ? 0 : 1;
                }

                public double[] derivativeSystem(double t, double[] y) {
                    double[] dy = new double[numEq];
                    //a
                    dy[0] = a1 * unitStep(D - Kd1) * unitStep(Kd2 - CLK) + a2 * unitStep(Kd3 - y[1]) - dt1 * y[0];
                    //ac
                    dy[1] = a1 * unitStep(Kd1 - D) * unitStep(Kd2 - CLK) + a2 * unitStep(Kd3 - y[0]) - dt1 * y[1];
                    //q
                    dy[2] = a3 * unitStep(y[0] - Kd4) * unitStep(CLK - Kd5) * unitStep(Kd7 - y[2]) + a4 * unitStep(Kd6 - y[3]) * unitStep(Kd7 - y[2]) - dt2 * y[2];
                    //qc
                    dy[3] = a3 * unitStep(y[1] - Kd4) * unitStep(CLK - Kd5) * unitStep(Kd7 - y[3]) + a4 * unitStep(Kd6 - y[2]) * unitStep(Kd7 - y[3]) - dt2 * y[3];

                    dy[0] *= t;
                    dy[1] *= t;
                    dy[2] *= t;
                    dy[3] *= t;

                    return dy;
                }

                // Set up external chemical level
                public void setExternalLevel(double d, double clk) {
                    D = d;
                    CLK = clk;
                }

                // Set up what the external chemical level is:
                public void setExternalQuorumLevel(double externalQuorumField) {
                    D = externalQuorumField;
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

                // Parameter Beta - ratio between mRNA and protein lifetimes
                public void generateBeta() {
                    // Garcia-Ojalvo paper part 1:
                    beta = 1.0 + 0.05 * r.nextGaussian();
                }

                public int getNumEq() {
                    return numEq;
                }
            }
        }


        /*********************************************************
         * Create the vector of all bacteria used in the simulation
         */
        final Vector<BSimDFlipFlopBacterium> bacteria = new Vector<BSimDFlipFlopBacterium>();

        // Add randomly positioned bacteria to the vector
        while (bacteria.size() < 200) {
            BSimDFlipFlopBacterium p = new BSimDFlipFlopBacterium(sim,
                    new Vector3d(Math.random() * sim.getBound().x,
                            Math.random() * sim.getBound().y,
                            Math.random() * sim.getBound().z), hField, iField, dField, qField, qcField);
            if (!p.intersection(bacteria)) bacteria.add(p);
        }


        /*********************************************************
         * Implement tick() on a BSimTicker and add the ticker to the simulation
         */
        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {

                // Update the bacteria at each time step
                for (BSimDFlipFlopBacterium p : bacteria) {
                    p.action();
                    p.updatePosition();
                }
                // Update the chemical field
                field.update();
                hField.update();
                iField.update();
                dField.update();
                qField.update();
                qcField.update();
            }
        });


        /*********************************************************
         * Implement draw(Graphics) on a BSimDrawer and add the
         * drawer to the simulation.
         *
         * Draw the particles such that they are yellow for low
         * levels of lacI mRNA and get more red as the level
         * increases.
         */
        BSimDrawer drawer = new BSimP3DDrawer(sim, 800, 600) {
            @Override
            public void scene(PGraphics3D p3d) {
                // Draw the chemical field a pretty pale blue colour
                draw(qField, new Color(Integer.parseInt("6899d3", 16)), 255 / 2);

                for (BSimDFlipFlopBacterium p : bacteria) {
                    // Oscillate from pale yellow to red.

                    int R = 255,
                            G = 255 - (int) (255 * p.y[2] * 1.0111),
                            B = 128 - (int) (128 * p.y[2] * 1.0111);
                    // Clamp these bad boys to [0, 255] to avoid errors
                    if (G < 0) G = 0;
                    if (B < 0) B = 0;
                    draw(p, new Color(R, G, B));
                }
            }
        };
        sim.setDrawer(drawer);


        /*********************************************************
         * Implement before(), during() and after() on BSimExporters
         * and add them to the simulation
         */
        // MOVIES
        BSimMovExporter movieExporter = new BSimMovExporter(sim, drawer, filePath + "D-flip-flop.mov");
        movieExporter.setSpeed(10);
        movieExporter.setDt(0.25);
        sim.addExporter(movieExporter);

        // IMAGES
        BSimPngExporter imageExporter = new BSimPngExporter(sim, drawer, filePath);
        imageExporter.setDt(10);
        sim.addExporter(imageExporter);


        /*********************************************************
         *  Create Loggers:
         *  - Simulation statistics
         *  - Time series of lacI mRNA concentration for every bacterium.
         *
         *  setDt() to reduce the amount of data.
         */
        BSimLogger stats_Logger = new BSimLogger(sim, filePath + "Settings.csv") {
            long tStart = 0;
            long tEnd = 0;

            @Override
            public void before() {
                super.before();
                tStart = Calendar.getInstance().getTimeInMillis();
                // Write parameters of the simulation
                write("Dt," + sim.getDt());
                write("Time (sec)," + sim.getSimulationTime());
                write("Diffusivity," + diffusivity);
                write("Decay rate," + decayRate);
                write("Cell wall diffusion," + cellWallDiffusivity);
                if (theInitialConditions == ICS_RANDOM) write("Initial Conditions, Random");
                else write("Initial Conditions, Uniform");

                String buffer = new String();
                buffer = "time(seconds),hFieldAvg,iFieldAvg,qFieldAvg,qcFieldAvg";
                write(buffer);
            }

            @Override
            public final void during() {
                String o = sim.getFormattedTime();

                int[] iBoxes = iField.getBoxes();
                double iConc = 0;
                for (int i = 0; i < iBoxes[0]; i++) {
                    for (int j = 0; j < iBoxes[1]; j++) {
                        iConc += iField.getConc(i, j, 0);
                    }
                }

                int[] hBoxes = hField.getBoxes();
                double hConc = 0;
                for (int i = 0; i < iBoxes[0]; i++) {
                    for (int j = 0; j < iBoxes[1]; j++) {
                        hConc += hField.getConc(i, j, 0);
                    }
                }

                int[] qBoxes = qField.getBoxes();
                double qConc = 0;
                for (int i = 0; i < qBoxes[0]; i++) {
                    for (int j = 0; j < qBoxes[1]; j++) {
                        qConc += qField.getConc(i, j, 0);
                    }
                }

                int[] qcBoxes = qcField.getBoxes();
                double qcConc = 0;
                for (int i = 0; i < qcBoxes[0]; i++) {
                    for (int j = 0; j < qcBoxes[1]; j++) {
                        qcConc += qcField.getConc(i, j, 0);
                    }
                }

                iConc /= iBoxes[0] * iBoxes[1];
                hConc /= hBoxes[0] * hBoxes[1];
                qConc /= qBoxes[0] * qBoxes[1];
                qcConc /= qcBoxes[0] * qcBoxes[1];
                write(o.replace(",", ".") + "," + String.format("%.0f", hConc) + "," + String.format("%.0f", iConc) + "," + String.format("%.0f", qConc) + "," + String.format("%.0f", qConc));

                if (sim.getSimulationTime() % 100 == 0) {
                    System.out.println("flushing");
                    try {
                        bufferedWriter.flush();
                    } catch (Exception e) {
                    }
                }
            }

            public void after() {
                // Elapsed time (real time)
                tEnd = Calendar.getInstance().getTimeInMillis();
                write("Elapsed time (sec)," + ((tEnd - tStart) / 1000.0));
                super.after();
            }
        };
        sim.addExporter(stats_Logger);


        // Print the level of lacI mRNA in all bacteria
        // WARNING: lots of data :)
        BSimLogger lacI_logger_ALL = new BSimLogger(sim, filePath + "lacI_ALL.csv") {
            @Override
            public void before() {
                super.before();
                write("time(seconds),lacI_mRNA");
            }

            @Override
            public void during() {
                String o = sim.getFormattedTime();
                String buffer = new String();
                for (int i = 0, n = bacteria.size(); i < n; i++) {
                    buffer = buffer + "," + bacteria.elementAt(i).y[2];
                }
                write(o + buffer);
            }
        };
        lacI_logger_ALL.setDt(1);            // Set export time step
        sim.addExporter(lacI_logger_ALL);


        // Print the level of internal AI in all bacteria
        // WARNING: lots of data :)
        BSimLogger AI_internal_logger_ALL = new BSimLogger(sim, filePath + "AI_internal_ALL.csv") {
            @Override
            public void before() {
                super.before();
                write("time(seconds),internal_AI");
            }

            @Override
            public void during() {
                String o = sim.getFormattedTime();
                String buffer = new String();
                for (int i = 0, n = bacteria.size(); i < n; i++) {
                    buffer = buffer + "," + bacteria.elementAt(i).y[3];
                }
                write(o + buffer);
            }
        };
        AI_internal_logger_ALL.setDt(1);    // Set export time step
        sim.addExporter(AI_internal_logger_ALL);


        /*********************************************************
         * Call sim.preview() to preview the scene
         * or sim.export() to set exporters working
         */
        //sim.preview();
        sim.export();
    }
}
