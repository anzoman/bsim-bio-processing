package BSimJohnsonCounter;

import BSimDFlipFlopKomac.SynchronousFlipFlop.ActivatorBacterium;
import BSimDFlipFlopKomac.SynchronousFlipFlop.RepressorBacterium;
import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.draw.BSimDrawer;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;
import bsim.export.BSimPngExporter;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import processing.core.PConstants;
import processing.core.PGraphics3D;

import java.awt.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

/**
 * Extended Chen oscillator in microfluidic chamber. Implemented synchronous master-slave D flip-flop
 * <p>
 * We adjust the diffusion, the boundary conditions, the density of cells, and the ratio of activators vs repressors.
 */
public class JohnsonCounterTest {

    @Parameter(names = "-export", description = "Enable export mode.")
    private boolean export = true;

    @Parameter(names = "-dim", arity = 3, description = "The dimensions (x, y, z) of simulation environment (um).")
    public List<Double> simDimensions = new ArrayList<>(Arrays.asList(100., 100., 1.));

    // Diffusion
    @Parameter(names = "-diff", arity = 1, description = "External diffusivity.")
    public double diffusivity = 80.;

    // External degradation
    @Parameter(names = "-mu_e", arity = 1, description = "External degradation.")
    public double mu_e = 0.0;

    // Boundaries
    @Parameter(names = "-fixedbounds", description = "Enable fixed boundaries. (If not, one boundary will be leaky as real uf chamber).")
    private boolean fixedBounds = false;

    @Parameter(names = "-pop", arity = 1, description = "Initial seed population (n_total).")
    public int initialPopulation = 100;

    // A:R ratio
    @Parameter(names = "-ratio", arity = 1, description = "Ratio of initial activator/repressor populations (proportion of activators).")
    public double populationRatio = 0.2;

    // Multipliers for the cell wall diffusion, and for the synthesis of QS molecules.
    @Parameter(names = "-qspars", arity = 4, description = "Multipliers for the quorum sensing parameters. [D_H, D_I, phi_H, phi_I].")
    public List<Double> qsPars = new ArrayList<>(Arrays.asList(new Double[]{1., 1., 1., 1.}));

    // Flip flop results
    SynchronousFlipFlopForJohnsonCounterTest.FlipFlopResult flipFlopResult1;
    SynchronousFlipFlopForJohnsonCounterTest.FlipFlopResult flipFlopResult2;
    SynchronousFlipFlopForJohnsonCounterTest.FlipFlopResult flipFlopResult3;

    // 3 flip flop path-result pairs
    ArrayList<FlipFlopPathResultPair> flipFlopPathResultPairs = new ArrayList<>();

    // Biochemical fields
    BSimChemicalField q1;
    BSimChemicalField q2;
    BSimChemicalField qc3;

    // Starting concentrations
    BSimChemicalField hFieldStart;
    BSimChemicalField iFieldStart;
    BSimChemicalField dFieldStart;
    BSimChemicalField qFieldStart;
    BSimChemicalField qcFieldStart;

    // Starting activators and repressors
    ArrayList<ActivatorBacterium> bacteriaActivatorsStart;
    ArrayList<RepressorBacterium> bacteriaRepressorsStart;

    public static void main(String[] args) {
        JohnsonCounterTest bsim_ex = new JohnsonCounterTest();
        new JCommander(bsim_ex, args);
        bsim_ex.run();
    }

    public void run() {
        int LOG_INTERVAL = 100;

        long simulationStartTime = System.nanoTime();
        double simX = simDimensions.get(0);
        double simY = simDimensions.get(1);
        double simZ = simDimensions.get(2);

        // create the simulation object
        BSim sim = new BSim();
        sim.setDt(0.25);                    // Simulation Timestep
        sim.setSimulationTime(43200);       // 36000 = 10 hours; 600 minutes.
        sim.setTimeFormat("0.00");            // Time Format for display
        sim.setBound(simX, simY, simZ);        // Simulation Boundaries

//        String simParameters = "" + BSimUtils.timeStamp() + "__dim_" + simX + "_" + simY + "_" + simZ
//                + "__ip_" + initialPopulation
//                + "__pr_" + populationRatio
//                + "__diff_" + diffusivity
//                + "__deg_" + mu_e
//                + "__qs_" + qsPars.get(0) + "_" + qsPars.get(1) + "_" + qsPars.get(2) + "_" + qsPars.get(3);
//
//        if (fixedBounds) {
//            simParameters += "__fixedBounds";
//        } else {
//            simParameters += "__leakyBounds";
//        }
//
//        // paths in which files will be generated
//        String filePath1 = BSimUtils.generateDirectoryPath("/Users/HP/Desktop/bsim-flip-flop/results/" + simParameters + "/flip-flop1" + "/");
//        String filePath2 = BSimUtils.generateDirectoryPath("/Users/HP/Desktop/bsim-flip-flop/results/" + simParameters + "/flip-flop2" + "/");
//        String filePath3 = BSimUtils.generateDirectoryPath("/Users/HP/Desktop/bsim-flip-flop/results/" + simParameters + "/flip-flop3" + "/");
//
//        flipFlopPathResultPairs.add(new FlipFlopPathResultPair(flipFlopResult1, filePath1));
//        flipFlopPathResultPairs.add(new FlipFlopPathResultPair(flipFlopResult2, filePath2));
//        flipFlopPathResultPairs.add(new FlipFlopPathResultPair(flipFlopResult3, filePath3));

        /*********************************************************
         * Set up the chemical fields
         */
        double external_diffusivity = diffusivity / 60.0;
        // 800/60 - repressor oscillates but activator levels out rapidly
        // 80/60 - more transients, only starts levelling at the end of the 10 hours

        // Boundaries are not periodic
        sim.setSolid(true, true, true);

        // Leaky on the bottom
        if (!fixedBounds) {
            sim.setLeaky(false, false, true, false, false, false);
            sim.setLeakyRate(0, 0, 0.1 / 60.0, 0, 0, 0);
        }

        double external_decay = mu_e / 60.0;

        qc3 = new BSimChemicalField(sim, new int[]{(int) simX, (int) simY, 1}, external_diffusivity, external_decay);

        hFieldStart = new BSimChemicalField(sim, new int[]{(int) simX, (int) simY, 1}, external_diffusivity, external_decay);
        iFieldStart = new BSimChemicalField(sim, new int[]{(int) simX, (int) simY, 1}, external_diffusivity, external_decay);
        dFieldStart = new BSimChemicalField(sim, new int[]{(int) simX, (int) simY, 1}, external_diffusivity, external_decay);
        qFieldStart = new BSimChemicalField(sim, new int[]{(int) simX, (int) simY, 1}, external_diffusivity, external_decay);
        qcFieldStart = new BSimChemicalField(sim, new int[]{(int) simX, (int) simY, 1}, external_diffusivity, external_decay);

        bacteriaActivatorsStart = new ArrayList<>();
        bacteriaRepressorsStart = new ArrayList<>();

        // ICs as in Chen paper (as in original DDEs)
        hFieldStart.setConc(10.0);
        iFieldStart.setConc(10.0);

        SynchronousFlipFlopForJohnsonCounterTest flipFlop1 = new SynchronousFlipFlopForJohnsonCounterTest(sim, qc3, LOG_INTERVAL);
        flipFlop1.create();

        SynchronousFlipFlopForJohnsonCounterTest flipFlop2 = new SynchronousFlipFlopForJohnsonCounterTest(sim, qc3, LOG_INTERVAL);
        flipFlop2.create();

        SynchronousFlipFlopForJohnsonCounterTest flipFlop3 = new SynchronousFlipFlopForJohnsonCounterTest(sim, qc3, LOG_INTERVAL);
        flipFlop3.create();


        /*********************************************************
         * Set up the ticker
         */

        BSimTicker ticker = new BSimTicker() {
            boolean toggled = true;

            @Override
            public void tick() {
                flipFlopResult1 = flipFlop1.animateFlipFlop(sim, qc3);
                q1 = flipFlopResult1.q;

                flipFlopResult2 = flipFlop2.animateFlipFlop(sim, q1);
                q2 = flipFlopResult2.q;

                flipFlopResult3 = flipFlop3.animateFlipFlop(sim, q2);
                qc3 = flipFlopResult3.qc;
            }
        };
        sim.setTicker(ticker);

        /*********************************************************
         * Set up the drawer
         */
        BSimDrawer drawer1 = new BSimP3DDrawer(sim, 800, 600) {
            /**
             * Draw the default cuboid boundary of the simulation as a partially transparent box
             * with a wireframe outline surrounding it.
             */
            @Override
            public void boundaries() {
                p3d.noFill();
                p3d.stroke(128, 128, 255);
                p3d.pushMatrix();
                p3d.translate((float) boundCentre.x, (float) boundCentre.y, (float) boundCentre.z);
                p3d.box((float) bound.x, (float) bound.y, (float) bound.z);
                p3d.popMatrix();
                p3d.noStroke();
            }

            /**
             * Original Chen/Bennet oscillator drawer function
             */
            @Override
            public void draw(Graphics2D g) {
                p3d.beginDraw();

                if (!cameraIsInitialised) {
                    // camera(eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ)
                    p3d.camera((float) bound.x * 0.5f, (float) bound.y * 0.5f,
                            // Set the Z offset to the largest of X/Y dimensions for a reasonable zoom-out distance:
                            simX > simY ? (float) simX : (float) simY,
//                            10,
                            (float) bound.x * 0.5f, (float) bound.y * 0.5f, 0,
                            0, 1, 0);
                    cameraIsInitialised = true;
                }

                p3d.textFont(font);
                p3d.textMode(PConstants.SCREEN);

                p3d.sphereDetail(10);
                p3d.noStroke();
                p3d.background(255, 255, 255);

                scene(p3d);
                boundaries();
                time();

                p3d.endDraw();
                g.drawImage(p3d.image, 0, 0, null);
            }

            /**
             * Draw the formatted simulation time to screen.
             */
            @Override
            public void time() {
                p3d.fill(0);
                p3d.text(sim.getFormattedTimeHours(), 50, 50);
            }

            /**
             * Original Chen/Bennet oscillator drawer function
             */
            @Override
            public void scene(PGraphics3D p3d) {
                p3d.ambientLight(128, 128, 128);
                p3d.directionalLight(128, 128, 128, 1, 1, -1);

                ArrayList<ActivatorBacterium> bacteriaActivators;
                ArrayList<RepressorBacterium> bacteriaRepressors;

                if (flipFlopResult1 != null) {
                    bacteriaActivators = flipFlopResult1.bacteriaActivators;
                    bacteriaRepressors = flipFlopResult1.bacteriaRepressors;
                } else {
                    bacteriaActivators = bacteriaActivatorsStart;
                    bacteriaRepressors = bacteriaRepressorsStart;
                }

                for (ActivatorBacterium b : bacteriaActivators) {
                    // Activators: blue 55,126,184
                    int R = 30 + (int) (25 * b.grn_state[4] / 3.5e4),
                            G = 30 + (int) (100 * b.grn_state[4] / 3.5e4),
                            B = 55 + (int) (135 * b.grn_state[4] / 3.5e4);
                    // Clamp these to [0, 255] to avoid errors
                    if (R < 0) R = 0;
                    if (R > 255) R = 255;
                    if (G < 0) G = 0;
                    if (G > 255) G = 255;
                    if (B < 0) B = 0;
                    if (B > 255) B = 255;
                    draw(b, new Color(R, G, B));
//                    draw(b, Color.blue);
                }

                for (RepressorBacterium b : bacteriaRepressors) {
                    // Repressors: green 77,175,74
                    int R = 30 + (int) (50 * b.grn_state[4] / 4.5e4),
                            G = 55 + (int) (120 * b.grn_state[4] / 4.5e4),
                            B = 30 + (int) (50 * b.grn_state[4] / 4.5e4);
                    // Clamp these to [0, 255] to avoid errors
                    if (R < 0) R = 0;
                    if (R > 255) R = 255;
                    if (G < 0) G = 0;
                    if (G > 255) G = 255;
                    if (B < 0) B = 0;
                    if (B > 255) B = 255;
                    draw(b, new Color(R, G, B));
//                    draw(b, Color.green);
                }

//                draw(reference_field, Color.BLUE, (float)(255/c));
            }
        };
        sim.setDrawer(drawer1);

        BSimDrawer drawer2 = new BSimP3DDrawer(sim, 800, 600) {
            /**
             * Draw the default cuboid boundary of the simulation as a partially transparent box
             * with a wireframe outline surrounding it.
             */
            @Override
            public void boundaries() {
                p3d.noFill();
                p3d.stroke(128, 128, 255);
                p3d.pushMatrix();
                p3d.translate((float) boundCentre.x, (float) boundCentre.y, (float) boundCentre.z);
                p3d.box((float) bound.x, (float) bound.y, (float) bound.z);
                p3d.popMatrix();
                p3d.noStroke();
            }

            /**
             * Original Chen/Bennet oscillator drawer function
             */
            @Override
            public void draw(Graphics2D g) {
                p3d.beginDraw();

                if (!cameraIsInitialised) {
                    // camera(eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ)
                    p3d.camera((float) bound.x * 0.5f, (float) bound.y * 0.5f,
                            // Set the Z offset to the largest of X/Y dimensions for a reasonable zoom-out distance:
                            simX > simY ? (float) simX : (float) simY,
//                            10,
                            (float) bound.x * 0.5f, (float) bound.y * 0.5f, 0,
                            0, 1, 0);
                    cameraIsInitialised = true;
                }

                p3d.textFont(font);
                p3d.textMode(PConstants.SCREEN);

                p3d.sphereDetail(10);
                p3d.noStroke();
                p3d.background(255, 255, 255);

                scene(p3d);
                boundaries();
                time();

                p3d.endDraw();
                g.drawImage(p3d.image, 0, 0, null);
            }

            /**
             * Draw the formatted simulation time to screen.
             */
            @Override
            public void time() {
                p3d.fill(0);
                p3d.text(sim.getFormattedTimeHours(), 50, 50);
            }

            /**
             * Original Chen/Bennet oscillator drawer function
             */
            @Override
            public void scene(PGraphics3D p3d) {
                p3d.ambientLight(128, 128, 128);
                p3d.directionalLight(128, 128, 128, 1, 1, -1);

                ArrayList<ActivatorBacterium> bacteriaActivators;
                ArrayList<RepressorBacterium> bacteriaRepressors;

                if (flipFlopResult2 != null) {
                    bacteriaActivators = flipFlopResult2.bacteriaActivators;
                    bacteriaRepressors = flipFlopResult2.bacteriaRepressors;
                } else {
                    bacteriaActivators = bacteriaActivatorsStart;
                    bacteriaRepressors = bacteriaRepressorsStart;
                }

                for (ActivatorBacterium b : bacteriaActivators) {
                    // Activators: blue 55,126,184
                    int R = 30 + (int) (25 * b.grn_state[4] / 3.5e4),
                            G = 30 + (int) (100 * b.grn_state[4] / 3.5e4),
                            B = 55 + (int) (135 * b.grn_state[4] / 3.5e4);
                    // Clamp these to [0, 255] to avoid errors
                    if (R < 0) R = 0;
                    if (R > 255) R = 255;
                    if (G < 0) G = 0;
                    if (G > 255) G = 255;
                    if (B < 0) B = 0;
                    if (B > 255) B = 255;
                    draw(b, new Color(R, G, B));
//                    draw(b, Color.blue);
                }

                for (RepressorBacterium b : bacteriaRepressors) {
                    // Repressors: green 77,175,74
                    int R = 30 + (int) (50 * b.grn_state[4] / 4.5e4),
                            G = 55 + (int) (120 * b.grn_state[4] / 4.5e4),
                            B = 30 + (int) (50 * b.grn_state[4] / 4.5e4);
                    // Clamp these to [0, 255] to avoid errors
                    if (R < 0) R = 0;
                    if (R > 255) R = 255;
                    if (G < 0) G = 0;
                    if (G > 255) G = 255;
                    if (B < 0) B = 0;
                    if (B > 255) B = 255;
                    draw(b, new Color(R, G, B));
//                    draw(b, Color.green);
                }

//                draw(reference_field, Color.BLUE, (float)(255/c));
            }
        };
        sim.setDrawer(drawer2);

        BSimDrawer drawer3 = new BSimP3DDrawer(sim, 800, 600) {
            /**
             * Draw the default cuboid boundary of the simulation as a partially transparent box
             * with a wireframe outline surrounding it.
             */
            @Override
            public void boundaries() {
                p3d.noFill();
                p3d.stroke(128, 128, 255);
                p3d.pushMatrix();
                p3d.translate((float) boundCentre.x, (float) boundCentre.y, (float) boundCentre.z);
                p3d.box((float) bound.x, (float) bound.y, (float) bound.z);
                p3d.popMatrix();
                p3d.noStroke();
            }

            /**
             * Original Chen/Bennet oscillator drawer function
             */
            @Override
            public void draw(Graphics2D g) {
                p3d.beginDraw();

                if (!cameraIsInitialised) {
                    // camera(eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ)
                    p3d.camera((float) bound.x * 0.5f, (float) bound.y * 0.5f,
                            // Set the Z offset to the largest of X/Y dimensions for a reasonable zoom-out distance:
                            simX > simY ? (float) simX : (float) simY,
//                            10,
                            (float) bound.x * 0.5f, (float) bound.y * 0.5f, 0,
                            0, 1, 0);
                    cameraIsInitialised = true;
                }

                p3d.textFont(font);
                p3d.textMode(PConstants.SCREEN);

                p3d.sphereDetail(10);
                p3d.noStroke();
                p3d.background(255, 255, 255);

                scene(p3d);
                boundaries();
                time();

                p3d.endDraw();
                g.drawImage(p3d.image, 0, 0, null);
            }

            /**
             * Draw the formatted simulation time to screen.
             */
            @Override
            public void time() {
                p3d.fill(0);
                p3d.text(sim.getFormattedTimeHours(), 50, 50);
            }

            /**
             * Original Chen/Bennet oscillator drawer function
             */
            @Override
            public void scene(PGraphics3D p3d) {
                p3d.ambientLight(128, 128, 128);
                p3d.directionalLight(128, 128, 128, 1, 1, -1);

                ArrayList<ActivatorBacterium> bacteriaActivators;
                ArrayList<RepressorBacterium> bacteriaRepressors;

                if (flipFlopResult3 != null) {
                    bacteriaActivators = flipFlopResult3.bacteriaActivators;
                    bacteriaRepressors = flipFlopResult3.bacteriaRepressors;
                } else {
                    bacteriaActivators = bacteriaActivatorsStart;
                    bacteriaRepressors = bacteriaRepressorsStart;
                }

                for (ActivatorBacterium b : bacteriaActivators) {
                    // Activators: blue 55,126,184
                    int R = 30 + (int) (25 * b.grn_state[4] / 3.5e4),
                            G = 30 + (int) (100 * b.grn_state[4] / 3.5e4),
                            B = 55 + (int) (135 * b.grn_state[4] / 3.5e4);
                    // Clamp these to [0, 255] to avoid errors
                    if (R < 0) R = 0;
                    if (R > 255) R = 255;
                    if (G < 0) G = 0;
                    if (G > 255) G = 255;
                    if (B < 0) B = 0;
                    if (B > 255) B = 255;
                    draw(b, new Color(R, G, B));
//                    draw(b, Color.blue);
                }

                for (RepressorBacterium b : bacteriaRepressors) {
                    // Repressors: green 77,175,74
                    int R = 30 + (int) (50 * b.grn_state[4] / 4.5e4),
                            G = 55 + (int) (120 * b.grn_state[4] / 4.5e4),
                            B = 30 + (int) (50 * b.grn_state[4] / 4.5e4);
                    // Clamp these to [0, 255] to avoid errors
                    if (R < 0) R = 0;
                    if (R > 255) R = 255;
                    if (G < 0) G = 0;
                    if (G > 255) G = 255;
                    if (B < 0) B = 0;
                    if (B > 255) B = 255;
                    draw(b, new Color(R, G, B));
//                    draw(b, Color.green);
                }

//                draw(reference_field, Color.BLUE, (float)(255/c));
            }
        };
        sim.setDrawer(drawer3);

        if (export) {
            String simParameters = "" + BSimUtils.timeStamp() + "__dim_" + simX + "_" + simY + "_" + simZ
                    + "__ip_" + initialPopulation
                    + "__pr_" + populationRatio
                    + "__diff_" + diffusivity
                    + "__deg_" + mu_e
                    + "__qs_" + qsPars.get(0) + "_" + qsPars.get(1) + "_" + qsPars.get(2) + "_" + qsPars.get(3);

            if (fixedBounds) {
                simParameters += "__fixedBounds";
            } else {
                simParameters += "__leakyBounds";
            }

            // paths in which files will be generated
            String filePath1 = BSimUtils.generateDirectoryPath("/Users/HP/Desktop/bsim-flip-flop/results/" + simParameters + "/flip-flop1" + "/");
            String filePath2 = BSimUtils.generateDirectoryPath("/Users/HP/Desktop/bsim-flip-flop/results/" + simParameters + "/flip-flop2" + "/");
            String filePath3 = BSimUtils.generateDirectoryPath("/Users/HP/Desktop/bsim-flip-flop/results/" + simParameters + "/flip-flop3" + "/");

            flipFlopPathResultPairs = new ArrayList<>();
            flipFlopPathResultPairs.add(new FlipFlopPathResultPair(flipFlopResult1, filePath1));
            flipFlopPathResultPairs.add(new FlipFlopPathResultPair(flipFlopResult2, filePath2));
            flipFlopPathResultPairs.add(new FlipFlopPathResultPair(flipFlopResult3, filePath3));


            /*********************************************************
             * Various properties of the simulation, for future reference.
             */
            BSimLogger metaLogger1 = new BSimLogger(sim, filePath1 + "simInfo.txt") {
                @Override
                public void before() {
                    super.before();
                    write("Simulation metadata.");
                    write("Chen/Bennett consortium oscillator system.");
                    write("Simulation dimensions: (" + simX + ", " + simY + ", " + simZ + ")");
                    write("Initial population: " + initialPopulation);
                    write("Ratio " + populationRatio);
                    write("Spatial signalling diffusivity: " + diffusivity);
                    write("Spatial degradation (mu_e): " + mu_e);

                    if (fixedBounds) {
                        write("Boundaries: fixed");
                    } else {
                        write("Boundaries: leaky");
                    }

                    write("Multiplier D_H: " + qsPars.get(0));
                    write("Multiplier D_I: " + qsPars.get(1));
                    write("Multiplier phi_H: " + qsPars.get(2));
                    write("Multiplier phi_I: " + qsPars.get(3));
                }

                @Override
                public void during() {

                }
            };
            metaLogger1.setDt(3600);            // Set export time step
            sim.addExporter(metaLogger1);

            BSimLogger metaLogger2 = new BSimLogger(sim, filePath2 + "simInfo.txt") {
                @Override
                public void before() {
                    super.before();
                    write("Simulation metadata.");
                    write("Chen/Bennett consortium oscillator system.");
                    write("Simulation dimensions: (" + simX + ", " + simY + ", " + simZ + ")");
                    write("Initial population: " + initialPopulation);
                    write("Ratio " + populationRatio);
                    write("Spatial signalling diffusivity: " + diffusivity);
                    write("Spatial degradation (mu_e): " + mu_e);

                    if (fixedBounds) {
                        write("Boundaries: fixed");
                    } else {
                        write("Boundaries: leaky");
                    }

                    write("Multiplier D_H: " + qsPars.get(0));
                    write("Multiplier D_I: " + qsPars.get(1));
                    write("Multiplier phi_H: " + qsPars.get(2));
                    write("Multiplier phi_I: " + qsPars.get(3));
                }

                @Override
                public void during() {

                }
            };
            metaLogger2.setDt(3600);            // Set export time step
            sim.addExporter(metaLogger2);

            BSimLogger metaLogger3 = new BSimLogger(sim, filePath3 + "simInfo.txt") {
                @Override
                public void before() {
                    super.before();
                    write("Simulation metadata.");
                    write("Chen/Bennett consortium oscillator system.");
                    write("Simulation dimensions: (" + simX + ", " + simY + ", " + simZ + ")");
                    write("Initial population: " + initialPopulation);
                    write("Ratio " + populationRatio);
                    write("Spatial signalling diffusivity: " + diffusivity);
                    write("Spatial degradation (mu_e): " + mu_e);

                    if (fixedBounds) {
                        write("Boundaries: fixed");
                    } else {
                        write("Boundaries: leaky");
                    }

                    write("Multiplier D_H: " + qsPars.get(0));
                    write("Multiplier D_I: " + qsPars.get(1));
                    write("Multiplier phi_H: " + qsPars.get(2));
                    write("Multiplier phi_I: " + qsPars.get(3));
                }

                @Override
                public void during() {

                }
            };
            metaLogger3.setDt(3600);            // Set export time step
            sim.addExporter(metaLogger3);

            BSimLogger dataLoggerConc1 = new BSimLogger(sim, filePath1 + "concentrations_average.csv") {
                @Override
                public void before() {
                    super.before();
                    String buffer = new String();
                    buffer = "time(seconds),activatory proteins(h),repressory proteins(i),Q proteins(q),Qc proteins(qc)";
                    write(buffer);
                }

                @Override
                public void during() {
                    String o = sim.getFormattedTime();
                    String buffer = new String();
                    //Contentration in the middle

                    BSimChemicalField iField;
                    BSimChemicalField hField;
                    BSimChemicalField qField;
                    BSimChemicalField qcField;

                    if (flipFlopResult1 != null) {
                        iField = flipFlopResult1.i;
                        hField = flipFlopResult1.h;
                        qField = flipFlopResult1.q;
                        qcField = flipFlopResult1.qc;
                    } else {
                        iField = iFieldStart;
                        hField = hFieldStart;
                        qField = qFieldStart;
                        qcField = qcFieldStart;
                    }

                    int[] i_boxes = iField.getBoxes();
                    //Average concentrations
                    double i_conc = 0;
                    for (int i = 0; i < i_boxes[0]; i++) {
                        for (int j = 0; j < i_boxes[1]; j++) {
                            i_conc += iField.getConc(i, j, 0);
                        }

                    }

                    int[] h_boxes = hField.getBoxes();
                    //Average concentrations
                    double h_conc = 0;
                    for (int i = 0; i < i_boxes[0]; i++) {
                        for (int j = 0; j < i_boxes[1]; j++) {
                            h_conc += hField.getConc(i, j, 0);
                        }

                    }

                    int[] q_boxes = qField.getBoxes();
                    //Average concentrations
                    double q_conc = 0;
                    for (int i = 0; i < q_boxes[0]; i++) {
                        for (int j = 0; j < q_boxes[1]; j++) {
                            q_conc += qField.getConc(i, j, 0);
                        }

                    }

                    int[] qc_boxes = qcField.getBoxes();
                    //Average concentrations
                    double qc_conc = 0;
                    for (int i = 0; i < qc_boxes[0]; i++) {
                        for (int j = 0; j < qc_boxes[1]; j++) {
                            qc_conc += qcField.getConc(i, j, 0);
                        }

                    }

                    i_conc /= i_boxes[0] * i_boxes[1];
                    h_conc /= h_boxes[0] * h_boxes[1];
                    q_conc /= q_boxes[0] * q_boxes[1];
                    qc_conc /= qc_boxes[0] * qc_boxes[1];
                    write(o + "," + h_conc + "," + i_conc + "," + q_conc + "," + qc_conc);

                    if (sim.getSimulationTime() % 100 == 0) {
                        System.out.println("flushing");
                        try {
                            bufferedWriter.flush();
                        } catch (Exception e) {
                        }
                    }
                }
            };
            dataLoggerConc1.setDt(30);
            sim.addExporter(dataLoggerConc1);

            BSimLogger dataLoggerConc2 = new BSimLogger(sim, filePath2 + "concentrations_average.csv") {
                @Override
                public void before() {
                    super.before();
                    String buffer = new String();
                    buffer = "time(seconds),activatory proteins(h),repressory proteins(i),Q proteins(q),Qc proteins(qc)";
                    write(buffer);
                }

                @Override
                public void during() {
                    String o = sim.getFormattedTime();
                    String buffer = new String();
                    //Contentration in the middle

                    BSimChemicalField iField;
                    BSimChemicalField hField;
                    BSimChemicalField qField;
                    BSimChemicalField qcField;

                    if (flipFlopResult2 != null) {
                        iField = flipFlopResult2.i;
                        hField = flipFlopResult2.h;
                        qField = flipFlopResult2.q;
                        qcField = flipFlopResult2.qc;
                    } else {
                        iField = iFieldStart;
                        hField = hFieldStart;
                        qField = qFieldStart;
                        qcField = qcFieldStart;
                    }

                    int[] i_boxes = iField.getBoxes();
                    //Average concentrations
                    double i_conc = 0;
                    for (int i = 0; i < i_boxes[0]; i++) {
                        for (int j = 0; j < i_boxes[1]; j++) {
                            i_conc += iField.getConc(i, j, 0);
                        }

                    }

                    int[] h_boxes = hField.getBoxes();
                    //Average concentrations
                    double h_conc = 0;
                    for (int i = 0; i < i_boxes[0]; i++) {
                        for (int j = 0; j < i_boxes[1]; j++) {
                            h_conc += hField.getConc(i, j, 0);
                        }

                    }

                    int[] q_boxes = qField.getBoxes();
                    //Average concentrations
                    double q_conc = 0;
                    for (int i = 0; i < q_boxes[0]; i++) {
                        for (int j = 0; j < q_boxes[1]; j++) {
                            q_conc += qField.getConc(i, j, 0);
                        }

                    }

                    int[] qc_boxes = qcField.getBoxes();
                    //Average concentrations
                    double qc_conc = 0;
                    for (int i = 0; i < qc_boxes[0]; i++) {
                        for (int j = 0; j < qc_boxes[1]; j++) {
                            qc_conc += qcField.getConc(i, j, 0);
                        }

                    }

                    i_conc /= i_boxes[0] * i_boxes[1];
                    h_conc /= h_boxes[0] * h_boxes[1];
                    q_conc /= q_boxes[0] * q_boxes[1];
                    qc_conc /= qc_boxes[0] * qc_boxes[1];
                    write(o + "," + h_conc + "," + i_conc + "," + q_conc + "," + qc_conc);

                    if (sim.getSimulationTime() % 100 == 0) {
                        System.out.println("flushing");
                        try {
                            bufferedWriter.flush();
                        } catch (Exception e) {
                        }
                    }
                }
            };
            dataLoggerConc2.setDt(30);
            sim.addExporter(dataLoggerConc2);

            BSimLogger dataLoggerConc3 = new BSimLogger(sim, filePath3 + "concentrations_average.csv") {
                @Override
                public void before() {
                    super.before();
                    String buffer = new String();
                    buffer = "time(seconds),activatory proteins(h),repressory proteins(i),Q proteins(q),Qc proteins(qc)";
                    write(buffer);
                }

                @Override
                public void during() {
                    String o = sim.getFormattedTime();
                    String buffer = new String();
                    //Contentration in the middle

                    BSimChemicalField iField;
                    BSimChemicalField hField;
                    BSimChemicalField qField;
                    BSimChemicalField qcField;

                    if (flipFlopResult3 != null) {
                        iField = flipFlopResult3.i;
                        hField = flipFlopResult3.h;
                        qField = flipFlopResult3.q;
                        qcField = flipFlopResult3.qc;
                    } else {
                        iField = iFieldStart;
                        hField = hFieldStart;
                        qField = qFieldStart;
                        qcField = qcFieldStart;
                    }

                    int[] i_boxes = iField.getBoxes();
                    //Average concentrations
                    double i_conc = 0;
                    for (int i = 0; i < i_boxes[0]; i++) {
                        for (int j = 0; j < i_boxes[1]; j++) {
                            i_conc += iField.getConc(i, j, 0);
                        }

                    }

                    int[] h_boxes = hField.getBoxes();
                    //Average concentrations
                    double h_conc = 0;
                    for (int i = 0; i < i_boxes[0]; i++) {
                        for (int j = 0; j < i_boxes[1]; j++) {
                            h_conc += hField.getConc(i, j, 0);
                        }

                    }

                    int[] q_boxes = qField.getBoxes();
                    //Average concentrations
                    double q_conc = 0;
                    for (int i = 0; i < q_boxes[0]; i++) {
                        for (int j = 0; j < q_boxes[1]; j++) {
                            q_conc += qField.getConc(i, j, 0);
                        }

                    }

                    int[] qc_boxes = qcField.getBoxes();
                    //Average concentrations
                    double qc_conc = 0;
                    for (int i = 0; i < qc_boxes[0]; i++) {
                        for (int j = 0; j < qc_boxes[1]; j++) {
                            qc_conc += qcField.getConc(i, j, 0);
                        }

                    }

                    i_conc /= i_boxes[0] * i_boxes[1];
                    h_conc /= h_boxes[0] * h_boxes[1];
                    q_conc /= q_boxes[0] * q_boxes[1];
                    qc_conc /= qc_boxes[0] * qc_boxes[1];
                    write(o + "," + h_conc + "," + i_conc + "," + q_conc + "," + qc_conc);

                    if (sim.getSimulationTime() % 100 == 0) {
                        System.out.println("flushing");
                        try {
                            bufferedWriter.flush();
                        } catch (Exception e) {
                        }
                    }
                }
            };
            dataLoggerConc3.setDt(30);
            sim.addExporter(dataLoggerConc3);

            BSimLogger posLogger1 = new BSimLogger(sim, filePath1 + "position.csv") {
                DecimalFormat formatter = new DecimalFormat("###.##", DecimalFormatSymbols.getInstance(Locale.ENGLISH));

                @Override
                public void before() {
                    super.before();
                    write("per Act; per Rep; id, p1x, p1y, p1z, p2x, p2y, p2z");
                }

                @Override
                public void during() {
                    String buffer = new String();

                    buffer += sim.getFormattedTime() + "\n";
                    write(buffer);

                    write("acts");

                    ArrayList<ActivatorBacterium> bacteriaActivators;
                    ArrayList<RepressorBacterium> bacteriaRepressors;

                    if (flipFlopResult1 != null) {
                        bacteriaActivators = flipFlopResult1.bacteriaActivators;
                        bacteriaRepressors = flipFlopResult1.bacteriaRepressors;
                    } else {
                        bacteriaActivators = bacteriaActivatorsStart;
                        bacteriaRepressors = bacteriaRepressorsStart;
                    }

                    buffer = "";
                    for (BSimCapsuleBacterium b : bacteriaActivators) {
                        buffer += b.id + "," + formatter.format(b.x1.x)
                                + "," + formatter.format(b.x1.y)
                                + "," + formatter.format(b.x1.z)
                                + "," + formatter.format(b.x2.x)
                                + "," + formatter.format(b.x2.y)
                                + "," + formatter.format(b.x2.z)
                                + "\n";
                    }

                    write(buffer);

                    write("reps");

                    buffer = "";
                    for (BSimCapsuleBacterium b : bacteriaRepressors) {
                        buffer += b.id + "," + formatter.format(b.x1.x)
                                + "," + formatter.format(b.x1.y)
                                + "," + formatter.format(b.x1.z)
                                + "," + formatter.format(b.x2.x)
                                + "," + formatter.format(b.x2.y)
                                + "," + formatter.format(b.x2.z)
                                + "\n";
                    }

                    write(buffer);

                }
            };
            posLogger1.setDt(30);            // Set export time step
            sim.addExporter(posLogger1);

            BSimLogger posLogger2 = new BSimLogger(sim, filePath2 + "position.csv") {
                DecimalFormat formatter = new DecimalFormat("###.##", DecimalFormatSymbols.getInstance(Locale.ENGLISH));

                @Override
                public void before() {
                    super.before();
                    write("per Act; per Rep; id, p1x, p1y, p1z, p2x, p2y, p2z");
                }

                @Override
                public void during() {
                    String buffer = new String();

                    buffer += sim.getFormattedTime() + "\n";
                    write(buffer);

                    write("acts");

                    ArrayList<ActivatorBacterium> bacteriaActivators;
                    ArrayList<RepressorBacterium> bacteriaRepressors;

                    if (flipFlopResult2 != null) {
                        bacteriaActivators = flipFlopResult2.bacteriaActivators;
                        bacteriaRepressors = flipFlopResult2.bacteriaRepressors;
                    } else {
                        bacteriaActivators = bacteriaActivatorsStart;
                        bacteriaRepressors = bacteriaRepressorsStart;
                    }

                    buffer = "";
                    for (BSimCapsuleBacterium b : bacteriaActivators) {
                        buffer += b.id + "," + formatter.format(b.x1.x)
                                + "," + formatter.format(b.x1.y)
                                + "," + formatter.format(b.x1.z)
                                + "," + formatter.format(b.x2.x)
                                + "," + formatter.format(b.x2.y)
                                + "," + formatter.format(b.x2.z)
                                + "\n";
                    }

                    write(buffer);

                    write("reps");

                    buffer = "";
                    for (BSimCapsuleBacterium b : bacteriaRepressors) {
                        buffer += b.id + "," + formatter.format(b.x1.x)
                                + "," + formatter.format(b.x1.y)
                                + "," + formatter.format(b.x1.z)
                                + "," + formatter.format(b.x2.x)
                                + "," + formatter.format(b.x2.y)
                                + "," + formatter.format(b.x2.z)
                                + "\n";
                    }

                    write(buffer);

                }
            };
            posLogger2.setDt(30);            // Set export time step
            sim.addExporter(posLogger2);

            BSimLogger posLogger3 = new BSimLogger(sim, filePath3 + "position.csv") {
                DecimalFormat formatter = new DecimalFormat("###.##", DecimalFormatSymbols.getInstance(Locale.ENGLISH));

                @Override
                public void before() {
                    super.before();
                    write("per Act; per Rep; id, p1x, p1y, p1z, p2x, p2y, p2z");
                }

                @Override
                public void during() {
                    String buffer = new String();

                    buffer += sim.getFormattedTime() + "\n";
                    write(buffer);

                    write("acts");

                    ArrayList<ActivatorBacterium> bacteriaActivators;
                    ArrayList<RepressorBacterium> bacteriaRepressors;

                    if (flipFlopResult3 != null) {
                        bacteriaActivators = flipFlopResult3.bacteriaActivators;
                        bacteriaRepressors = flipFlopResult3.bacteriaRepressors;
                    } else {
                        bacteriaActivators = bacteriaActivatorsStart;
                        bacteriaRepressors = bacteriaRepressorsStart;
                    }

                    buffer = "";
                    for (BSimCapsuleBacterium b : bacteriaActivators) {
                        buffer += b.id + "," + formatter.format(b.x1.x)
                                + "," + formatter.format(b.x1.y)
                                + "," + formatter.format(b.x1.z)
                                + "," + formatter.format(b.x2.x)
                                + "," + formatter.format(b.x2.y)
                                + "," + formatter.format(b.x2.z)
                                + "\n";
                    }

                    write(buffer);

                    write("reps");

                    buffer = "";
                    for (BSimCapsuleBacterium b : bacteriaRepressors) {
                        buffer += b.id + "," + formatter.format(b.x1.x)
                                + "," + formatter.format(b.x1.y)
                                + "," + formatter.format(b.x1.z)
                                + "," + formatter.format(b.x2.x)
                                + "," + formatter.format(b.x2.y)
                                + "," + formatter.format(b.x2.z)
                                + "\n";
                    }

                    write(buffer);

                }
            };
            posLogger3.setDt(30);            // Set export time step
            sim.addExporter(posLogger3);

            /**
             * Export a rendered image file
             */
            BSimPngExporter imageExporter1 = new BSimPngExporter(sim, drawer1, filePath1);
            System.out.println("Exporting pngs to:");
            System.out.println(filePath1);
            imageExporter1.setDt(30);
            sim.addExporter(imageExporter1);

            BSimPngExporter imageExporter2 = new BSimPngExporter(sim, drawer2, filePath2);
            System.out.println("Exporting pngs to:");
            System.out.println(filePath2);
            imageExporter2.setDt(30);
            sim.addExporter(imageExporter2);

            BSimPngExporter imageExporter3 = new BSimPngExporter(sim, drawer3, filePath3);
            System.out.println("Exporting pngs to:");
            System.out.println(filePath3);
            imageExporter3.setDt(30);
            sim.addExporter(imageExporter3);

            sim.export();

        } else {
            sim.preview();
        }

        long simulationEndTime = System.nanoTime();

        System.out.println("Total simulation time: " + (simulationEndTime - simulationStartTime) / 1e9 + " sec.");
    }

    public static class FlipFlopPathResultPair {
        SynchronousFlipFlopForJohnsonCounterTest.FlipFlopResult flipFlopResult;
        String filePath;

        public FlipFlopPathResultPair(SynchronousFlipFlopForJohnsonCounterTest.FlipFlopResult flipFlopResult, String filePath) {
            this.flipFlopResult = flipFlopResult;
            this.filePath = filePath;
        }
    }
}
