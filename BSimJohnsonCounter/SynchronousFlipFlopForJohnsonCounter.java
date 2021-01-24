package BSimJohnsonCounter;

import BSimDFlipFlopKomac.SynchronousFlipFlop.ActivatorBacterium;
import BSimDFlipFlopKomac.SynchronousFlipFlop.BSimDBacterium;
import BSimDFlipFlopKomac.SynchronousFlipFlop.ChenParameters;
import BSimDFlipFlopKomac.SynchronousFlipFlop.RepressorBacterium;
import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.capsule.Mover;
import bsim.capsule.RelaxationMoverGrid;
import bsim.draw.BSimDrawer;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;
import bsim.export.BSimPngExporter;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import processing.core.PConstants;
import processing.core.PGraphics3D;

import javax.vecmath.Vector3d;
import java.awt.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.*;

/**
 * Extended Chen oscillator in microfluidic chamber. Implemented synchronous master-slave D flip-flop
 *
 * We adjust the diffusion, the boundary conditions, the density of cells, and the ratio of activators vs repressors.
 */
public class SynchronousFlipFlopForJohnsonCounter {

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
    public List<Double> qsPars = new ArrayList<>(Arrays.asList(new Double[] {1., 1., 1., 1.}));


    /**
     * Whether to enable growth
     */
    private static final boolean WITH_GROWTH = false;


    public static void main(String[] args) {
        SynchronousFlipFlopForJohnsonCounter bsim_ex = new SynchronousFlipFlopForJohnsonCounter();

        new JCommander(bsim_ex, args);

        //bsim_ex.run();
    }

    public FlipFlopResult run(BSim sim, BSimChemicalField d) {
        double simX = simDimensions.get(0);
        double simY = simDimensions.get(1);
        double simZ = simDimensions.get(2);

        int nActivatorStart = (int)Math.round(populationRatio*initialPopulation);
        int nRepressorStart = (int)Math.round((1 - populationRatio)*initialPopulation);
        int nD1Bacterium = initialPopulation;

        long simulationStartTime = System.nanoTime();

        double def_D_H = ChenParameters.p.get("D_H");
        ChenParameters.p.put("D_H", def_D_H*qsPars.get(0));
        double def_D_I = ChenParameters.p.get("D_I");
        ChenParameters.p.put("D_I", def_D_I*qsPars.get(1));

        double def_phi_H = ChenParameters.p.get("phi_H");
        ChenParameters.p.put("phi_H", def_phi_H*qsPars.get(2));
        double def_phi_I = ChenParameters.p.get("phi_I");
        ChenParameters.p.put("phi_I", def_phi_I*qsPars.get(3));


        /*********************************************************
         * Set up the chemical fields
         */
        double external_diffusivity = diffusivity/60.0;
        // 800/60 - repressor oscillates but activator levels out rapidly
        // 80/60 - more transients, only starts levelling at the end of the 10 hours

        // Boundaries are not periodic
        sim.setSolid(true, true, true);

        // Leaky on the bottom
        if(!fixedBounds) {
            sim.setLeaky(false, false, true, false, false, false);
            sim.setLeakyRate(0, 0, 0.1/60.0, 0, 0, 0);
        }

        double external_decay = mu_e/60.0;

        BSimChemicalField h_e_field  = new BSimChemicalField(sim, new int[] {(int) simX, (int)simY, 1}, external_diffusivity, external_decay);
        BSimChemicalField i_e_field  = new BSimChemicalField(sim, new int[] {(int) simX, (int)simY, 1}, external_diffusivity, external_decay);

        BSimChemicalField d_e_field  = d;
        BSimChemicalField q_e_field  = new BSimChemicalField(sim, new int[] {(int) simX, (int)simY, 1}, external_diffusivity, external_decay);
        BSimChemicalField qc_e_field = new BSimChemicalField(sim, new int[] {(int) simX, (int)simY, 1}, external_diffusivity, external_decay);

        // ICs as in Chen paper (as in original DDEs)
        h_e_field.setConc(10.0);
        i_e_field.setConc(10.0);


        /*********************************************************
         * Create the bacteria
         */
        // Separate lists of bacteria in case we want to manipulate the species individually
        final ArrayList<ActivatorBacterium> bacteriaActivators = new ArrayList();
        final ArrayList<RepressorBacterium> bacteriaRepressors = new ArrayList();
        final ArrayList<BSimDBacterium> bacteriaD = new ArrayList();

        // Track all of the bacteria in the simulation, for use of common methods etc
        final ArrayList<BSimCapsuleBacterium> bacteriaAll = new ArrayList();


        Random bacRng = new Random();

        generator:
        while(bacteriaActivators.size() < nActivatorStart) {
            double bL = 1. + 0.1*(bacRng.nextDouble() - 0.5);
            double angle = bacRng.nextDouble()*2*Math.PI;

            Vector3d pos = new Vector3d(1.1 + bacRng.nextDouble()*(sim.getBound().x - 2.2), 1.1 + bacRng.nextDouble()*(sim.getBound().y - 2.2), bacRng.nextDouble()*0.1*(simZ - 0.1)/2.0);
            // Test intersection

            Vector3d distance = new Vector3d(0,0,0);

            for(BSimCapsuleBacterium otherBac : bacteriaAll){
                distance.sub(otherBac.position, pos);
                if(distance.lengthSquared() < 4.5){
                    continue generator;
                }
            }

            double[] ICs = {10, 1, 10, 10, 10, 10, 10, 0};

            ActivatorBacterium bac = new ActivatorBacterium(sim,
                    new Vector3d(pos.x - bL*Math.sin(angle), pos.y - bL*Math.cos(angle), pos.z),
                    new Vector3d(bL*Math.sin(angle) + pos.x, bL*Math.cos(angle) + pos.y, pos.z),
                    h_e_field, i_e_field, ICs);

            bac.L = bL;

            bacteriaActivators.add(bac);
            bacteriaAll.add(bac);
        }

        generator:
        while(bacteriaRepressors.size() < nRepressorStart) {
            double bL = 1. + 0.1*(bacRng.nextDouble() - 0.5);
            double angle = bacRng.nextDouble()*2*Math.PI;

            Vector3d pos = new Vector3d(1.1 + bacRng.nextDouble()*(sim.getBound().x - 2.2), 1.1 + bacRng.nextDouble()*(sim.getBound().y - 2.2), simZ/2.0);
            // Test intersection

            Vector3d distance = new Vector3d(0,0,0);

            for(BSimCapsuleBacterium otherBac : bacteriaAll){
                distance.sub(otherBac.position, pos);
                if(distance.lengthSquared() < 4.5){
                    continue generator;
                }
            }

            double[] ICs = {10, 1, 10, 10, 10, 10, 10, 0};

            RepressorBacterium bac = new RepressorBacterium (sim,
                    new Vector3d(pos.x - bL*Math.sin(angle), pos.y - bL*Math.cos(angle), pos.z),
                    new Vector3d(bL*Math.sin(angle) + pos.x, bL*Math.cos(angle) + pos.y, pos.z),
                    h_e_field, i_e_field, ICs);

            bac.L = bL;

            bacteriaRepressors.add(bac);
            bacteriaAll.add(bac);
        }


        generator:
        while(bacteriaD.size() < nD1Bacterium) {
            double bL = 1. + 0.1*(bacRng.nextDouble() - 0.5);
            double angle = bacRng.nextDouble()*2*Math.PI;

            Vector3d pos = new Vector3d(1.1 + bacRng.nextDouble()*(sim.getBound().x - 2.2), 1.1 + bacRng.nextDouble()*(sim.getBound().y - 2.2), simZ/2.0);
            // Test intersection

            Vector3d distance = new Vector3d(0,0,0);

            for(BSimCapsuleBacterium otherBac : bacteriaAll){
                distance.sub(otherBac.position, pos);
                if(distance.lengthSquared() < 4.5){
                    continue generator;
                }
            }

            BSimDBacterium bac = new BSimDBacterium (sim,
                    new Vector3d(pos.x - bL*Math.sin(angle), pos.y - bL*Math.cos(angle), pos.z),
                    new Vector3d(bL*Math.sin(angle) + pos.x, bL*Math.cos(angle) + pos.y, pos.z),
                    h_e_field, i_e_field, d_e_field, q_e_field, qc_e_field);

            bac.L = bL;

            bacteriaD.add(bac);
            bacteriaAll.add(bac);
        }

        // Set up stuff for growth.
        final ArrayList<ActivatorBacterium> act_born = new ArrayList();
        final ArrayList<ActivatorBacterium> act_dead = new ArrayList();

        final ArrayList<RepressorBacterium> rep_born = new ArrayList();
        final ArrayList<RepressorBacterium> rep_dead = new ArrayList();

        final ArrayList<BSimDBacterium> d_born = new ArrayList();
        final ArrayList<BSimDBacterium> d_dead = new ArrayList();

        final Mover mover;


        mover = new RelaxationMoverGrid(bacteriaAll, sim);

        /*********************************************************
         * Set up the ticker
         */
        int LOG_INTERVAL = 100;

        BSimTicker ticker = new BSimTicker() {
            boolean toggled = true;
            @Override
            public void tick() {
                // ********************************************** Action
                long startTimeAction = System.nanoTime();


                for(BSimCapsuleBacterium b : bacteriaAll) {
                    b.action();
                }

                long endTimeAction = System.nanoTime();
                if((sim.getTimestep() % LOG_INTERVAL) == 0) {
                    System.out.println("Action update for " + bacteriaAll.size() + " bacteria took " + (endTimeAction - startTimeAction)/1e6 + " ms.");
                }

                // ********************************************** Chemical fields
                startTimeAction = System.nanoTime();

                h_e_field.update();
                i_e_field.update();
                d_e_field.update();
                q_e_field.update();
                qc_e_field.update();

                endTimeAction = System.nanoTime();
                if((sim.getTimestep() % LOG_INTERVAL) == 0) {
                    System.out.println("Chemical field update took " + (endTimeAction - startTimeAction)/1e6 + " ms.");
                }

                // ********************************************** Growth related activities if enabled.
                if(WITH_GROWTH) {

                    // ********************************************** Growth and division
                    startTimeAction = System.nanoTime();

                    for (ActivatorBacterium b : bacteriaActivators) {
                        b.grow();

                        // Divide if grown past threshold
                        if (b.L > b.L_th) {
                            act_born.add(b.divide());
                        }
                    }
                    bacteriaActivators.addAll(act_born);
                    bacteriaAll.addAll(act_born);
                    act_born.clear();

                    for (RepressorBacterium b : bacteriaRepressors) {
                        b.grow();

                        // Divide if grown past threshold
                        if (b.L > b.L_th) {
                            rep_born.add(b.divide());
                        }
                    }
                    bacteriaRepressors.addAll(rep_born);
                    bacteriaAll.addAll(rep_born);
                    rep_born.clear();

                    for (BSimDBacterium b : bacteriaD) {
                        b.grow();

                        // Divide if grown past threshold
                        if (b.L > b.L_th) {
                            d_born.add(b.divide());
                        }
                    }
                    bacteriaD.addAll(d_born);
                    bacteriaAll.addAll(d_born);
                    d_born.clear();

                    endTimeAction = System.nanoTime();
                    if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                        System.out.println("Growth and division took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
                    }

                    // ********************************************** Neighbour interactions
                    startTimeAction = System.nanoTime();

                    mover.move();

                    endTimeAction = System.nanoTime();
                    if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                        System.out.println("Wall and neighbour interactions took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
                    }

                    // ********************************************** Boundaries/removal
                    startTimeAction = System.nanoTime();
                    // Removal
                    for (ActivatorBacterium b : bacteriaActivators) {
                        // remove if past any boundary
                        if(b.position.x < 0 || b.position.x > sim.getBound().x || b.position.y < 0 || b.position.y > sim.getBound().y || b.position.z < 0 || b.position.z > sim.getBound().z){
                            act_dead.add(b);
                        }
                    }
                    bacteriaActivators.removeAll(act_dead);
                    bacteriaAll.removeAll(act_dead);
                    act_dead.clear();

                    // Removal
                    for (RepressorBacterium b : bacteriaRepressors) {
                        // remove if past the boundary
                        if(b.position.x < 0 || b.position.x > sim.getBound().x || b.position.y < 0 || b.position.y > sim.getBound().y || b.position.z < 0 || b.position.z > sim.getBound().z){
                            rep_dead.add(b);
                        }
                    }
                    bacteriaRepressors.removeAll(rep_dead);
                    bacteriaAll.removeAll(rep_dead);
                    rep_dead.clear();

                    for (BSimDBacterium b : bacteriaD) {
                        // remove if past the boundary
                        if(b.position.x < 0 || b.position.x > sim.getBound().x || b.position.y < 0 || b.position.y > sim.getBound().y || b.position.z < 0 || b.position.z > sim.getBound().z){
                            d_dead.add(b);
                        }
                    }
                    bacteriaD.removeAll(d_dead);
                    bacteriaAll.removeAll(d_dead);
                    d_dead.clear();

                    endTimeAction = System.nanoTime();
                    if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                        System.out.println("Death and removal took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
                    }
                }
            }
        };
        sim.setTicker(ticker);

        /*********************************************************
         * Set up the drawer
         */
        BSimDrawer drawer = new BSimP3DDrawer(sim, 800, 600) {
            /**
             * Draw the default cuboid boundary of the simulation as a partially transparent box
             * with a wireframe outline surrounding it.
             */
            @Override
            public void boundaries() {
                p3d.noFill();
                p3d.stroke(128, 128, 255);
                p3d.pushMatrix();
                p3d.translate((float)boundCentre.x,(float)boundCentre.y,(float)boundCentre.z);
                p3d.box((float)bound.x, (float)bound.y, (float)bound.z);
                p3d.popMatrix();
                p3d.noStroke();
            }

            /**
             * Original Chen/Bennet oscillator drawer function
             */
            @Override
            public void draw(Graphics2D g) {
                p3d.beginDraw();

                if(!cameraIsInitialised){
                    // camera(eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ)
                    p3d.camera((float)bound.x*0.5f, (float)bound.y*0.5f,
                            // Set the Z offset to the largest of X/Y dimensions for a reasonable zoom-out distance:
                            simX > simY ? (float)simX : (float)simY,
//                            10,
                            (float)bound.x*0.5f, (float)bound.y*0.5f, 0,
                            0,1,0);
                    cameraIsInitialised = true;
                }

                p3d.textFont(font);
                p3d.textMode(PConstants.SCREEN);

                p3d.sphereDetail(10);
                p3d.noStroke();
                p3d.background(255, 255,255);

                scene(p3d);
                boundaries();
                time();

                p3d.endDraw();
                g.drawImage(p3d.image, 0,0, null);
            }

            /**
             * Draw the formatted simulation time to screen.
             */
            @Override
            public void time() {
                p3d.fill(0);
                p3d.text(sim.getFormattedTimeHours(), 50, 50);
//                p3d.text(sim.getFormattedTime(), 50, 50);
            }

            /**
             * Original Chen/Bennet oscillator drawer function
             */
            @Override
            public void scene(PGraphics3D p3d) {
                p3d.ambientLight(128, 128, 128);
                p3d.directionalLight(128, 128, 128, 1, 1, -1);

//                draw(bac_act, new Color(55, 126, 184));
//                draw(bac_rep, new Color(228, 26, 28));

                for(ActivatorBacterium b : bacteriaActivators) {
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

                for(RepressorBacterium b : bacteriaRepressors) {
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
        sim.setDrawer(drawer);


        if(export) {
            String simParameters = "" + BSimUtils.timeStamp() + "__dim_" + simX + "_" + simY + "_" + simZ
                                            + "__ip_" + initialPopulation
                                            + "__pr_" + populationRatio
                                            + "__diff_" + diffusivity
                                            + "__deg_" + mu_e
                                            + "__qs_" + qsPars.get(0) + "_" + qsPars.get(1) + "_" + qsPars.get(2) + "_" + qsPars.get(3);

            if(fixedBounds){
                simParameters += "__fixedBounds";
            } else {
                simParameters += "__leakyBounds";
            }

            // !path in which files will be generated!
            String filePath = BSimUtils.generateDirectoryPath("/Users/HP/Desktop/bsim-flip-flop/tmp-results/" + simParameters + "/");


            /*********************************************************
             * Various properties of the simulation, for future reference.
             */
            BSimLogger metaLogger = new BSimLogger(sim, filePath + "simInfo.txt") {
                @Override
                public void before() {
                    super.before();
                    write("Simulation metadata.");
                    write("Chen/Bennett consortium oscillator system.");
                    write("Simulation dimensions: (" + simX + ", " + simY + ", " + simZ + ")");
                    write("Initial population: "+ initialPopulation);
                    write("Ratio " + populationRatio);
                    write("Spatial signalling diffusivity: " + diffusivity);
                    write("Spatial degradation (mu_e): " + mu_e);

                    if(fixedBounds){
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
            metaLogger.setDt(3600);			// Set export time step
            sim.addExporter(metaLogger);


            BSimLogger dataLoggerConc = new BSimLogger(sim, filePath + "Concentrations_average.csv") {
                @Override
                public void before() {
                    super.before();
                    String buffer = new String();
                    buffer = "time(seconds),h_e_field_avg,i_e_field_avg,q_e_field_avg,qc_e_field_avg";
                    write(buffer);
                }
                @Override
                public void during() {
                    String o = sim.getFormattedTime();
                    String buffer = new String();
                    //Contentration in the middle

                    int[] i_boxes = i_e_field.getBoxes();
                    //Average concentrations
                    double i_conc = 0;
                    for(int i = 0; i < i_boxes[0]; i++) {
                        for(int j = 0; j < i_boxes[1]; j++) {
                            i_conc += i_e_field.getConc(i, j, 0);
                        }

                    }

                    int[] h_boxes = h_e_field.getBoxes();
                    //Average concentrations
                    double h_conc = 0;
                    for(int i = 0; i < i_boxes[0]; i++) {
                        for(int j = 0; j < i_boxes[1]; j++) {
                            h_conc += h_e_field.getConc(i, j, 0);
                        }

                    }

                    int[] q_boxes = q_e_field.getBoxes();
                    //Average concentrations
                    double q_conc = 0;
                    for(int i = 0; i < q_boxes[0]; i++) {
                        for(int j = 0; j < q_boxes[1]; j++) {
                            q_conc += q_e_field.getConc(i, j, 0);
                        }

                    }

                    int[] qc_boxes = qc_e_field.getBoxes();
                    //Average concentrations
                    double qc_conc = 0;
                    for(int i = 0; i < qc_boxes[0]; i++) {
                        for(int j = 0; j < qc_boxes[1]; j++) {
                            qc_conc += qc_e_field.getConc(i, j, 0);
                        }

                    }

                    i_conc /= i_boxes[0] * i_boxes[1];
                    h_conc /= h_boxes[0] * h_boxes[1];
                    q_conc /= q_boxes[0] * q_boxes[1];
                    qc_conc /= qc_boxes[0] * qc_boxes[1];
                    write(o + "," + h_conc+ "," + i_conc+ "," + q_conc+ "," + qc_conc);

                    if(sim.getSimulationTime()%100 == 0){
                        System.out.println("flushing");
                        try {
                            bufferedWriter.flush();
                        } catch(Exception e){}
                    }
                }
            };
            dataLoggerConc.setDt(30);
            sim.addExporter(dataLoggerConc);


            BSimLogger posLogger = new BSimLogger(sim, filePath + "position.csv") {
                DecimalFormat formatter = new DecimalFormat("###.##", DecimalFormatSymbols.getInstance( Locale.ENGLISH ));

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

                    buffer = "";
                    for(BSimCapsuleBacterium b : bacteriaActivators) {
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
                    for(BSimCapsuleBacterium b : bacteriaRepressors) {
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
            posLogger.setDt(30);			// Set export time step
            sim.addExporter(posLogger);

            /**
             * Export a rendered image file
             */
            BSimPngExporter imageExporter = new BSimPngExporter(sim, drawer, filePath);
            System.out.println("Exporting pngs to:");
            System.out.println(filePath);
            imageExporter.setDt(30);
            sim.addExporter(imageExporter);

            sim.export();


        } else {
            sim.preview();
        }

        long simulationEndTime = System.nanoTime();

        System.out.println("Total simulation time: " + (simulationEndTime - simulationStartTime)/1e9 + " sec.");
        return new FlipFlopResult(q_e_field, qc_e_field);
    }

    public static class FlipFlopResult {
        BSimChemicalField q;
        BSimChemicalField qc;

        public FlipFlopResult(BSimChemicalField q, BSimChemicalField qc) {
            this.q = q;
            this.qc = qc;
        }
    }
}
