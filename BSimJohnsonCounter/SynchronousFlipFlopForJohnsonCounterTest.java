package BSimJohnsonCounter;

import BSimDFlipFlopKomac.SynchronousFlipFlop.ActivatorBacterium;
import BSimDFlipFlopKomac.SynchronousFlipFlop.BSimDBacterium;
import BSimDFlipFlopKomac.SynchronousFlipFlop.ChenParameters;
import BSimDFlipFlopKomac.SynchronousFlipFlop.RepressorBacterium;
import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.capsule.Mover;
import bsim.capsule.RelaxationMoverGrid;
import com.beust.jcommander.Parameter;

import javax.vecmath.Vector3d;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Extended Chen oscillator in microfluidic chamber. Implemented synchronous master-slave D flip-flop
 * <p>
 * We adjust the diffusion, the boundary conditions, the density of cells, and the ratio of activators vs repressors.
 */
public class SynchronousFlipFlopForJohnsonCounterTest {

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


    /**
     * Whether to enable growth
     */
    private static final boolean WITH_GROWTH = false;

    public BSim sim;
    public BSimChemicalField d;
    int LOG_INTERVAL;

    public SynchronousFlipFlopForJohnsonCounterTest(BSim sim, BSimChemicalField d, int logInterval) {
        this.sim = sim;
        this.d = d;
        this.LOG_INTERVAL = logInterval;
        this.mover = new RelaxationMoverGrid(bacteriaAll, sim);
    }

    double simX;
    double simY;
    double simZ;

    int nActivatorStart;
    int nRepressorStart;
    int nD1Bacterium;

    long simulationStartTime;

    double def_D_H;
    double def_D_I;
    double def_phi_H;
    double def_phi_I;

    double external_diffusivity;
    double external_decay = mu_e;

    BSimChemicalField h_e_field;
    BSimChemicalField i_e_field;
    BSimChemicalField d_e_field;
    BSimChemicalField q_e_field;
    BSimChemicalField qc_e_field;

    // Separate lists of bacteria in case we want to manipulate the species individually
    final ArrayList<ActivatorBacterium> bacteriaActivators = new ArrayList();
    final ArrayList<RepressorBacterium> bacteriaRepressors = new ArrayList();
    final ArrayList<BSimDBacterium> bacteriaD = new ArrayList();

    // Track all of the bacteria in the simulation, for use of common methods etc
    final ArrayList<BSimCapsuleBacterium> bacteriaAll = new ArrayList();

    // Set up stuff for growth.
    final ArrayList<ActivatorBacterium> act_born = new ArrayList();
    final ArrayList<ActivatorBacterium> act_dead = new ArrayList();
    final ArrayList<RepressorBacterium> rep_born = new ArrayList();
    final ArrayList<RepressorBacterium> rep_dead = new ArrayList();
    final ArrayList<BSimDBacterium> d_born = new ArrayList();
    final ArrayList<BSimDBacterium> d_dead = new ArrayList();
    final Mover mover;

    public void create() {
        simX = simDimensions.get(0);
        simY = simDimensions.get(1);
        simZ = simDimensions.get(2);

        nActivatorStart = (int) Math.round(populationRatio * initialPopulation);
        nRepressorStart = (int) Math.round((1 - populationRatio) * initialPopulation);
        nD1Bacterium = initialPopulation;

        simulationStartTime = System.nanoTime();

        def_D_H = ChenParameters.p.get("D_H");
        ChenParameters.p.put("D_H", def_D_H * qsPars.get(0));
        def_D_I = ChenParameters.p.get("D_I");
        ChenParameters.p.put("D_I", def_D_I * qsPars.get(1));

        def_phi_H = ChenParameters.p.get("phi_H");
        ChenParameters.p.put("phi_H", def_phi_H * qsPars.get(2));
        def_phi_I = ChenParameters.p.get("phi_I");
        ChenParameters.p.put("phi_I", def_phi_I * qsPars.get(3));


        /*********************************************************
         * Set up the chemical fields
         */
        external_diffusivity = diffusivity / 60.0;
        // 800/60 - repressor oscillates but activator levels out rapidly
        // 80/60 - more transients, only starts levelling at the end of the 10 hours

        // Boundaries are not periodic
        sim.setSolid(true, true, true);

        // Leaky on the bottom
        if (!fixedBounds) {
            sim.setLeaky(false, false, true, false, false, false);
            sim.setLeakyRate(0, 0, 0.1 / 60.0, 0, 0, 0);
        }

        external_decay = mu_e / 60.0;

        h_e_field = new BSimChemicalField(sim, new int[]{(int) simX, (int) simY, 1}, external_diffusivity, external_decay);
        i_e_field = new BSimChemicalField(sim, new int[]{(int) simX, (int) simY, 1}, external_diffusivity, external_decay);

        d_e_field = d;
        q_e_field = new BSimChemicalField(sim, new int[]{(int) simX, (int) simY, 1}, external_diffusivity, external_decay);
        qc_e_field = new BSimChemicalField(sim, new int[]{(int) simX, (int) simY, 1}, external_diffusivity, external_decay);

        // ICs as in Chen paper (as in original DDEs)
        h_e_field.setConc(10.0);
        i_e_field.setConc(10.0);


        /*********************************************************
         * Create the bacteria
         */

        Random bacRng = new Random();

        generator:
        while (bacteriaActivators.size() < nActivatorStart) {
            double bL = 1. + 0.1 * (bacRng.nextDouble() - 0.5);
            double angle = bacRng.nextDouble() * 2 * Math.PI;

            Vector3d pos = new Vector3d(1.1 + bacRng.nextDouble() * (sim.getBound().x - 2.2), 1.1 + bacRng.nextDouble() * (sim.getBound().y - 2.2), bacRng.nextDouble() * 0.1 * (simZ - 0.1) / 2.0);
            // Test intersection

            Vector3d distance = new Vector3d(0, 0, 0);

            for (BSimCapsuleBacterium otherBac : bacteriaAll) {
                distance.sub(otherBac.position, pos);
                if (distance.lengthSquared() < 4.5) {
                    continue generator;
                }
            }

            double[] ICs = {10, 1, 10, 10, 10, 10, 10, 0};

            ActivatorBacterium bac = new ActivatorBacterium(sim,
                    new Vector3d(pos.x - bL * Math.sin(angle), pos.y - bL * Math.cos(angle), pos.z),
                    new Vector3d(bL * Math.sin(angle) + pos.x, bL * Math.cos(angle) + pos.y, pos.z),
                    h_e_field, i_e_field, ICs);

            bac.L = bL;

            bacteriaActivators.add(bac);
            bacteriaAll.add(bac);
        }

        generator:
        while (bacteriaRepressors.size() < nRepressorStart) {
            double bL = 1. + 0.1 * (bacRng.nextDouble() - 0.5);
            double angle = bacRng.nextDouble() * 2 * Math.PI;

            Vector3d pos = new Vector3d(1.1 + bacRng.nextDouble() * (sim.getBound().x - 2.2), 1.1 + bacRng.nextDouble() * (sim.getBound().y - 2.2), simZ / 2.0);
            // Test intersection

            Vector3d distance = new Vector3d(0, 0, 0);

            for (BSimCapsuleBacterium otherBac : bacteriaAll) {
                distance.sub(otherBac.position, pos);
                if (distance.lengthSquared() < 4.5) {
                    continue generator;
                }
            }

            double[] ICs = {10, 1, 10, 10, 10, 10, 10, 0};

            RepressorBacterium bac = new RepressorBacterium(sim,
                    new Vector3d(pos.x - bL * Math.sin(angle), pos.y - bL * Math.cos(angle), pos.z),
                    new Vector3d(bL * Math.sin(angle) + pos.x, bL * Math.cos(angle) + pos.y, pos.z),
                    h_e_field, i_e_field, ICs);

            bac.L = bL;

            bacteriaRepressors.add(bac);
            bacteriaAll.add(bac);
        }


        generator:
        while (bacteriaD.size() < nD1Bacterium) {
            double bL = 1. + 0.1 * (bacRng.nextDouble() - 0.5);
            double angle = bacRng.nextDouble() * 2 * Math.PI;

            Vector3d pos = new Vector3d(1.1 + bacRng.nextDouble() * (sim.getBound().x - 2.2), 1.1 + bacRng.nextDouble() * (sim.getBound().y - 2.2), simZ / 2.0);
            // Test intersection

            Vector3d distance = new Vector3d(0, 0, 0);

            for (BSimCapsuleBacterium otherBac : bacteriaAll) {
                distance.sub(otherBac.position, pos);
                if (distance.lengthSquared() < 4.5) {
                    continue generator;
                }
            }

            BSimDBacterium bac = new BSimDBacterium(sim,
                    new Vector3d(pos.x - bL * Math.sin(angle), pos.y - bL * Math.cos(angle), pos.z),
                    new Vector3d(bL * Math.sin(angle) + pos.x, bL * Math.cos(angle) + pos.y, pos.z),
                    h_e_field, i_e_field, d_e_field, q_e_field, qc_e_field);

            bac.L = bL;

            bacteriaD.add(bac);
            bacteriaAll.add(bac);
        }
    }

    public FlipFlopResult animateFlipFlop(BSim sim, BSimChemicalField d) {
        d_e_field = d;

        // ********************************************** Action
        long startTimeAction = System.nanoTime();

        for (BSimCapsuleBacterium b : bacteriaAll) {
            b.action();
        }

        long endTimeAction = System.nanoTime();
        if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
            System.out.println("Action update for " + bacteriaAll.size() + " bacteria took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
        }

        // ********************************************** Chemical fields
        startTimeAction = System.nanoTime();

        h_e_field.update();
        i_e_field.update();
        d_e_field.update();
        q_e_field.update();
        qc_e_field.update();

        endTimeAction = System.nanoTime();
        if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
            System.out.println("Chemical field update took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
        }

        // ********************************************** Growth related activities if enabled.
        if (WITH_GROWTH) {

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
                if (b.position.x < 0 || b.position.x > sim.getBound().x || b.position.y < 0 || b.position.y > sim.getBound().y || b.position.z < 0 || b.position.z > sim.getBound().z) {
                    act_dead.add(b);
                }
            }
            bacteriaActivators.removeAll(act_dead);
            bacteriaAll.removeAll(act_dead);
            act_dead.clear();

            // Removal
            for (RepressorBacterium b : bacteriaRepressors) {
                // remove if past the boundary
                if (b.position.x < 0 || b.position.x > sim.getBound().x || b.position.y < 0 || b.position.y > sim.getBound().y || b.position.z < 0 || b.position.z > sim.getBound().z) {
                    rep_dead.add(b);
                }
            }
            bacteriaRepressors.removeAll(rep_dead);
            bacteriaAll.removeAll(rep_dead);
            rep_dead.clear();

            for (BSimDBacterium b : bacteriaD) {
                // remove if past the boundary
                if (b.position.x < 0 || b.position.x > sim.getBound().x || b.position.y < 0 || b.position.y > sim.getBound().y || b.position.z < 0 || b.position.z > sim.getBound().z) {
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
        return new FlipFlopResult(i_e_field, h_e_field, q_e_field, qc_e_field, bacteriaActivators, bacteriaRepressors);
    }

    public static class FlipFlopResult {
        BSimChemicalField i;
        BSimChemicalField h;
        BSimChemicalField q;
        BSimChemicalField qc;
        ArrayList<ActivatorBacterium> bacteriaActivators;
        ArrayList<RepressorBacterium> bacteriaRepressors;

        public FlipFlopResult(BSimChemicalField i, BSimChemicalField h, BSimChemicalField q, BSimChemicalField qc,
                              ArrayList<ActivatorBacterium> bacteriaActivators, ArrayList<RepressorBacterium> bacteriaRepressors) {
            this.i = h;
            this.h = h;
            this.q = q;
            this.qc = qc;
            this.bacteriaActivators = bacteriaActivators;
            this.bacteriaRepressors = bacteriaRepressors;
        }
    }
}
