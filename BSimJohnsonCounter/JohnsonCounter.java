package BSimJohnsonCounter;

import bsim.BSim;
import bsim.BSimChemicalField;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import java.util.List;
import java.util.*;

/**
 * Extended Chen oscillator in microfluidic chamber. Implemented synchronous master-slave D flip-flop
 *
 * We adjust the diffusion, the boundary conditions, the density of cells, and the ratio of activators vs repressors.
 */
public class JohnsonCounter {

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

    public static void main(String[] args) {
        JohnsonCounter bsim_ex = new JohnsonCounter();
        new JCommander(bsim_ex, args);
        bsim_ex.run();
    }

    public void run() {
        long simulationStartTime = System.nanoTime();
        double simX = simDimensions.get(0);
        double simY = simDimensions.get(1);
        double simZ = simDimensions.get(2);

        // create the simulation object
        BSim sim = new BSim();
        sim.setDt(0.25);				    // Simulation Timestep
        sim.setSimulationTime(43200);       // 36000 = 10 hours; 600 minutes.
        sim.setTimeFormat("0.00");		    // Time Format for display
        sim.setBound(simX, simY, simZ);		// Simulation Boundaries

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

        BSimChemicalField qc3 = new BSimChemicalField(sim, new int[] {(int) simX, (int)simY, 1}, external_diffusivity, external_decay);

        while (true) {
            SynchronousFlipFlopForJohnsonCounter flipFlop1 = new SynchronousFlipFlopForJohnsonCounter();
            SynchronousFlipFlopForJohnsonCounter.FlipFlopResult r1 = flipFlop1.run(sim, qc3);
            BSimChemicalField q1 = r1.q;
            BSimChemicalField qc1 = r1.qc;

            SynchronousFlipFlopForJohnsonCounter flipFlop2 = new SynchronousFlipFlopForJohnsonCounter();
            SynchronousFlipFlopForJohnsonCounter.FlipFlopResult r2 = flipFlop2.run(sim, q1);
            BSimChemicalField q2 = r2.q;
            BSimChemicalField qc2 = r2.qc;

            SynchronousFlipFlopForJohnsonCounter flipFlop3 = new SynchronousFlipFlopForJohnsonCounter();
            SynchronousFlipFlopForJohnsonCounter.FlipFlopResult r3 = flipFlop3.run(sim, q2);
            BSimChemicalField q3 = r3.q;
            qc3 = r3.qc;
        }
    }
}
