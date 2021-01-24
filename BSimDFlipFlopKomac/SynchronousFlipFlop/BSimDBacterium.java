package BSimDFlipFlopKomac.SynchronousFlipFlop;

import BSimDFlipFlopKomac.SynchronousFlipFlop.grn.Dsys;
import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.ode.BSimOdeSolver;

import javax.vecmath.Vector3d;

public class BSimDBacterium extends BSimCapsuleBacterium {
    protected Dsys odesys;	// Instance of ODE system
    protected double[] y, yNew;				// Local values of ODE variables
    final double cellWallDiffusivity = 2.0; 		// Cell wall diffusivity, taken from other implementations using BSimCapsuleBacterium
    BSimChemicalField _h_field;
    BSimChemicalField _i_field;
    BSimChemicalField _d_field;
    BSimChemicalField _q_field;
    BSimChemicalField _qc_field;


    public BSimDBacterium(BSim sim, Vector3d position, Vector3d position2, BSimChemicalField h_field, BSimChemicalField i_field, BSimChemicalField d_field, BSimChemicalField q_field, BSimChemicalField qc_field){
        super(sim, position, position2);
        // Setting up fields
        this._h_field  = h_field;
        this._i_field  = i_field;
        this._d_field  = d_field;
        this._q_field  = q_field;
        this._qc_field = qc_field;

        // Create the parameters and initial conditions for the ODE system
        odesys = new Dsys();
        y = odesys.getICs();
    }

    /*
     * Action each time step
     */
    @Override
    public void action() {

        // Movement
        super.action();

        // Chemical fields concentrations
        double externalChemD;	// External data chem. field
        double externalChemCLK;	// External clock chem. field
        double externalChemQ;	// External q chem. field
        double externalChemQc;	// External qc chem. field
        double deltaChemQ;		// Change in q chemical quantity
        double deltaChemQc;		// Change in qc chemical quantity

        // external chemical level at position of the bacterium
        externalChemD   = _d_field.getConc(position);
        externalChemCLK = _i_field.getConc(position);
        externalChemQ   = _q_field.getConc(position);
        externalChemQc  = _qc_field.getConc(position);

        // Get the external chemical field level for the GRN ode system later on:
        /* Qc for reverse! */
        odesys.setExternalLevel(externalChemQc, externalChemCLK);

        // re-scaled time units
        yNew = BSimOdeSolver.rungeKutta45(odesys, sim.getTime()/60, y, sim.getDt()/60);
        y = yNew;

        // Adjust the external chemical field
        deltaChemQ  = externalChemQ  - y[2];
        deltaChemQc = externalChemQc - y[3];

        // Changing external concentration
        _q_field.addQuantity(position, cellWallDiffusivity*(-deltaChemQ));
        _qc_field.addQuantity(position, cellWallDiffusivity*(-deltaChemQc));
    }

    @Override
    public BSimDBacterium divide() {

        System.out.println("D bacterium " + this.id + " is dividing...");

        Vector3d u = new Vector3d(); u.sub(this.x2, this.x1);

        // Uniform pertrubation
        double divPert = 0.1*L_max*(rng.nextDouble() - 0.5);

        double L_actual = u.length();

        double L1 = L_actual*0.5*(1 + divPert) - radius;
        double L2 = L_actual*0.5*(1 - divPert) - radius;

        Vector3d x2_new = new Vector3d();
        x2_new.scaleAdd(L1/L_actual, u, this.x1);
        x2_new.add(new Vector3d(0.05*L_initial*(rng.nextDouble() - 0.5),
                0.05*L_initial*(rng.nextDouble() - 0.5),
                0.05*L_initial*(rng.nextDouble() - 0.5)));

        Vector3d x1_child = new Vector3d();
        x1_child.scaleAdd(-(L2/L_actual), u, this.x2);
        x1_child.add(new Vector3d(0.05*L_initial*(rng.nextDouble() - 0.5),
                0.05*L_initial*(rng.nextDouble() - 0.5),
                0.05*L_initial*(rng.nextDouble() - 0.5)));

        // ICs must be slightly perturbed for both the mother and the daughter.
        // the current state of the mother cell

        // The new states that will be slightly perturbed based on division
        double[] new_state = new double[this.odesys.getNumEq()];
        double[] child_state = new double[this.odesys.getNumEq()];

        // Store the current state of the ODE system:
        System.arraycopy(this.y, 0, new_state, 0, this.odesys.getNumEq());
        System.arraycopy(this.y, 0, child_state, 0, this.odesys.getNumEq());

        // Iterate to length-1 as we don't care about time
        for(int i = 0; i < y.length; i++){
            double pert_state = 0.1*rng.nextGaussian()*this.y[i];

            new_state[i] = new_state[i] + pert_state;
            child_state[i] = child_state[i] - pert_state;
        }


        // Set the new GRN for this cell.
        Dsys new_odesys = new Dsys();
        this.odesys = new_odesys;


        BSimDBacterium child = new BSimDBacterium(sim, x1_child, new Vector3d(this.x2), this._h_field, this._i_field, this._d_field, this._q_field, this._qc_field);
        this.initialise(L1, this.x1, x2_new);

        child.L = L2;

        System.out.println("Child ID id " + child.id);

        return child;
    }


}
