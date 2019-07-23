package optimisation;

/**
 * An abstract function to define the residual function used for optimisation
 *
 * @author Ben Hui
 * @version 20131122
 */
public abstract class AbstractResidualFunc {

    public abstract double[] generateResidual(double[] param);

    // Preset seed for simulation
    protected long[] preset_Seed = null;

    public long[] getPreset_Seed() {
        return preset_Seed;
    }

    public void setPreset_Seed(long[] preset_Seed) {
        this.preset_Seed = preset_Seed;
    }

}
