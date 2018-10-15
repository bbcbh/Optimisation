package transform;

/**
 * An abstract class define the object for the transformation of constraint and unconstraint parameters
 *
 * @author Ben Hui
 */
public abstract class ParameterConstraintTransform {

    protected final double[] limits;
    protected final double diff;

    public static final int LIMIT_L = 0;
    public static final int LIMIT_U = 1;

    public ParameterConstraintTransform(double[] limits) {
        this.limits = limits;
        diff = limits[LIMIT_U] - limits[LIMIT_L];
    }

    public abstract double toContrainted(double x); // To be called within optimisation procedure getResidual

    public abstract double toUncontrainted(double p); // To be called when defining initial value

}
