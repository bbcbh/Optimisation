
package transform;

/**
 * Parameter transformation of constraint and unconstraint parameters using 
 * simple linear approach
 * 
 * @author Ben Hui
 */
public class ParameterConstraintTransformLinear extends ParameterConstraintTransform{

    public ParameterConstraintTransformLinear(double[] limits) {
        super(limits);
    }

    @Override
    public double toContrainted(double x) {
        // To be called within optimisation procedure getResidual
        return limits[LIMIT_L] + diff * x;       
    }

    @Override
    public double toUncontrainted(double p) {
        // To be called when defining initial value
        return (p - limits[LIMIT_L]) / diff;
    }
    
}
