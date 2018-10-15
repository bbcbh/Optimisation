package transform;

/**
 * Parameter transformation of constraint and unconstraint parameters Based on proposal in R Fletcher, Practical methods of optimization, 2nd ed. Chichester ; New York: Wiley, 1987
 *
 * @author Ben Hui
 */
public class ParameterConstraintTransformSineCurve extends ParameterConstraintTransform {
    
    public ParameterConstraintTransformSineCurve(double[] limits) {
        super(limits);
    }
    
    public double[] getLimits() {
        return limits;
    }
    
    @Override
    public double toContrainted(double x) {
        // To be called within optimisation procedure getResidual
        return limits[LIMIT_L] + diff * Math.sin(x) * Math.sin(x);        
    }
    
    @Override
    public double toUncontrainted(double p) {
        // To be called when defining initial value
        return Math.asin(Math.sqrt((p - limits[LIMIT_L]) / diff));        
    }    
    
}
