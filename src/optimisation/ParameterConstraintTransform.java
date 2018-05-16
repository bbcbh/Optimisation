
package optimisation;

/**
 * Parameter transformation of constraint and unconstraint parameters
 * Based on proposal in
 * R Fletcher, Practical methods of optimization, 2nd ed. Chichester ; New York: Wiley, 1987
 * 
 * @author Ben Hui
 */
public class ParameterConstraintTransform {
    
    private final double[] limits;
    private final double diff;
    
    public static final int LIMIT_L = 0;
    public static final int LIMIT_U = 1;
    
   

    public ParameterConstraintTransform(double[] limits) {
        this.limits = limits;
        diff = limits[LIMIT_U] - limits[LIMIT_L];
    }

    public double[] getLimits() {
        return limits;
    }
    
    public double toContrainted(double x){
        // To be called within optimisation procedure getResidual
        return limits[LIMIT_L] + diff * Math.sin(x) * Math.sin(x);                
    }
    
    public double toUncontrainted(double p){
        // To be called when defining initial value
        return Math.asin(Math.sqrt((p - limits[LIMIT_L]) / diff));        
    }        
    
}
