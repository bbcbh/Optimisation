
package optimisation;

/**
 * An abstract function to define the residual function used for optimisation 
 * @author Ben Hui
 * @version 20131122
 */
public abstract class AbstractResidualFunc {          
    public abstract double[] generateResidual(double[] param);    
}
