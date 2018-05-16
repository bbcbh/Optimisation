
package optimisation;

/**
 * Dx function (e.g. curvilinear) for line search 
 * @author Ben Hui
 */
public abstract class AbstractDxFunc {    
    public abstract double[] getDx(double[] alpha);            
}
