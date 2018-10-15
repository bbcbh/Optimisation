package test;

import optimisation.AbstractDxFunc;
import optimisation.AbstractResidualFunc;
import transform.ParameterConstraintTransform;
import optimisation.LineSearchOptimisier;

/**
 * Optimisation test using Rosenbrock banana function, with minimum at (1,1)
 * See <a href="http://www.mathworks.com.au/help/matlab/ref/fminsearch.html> Matlab </a> for an example
 * @author Ben Hui
 */
public class Test_Opt_Rosenbrock {

    public static void main(String[] arg) {

        // Rosenbrock banana function 
        AbstractResidualFunc rosenbrock = new AbstractResidualFunc() {
            @Override
            public double[] generateResidual(double[] param) {
                double[] r = new double[1];  
                r[0] = 100 * Math.pow(param[1] - param[0] * param[0],2) + (1-param[0])*(1-param[0]);           
                return r;
            }
        };
        
        double[] x0 = {-1.2,1};
        
        LineSearchOptimisier opt = new LineSearchOptimisier(rosenbrock);
        opt.setResOptions(false, LineSearchOptimisier.RES_OPTIONS_OBJ);
        opt.setResOptions(false, LineSearchOptimisier.RES_OPTIONS_CSV);                
        
        opt.setP0(x0, new ParameterConstraintTransform[]{
            new transform.ParameterConstraintTransformSineCurve(new double[]{-1.5,1.5}),
            new transform.ParameterConstraintTransformSineCurve(new double[]{-1.5,1.5}),            
        });                               
        
        opt.initialise();                                    
        opt.optimise();
    }
}
