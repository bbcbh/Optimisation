package optimisation;

import transform.ParameterConstraintTransform;
import java.util.Arrays;
import java.util.HashMap;

/**
 * An abstract version of Parameter optimiser
 *
 * @author Ben Hui
 * @version 20140219
 *
 */
public abstract class AbstractParameterOptimiser {

    public static final double ZERO_TOL = 1E-16;
    public static final int RES_OPTIONS_PRINT = 0;
    public static final int RES_OPTIONS_CSV = RES_OPTIONS_PRINT + 1;
    public static final int RES_OPTIONS_OBJ = RES_OPTIONS_CSV + 1;
    protected static final int RES_OPTIONS_LEN = RES_OPTIONS_OBJ + 1;
    protected final boolean[] res_options;

    public static final int FILEPATH_CSV = 0;
    public static final int FILEPATH_OBJ = FILEPATH_CSV + 1;
    protected static final int FILEPATH_LENGTH = FILEPATH_OBJ + 1;
    protected final String[] filepaths = {"optRes.csv", "optRes.obj"};

    protected double[] X0; // Dimension n - if has constraint this will be transformed parameter
    protected double[] R0; // Dimension m
    protected ParameterConstraintTransform[] constraints; // Dimension n        

    protected double baseCost;

    protected final HashMap<String, double[]> pastResults = new HashMap<>();

    protected boolean optStopped = false;

    protected AbstractResidualFunc func;

    public AbstractParameterOptimiser(AbstractResidualFunc func) {
        this.func = func;
        res_options = new boolean[RES_OPTIONS_LEN];
        Arrays.fill(res_options, true);
    }

    public boolean isOptStopped() {
        return optStopped;
    }

    public void setOptStopped(boolean optStopped) {
        this.optStopped = optStopped;
    }

    public void setResOptions(boolean b, int index) {
        if (index < res_options.length) {
            res_options[index] = b;
        }
    }

    public void setFilename(String filename, int index) {
        if (index < filepaths.length) {
            filepaths[index] = filename;
        }
    }

    protected abstract void handleResults();

    public abstract void optimise();

    public abstract void initialise();

    public ParameterConstraintTransform[] getConstraints() {
        return constraints;
    }

    public double[] getR0() {
        return R0;
    }

    public void setR0(double[] R0) {
        this.R0 = R0;
    }


    public double[] getX0() {
        return X0;
    }

    public void setX0(double[] X0) {
        this.X0 = X0;
        constraints = new ParameterConstraintTransform[X0.length];
    }

    public void setP0(double[] p0, ParameterConstraintTransform[] constraints) {
        this.X0 = new double[p0.length];
        this.constraints = new ParameterConstraintTransform[p0.length];

        for (int i = 0; i < p0.length; i++) {
            this.X0[i] = p0[i];
            this.constraints[i] = constraints[i];
            if (this.constraints[i] != null) {
                this.X0[i] = constraints[i].toUncontrainted(this.X0[i]);
            }
        }
    }

    public void clearPastResults() {
        pastResults.clear();
    }

    public void clearAll() {
    }

    public double[] getResidual(double[] x) {

        double[] p = Arrays.copyOf(x, x.length);

        // Rounding 
        if (ZERO_TOL > 0) {
            for (int i = 0; i < p.length; i++) {
                p[i] = Math.round(x[i] / ZERO_TOL) * ZERO_TOL;
            }
        }

        String xStr = Arrays.toString(p);
        if (pastResults.containsKey(xStr)) {
            return pastResults.get(xStr);
        } else {
            for (int i = 0; i < p.length; i++) {
                if (constraints[i] != null) {
                    p[i] = constraints[i].toContrainted(p[i]);
                }
            }

            double[] r = func.generateResidual(p);
            pastResults.put(xStr, r);
            return r;
        }
    }

    public double[] convertParameter(double[] p, boolean toConstrainted) {
        double[] res = Arrays.copyOf(p, p.length);
        if (constraints != null) {
            for (int i = 0; i < res.length; i++) {
                if (constraints[i] != null) {
                    if (toConstrainted) {
                        res[i] = constraints[i].toContrainted(res[i]);
                    } else {
                        res[i] = constraints[i].toUncontrainted(res[i]);
                    }
                }
            }
        }
        return res;
    }
    
    // General method for getting and setting parameter
    public Object getParameter(int paramId){
        return null;        
    }
    
    public Object setParameter(int paramId, Object newParm){
        return null;
    }
    

}
