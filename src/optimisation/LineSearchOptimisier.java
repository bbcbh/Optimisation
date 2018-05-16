package optimisation;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import util.AppendableObjOutstreamFactory;

/**
 * Running of parameter optimisation algorithm
 *
 * @author Ben Hui
 * @version 20140219
 *
 * History: 20140214: Add use line search parameter 20140219: Restructure to use abstract class
 */
public class LineSearchOptimisier extends AbstractParameterOptimiser {

    private final double JACOBIAN_DIFF = 1E-2; // As a ratio of x0          
    private double[][] jacobian; // J, dimiension m x n 
    //private final Algebra MatOp = new Algebra(ZERO_TOL);
    private RealMatrix a;
    private RealMatrix H;
    private boolean usingLineSearch = false;

    public boolean isUsingLineSearch() {
        return usingLineSearch;
    }

    public void setUsingLineSearch(boolean usingLineSearch) {
        this.usingLineSearch = usingLineSearch;
    }

    public LineSearchOptimisier(AbstractResidualFunc func) {
        super(func);
    }

    @Override
    protected void handleResults() {
        // For now just print the CSV file
        final double[] x_wri = Arrays.copyOf(X0, X0.length);
        final double[] r_wri = Arrays.copyOf(R0, R0.length);
        final boolean[] resOpt = Arrays.copyOf(res_options, RES_OPTIONS_LEN);
        final double cost_wri = baseCost;

        // Asynchronous IO 
        Runnable handelResultsThread = new Runnable() {
            @Override
            public void run() {
                double[] p_wri = null;
                if (constraints != null) {
                    p_wri = Arrays.copyOf(X0, X0.length);
                    for (int i = 0; i < p_wri.length; i++) {
                        if (constraints[i] != null) {
                            p_wri[i] = constraints[i].toContrainted(x_wri[i]);
                        } else {
                            p_wri[i] = x_wri[i];
                        }
                    }
                }

                StringBuilder str = new StringBuilder();
                for (int i = 0; i < x_wri.length; i++) {
                    if (str.length() > 0) {
                        str.append(',');
                    }
                    str.append(x_wri[i]);
                }
                if (p_wri != null) {
                    for (int i = 0; i < x_wri.length; i++) {
                        str.append(',');
                        str.append(p_wri[i]);
                    }
                }

                for (int i = 0; i < r_wri.length; i++) {
                    str.append(',');
                    str.append(r_wri[i]);
                }

                str.append(',');
                str.append(cost_wri);

                if (resOpt[RES_OPTIONS_PRINT]) {
                    System.out.println(str.toString());
                }

                // For CSV file        
                if (resOpt[RES_OPTIONS_CSV]) {
                    PrintWriter wri;
                    File CSV_File = new File(filepaths[FILEPATH_CSV]);

                    try {
                        wri = new PrintWriter(new FileWriter(CSV_File, true));
                        wri.println(str.toString());
                        wri.close();
                    } catch (IOException ex) {
                        System.err.println("Error in printing CSV to " + CSV_File.getAbsolutePath());
                    }
                }

                // For ResultObj
                if (resOpt[RES_OPTIONS_OBJ]) {
                    ObjectOutputStream objStr;
                    File OBJ_File = new File(filepaths[FILEPATH_OBJ]);

                    try {
                        objStr = AppendableObjOutstreamFactory.generateFromFile(OBJ_File);
                        objStr.writeObject(x_wri);
                        objStr.writeObject(p_wri);
                        objStr.writeObject(r_wri);
                        objStr.writeObject(cost_wri);
                        objStr.close();
                    } catch (IOException ex) {
                        System.err.println("Error in printing object file to " + OBJ_File.getAbsolutePath());
                    }
                }
            }
        };

        ExecutorService executor = Executors.newFixedThreadPool(1);
        executor.submit(handelResultsThread);
        executor.shutdown();
        //handelResultsThread.run();
    }

    @Override
    public void optimise() {
        initialise();
        baseCost = 0;
        for (int i = 0; i < R0.length; i++) {
            baseCost += R0[i] * R0[i];
        }

        handleResults();

        while (!optStopped) {
            // Generate Jacobian if need
            initialise();
            // Generate gradient and Hessian      
            generateMatrices();
            // Calculate curvlinear path 
            AbstractDxFunc curvlinFunc = generatedDxFunc();
            // Determine best Dx
            double[] dx;

            if (usingLineSearch) {
                dx = dxLineSearch(curvlinFunc);
            } else {
                dx = curvlinFunc.getDx(new double[]{1e-9});
            }

            double[] x_org = Arrays.copyOf(X0, X0.length);

            clearAll();
            X0 = new double[x_org.length];

            boolean allSmall = true;

            for (int i = 0; i < X0.length; i++) {
                X0[i] = x_org[i] + dx[i];
                allSmall &= Math.abs(dx[i]) <= ZERO_TOL;
            }
            R0 = getResidual(X0);

            baseCost = 0;
            for (int i = 0; i < R0.length; i++) {
                baseCost += R0[i] * R0[i];
            }
            handleResults();

            optStopped |= allSmall;
            optStopped |= baseCost <= ZERO_TOL;

        }

    }

    // Return best DX based on line search           
    // <editor-fold defaultstate="collapsed" desc="public double[] dxLineSearch(AbstractDxFunc dxFunc)">
    public double[] dxLineSearch(AbstractDxFunc dxFunc) {

        // See Brent's Method for line seach, see  Numberical recipes in C 
        final double GOLDEN_RATIO = (3 - Math.sqrt(5)) / 2.0;
        final double DEFAULT_MAX_ALPHA = 1e9;
        final double DEFAULT_ALPHA_RANGE_MIN = 1e-8;
        final int DEFAULT_MAX_ITERNATION = 25;

        double[] a_interval;
        double interval_size;

        /*
         * Used in Brent's Method with
         * 0 - > Third least function value(or Second least in previous cycle) .
         * 1 - > Second least function value
         * 2 - > Least function value so far
         * 3 - > Recent function
         */
        final int I_3_BEST = 0;
        final int I_2_BEST = 1;
        final int I_BEST = 2;
        final int I_RECENT = 3;

        double[] i_store = new double[4];
        double[] i_cost = new double[4];
        double[][] i_dx = new double[4][];

        double tolerance_interval;

        double dist_move_pre = 0;
        double dist_move = 0;
        double dist_move_temp;
        double r, q, p; // for trial parabolic fit

        int k; // # Iternation

        boolean[] bCond; // Break condition

        // Initalisation
        a_interval = new double[3];
        a_interval[0] = 0; // L(0) = zero vector - > inital value
        a_interval[1] = DEFAULT_MAX_ALPHA / 2;
        a_interval[2] = DEFAULT_MAX_ALPHA;
        k = 0;

        // L(0) = zero vector -> dx = 0 - > Initial value        
        for (int i = 0; i < i_dx.length; i++) {
            i_store[i] = 0;
            if (i != I_RECENT) {
                i_cost[i] = baseCost;
            } else {
                i_cost[i] = Double.POSITIVE_INFINITY;
            }
            i_dx[i] = new double[X0.length];
        }

        while ((k < DEFAULT_MAX_ITERNATION) && i_store[I_BEST] < DEFAULT_MAX_ALPHA) {
            a_interval[1] = 0.5 * (a_interval[2] + a_interval[0]);
            tolerance_interval = 2 * (DEFAULT_ALPHA_RANGE_MIN * Math.abs(i_store[I_BEST]) + ZERO_TOL);
            interval_size = a_interval[2] - a_interval[0];

            // Breaking conditions
            bCond = new boolean[2];
            bCond[0] = (i_store[I_BEST] > DEFAULT_MAX_ALPHA);
            bCond[1] = Math.abs(i_store[I_BEST] - a_interval[1])
                    <= (tolerance_interval - 0.5 * interval_size);

            for (int b = 0; b < bCond.length; b++) {
                if (bCond[b]) {
                    break;
                }
            }

            if (Math.abs(dist_move_pre) > tolerance_interval / 2) {
                // Parabolic fit

                // alpha_next = alpha  - 0.5 * ( p / (q-r))
                r = (i_store[I_BEST] - i_store[I_2_BEST]) * (i_cost[I_BEST] - i_cost[I_3_BEST]);
                q = (i_store[I_BEST] - i_store[I_3_BEST]) * (i_cost[I_BEST] - i_cost[I_2_BEST]);
                p = (i_store[I_BEST] - i_store[I_3_BEST]) * q - (i_store[I_BEST] - i_store[I_2_BEST]) * r;

                q = 2 * (q - r);
                if (q > 0) {
                    p = -p;
                }
                q = Math.abs(q);

                dist_move_temp = dist_move_pre;
                dist_move_pre = dist_move;

                // Check accepatabity of the parabolic fit
                if (Math.abs(p) >= Math.abs(0.5 * q * dist_move_temp) // Movement less than half of dist_move_pre
                        || p <= q * (a_interval[0] - i_store[2]) // Fall outside parabola on left
                        || p >= q * (a_interval[2] - i_store[2]) // Fall outside parabola on right
                        ) {
                    // Take Golden section of the larger segment
                    if (i_store[I_BEST] >= a_interval[1]) {
                        dist_move_pre = a_interval[0] - i_store[I_BEST];
                    } else {
                        dist_move_pre = a_interval[2] - i_store[I_BEST];
                    }
                    dist_move = GOLDEN_RATIO * dist_move_pre;
                } else {
                    dist_move = p / q;
                    i_store[I_RECENT] = i_store[I_BEST] + dist_move;

                    // Check if movement is valid
                    if (i_store[I_RECENT] - a_interval[0] < tolerance_interval
                            || a_interval[2] - i_store[I_RECENT] < tolerance_interval) {
                        dist_move = ((a_interval[1] - i_store[I_BEST]) >= 0 ? tolerance_interval : -tolerance_interval) / 2;
                    }
                }
            } else {
                // Take Golden section of the larger segment                
                dist_move_pre = (i_store[I_BEST] >= a_interval[1])
                        ? a_interval[0] - i_store[I_BEST] : a_interval[2] - i_store[I_BEST];
                dist_move = GOLDEN_RATIO * dist_move_pre;

            }
            // Check if movemment is valid.
            // i.e. If move outside the tolerance interval
            if (Math.abs(dist_move) >= tolerance_interval / 2) {
                i_store[I_RECENT] = i_store[I_BEST] + dist_move;
            } else {
                // Use half tolerance interval as range.                
                i_store[I_RECENT] = i_store[2]
                        + ((dist_move > 0) ? tolerance_interval / 2 : -tolerance_interval / 2);
            }

            // Cost evalation for the newest point
            i_dx[I_RECENT] = dxFunc.getDx(new double[]{i_store[I_RECENT]});

            double[] newX = Arrays.copyOf(X0, X0.length);
            for (int i = 0; i < newX.length; i++) {
                newX[i] += i_dx[I_RECENT][i];
            }
            double[] newR = getResidual(newX);
            i_cost[I_RECENT] = 0;
            for (int i = 0; i < newR.length; i++) {
                i_cost[I_RECENT] += newR[i] * newR[i];
            }

            if (i_cost[I_RECENT] <= i_cost[I_BEST]) {
                // A newest point is minimum value
                if (i_store[I_RECENT] >= i_store[I_BEST]) {
                    a_interval[0] = i_store[I_BEST];
                } else {
                    a_interval[2] = i_store[I_BEST];
                }
                // Shift bracket
                for (int i = 0; i < i_store.length - 1; i++) {
                    i_store[i] = i_store[i + 1];
                    i_cost[i] = i_cost[i + 1];
                    i_dx[i] = i_dx[i + 1];
                }
            } else {
                // New point non-minimum - shift bracket range to new point
                if (i_store[I_RECENT] < i_store[I_BEST]) {
                    a_interval[0] = i_store[I_RECENT];
                } else {
                    a_interval[2] = i_store[I_RECENT];
                }
                // Shift bracket
                if (i_cost[I_RECENT] <= i_cost[I_2_BEST]
                        || i_store[I_2_BEST] == i_store[I_BEST]) {
                    // Newest point is second-least
                    i_store[I_3_BEST] = i_store[I_2_BEST];
                    i_store[I_2_BEST] = i_store[I_RECENT];
                    i_cost[I_3_BEST] = i_cost[I_2_BEST];
                    i_cost[I_2_BEST] = i_cost[I_RECENT];
                    i_dx[I_3_BEST] = i_dx[I_2_BEST];
                    i_dx[I_2_BEST] = i_dx[I_RECENT];

                } else if (i_cost[I_RECENT] <= i_cost[I_3_BEST]
                        || i_store[I_3_BEST] == i_store[I_BEST]
                        || i_store[I_3_BEST] == i_store[I_2_BEST]) {

                    // Newest point point is third least 
                    i_store[I_3_BEST] = i_store[I_RECENT];
                    i_cost[I_3_BEST] = i_cost[I_RECENT];
                    i_dx[I_3_BEST] = i_dx[I_RECENT];
                } else {
                    // Newest point is greatest                    
                }
            }
            k++;
        }
        return i_dx[I_BEST];
    }

    // </editor-fold>
    public AbstractDxFunc generatedDxFunc() {
        if (a == null || H == null) {
            generateMatrices();
        }
        
        EigenDecomposition eig = new EigenDecomposition(H);
        
        final RealMatrix V = eig.getV();
        final RealMatrix invV = eig.getVT(); // V_T = invT for symmetric matrix H
        final RealMatrix D = eig.getD();

        // Using curvilinear path
        AbstractDxFunc dxFunc = new AbstractDxFunc() {
            @Override
            public double[] getDx(double[] alpha) {
                RealMatrix dx;
                RealMatrix M =  D.copy();  //new DenseDoubleMatrix2D(D.toArray());
                for (int d = 0; d < D.getRowDimension(); d++) {
                    if (M.getEntry(d, d) == 0) {
                        M.setEntry(d, d, -alpha[0]);
                    } else {
                        M.setEntry(d, d, (Math.exp(-M.getEntry(d, d) * alpha[0]) - 1) / M.getEntry(d, d));
                    }
                }
                dx =  V.multiply(M).multiply(invV).multiply(a);    //MatOp.mult(MatOp.mult(MatOp.mult(V, M), invV), a);                                
                
                return dx.getColumn(0);
            }
        };
        return dxFunc;

    }

    private void generateMatrices() {
        if (R0 == null || jacobian == null) {
            initialise();
        }

        RealMatrix J = new Array2DRowRealMatrix(jacobian);
        RealMatrix JT = J.transpose();
        RealMatrix r0 = new Array2DRowRealMatrix(R0);

        // a = 2 * JT *r0  , H = 2 * JT * J                  
        a =  JT.scalarMultiply(2).multiply(r0);   
        H =  JT.scalarMultiply(2).multiply(J);    

    }

    @Override
    public void initialise() {
        if (X0 != null) {

            if (R0 == null) {
                R0 = getResidual(X0);
            }
            if (jacobian == null) {
                int n = X0.length;
                int m = R0.length;

                jacobian = new double[m][n];

                for (int j = 0; j < n; j++) {
                    double[] x1 = Arrays.copyOf(X0, X0.length);

                    // Jacobian different
                    double dxj;
                    //dxj = x1[j] * JACOBIAN_DIFF;
                    dxj = JACOBIAN_DIFF;
                    x1[j] = x1[j] + dxj;

                    // f(x1)
                    double[] r1 = getResidual(x1);
                    for (int i = 0; i < m; i++) {
                        // At this stage, no weight
                        jacobian[i][j] = (r1[i] - R0[i]) / dxj;
                    }
                }
            }

            generateMatrices();
        } else {
            throw new NullPointerException("X0 not defined");
        }
    }

    @Override
    public void clearAll() {
        R0 = null;
        jacobian = null;
        a = null;
        H = null;
    }

    public double[][] getJacobian() {
        return jacobian;
    }

    public void setJacobian(double[][] jacobian) {
        this.jacobian = jacobian;
    }
}
