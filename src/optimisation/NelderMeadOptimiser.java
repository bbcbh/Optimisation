package optimisation;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import static optimisation.AbstractParameterOptimiser.FILEPATH_CSV;
import static optimisation.AbstractParameterOptimiser.FILEPATH_OBJ;
import static optimisation.AbstractParameterOptimiser.RES_OPTIONS_CSV;
import static optimisation.AbstractParameterOptimiser.RES_OPTIONS_OBJ;
import static optimisation.AbstractParameterOptimiser.RES_OPTIONS_PRINT;
import util.AppendableObjOutstreamFactory;

/**
 * Nelder-Mead optimisatier
 *
 * See:
 * <pre>
 * Lagariai JC, Reedi JA, Wright MH, Wright PE.
 * Convergence properties of the Nelder-Mead simplex method in low dimensions.
 * Siam J Optimiz. 1998 Dec 21;9(1):112-47.
 * </pre>
 *
 * @author Ben Hui
 * @version 20190121
 *
 * <p>
 * History</p>
 * <p>
 * 20140226: Bug Fix: Sorting under same cost </p>
 * <p>
 * 20140225: Added stopping condition.</p>  
 * <p>
 * 20190121: Added printing of best result to output.</p>
 *
 */
public class NelderMeadOptimiser extends AbstractParameterOptimiser {

    public static final int FILEPATH_SIMPLEX = AbstractParameterOptimiser.FILEPATH_LENGTH + 1;
    protected final String[] NM_filepaths = {null};

    private double[][] simplex_x; // unconstrainted
    //private double[][] simplex_p; // constrainted (or "real")
    private double[][] simplex_r;
    private double[] simplex_cost;

    private final double rho = 1;
    private final double chi = 2;
    private final double gamma = 0.5;
    private final double signma = 0.5;
    private boolean preFilled = false;

    private double D_X0 = 0.005;

    public double getD_X0() {
        return D_X0;
    }

    public void setD_X0(double D_X0) {
        this.D_X0 = D_X0;
    }

    public NelderMeadOptimiser(AbstractResidualFunc func) {
        super(func);
        baseCost = Double.POSITIVE_INFINITY;
    }

    @Override
    public void setFilename(String filename, int index) {
        if (index < AbstractParameterOptimiser.FILEPATH_LENGTH) {
            super.setFilename(filename, index);
        } else {
            NM_filepaths[index - AbstractParameterOptimiser.FILEPATH_LENGTH - 1] = filename;
        }
    }

    @Override
    protected void handleResults() {

        final boolean newMin = (baseCost > simplex_cost[0]);

        X0 = Arrays.copyOf(simplex_x[0], simplex_x[0].length);
        R0 = Arrays.copyOf(simplex_r[0], simplex_r[0].length);
        baseCost = simplex_cost[0];

        final double[] bestX = Arrays.copyOf(simplex_x[0], simplex_x[0].length);
        final double[] bestP = convertParameter(simplex_x[0], true);
        final double[] bestR = Arrays.copyOf(simplex_r[0], simplex_r[0].length);
        final double bestCost = simplex_cost[0];
        final boolean[] resOpt = Arrays.copyOf(res_options, RES_OPTIONS_LEN);

        final double[] sCost = Arrays.copyOf(simplex_cost, simplex_cost.length);
        final double[][] sX = new double[simplex_x.length][];
        final double[][] sR = new double[simplex_r.length][];
        for (int i = 0; i < sR.length; i++) {
            sX[i] = Arrays.copyOf(simplex_x[i], simplex_x[i].length);
            sR[i] = Arrays.copyOf(simplex_r[i], simplex_r[i].length);
        }

        // Asynchronous IO 
        Runnable handelResultsThread = new Runnable() {
            @Override
            public void run() {

                java.text.DateFormat df = java.text.DateFormat.getDateTimeInstance();
                if (newMin) {
                    StringBuilder str = new StringBuilder();
                    for (int i = 0; i < bestX.length; i++) {
                        if (str.length() > 0) {
                            str.append(',');
                        }
                        str.append(bestX[i]);
                    }
                    for (int i = 0; i < bestP.length; i++) {
                        str.append(',');
                        str.append(bestP[i]);
                    }
                    for (int i = 0; i < bestR.length; i++) {
                        str.append(',');
                        str.append(bestR[i]);
                    }
                    str.append(',');
                    str.append(bestCost);
                    if (resOpt[RES_OPTIONS_PRINT]) {
                        System.out.println(str.toString());
                    }
                    // For CSV file   
                    if (resOpt[RES_OPTIONS_CSV] && filepaths[FILEPATH_CSV] != null) {
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
                        File OBJ_File;
                        if (filepaths[FILEPATH_OBJ] != null) {
                            OBJ_File = new File(filepaths[FILEPATH_OBJ]);

                            try {
                                objStr = AppendableObjOutstreamFactory.generateFromFile(OBJ_File);
                                objStr.writeObject(bestX);
                                objStr.writeObject(bestP);
                                objStr.writeObject(bestR);
                                objStr.writeObject(bestCost);
                                objStr.close();
                            } catch (IOException ex) {
                                System.err.println("Error in printing object file to " + OBJ_File.getAbsolutePath());
                            }
                        }

                    }

                    System.out.println("Last result update = " + df.format(new java.util.Date()) + 
                            " P = " + Arrays.toString(bestP) + " Cost = " + bestCost );

                }
                // Simplex
                if (NM_filepaths[FILEPATH_SIMPLEX - AbstractParameterOptimiser.FILEPATH_LENGTH - 1] != null) {
                    ObjectOutputStream objStr;
                    File OBJ_File = new File(NM_filepaths[FILEPATH_SIMPLEX - AbstractParameterOptimiser.FILEPATH_LENGTH - 1]);
                    try {
                        objStr = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(OBJ_File))); // Only store the last simplex
                        objStr.writeObject(sX);
                        objStr.writeObject(sR);
                        objStr.writeObject(sCost);
                        objStr.close();
                    } catch (IOException ex) {
                        System.err.println("Error in printing object file to " + OBJ_File.getAbsolutePath());
                    }

                    System.out.println("Last simplex update = " + df.format(new java.util.Date()));
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

        int n = simplex_x.length - 1; // Second last entry is simplex_x[n-1]

        while (!optStopped) {
            handleResults();
            double pre_step_max = simplex_cost[simplex_cost.length - 1];

            double[] x_bar_sum, x_r; // centroid, refection points
            double[] f_r;
            double cost_r;

            x_bar_sum = new double[simplex_x[0].length];
            x_r = new double[simplex_x[0].length];

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < simplex_x[i].length; j++) {
                    x_bar_sum[j] += simplex_x[i][j];
                }
            }

            for (int j = 0; j < x_r.length; j++) {
                x_r[j] = ((1 + rho) * x_bar_sum[j] / n) - rho * (simplex_x[n][j]);
            }

            f_r = getResidual(x_r);
            cost_r = sumOfSqCost(f_r);

            if (simplex_cost[0] <= cost_r && cost_r < simplex_cost[n - 1]) {

                acceptSimplexPoint(cost_r, x_r, f_r);

            } else if (cost_r < simplex_cost[0]) {

                // Expansion                
                double[] x_e, f_e;
                double cost_e;

                x_e = new double[simplex_x[0].length];

                for (int j = 0; j < x_e.length; j++) {
                    x_e[j] = ((1 + rho * chi) * x_bar_sum[j] / n) - rho * chi * (simplex_x[n][j]);
                }

                f_e = getResidual(x_e);
                cost_e = sumOfSqCost(f_e);

                if (cost_e < cost_r) {
                    acceptSimplexPoint(cost_e, x_e, f_e);
                } else {
                    acceptSimplexPoint(cost_r, x_r, f_r);
                }

            } else if (cost_r >= simplex_cost[n - 1]) {
                // Contraction        
                double[] x_c, f_c;
                double cost_c;
                x_c = new double[simplex_x[0].length];

                if (simplex_cost[n - 1] <= cost_r && cost_r < simplex_cost[n]) {
                    // Outside contraction
                    for (int j = 0; j < x_c.length; j++) {
                        x_c[j] = ((1 + rho * gamma) * x_bar_sum[j] / n) - rho * gamma * (simplex_x[n][j]);
                    }

                    f_c = getResidual(x_c);
                    cost_c = sumOfSqCost(f_c);

                    if (cost_c <= cost_r) {
                        acceptSimplexPoint(cost_c, x_c, f_c);
                    } else {
                        shrinkStep();
                    }

                } else {
                    // Inside contraction
                    for (int j = 0; j < x_c.length; j++) {
                        x_c[j] = ((1 - gamma) * x_bar_sum[j] / n) + gamma * (simplex_x[n][j]);
                    }

                    f_c = getResidual(x_c);
                    cost_c = sumOfSqCost(f_c);

                    if (cost_c < simplex_cost[n]) {
                        acceptSimplexPoint(cost_c, x_c, f_c);
                    } else {
                        shrinkStep();
                    }

                }

            }

            // Check for stopping condition
            optStopped = simplex_cost[0] <= ZERO_TOL;                       
            
            if (optStopped) {
                System.out.println("Optimisaiton stopped due to min cost ("
                        + simplex_cost[0] + ") smaller than " + ZERO_TOL);               
                
                
            }else{   
                /*
                optStopped = true;
                for(int i = 1; i < simplex_cost.length; i++){
                    optStopped &= Math.abs(simplex_cost[i]- simplex_cost[0]) <= ZERO_TOL;
                }                                
                
                if(optStopped){
                    System.out.println("Optimisaiton stopped due to max - min cost ("
                        + Math.abs(simplex_cost[simplex_cost.length - 1] - simplex_cost[0])
                        + ") smaller than " + ZERO_TOL);
                }
                */
                
            }

        }

    }

    private void shrinkStep() {
        // Shrink / reduction
        for (int i = 1; i < simplex_cost.length; i++) {
            for (int ss = 0; ss < simplex_x[i].length; ss++) {
                simplex_x[i][ss] = simplex_x[0][ss] + signma * (simplex_x[i][ss] - simplex_x[0][ss]);
            }
            simplex_r[i] = getResidual(simplex_x[i]);
            simplex_cost[i] = sumOfSqCost(simplex_r[i]);
        }
        sortSimplexFull();
    }

    private void acceptSimplexPoint(double cost, double[] x, double[] f) {
        // Accept refelction
        int index_r = Arrays.binarySearch(simplex_cost, cost);

        if (index_r < 0) {
            index_r = -(index_r + 1); // insert_pt
        }

        // Shift all to left by one
        for (int r = simplex_cost.length - 2; r >= index_r; r--) {
            simplex_x[r + 1] = simplex_x[r];
            simplex_r[r + 1] = simplex_r[r];
            simplex_cost[r + 1] = simplex_cost[r];
        }
        simplex_x[index_r] = x;
        simplex_r[index_r] = f;
        simplex_cost[index_r] = cost;

    }

    private static double sumOfSqCost(double[] r) {
        double res = 0;
        for (int i = 0; i < r.length; i++) {
            res += r[i] * r[i];
        }
        return res;
    }

    public double[] setPreDefineSimplex(double[] x0, double[] r0, int index) {
        if (X0 != null) {
            if (simplex_cost == null) {
                simplex_cost = new double[X0.length + 1];
                simplex_r = new double[X0.length + 1][];
                simplex_x = new double[X0.length + 1][X0.length];
            }
            simplex_x[index] = x0;
            if (r0 == null) {
                simplex_r[index] = getResidual(x0);
            } else {
                simplex_r[index] = r0;
            }
            simplex_cost[index] = sumOfSqCost(simplex_r[index]);
            preFilled = true;
            return simplex_r[index];

        } else {
            throw new NullPointerException("X0 not defined");
        }

    }

    @Override
    public void initialise() {
        if (X0 != null) {
            System.out.println("X0 = " + Arrays.toString(X0));
            if (R0 == null) {
                R0 = getResidual(X0);
                System.out.println("R0 = " + Arrays.toString(R0));
                baseCost = sumOfSqCost(R0);
            } else if (baseCost == Double.POSITIVE_INFINITY) {
                baseCost = sumOfSqCost(R0);
            }

            // Intialise simplex
            if (!preFilled) {
                simplex_cost = new double[X0.length + 1];
                simplex_r = new double[X0.length + 1][];
                simplex_x = new double[X0.length + 1][X0.length];
                Arrays.fill(simplex_cost, -1);
            }

            simplex_cost[0] = baseCost;
            simplex_x[0] = Arrays.copyOf(X0, X0.length);
            simplex_r[0] = Arrays.copyOf(R0, R0.length);

            for (int i = 1; i < simplex_cost.length; i++) {
                if (simplex_cost[i] <= 0) { // Skip those already defined
                    simplex_x[i] = Arrays.copyOf(simplex_x[0], X0.length);
                    simplex_x[i][i - 1] += D_X0;
                    System.out.println("Setting simplex vertex #" + i + " with " + Arrays.toString(simplex_x[i]));
                }
                if (simplex_r[i] == null) {
                    System.out.print("Residual of simplex vertex #" + i + " = ");
                    simplex_r[i] = getResidual(simplex_x[i]);
                    simplex_cost[i] = sumOfSqCost(simplex_r[i]);
                    System.out.println(Arrays.toString(simplex_r[i]));
                }

            }

            sortSimplexFull();

        } else {
            throw new NullPointerException("X0 not defined");
        }

    }

    private void sortSimplexFull() {
        double[] sortedCost = Arrays.copyOf(simplex_cost, simplex_cost.length);
        Arrays.sort(sortedCost);

        double[] sort_simplex_cost = new double[simplex_cost.length];
        double[][] sort_simplex_r = new double[simplex_cost.length][];
        double[][] sort_simplex_x = new double[simplex_cost.length][X0.length];

        Arrays.fill(sort_simplex_cost, -1);

        for (int i = 0; i < simplex_cost.length; i++) {
            int bI = 0;
            while ((sortedCost[bI] < simplex_cost[i] 
                    ||  sort_simplex_cost[bI] >= 0)          // In case of equal, skip to next
                    && (bI + 1 < sort_simplex_cost.length)) { 
                bI++;
            }
            sort_simplex_cost[bI] = simplex_cost[i];
            sort_simplex_r[bI] = simplex_r[i];
            sort_simplex_x[bI] = simplex_x[i];

        }
        simplex_cost = sort_simplex_cost;
        simplex_r = sort_simplex_r;
        simplex_x = sort_simplex_x;

    }

}
