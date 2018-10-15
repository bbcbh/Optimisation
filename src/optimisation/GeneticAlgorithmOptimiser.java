package optimisation;

import transform.ParameterConstraintTransform;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import util.AppendableObjOutstreamFactory;

/**
 *
 * @author Ben Hui
 */
public class GeneticAlgorithmOptimiser extends AbstractParameterOptimiser {

    public static final int FILEPATH_GA_POP = AbstractParameterOptimiser.FILEPATH_LENGTH + 1;
    public static final int PARAM_GA_OPT_POP_FILE = 0;
    public static final int PARAM_GA_OPT_POP_SIZE = PARAM_GA_OPT_POP_FILE + 1;
    public static final int PARAM_GA_OPT_POP_CROSSOVER = PARAM_GA_OPT_POP_SIZE + 1;
    public static final int PARAM_GA_OPT_USE_PARALLEL = PARAM_GA_OPT_POP_CROSSOVER + 1;
    public static final int PARAM_GA_OPT_RNG = PARAM_GA_OPT_USE_PARALLEL + 1;
    public static final int PARAM_GA_OPT_NUM_PARAM = PARAM_GA_OPT_RNG + 1;
    public static final int PARAM_GA_OPT_NUM_R0 = PARAM_GA_OPT_NUM_PARAM + 1;
    public static final int PARAM_GA_OPT_LAST_EXPORT_SIM_COUNTER = PARAM_GA_OPT_NUM_R0 + 1;
    public static final int PARAM_GA_OPT_LAST_COMPLETED_SIM_COUNTER = PARAM_GA_OPT_LAST_EXPORT_SIM_COUNTER + 1;

    protected Object[] paramList = {
        //0: PARAM_GA_OPT_POP_FILE
        new File("GA_POP.obj"),
        //1: PARAM_GA_OPT_POP_SIZE
        new Integer(1000),
        //2: PARAM_GA_OPT_POP_CROSSOVER
        // For example, 0.4 means keep 40% of simulation run, and 40% will be generated, and 20% from random
        new Float(0.40f),
        //3: PARAM_GA_OPT_USE_PARALLEL
        new Integer(Runtime.getRuntime().availableProcessors()),
        //4: PARAM_GA_OPT_RNG
        new MersenneTwister(2251912970037127827l),
        //5: PARAM_GA_OPT_NUM_PARAM
        new Integer(2),
        //6: PARAM_GA_OPT_NUM_R0
        new Integer(1),
        //7: PARAM_GA_OPT_LAST_EXPORT_SIM_COUNTER
        new Integer(0),
        //8: PARAM_GA_OPT_LAST_COMPLETED_SIM_COUNTER
        new Integer(0),};

    protected double[][] GA_POP = null; // Format: {residue, x0,.., r0}
    protected double[] bestRow = null;

    private final Comparator<double[]> GA_POP_COMPARATOR = new Comparator<double[]>() {
        @Override
        public int compare(double[] t, double[] t1) {
            // Only compare the first entry
            return Double.compare(t[0], t1[0]);
        }
    };

    public GeneticAlgorithmOptimiser(AbstractResidualFunc func) {
        super(func);
    }

    // For backward compatability
    @Override
    public void setFilename(String filename, int index) {
        if (index < AbstractParameterOptimiser.FILEPATH_LENGTH) {
            super.setFilename(filename, index);
        } else if (index == FILEPATH_GA_POP) {
            paramList[PARAM_GA_OPT_POP_FILE] = new File(filename);
        }
    }

    @Override
    protected void handleResults() {
        sortGA_POP();
        double bestSoFar = GA_POP[0][0];

        int NUM_PARAM = (Integer) paramList[PARAM_GA_OPT_NUM_PARAM];

        if (bestRow == null || bestSoFar < bestRow[0]) {
            bestRow = Arrays.copyOf(GA_POP[0], GA_POP[0].length);
            final double[] bestX = Arrays.copyOfRange(bestRow, 1, 1 + NUM_PARAM);
            final double[] bestP = convertParameter(bestX, true);
            final double[] bestR = Arrays.copyOfRange(bestRow, 1 + NUM_PARAM, bestRow.length);
            final double bestCost = GA_POP[0][0];
            final boolean[] resOpt = Arrays.copyOf(res_options, RES_OPTIONS_LEN);

            // Asynchronous IO 
            Runnable handelResultsThread = new Runnable() {
                @Override
                public void run() {
                    java.text.DateFormat df = java.text.DateFormat.getDateTimeInstance();

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
                            ex.printStackTrace(System.err);
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
                                ex.printStackTrace(System.err);
                            }
                        }

                    }
                    System.out.println("Last result update = " + df.format(new java.util.Date()));

                }
            };

            ExecutorService executor = Executors.newFixedThreadPool(1);
            executor.submit(handelResultsThread);
            executor.shutdown();
        }
        perform_GA();

    }

    protected void perform_GA() {

        java.text.DateFormat df = java.text.DateFormat.getDateTimeInstance();
        System.out.println("Perform GA at = " + df.format(new java.util.Date()));

        RandomGenerator rng = (RandomGenerator) paramList[PARAM_GA_OPT_RNG];

        // Selection
        int numToKeep = Math.max(2, Math.round(GA_POP.length * ((Float) paramList[PARAM_GA_OPT_POP_CROSSOVER])));

        int numNewGen = 0;
        for (int row = numToKeep; row < GA_POP.length; row++) {
            Arrays.fill(GA_POP[row], Double.NaN);
            for (int i = 1; i < (Integer) paramList[PARAM_GA_OPT_NUM_PARAM] + 1; i++) {
                if (numNewGen < numToKeep) {
                    // Cross Over                                          
                    double p1 = GA_POP[rng.nextInt(numToKeep)][i];
                    double p2 = GA_POP[rng.nextInt(numToKeep)][i];
                    GA_POP[row][i] = (p1 + p2) / 2;
                } else {
                    // Heuristics / mutation
                    GA_POP[row][i] = rng.nextDouble() * Math.PI / 2;
                }
            }

            if (numNewGen < numToKeep) {
                numNewGen++;
            }

        }
    }

    @Override
    public void optimise() {

        ExecutorService executor = null;

        int numParam = (Integer) paramList[PARAM_GA_OPT_NUM_PARAM];

        while (!isOptStopped()) {

            double bestSoFar = GA_POP[0][0];
            double worstSoFar = Double.NaN;

            int lastRowPt = GA_POP.length - 1;

            while (Double.isNaN(worstSoFar) && lastRowPt > 0) {
                worstSoFar = GA_POP[lastRowPt][0];
                lastRowPt--;
            }

            if (bestSoFar <= ZERO_TOL
                    || ((worstSoFar - bestSoFar) <= ZERO_TOL)) {
                setOptStopped(true);
                //exportGAPopTimer.cancel();
            }

            if (!isOptStopped()) {

                int numThreadAtTime = ((Integer) paramList[PARAM_GA_OPT_USE_PARALLEL]);
                int numThreadInPool = 0;
                Future<double[]>[] r0Collection = new Future[GA_POP.length];

                if (numThreadAtTime > 1) {
                    executor = Executors.newFixedThreadPool(((Integer) paramList[PARAM_GA_OPT_USE_PARALLEL]));
                }

                for (int r = 0; r < GA_POP.length; r++) {
                    double[] GA_POP_ROW = GA_POP[r];
                    if (Double.isNaN(GA_POP_ROW[0])) {
                        // Need to generate new result
                        if (numThreadAtTime <= 1) {
                            double[] r0 = getResidual(Arrays.copyOfRange(GA_POP_ROW, 1, 1 + numParam));
                            GA_POP_ROW[0] = sumOfSqCost(r0);
                            System.arraycopy(r0, 0, GA_POP_ROW, 1 + numParam, r0.length);
                            paramList[PARAM_GA_OPT_LAST_COMPLETED_SIM_COUNTER]
                                    = ((Integer) paramList[PARAM_GA_OPT_LAST_COMPLETED_SIM_COUNTER]) + 1;
                            exportGAPop();

                        } else {
                            if (executor == null) {
                                executor = Executors.newFixedThreadPool(numThreadAtTime);
                                numThreadInPool = 0;
                            }

                            Callable_CalculateResidue runnable = new Callable_CalculateResidue(GA_POP_ROW, numParam);
                            r0Collection[r] = executor.submit(runnable);

                            numThreadInPool++;

                            if (numThreadInPool == numThreadAtTime) {
                                runExecutor(executor, r0Collection, numParam);
                                executor = null;
                            }
                        }
                    }
                }
                // Run the ending thread pool if needed                
                if (executor != null) {
                    runExecutor(executor, r0Collection, numParam);
                    executor = null;
                }
                handleResults();
            }

        }

    }

    private void runExecutor(ExecutorService executor,
            Future<double[]>[] r0Collection,
            int numParam) {
        try {
            executor.shutdown();
            if (!executor.awaitTermination(2, TimeUnit.DAYS)) {
                System.out.println("Inf Thread time-out!");
            }
            for (int updateGAPopPt = 0; updateGAPopPt < r0Collection.length; updateGAPopPt++) {
                if (r0Collection[updateGAPopPt] != null
                        && Double.isNaN(GA_POP[updateGAPopPt][0])) {
                    double[] r0 = r0Collection[updateGAPopPt].get();
                    GA_POP[updateGAPopPt][0] = sumOfSqCost(r0);
                    System.arraycopy(r0, 0, GA_POP[updateGAPopPt], 1 + numParam, r0.length);
                    paramList[PARAM_GA_OPT_LAST_COMPLETED_SIM_COUNTER]
                            = ((Integer) paramList[PARAM_GA_OPT_LAST_COMPLETED_SIM_COUNTER]) + 1;

                }
            }

        } catch (InterruptedException | ExecutionException ex) {
            ex.printStackTrace(System.err);
        }
        exportGAPop();
    }

    private class Callable_CalculateResidue implements Callable<double[]> {

        double[] paramX0;
        double[] rowEnt;
        int numParam;

        public Callable_CalculateResidue(double[] rowEnt, int numParam) {
            this.rowEnt = rowEnt;
            this.numParam = numParam;
            this.paramX0 = Arrays.copyOfRange(rowEnt, 1, 1 + numParam);

        }

        @Override
        public double[] call() throws Exception {
            double[] r0 = getResidual(paramX0);
            return r0;
        }

    }

    @Override
    public void setP0(double[] p0, ParameterConstraintTransform[] constraints) {
        super.setP0(p0, constraints);
        paramList[PARAM_GA_OPT_NUM_PARAM] = new Integer(p0.length);

    }

    @Override
    public void setR0(double[] R0) {
        super.setR0(R0);
        paramList[PARAM_GA_OPT_NUM_R0] = new Integer(R0.length);
    }

    @Override
    public void initialise() {
        if (paramList[PARAM_GA_OPT_POP_FILE] != null) {
            File f = (File) paramList[PARAM_GA_OPT_POP_FILE];
            if (f.exists()) {
                try {
                    System.out.println(this.getClass().getName() + ": Reading pervious GA population from " + f.getAbsolutePath());
                    ObjectInputStream inS = new ObjectInputStream(new FileInputStream(f));
                    GA_POP = (double[][]) inS.readObject();
                    inS.close();
                } catch (IOException | ClassNotFoundException ex) {
                    ex.printStackTrace(System.err);
                }
                sortGA_POP();
                int numValidEnt = 0;

                for (double[] GA_POP_ROW : GA_POP) {
                    if (!Double.isNaN(GA_POP_ROW[0])) {
                        numValidEnt++;
                    }
                }
                paramList[PARAM_GA_OPT_LAST_COMPLETED_SIM_COUNTER] = numValidEnt;
                paramList[PARAM_GA_OPT_LAST_EXPORT_SIM_COUNTER] = numValidEnt;

                f.renameTo(new File(f.getAbsoluteFile()  
                        + "_" + Integer.toString(numValidEnt) 
                        + "_" + Long.toString(System.currentTimeMillis())));

                System.out.println(this.getClass().getName() + ": Number of valid entries = " + numValidEnt);

            }
        }

        if (GA_POP == null) {
            System.out.println(this.getClass().getName() + ": Generating new GA population");
            double[] x0 = getX0(); // x0 - transform parameter
            double[] r0 = getR0();
            paramList[PARAM_GA_OPT_NUM_PARAM] = x0.length;
            paramList[PARAM_GA_OPT_NUM_R0] = r0.length;
            RandomGenerator rng = (RandomGenerator) paramList[PARAM_GA_OPT_RNG];

            GA_POP = new double[((Integer) paramList[PARAM_GA_OPT_POP_SIZE])][1 + x0.length + r0.length];

            // Always include the first one               
            GA_POP[0][0] = Double.NaN;
            System.arraycopy(x0, 0, GA_POP[0], 1, x0.length);
            System.arraycopy(r0, 0, GA_POP[0], 1 + x0.length, r0.length);

            //Populating GA_Population 
            for (int i = 1; i < GA_POP.length; i++) {
                GA_POP[i] = new double[1 + x0.length + +r0.length];
                GA_POP[i][0] = Double.NaN;
                for (int j = 0; j < (Integer) paramList[PARAM_GA_OPT_NUM_PARAM]; j++) {
                    GA_POP[i][j + 1] = rng.nextDouble() * Math.PI / 2; // j+ 1 as the first term is residue
                }
                for (int j = (Integer) paramList[PARAM_GA_OPT_NUM_PARAM] + 1; j < GA_POP[i].length; j++) {
                    GA_POP[i][j] = Double.NaN;
                }

            }

            GA_POP[0][0] = sumOfSqCost(r0);
            sortGA_POP();

        }

        System.out.println("PARAM_GA_OPT_POP_SIZE = " + GA_POP.length);
        System.out.println("PARAM_GA_OPT_PROP_CROSSOVER = " + paramList[PARAM_GA_OPT_POP_CROSSOVER].toString());
        System.out.println("PARAM_GA_USE_PARALLEL = " + paramList[PARAM_GA_OPT_USE_PARALLEL].toString());

    }

    @Override
    public Object setParameter(int paramId, Object newParm) {
        Object oldParam = paramList[paramId];
        paramList[paramId] = newParm;
        return oldParam;
    }

    @Override
    public Object getParameter(int paramId) {
        return paramList[paramId];
    }

    protected void sortGA_POP() {
        Arrays.sort(GA_POP, GA_POP_COMPARATOR);
    }

    protected static double sumOfSqCost(double[] r) {
        double res = 0;
        for (int i = 0; i < r.length; i++) {
            if (Double.isNaN(r[i])) {
                return Double.NaN;
            }
            res += r[i] * r[i];
        }
        return res;
    }

    protected void exportGAPop() {

        int numSim = (Integer) paramList[PARAM_GA_OPT_LAST_COMPLETED_SIM_COUNTER];
        int numExport = (Integer) paramList[PARAM_GA_OPT_LAST_EXPORT_SIM_COUNTER];

        if (numExport < numSim) {

            final double[][] GA_Pop_copy = Arrays.copyOf(GA_POP, GA_POP.length);
            final File tarFile = (File) paramList[PARAM_GA_OPT_POP_FILE];

            Runnable exportPopThread = new Runnable() {
                @Override
                public void run() {
                    ObjectOutputStream objOut = null;
                    try {
                        objOut = new ObjectOutputStream(new FileOutputStream(tarFile));
                        objOut.writeObject(GA_Pop_copy);
                        objOut.close();
                    } catch (IOException ex) {
                        ex.printStackTrace(System.err);
                    }
                }
            };

            ExecutorService executor = Executors.newFixedThreadPool(1);
            executor.submit(exportPopThread);
            executor.shutdown();
            paramList[PARAM_GA_OPT_LAST_EXPORT_SIM_COUNTER] = numSim;

        }

    }

}
