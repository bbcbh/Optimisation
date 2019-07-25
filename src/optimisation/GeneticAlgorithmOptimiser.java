package optimisation;

import java.io.BufferedReader;
import transform.ParameterConstraintTransform;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
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
import java.util.regex.Matcher;
import java.util.regex.Pattern;
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
    public static final int PARAM_GA_RO_TOL = PARAM_GA_OPT_LAST_COMPLETED_SIM_COUNTER + 1;
    public static final int PARAM_GA_NUM_SEED_PER_GA_POP_ENT = PARAM_GA_RO_TOL + 1;

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
        new Integer(0),
        //9: PARAM_GA_RO_TOL
        // if used (i.e. not null), then it will have format of [[r0_min, r1_min...], [[r0_max, r1_max...]]
        null,
        //10: PARAM_GA_NUM_SEED_PER_GA_POP_ENT
        // if used, will store preset seed in GA_POP as well
        new Integer(0),};

    protected Number[][] GA_POP = null; // Format: {residue, x0,.., r0,..., seed_0,....}
    protected Number[] bestRow = null;
    double[][] p0_collection = null;

    private final Comparator<Number[]> GA_POP_COMPARATOR = new Comparator<Number[]>() {
        @Override
        public int compare(Number[] t, Number[] t1) {
            // Only compare the first entry
            return Double.compare((double) t[0], (double) t1[0]);
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

    private double[] numberToDoubleArray(Number[] ent) {
        double[] res = new double[ent.length];

        for (int i = 0; i < res.length; i++) {
            res[i] = ent[i].doubleValue();
        }
        return res;
    }

    private Number[] doubleToNumberArray(double[] ent) {
        Number[] res = new Number[ent.length];
        for (int i = 0; i < res.length; i++) {
            res[i] = ent[i];
        }
        return res;
    }

    @Override
    protected void handleResults() {
        sortGA_POP();
        double bestSoFar = (double) GA_POP[0][0];

        int NUM_PARAM = (Integer) paramList[PARAM_GA_OPT_NUM_PARAM];
        int NUM_R0 = (int) paramList[PARAM_GA_OPT_NUM_R0];

        if (bestRow == null || bestSoFar < (double) bestRow[0]) {
            bestRow = Arrays.copyOf(GA_POP[0], GA_POP[0].length);
            final double[] bestX = numberToDoubleArray(Arrays.copyOfRange(bestRow, 1, 1 + NUM_PARAM));
            final double[] bestP = convertParameter(bestX, true);
            final double[] bestR = numberToDoubleArray(Arrays.copyOfRange(bestRow, 1 + NUM_PARAM, 1 + NUM_PARAM + NUM_R0));
            final double bestCost = (double) GA_POP[0][0];
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

        // Check if all results falls into tol range
        if (paramList[PARAM_GA_RO_TOL] != null) {
            double[][] r0_tol = (double[][]) paramList[PARAM_GA_RO_TOL];
            boolean r0_all_in_range = true;

            checkLoop:
            for (Number[] entPop : GA_POP) {
                double[] r0_Ent = numberToDoubleArray(Arrays.copyOfRange(entPop, 1 + NUM_PARAM, 1 + NUM_PARAM + NUM_R0));

                for (int i = 0; i < r0_Ent.length; i++) {
                    if (r0_tol[0] != null) {
                        r0_all_in_range &= r0_tol[0][i] <= r0_Ent[i];
                    }
                    if (r0_tol[1] != null) {
                        r0_all_in_range &= r0_Ent[i] <= r0_tol[1][i];
                    }
                    if (!r0_all_in_range) {
                        break checkLoop;
                    }
                }
            }

            if (r0_all_in_range) {
                java.text.DateFormat df = java.text.DateFormat.getDateTimeInstance();
                System.out.println("GA stopped at " + df.format(new java.util.Date()) + " as all r0 in population falls within r0_tol");

                exportGAPop();
                File GA_pop_file = (File) paramList[PARAM_GA_OPT_POP_FILE];

                try {
                    printNumberArraysToCSV(GA_POP, new File(GA_pop_file.getParent(), GA_pop_file.getName() + "_"
                            + Long.toString(System.currentTimeMillis()) + ".csv"));
                } catch (FileNotFoundException ex) {
                    ex.printStackTrace(System.err);
                }

                setOptStopped(true);
            }

        }

        if (!isOptStopped()) {
            perform_GA();
        }

    }

    protected void perform_GA() {

        sortGA_POP();
        java.text.DateFormat df = java.text.DateFormat.getDateTimeInstance();
        System.out.println("Perform GA at = " + df.format(new java.util.Date()));

        File gaPop = (File) paramList[PARAM_GA_OPT_POP_FILE];

        File exportCSVFile = new File(gaPop.getParentFile(), gaPop.getName()
                + "_" + Long.toString(System.currentTimeMillis()) + ".csv");

        try {
            GeneticAlgorithmOptimiser.printNumberArraysToCSV(GA_POP, exportCSVFile);
            System.out.println("GA Pop at " + gaPop.getAbsolutePath() + " exported as " + exportCSVFile.getAbsolutePath());
        } catch (FileNotFoundException ex) {
            ex.printStackTrace(System.err);
        }

        RandomGenerator rng = (RandomGenerator) paramList[PARAM_GA_OPT_RNG];

        // Selection
        int numToKeep = Math.max(2, Math.round(GA_POP.length * ((Float) paramList[PARAM_GA_OPT_POP_CROSSOVER])));

        int numNewGen = 0;
        for (int row = numToKeep; row < GA_POP.length; row++) {
            Arrays.fill(GA_POP[row], Double.NaN);
            for (int i = 1; i < (Integer) paramList[PARAM_GA_OPT_NUM_PARAM] + 1; i++) {
                if (numNewGen < numToKeep) {
                    // Cross Over                                          
                    double p1 = (double) GA_POP[rng.nextInt(numToKeep)][i];
                    double p2 = (double) GA_POP[rng.nextInt(numToKeep)][i];
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

            double bestSoFar = (double) GA_POP[0][0];
            double worstSoFar = Double.NaN;

            int lastRowPt = GA_POP.length - 1;

            while (Double.isNaN(worstSoFar) && lastRowPt > 0) {
                worstSoFar = (double) GA_POP[lastRowPt][0];
                lastRowPt--;
            }

            if (bestSoFar <= ZERO_TOL
                    || ((worstSoFar - bestSoFar) <= ZERO_TOL)) {
                setOptStopped(true);
            }

            if (!isOptStopped()) {

                int numThreadAtTime = ((Integer) paramList[PARAM_GA_OPT_USE_PARALLEL]);
                int numThreadInPool = 0;
                Future<Object[]>[] r0Collection = new Future[GA_POP.length];

                if (numThreadAtTime > 1) {
                    executor = Executors.newFixedThreadPool(((Integer) paramList[PARAM_GA_OPT_USE_PARALLEL]));
                }

                for (int r = 0; r < GA_POP.length; r++) {
                    Number[] GA_POP_ROW = GA_POP[r];
                    if (Double.isNaN((double) GA_POP_ROW[0])) {
                        // Need to generate new result
                        if (numThreadAtTime <= 1) {
                            Object[] resAndSeed = getResidualAtSeed(numberToDoubleArray(Arrays.copyOfRange(GA_POP_ROW, 1, 1 + numParam)));

                            double[] r0 = (double[]) resAndSeed[0];
                            long[] presetSeed = (long[]) resAndSeed[1];
                            GA_POP_ROW[0] = sumOfSqCost(r0);
                            System.arraycopy(doubleToNumberArray(r0), 0, GA_POP_ROW, 1 + numParam, r0.length);

                            int offset = 1 + numParam + r0.length;
                            for (int i = 0; i < presetSeed.length; i++) {
                                GA_POP_ROW[offset + i] = presetSeed[i];
                            }

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
            Future<Object[]>[] r0Collection,
            int numParam) {
        try {
            executor.shutdown();
            if (!executor.awaitTermination(2, TimeUnit.DAYS)) {
                System.out.println("Inf Thread time-out!");
            }
            for (int updateGAPopPt = 0; updateGAPopPt < r0Collection.length; updateGAPopPt++) {
                if (r0Collection[updateGAPopPt] != null
                        && Double.isNaN((double) GA_POP[updateGAPopPt][0])) {

                    Object[] resAndSeed = r0Collection[updateGAPopPt].get();
                    double[] r0 = (double[]) resAndSeed[0];
                    long[] presetSeed = (long[]) resAndSeed[1];

                    GA_POP[updateGAPopPt][0] = sumOfSqCost(r0);
                    System.arraycopy(doubleToNumberArray(r0), 0, GA_POP[updateGAPopPt], 1 + numParam, r0.length);
                    int offset = 1 + numParam + r0.length;

                    for (int i = 0; i < presetSeed.length; i++) {
                        GA_POP[updateGAPopPt][offset + i] = presetSeed[i];
                    }

                    paramList[PARAM_GA_OPT_LAST_COMPLETED_SIM_COUNTER]
                            = ((Integer) paramList[PARAM_GA_OPT_LAST_COMPLETED_SIM_COUNTER]) + 1;

                }
            }

        } catch (InterruptedException | ExecutionException ex) {
            ex.printStackTrace(System.err);
        }
        exportGAPop();
    }

    private class Callable_CalculateResidue implements Callable<Object[]> {

        double[] paramX0;
        Number[] rowEnt;
        int numParam;

        public Callable_CalculateResidue(Number[] rowEnt, int numParam) {
            this.rowEnt = rowEnt;
            this.numParam = numParam;
            this.paramX0 = numberToDoubleArray(Arrays.copyOfRange(rowEnt, 1, 1 + numParam));

        }

        @Override
        public Object[] call() throws Exception {
            return getResidualAtSeed(paramX0);
        }

    }

    @Override
    public void setP0(double[] p0, ParameterConstraintTransform[] constraints) {

        int numP0 = p0.length / constraints.length;
        p0_collection = new double[numP0][constraints.length];

        for (int i = 0; i < numP0; i++) {
            System.arraycopy(p0, i * constraints.length, p0_collection[i], 0, constraints.length);
        }

        super.setP0(Arrays.copyOf(p0, constraints.length), constraints);
        paramList[PARAM_GA_OPT_NUM_PARAM] = new Integer(constraints.length);

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
                    GA_POP = (Number[][]) inS.readObject();
                    inS.close();
                } catch (IOException | ClassNotFoundException ex) {
                    ex.printStackTrace(System.err);
                }
                sortGA_POP();

                int numValidEnt = 0;
                for (Number[] GA_POP_ROW : GA_POP) {
                    if (!Double.isNaN((double) GA_POP_ROW[0])) {
                        numValidEnt++;
                    }
                }
                paramList[PARAM_GA_OPT_LAST_COMPLETED_SIM_COUNTER] = numValidEnt;
                paramList[PARAM_GA_OPT_LAST_EXPORT_SIM_COUNTER] = numValidEnt;

                f.renameTo(new File(f.getAbsoluteFile()
                        + "_" + Long.toString(System.currentTimeMillis())));

                System.out.println(this.getClass().getName() + ": Number of valid entries = " + numValidEnt);

            }

        }

        if (GA_POP == null) {

            double[] x0 = getX0(); // x0 - transform parameter
            double[] r0 = getR0();

            paramList[PARAM_GA_OPT_NUM_PARAM] = x0.length;
            paramList[PARAM_GA_OPT_NUM_R0] = r0.length;

            int numSeedPerEnt = (Integer) paramList[PARAM_GA_NUM_SEED_PER_GA_POP_ENT];
            RandomGenerator rng = (RandomGenerator) paramList[PARAM_GA_OPT_RNG];

            GA_POP = new Number[((Integer) paramList[PARAM_GA_OPT_POP_SIZE])][1 + x0.length + r0.length + numSeedPerEnt];
            
            
            boolean genGA = true;

            if (paramList[PARAM_GA_OPT_POP_FILE] != null) {
                File f = (File) paramList[PARAM_GA_OPT_POP_FILE];
                // Test for pre-exisiting population file            
                final Pattern Pattern_PRE_GA_POP_CSV = Pattern.compile(f.getName() + "_(\\d+).csv");

                File[] PRE_GA_POP_CSV_FILES = f.getParentFile().listFiles(new FileFilter() {
                    @Override
                    public boolean accept(File file) {
                        return Pattern_PRE_GA_POP_CSV.matcher(file.getName()).matches();

                    }
                });

                if (PRE_GA_POP_CSV_FILES.length > 0) {

                    Arrays.sort(PRE_GA_POP_CSV_FILES, new Comparator<File>() {
                        @Override
                        public int compare(File t, File t1) {
                            Matcher m = Pattern_PRE_GA_POP_CSV.matcher(t.getName());
                            Matcher m1 = Pattern_PRE_GA_POP_CSV.matcher(t1.getName());

                            m.matches();
                            m1.matches();

                            return Long.compare(Long.parseLong(m.group(1)),
                                    Long.parseLong(m1.group(1)));

                        }
                    });

                    System.out.println(this.getClass().getName() + ": Number of GA_POP CSV = " + PRE_GA_POP_CSV_FILES.length
                            + ". Importing the latest one (" + PRE_GA_POP_CSV_FILES[PRE_GA_POP_CSV_FILES.length - 1].getAbsolutePath()
                            + ") as GA_POP.");

                    File selCSV = PRE_GA_POP_CSV_FILES[PRE_GA_POP_CSV_FILES.length - 1];

                    try (BufferedReader reader = new BufferedReader(new FileReader(selCSV))) {
                        String line;

                        int r = 0;
                        while ((line = reader.readLine()) != null) {

                            if (r > GA_POP.length) {
                                System.out.println(this.getClass().getName() + ": Number of row in "
                                        + selCSV.getAbsolutePath()
                                        + " > number of row in GA_POP. Row #" + r + " ignored");
                                break;
                            }

                            String[] ent = line.split(",");
                            for (int i = 0; i < ent.length; i++) {
                                if (i > GA_POP[r].length) {
                                    System.out.println(this.getClass().getName() + ": Number of col in "
                                            + selCSV.getAbsolutePath()
                                            + " > number of col in GA_POP. Col #" + i + " ignored");
                                } else {
                                    if (i < 1 + x0.length + r0.length) {
                                        GA_POP[r][i] = new Double(ent[i]);
                                    } else {
                                        GA_POP[r][i] = new Long(ent[i]);
                                    }
                                }
                            }
                            r++;
                        }
                        
                        genGA = false;

                    } catch (IOException ex) {
                        ex.printStackTrace(System.err);                        
                    }

                }

            }

            if (genGA) {
                
                System.out.println(this.getClass().getName() + ": Generating new GA population");

                // Always include the first one               
                GA_POP[0][0] = Double.NaN;
                System.arraycopy(doubleToNumberArray(x0), 0, GA_POP[0], 1, x0.length);
                System.arraycopy(doubleToNumberArray(r0), 0, GA_POP[0], 1 + x0.length, r0.length);

                //Populating GA_Population 
                for (int i = 1; i < GA_POP.length; i++) {
                    GA_POP[i] = new Number[1 + x0.length + +r0.length + numSeedPerEnt];
                    GA_POP[i][0] = Double.NaN;
                    for (int j = 0; j < (Integer) paramList[PARAM_GA_OPT_NUM_PARAM]; j++) {
                        if (i < p0_collection.length) {
                            GA_POP[i][j + 1] = constraints[j].toUncontrainted(p0_collection[i][j]);

                        } else {
                            GA_POP[i][j + 1] = rng.nextDouble() * Math.PI / 2; // j+ 1 as the first term is residue
                        }
                    }
                    for (int j = (Integer) paramList[PARAM_GA_OPT_NUM_PARAM] + 1; j < GA_POP[i].length; j++) {
                        GA_POP[i][j] = Double.NaN;
                    }

                }
                GA_POP[0][0] = sumOfSqCost(r0);
            }

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

    public static void printNumberArraysToCSV(Number[][] entArr, File csvfile) throws FileNotFoundException {
        try (PrintWriter pri = new PrintWriter(csvfile)) {
            for (Number[] ent : entArr) {
                StringBuilder line = new StringBuilder();
                for (Number row : ent) {
                    if (line.length() != 0) {
                        line.append(',');
                    }
                    line.append(row);
                }
                pri.println(line.toString());
            }
        }
    }

    protected void exportGAPop() {

        final int numSim = (Integer) paramList[PARAM_GA_OPT_LAST_COMPLETED_SIM_COUNTER];
        final int numExport = (Integer) paramList[PARAM_GA_OPT_LAST_EXPORT_SIM_COUNTER];

        //System.out.println("Export GA Pop called. Number of completed sim = " + numSim + " Number of exported sim = " +  numExport);
        if (numExport < numSim) {

            final Number[][] GA_Pop_copy = Arrays.copyOf(GA_POP, GA_POP.length);
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

                    int numValidEnt = 0;
                    for (Number[] GA_POP_ROW : GA_Pop_copy) {
                        if (!Double.isNaN((double) GA_POP_ROW[0])) {
                            numValidEnt++;
                        }
                    }

                    System.out.println(this.getClass().getName() + ": Pop Exported at " + tarFile.getAbsolutePath()
                            + " Sim #" + numSim);
                }
            };

            ExecutorService executor = Executors.newFixedThreadPool(1);
            executor.submit(exportPopThread);
            executor.shutdown();

            paramList[PARAM_GA_OPT_LAST_EXPORT_SIM_COUNTER] = numSim;

        }

    }

    public Object[] getResidualAtSeed(double[] x) {
        int numSeedPerEnt = (Integer) paramList[PARAM_GA_NUM_SEED_PER_GA_POP_ENT];
        long[] presetSeed = new long[numSeedPerEnt];

        if (numSeedPerEnt > 0) {
            RandomGenerator rng = (RandomGenerator) paramList[PARAM_GA_OPT_RNG];
            for (int i = 0; i < presetSeed.length; i++) {
                presetSeed[i] = rng.nextLong();
            }
            func.setPreset_Seed(presetSeed);
        }

        Object[] resAndSeed = new Object[2];

        resAndSeed[0] = super.getResidual(x);
        resAndSeed[1] = presetSeed;

        return resAndSeed;
    }

}
