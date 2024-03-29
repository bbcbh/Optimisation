package test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.ObjectInputStream;
import java.util.Arrays;
import optimisation.GeneticAlgorithmOptimiser;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import transform.ParameterConstraintTransform;

public class Test_GA_Pop_Analysis {
    
    public static void main(String[] arg) {
        final String BASE_DIR_STR
                = "C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\Test\\Optim_Preval_Range_Fit";
                //"C:\\Users\\Bhui\\Desktop\\FTP\\MSM\\Optim_7_Range_GA";
        
        final File BASE_DIR = new File(BASE_DIR_STR);
        
        File[] filelist = new File[]{
            new File(BASE_DIR, "GA_POP.obj")
        };
        
        File maxEntryFile = null;
        int maxEntryNumEntry = 0;
        
        File lastEntryFile = null;
        int lastEntryFileNumEntry = 0;
        long lastEntryDateStamp = 0;
        
        File currentEntryFile = null;
        int currentEntryNumEntry = 0;
        
        Arrays.sort(filelist);
        
        for (File gaPop : filelist) {
            Number[][] GA_POP = null;
            
            String[] splitFilename = gaPop.getName().split("_");
            
            try {
                ObjectInputStream inS = new ObjectInputStream(new FileInputStream(gaPop));
                GA_POP = (Number[][]) inS.readObject();
                inS.close();
            } catch (Exception ex) {
                System.out.println("Error in reading " + gaPop.getAbsolutePath() + ". It might not be a valid GA_POP_file ");
                
            }
            if (GA_POP != null) {
                int numValidEnt = 0;
                
                int numCol = GA_POP[0].length;
                
                for (Number[] GA_POP_ROW : GA_POP) {
                    if (!Double.isNaN((double) GA_POP_ROW[0])) {
                        numValidEnt++;
                    }
                }
                System.out.println("File = " + gaPop.getAbsolutePath() + ": # Entry " + numValidEnt + ", # Col = " + numCol);
                
                if (filelist.length > 1 && maxEntryNumEntry <= numValidEnt) {
                    maxEntryNumEntry = numValidEnt;
                    maxEntryFile = gaPop;
                }
                
                if (splitFilename.length == 2) {
                    currentEntryFile = gaPop;
                    currentEntryNumEntry = numValidEnt;
                } else {
                    long dateStamp = Long.parseLong(splitFilename[splitFilename.length - 1]);
                    if (dateStamp > lastEntryDateStamp) {
                        lastEntryFile = gaPop;
                        lastEntryFileNumEntry = numValidEnt;
                        lastEntryDateStamp = dateStamp;
                    }
                }
                
            }
            
            File exportCSVFile = new File(gaPop.getParentFile(), gaPop.getName()
                    + "_Analysis_" + Long.toString(System.currentTimeMillis()) + ".csv");
            
            try {
                GeneticAlgorithmOptimiser.printNumberArraysToCSV(GA_POP, exportCSVFile);                
                System.out.println("GA Pop at " + gaPop.getAbsolutePath() + " exported as " + exportCSVFile.getAbsolutePath());
            } catch (FileNotFoundException ex) {
                ex.printStackTrace(System.err);
            }
        }
        
        System.out.println();
        if (currentEntryFile != null) {
            System.out.println("Current Pop file: " + currentEntryFile.getAbsolutePath() + " with " + currentEntryNumEntry + " entries");
            extractPopStat(currentEntryFile, currentEntryNumEntry, BASE_DIR);
            System.out.println();
        }
        
        if (maxEntryFile != null) {
            System.out.println("Pop file with max number of entries: " + maxEntryFile.getAbsolutePath() + " with " + maxEntryNumEntry + " entries");
            extractPopStat(maxEntryFile, maxEntryNumEntry, BASE_DIR);
            System.out.println();
        }
        
        if (lastEntryFile != null) {
            System.out.println("Pop file with latest timestamp: " + lastEntryFile.getAbsolutePath() + " with " + lastEntryFileNumEntry + " entries");
            extractPopStat(lastEntryFile, lastEntryFileNumEntry, BASE_DIR);
            System.out.println();
        }
        
    }
    
    private static void extractPopStat(File maxEntryFile, int maxEntry, final File BASE_DIR) throws MathIllegalArgumentException {
        
        System.out.println("Analysis GA Pop at " + maxEntryFile.getAbsolutePath());
        
        Number[][] GA_POP = null;
        
        try {
            ObjectInputStream inS = new ObjectInputStream(new FileInputStream(maxEntryFile));
            GA_POP = (Number[][]) inS.readObject();
            inS.close();
        } catch (Exception ex) {
            System.out.println("Error in reading " + maxEntryFile.getAbsolutePath() + " for extractPopStat.");
        }
        
        if (GA_POP != null) {
            File contraintFile = new File(BASE_DIR, "ParamConstriants.csv");
            
            ParameterConstraintTransform[] constraints = null;
            try (BufferedReader constraintReader = new BufferedReader(new FileReader(contraintFile))) {
                int numParam = 0;
                while (constraintReader.readLine() != null) {
                    numParam++;
                }
                
                String line;
                int p = 0;
                
                constraints = new ParameterConstraintTransform[numParam];
                
                BufferedReader constraintReader2 = new BufferedReader(new FileReader(contraintFile));
                
                while ((line = constraintReader2.readLine()) != null) {
                    String[] ent = line.split(",");
                    constraints[p] = new transform.ParameterConstraintTransformSineCurve(new double[]{
                        Double.parseDouble(ent[0]), Double.parseDouble(ent[1])});
                    p++;
                }
                constraintReader2.close();
            } catch (Exception ex) {
                System.out.println("Error in reading " + contraintFile.getAbsolutePath() + ".");
            }
            
            if (constraints != null) {
                System.out.println("Number of parameters (from " + contraintFile.getName() + ") = " + constraints.length);
                
                int bestRow = -1;
                double bestRes = Double.POSITIVE_INFINITY;
                
                double[][] entryByParam = new double[constraints.length][maxEntry];
                for (int r = 0; r < maxEntry; r++) {
                    
                    if ((double) GA_POP[r][0] < bestRes) {
                        bestRow = r;
                        bestRes = (double) GA_POP[r][0];
                    }
                    
                    for (int p = 0; p < constraints.length; p++) {
                        entryByParam[p][r] = (double) GA_POP[r][p + 1];
                    }
                }
                
                if (bestRow >= 0) {
                    System.out.println("Best row at #" + Integer.toString(bestRow));
                }
                
                Percentile percentile = new Percentile();
                
                for (int p = 0; p < constraints.length; p++) {
                    
                    Arrays.sort(entryByParam[p]);
                    System.out.println("Param #" + (p + 1) + ":");
                    
                    if (bestRow >= 0) {
                        System.out.println("Best   = " + Double.toString(constraints[p].toContrainted((double) GA_POP[bestRow][p + 1])));
                    }
                    System.out.println("Range  = ["
                            + Double.toString(constraints[p].toContrainted(entryByParam[p][0]))
                            + ","
                            + Double.toString(constraints[p].toContrainted(entryByParam[p][maxEntry - 1]))
                            + "]");
                    System.out.println("Median = "
                            + ((entryByParam[p].length % 2 == 0)
                                    ? Double.toString(constraints[p].toContrainted(
                                            (entryByParam[p][entryByParam[p].length / 2]
                                            + entryByParam[p][entryByParam[p].length / 2 - 1]) / 2))
                                    : Double.toString(constraints[p].toContrainted(entryByParam[p][entryByParam[p].length / 2]))));
                    System.out.println("25%    = "
                            + Double.toString(constraints[p].toContrainted(percentile.evaluate(entryByParam[p], 25))));
                    
                    System.out.println("75%    = "
                            + Double.toString(constraints[p].toContrainted(percentile.evaluate(entryByParam[p], 75))));
                    
                }
                
            }
            
        }
    }
    
}
