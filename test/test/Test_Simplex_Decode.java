/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.Arrays;
import transform.ParameterConstraintTransform;

/**
 *
 * @author Bhui
 */
public class Test_Simplex_Decode {

    public static void main(String[] arg) throws FileNotFoundException, IOException, ClassNotFoundException {
        String filepath = "C:\\Users\\Bhui\\Desktop\\FTP\\RMP\\OptResults\\ParamSimplex.obj_1520156438339";
        String constraintFilePath = "C:\\Users\\Bhui\\Desktop\\FTP\\RMP\\OptResults\\ParamConstriants.csv";

        File simplexFile = new File(filepath);
        File constraintFile = new File(constraintFilePath);

        System.out.print("Reading simplex from " + simplexFile.getAbsolutePath() + "....");

        double[][] sX = null;
        double[][] sR = null;

        try (ObjectInputStream objStr = new ObjectInputStream(new FileInputStream(simplexFile))) {
            sX = (double[][]) objStr.readObject();
            sR = (double[][]) objStr.readObject();
        }

        System.out.println(" done");

        if (sX != null) {
            System.out.println("SX:");
            for (double[] r : sX) {
                System.out.println(Arrays.toString(r));
            }

            if (constraintFile.isFile()) {
                ParameterConstraintTransform[] constraints;
                try (BufferedReader constraintReader = new BufferedReader(new FileReader(constraintFile))) {
                    int lnNum = 0;
                    String line;
                    while (constraintReader.readLine() != null) {
                        lnNum++;
                    }
                    constraints = new ParameterConstraintTransform[lnNum];
                    lnNum = 0;
                    BufferedReader constraintReader2 = new BufferedReader(new FileReader(constraintFile));

                    while ((line = constraintReader2.readLine()) != null) {
                        String[] ent = line.split(",");
                        constraints[lnNum] = new transform.ParameterConstraintTransformSineCurve(new double[]{
                            Double.parseDouble(ent[0]), Double.parseDouble(ent[1])});
                        lnNum++;
                    }
                }

                System.out.println("SX (raw parameter):");
                
                for (double[] r : sX) {
                    StringBuilder s = new StringBuilder();   
                    int i = 0;
                    for(double re : r){                        
                        if(s.length() == 0){
                            s.append('[');
                        }else{
                            s.append(", ");                            
                        }
                        s.append(Double.toString(constraints[i].toContrainted(re)));                        
                        i++;
                    }
                    s.append(']');
                    System.out.println(s.toString());
                }
            }

        }

        if (sR != null) {
            System.out.println("SR:");
            for (double[] r : sR) {
                System.out.println(Arrays.toString(r));
            }
        }

    }

}
