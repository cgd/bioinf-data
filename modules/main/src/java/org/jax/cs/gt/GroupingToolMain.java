/*
* Copyright (c) 2011 The Jackson Laboratory
*
* This is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This software is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this software. If not, see <http://www.gnu.org/licenses/>.
*/

package org.jax.cs.gt;

import java.io.File;
import ch.systemsx.cisd.hdf5.HDF5FactoryProvider;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import ncsa.hdf.hdf5lib.exceptions.*;

/**
 * Group MooseDB Array Genotyping SNP probe intensities and format for Alchemy/MDG
 *
 * @author gbeane
 */
public class GroupingToolMain {

    private static boolean verbose = false;
    private static boolean stripBaseFromID = false;
    private static File inputFile;
    private static File outputFile;
    private static String[] remainingArgs;
    private static Map<String, snpProbeGrouping> snpProbeGroups;

   /**
    * parse the command line arguments and set some of our private member variables
    * @param args the command line arguments
    */
   private static void parseCommandLineArgs(String[] args)
   {

       // users can set a property called verbose to true (case insensitive)
       // to increase how much information is written to the console.
       String verboseProperty = System.getProperty("verbose");

       if (verboseProperty != null && verboseProperty.equalsIgnoreCase("true")) {
           verbose = true;
       }

       // setting the property to stripBaseFromID  will cause the program to
       // strip the last character from the "probesetid" annotations in the HDF5 file
       String stripBaseProperty = System.getProperty("stripBaseFromID");
       
       if (stripBaseProperty != null && stripBaseProperty.equalsIgnoreCase("true")) {
           stripBaseFromID = true;
       }


       // get the two required arguments, the input file name and the output file name
       if (args.length == 2) {
           inputFile = new File(args[0]);
           outputFile = new File(args[1]);
       } else {
           throw new RuntimeException("TODO PRINT USAGE MESSAGE");
       }
   }

   /**
    *
    * @param reader the HDF5 reader for the input file
    * @param numProbes number of probes in the intensity matrix
    * @param sampleIndex index of the sample we want to read
    * @return an array of floats containing the probe intensities for this sample
    */
   private static float[] readSampleIntensities(IHDF5Reader reader, int numProbes, int sampleIndex)
   {
       float block[][] = reader.readFloatMatrixBlock(
               "/intensities", numProbes, 1, 0, sampleIndex);

       float intensities[] = new float[numProbes];
       for (int row = 0; row < numProbes; row++) {
           intensities[row] = block[row][0];
       }

       return intensities;
   }

   /**
    * take a Probe Set ID and trim a base off the end to generate the SNP ID
    * @param id the ID string as contained in the HDF5 file
    * @return id with the last character trimmed off
    */
   private static String prepareProbeSetID(String id)
   {
       return id.substring(0, id.length() - 1);
   }

   /**
    * log normalize an array of intensities
    * @param intensities an array of intensities to normalize, this array is modified
    */
   private static void logNormalize(float intensities[])
   {
       for (int i = 0; i < intensities.length; i++) {
           // probeSet intensities[i] to log2(intensities[i])
           intensities[i] = (float)Math.log(intensities[i]) / (float)Math.log(2);
       }
   }

   /**
    * apply corrections to probe intensities, assumes intensities.length == correction.lenght,
    * which is validated prior to calling
    * @param intensities probe intensity values, this array is modified
    * @param correction correction values for each probe
    */
   private static void applyCorrection(float intensities[], float correction[])
   {
       for (int i = 0; i < intensities.length; i++) {
           intensities[i] += correction[i];
       }
   }

    /**
     * the main entry point
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        long rowCount;
        long columnCount;

        parseCommandLineArgs(args);


        snpProbeGroups = new HashMap<String, snpProbeGrouping>();

        if (verbose) {
            System.out.println("running in verbose mode...");
        }

        System.out.println("Input HDF5 file is " + inputFile);

        //try to open and parse hdf5 file
        IHDF5Reader reader;
        try {
            reader =  HDF5FactoryProvider.get().openForReading(inputFile);

            // read in our annotations
            String sampleNames[] = reader.readStringArray("/annotations/samples/names");
            String probeSetIDs[] = reader.readStringArray("/annotations/probes/probesetid");
            String allele[] = reader.readStringArray("/annotations/probes/allele");
            float correction[] = reader.readFloatArray("/annotations/probes/correction");


            if (verbose) {
                System.out.println(" read " + sampleNames.length + " sample ID annotations");
                System.out.println(" read " + probeSetIDs.length + " probe set ID annotations");
                System.out.println(" read " + allele.length + " allele annotations");
                System.out.println(" read " + correction.length + " correction annotations");
            }

            long[] dimensions = reader.getDataSetInformation("/intensities").getDimensions();

            if(dimensions.length == 2) {
                rowCount = dimensions[0];
                columnCount = dimensions[1];
                if (verbose) {
                    System.out.println(" intensities matrix is " + rowCount + " by " + columnCount);
                }
            } else {
                throw new IllegalArgumentException(
                    "The intensities matrix is expected to have two dimensions but " +
                    "instead it has " + dimensions.length);
            }

            // do some sanity checking
            if (columnCount != sampleNames.length) {
                throw new IllegalArgumentException(
                        "The number of sample name annotations must match the number of " +
                        "columns in the intensities matrix.");
            }

            if (rowCount != probeSetIDs.length) {
                throw new IllegalArgumentException(
                        "The number of probe set ID annotations must match the " +
                        "number of rows in the intensities matrix.");
            }

            if (probeSetIDs.length != allele.length 
                    || probeSetIDs.length != correction.length) {
                throw new IllegalArgumentException(
                        "Probe annotation arrays must all be the same length.");
            }



            if (stripBaseFromID) {

                if (verbose) {
                    System.out.println("Stripping base from Probe Set ID");
                }

                for (int i = 0; i < probeSetIDs.length; i++) {
                    probeSetIDs[i] = prepareProbeSetID(probeSetIDs[i]);
                }
            }

            
            // Open file for output
            BufferedWriter output = new BufferedWriter(new FileWriter(outputFile));

            for (int i = 0; i < sampleNames.length; i++) {

                System.out.println("grouping intensities for sample " + sampleNames[i] + "...");
                

                float intensities[] = readSampleIntensities(reader, (int)rowCount, i);

                applyCorrection(intensities, correction);

                if (verbose) {
                    System.out.println("  log normalizing intensities...");
                }
                logNormalize(intensities);

                // put each intensity value into the right probeset
                for (int j = 0; j < intensities.length; j++) {
                    snpProbeGrouping group = snpProbeGroups.get(probeSetIDs[j]);
                    if (group == null) {
                        group = new snpProbeGrouping(probeSetIDs[j]);
                        snpProbeGroups.put(probeSetIDs[j], group);
                    }
                    if (allele[j].equalsIgnoreCase("A")) {
                        group.insertAIntensity(intensities[j]);
                    } else {
                        group.insertBIntensity(intensities[j]);
                    }                  
                }

                if (verbose) {
                    System.out.println("  grouped " + intensities.length
                            + " intensities into " + snpProbeGroups.size()
                            + " SNPs");
                }

                Collection probeSetCollection = snpProbeGroups.values();
                Iterator itr = probeSetCollection.iterator();

                while (itr.hasNext()) {
                    snpProbeGrouping probeSet = (snpProbeGrouping)itr.next();

                    //sample|SNPID|A|B
                    String line = sampleNames[i] + "\t" + probeSet.getID()
                            + "\t" + snpProbeGrouping.mean(probeSet.getAIntensities())
                            + "\t" + snpProbeGrouping.mean(probeSet.getBIntensities()) + "\n";

                    output.write(line);

                }
            }

            output.close();
        }
        catch (HDF5JavaException e) {
            System.err.println("Error opening input: " + e.getMessage());
            System.exit(1);
        }
        catch (IllegalArgumentException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
        catch (IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
 
    }

}
