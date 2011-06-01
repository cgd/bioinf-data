/*
 * Copyright (c) 2010 The Jackson Laboratory
 * 
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.jax.emma;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jax.bioinfdata.genocall.GenotypeCallMatrix;
import org.jax.bioinfdata.genocall.GenotypeCallMatrix.AlleleCallCode;
import org.jax.util.io.IllegalFormatException;
import org.jax.util.nativeutil.NativeLibraryUtilities;

/**
 * For performing an association test using native EMMA library
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class EMMAAssociationTest
{
    private static final Logger LOG = Logger.getLogger(EMMAAssociationTest.class.getName());
    
    /**
     * You must have at least this many strains in common between the genotype
     * and phenotype data before an association test is valid.
     */
    public static final int MIN_STRAIN_COUNT = 3;
    
    static
    {
        NativeLibraryUtilities.loadNativeLibrary("emma");
    }
    
    /**
     * Emma scan on the given genotype matrix
     * @param genoMat
     *          the genotype matrix
     * @param phenoFileName
     *          the phenotype file
     * @param phenotype
     *          the name of the phenotype to scan (you can use null here only
     *          if there is only a single phenotype in the phenotype file)
     * @param sexToScan
     *          the sex type used to filter phenotype data
     * @return
     *          the flattened kinship matrix
     * @throws IllegalFormatException
     *          if there is a problem with how with how the file is formatted
     * @throws IOException
     *          if there is a problem with file IO while reading the flat file
     */
    public double[] emmaScan(
            GenotypeCallMatrix genoMat,
            String phenoFileName,
            String phenotype,
            SexFilter sexToScan) throws IllegalFormatException, IOException
    {
        String[] headerStrains = genoMat.getSampleIds();
        
        // now get the strains in common with phenotype data
        MPDIndividualStrainPhenotypeParser phenoParser = new MPDIndividualStrainPhenotypeParser();
        FileInputStream phenoIn;
        if(phenotype == null || phenotype.length() == 0)
        {
            phenoIn = new FileInputStream(phenoFileName);
            Set<String> phenos = phenoParser.parseAvailablePhenotypes(phenoIn);
            phenoIn.close();
            
            if(phenos.size() != 1)
            {
                throw new IllegalFormatException();
            }
            else
            {
                phenotype = phenos.iterator().next();
            }
        }
        phenoIn = new FileInputStream(phenoFileName);
        Set<String> phenoStrains = phenoParser.parseAvailableStrainNames(
                phenotype,
                phenoIn,
                sexToScan);
        phenoIn.close();
        
        Set<String> commonStrainSet = new HashSet<String>(phenoStrains);
        commonStrainSet.retainAll(Arrays.asList(headerStrains));
        int strainCount = commonStrainSet.size();
        if(strainCount < MIN_STRAIN_COUNT)
        {
            throw new IllegalArgumentException(
                    "In order to perform a haplotype test there must be at least " +
                    MIN_STRAIN_COUNT + " strains in common between the phenotype " +
                    "and genotype data but the number of strains found in common " +
                    "is " + strainCount);
        }
        String[] commonStrainArray = commonStrainSet.toArray(new String[0]);
        Arrays.sort(commonStrainArray);
        int[] commonStrainIndices = new int[strainCount];
        
        {
            List<String> headerList = Arrays.asList(headerStrains);
            for(int i = 0; i < strainCount; i++)
            {
                commonStrainIndices[i] = headerList.indexOf(commonStrainArray[i]);
            }
        }
        
        // read the phenotype data
        if(LOG.isLoggable(Level.FINE))
        {
            LOG.fine("the phenotype is: " + phenotype);
        }
        
        phenoIn = new FileInputStream(phenoFileName);
        Map<String, List<Double>> phenoData = phenoParser.parsePhenotypesFromStream(
                phenotype,
                phenoIn,
                sexToScan,
                commonStrainSet);
        
        if(LOG.isLoggable(Level.FINE))
        {
            LOG.fine("the # of phenotypes is: " + phenoData.size());
        }
        
        double[] phenotypeMeans = new double[phenoData.size()];
        for(int i = 0; i < strainCount; i++)
        {
            phenotypeMeans[i] = 0.0;
            List<Double> currData = phenoData.get(commonStrainArray[i]);
            for(Double measure: currData)
            {
                phenotypeMeans[i] += measure;
            }
            phenotypeMeans[i] /= currData.size();
        }
        
        // flatten the genotype matrix and convert to double
        byte[][] callMatrix = genoMat.getCallMatrix();
        double[] flatCallValues = new double[strainCount * callMatrix.length];
        for(int rowIndex = 0; rowIndex < callMatrix.length; rowIndex++)
        {
            for(int strainIndex = 0; strainIndex < strainCount; strainIndex++)
            {
                int flatIndex = rowIndex * strainCount + strainIndex;
                flatCallValues[flatIndex] = toCallValue(AlleleCallCode.numCodeToEnum(
                        callMatrix[rowIndex][commonStrainIndices[strainIndex]]));
            }
        }
        
        // calculate kinship matrix and do the scan
        double[] kinship = calculateKinship(strainCount, flatCallValues);
        return emmaScan(strainCount, phenotypeMeans, flatCallValues, kinship);
    }
    
    private double toCallValue(AlleleCallCode callCode)
    {
        switch(callCode)
        {
            case ACall: return 1.0;
            case BCall: return 0.0;
            case HCall: return 0.5;
            case NCall: return Double.NaN;
            default: throw new IllegalArgumentException(
                    "unexpected call code: " + callCode);
        }
    }
    
    /**
     * Native function for performing an EMMA scan
     * @param strainCount   the number of strains
     * @param phenos        the (flattened) phenotype matrix
     * @param genos         the (flattened) genotype matrix
     * @param kinship       the (flattened) kinship matrix. this can be null
     *                      in which case it will be estimated from the given
     *                      genotypes
     * @return              the resulting (flattened) pvalue matrix
     */
    private static native double[] emmaScan(
            int strainCount,
            double[] phenos,
            double[] genos,
            double[] kinship);
    
    /**
     * Native function for performing an EMMA scan
     * @param strainCount   the number of strains
     * @param genos         the (flattened) genotype matrix
     * @return              the resulting (flattened) kinship matrix
     */
    private static native double[] calculateKinship(
            int strainCount,
            double[] genos);
}
