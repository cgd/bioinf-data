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
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jax.util.io.CommonFlatFileFormat;
import org.jax.util.io.FlatFileReader;
import org.jax.util.io.IllegalFormatException;
import org.jax.util.nativeutil.NativeLibraryUtilities;

/**
 * For performing an association test using native EMMA library
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class EMMAAssociationTest
{
    private static final Logger LOG = Logger.getLogger(EMMAAssociationTest.class.getName());
    
    private enum AlleleCallCode {ACall, BCall, HCall, NCall}
    
    static
    {
        NativeLibraryUtilities.loadNativeLibrary("emma");
    }
    
    /**
     * Emma scan using flat files
     * @param genoFileName
     *          the genotype file
     * @param aAlleleColumn
     *          the column index for the A allele
     * @param bAlleleColumn
     *          the column index for the B allele
     * @param firstGenotypeColumn
     *          the 1st genotype column
     * @param lastGenotypeColumnExclusive
     *          the index after the last genotype column. You can use -1 to indicate that
     *          all of the remaining columns after firstGenotypeColumn are
     *          genotype columns
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
            String genoFileName,
            int aAlleleColumn,
            int bAlleleColumn,
            int firstGenotypeColumn,
            int lastGenotypeColumnExclusive,
            String phenoFileName,
            String phenotype,
            SexFilter sexToScan) throws IllegalFormatException, IOException
    {
        // start with the geno headers
        FlatFileReader genoFFR = new FlatFileReader(
                new FileReader(genoFileName),
                CommonFlatFileFormat.CSV_UNIX);
        String[] currRow = genoFFR.readRow();
        if(currRow == null)
        {
            throw new IllegalFormatException("Failed to read the header");
        }
        
        if(lastGenotypeColumnExclusive == -1)
        {
            lastGenotypeColumnExclusive = currRow.length;
        }
        String[] headerStrains = new String[lastGenotypeColumnExclusive - firstGenotypeColumn];
        for(int i = 0; i < headerStrains.length; i++)
        {
            headerStrains[i] = currRow[i + firstGenotypeColumn];
        }
        
        // now get the strains in common with phenotype data
        MPDIndividualStrainPhenotypeParser phenoParser = new MPDIndividualStrainPhenotypeParser();
        FileInputStream phenoIn = new FileInputStream(phenoFileName);
        Set<String> phenoStrains = phenoParser.parseAvailableStrainNames(phenoIn);
        phenoIn.close();
        
        Set<String> commonStrainSet = new HashSet<String>(phenoStrains);
        commonStrainSet.retainAll(Arrays.asList(headerStrains));
        int strainCount = commonStrainSet.size();
        String[] commonStrainArray = commonStrainSet.toArray(new String[0]);
        Arrays.sort(commonStrainArray);
        int[] commonStrainIndices = new int[strainCount];
        
        {
            List<String> headerRowList = Arrays.asList(currRow);
            for(int i = 0; i < strainCount; i++)
            {
                commonStrainIndices[i] = headerRowList.indexOf(commonStrainArray[i]);
            }
        }
        
        // TODO need failure if there are fewer than 3 strains in common
        
        // read the phenotype data
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
        
        // read the genotype data
        List<double[]> callValues = new LinkedList<double[]>();
        while((currRow = genoFFR.readRow()) != null)
        {
            double[] currSnpGenos = new double[strainCount];
            for(int strainIndex = 0; strainIndex < strainCount; strainIndex++)
            {
                currSnpGenos[strainIndex] = toCallValue(
                        currRow[aAlleleColumn],
                        currRow[bAlleleColumn],
                        currRow[commonStrainIndices[strainIndex]]);
            }
            callValues.add(currSnpGenos);
        }
        
        // flatten the genotype matrix
        double[] flatCallValues = new double[strainCount * callValues.size()];
        Iterator<double[]> iter = callValues.iterator();
        for(int rowIndex = 0; iter.hasNext(); rowIndex++)
        {
            double[] currSnpGenos = iter.next();
            for(int strainIndex = 0; strainIndex < currSnpGenos.length; strainIndex++)
            {
                int flatIndex = rowIndex * strainCount + strainIndex;
                flatCallValues[flatIndex] = currSnpGenos[strainIndex];
            }
        }
        
        // calculate kinship matrix and do the scan
        double[] kinship = calculateKinship(strainCount, flatCallValues);
        return emmaScan(strainCount, phenotypeMeans, flatCallValues, kinship);
    }
    
    private double toCallValue(String aAllele, String bAllele, String genoCall) throws IllegalFormatException
    {
        switch(toCallCode(aAllele, bAllele, genoCall))
        {
            case ACall: return 1.0;
            case BCall: return 0.0;
            case HCall: return 0.5;
            case NCall: return Double.NaN;
            default:    throw new IllegalStateException("this should never happen");
        }
    }
    
    private AlleleCallCode toCallCode(String aAllele, String bAllele, String genoCall) throws IllegalFormatException
    {
        aAllele = aAllele.toUpperCase();
        bAllele = bAllele.toUpperCase();
        genoCall = genoCall.toUpperCase();
        if(genoCall.equals(aAllele))
        {
            return AlleleCallCode.ACall;
        }
        else if(genoCall.equals(bAllele))
        {
            return AlleleCallCode.BCall;
        }
        else if(genoCall.equals("H") || genoCall.equals("HH"))
        {
            return AlleleCallCode.HCall;
        }
        else if(genoCall.length() == 0 || genoCall.equals("N") || genoCall.equals("-") || genoCall.equals("NN"))
        {
            return AlleleCallCode.NCall;
        }
        else
        {
            //throw new IllegalFormatException("Unknown genotype: " + genoCall);
            return AlleleCallCode.NCall;
        }
    }
    
//    /**
//     * Calculate the kinship values
//     * @param genoData
//     *          the genotype data to base it on
//     * @param strains
//     *          a map of strain names
//     * @return
//     *          the kinship
//     * @throws
//     *          IOException
//     */
//    public double[] calculateKinship(
//            GenomeDataSource genoData,
//            Set<String> strains)
//    throws IOException
//    {
//        strains = new HashSet<String>(strains);
//        strains.retainAll(genoData.getAvailableStrains());
//        String[] commonStrains = strains.toArray(new String[0]);
//        Arrays.sort(commonStrains);
//        int strainCount = commonStrains.length;
//        
//        int snpCount = 0;
//        for(ChromosomeDataSource currChr : genoData.getChromosomeDataSources().values())
//        {
//            snpCount += currChr.getSnpPositionInputStream().getSnpCount();
//        }
//        
//        double[] genos = new double[snpCount * strainCount];
//        int currStartIndex = 0;
//        for(ChromosomeDataSource currChr : genoData.getChromosomeDataSources().values())
//        {
//            SdpInputStream sdpStream = currChr.getSdpInputStream(commonStrains);
//            int currSnpCount = (int)currChr.getSnpPositionInputStream().getSnpCount();
//            for(int i = 0; i < currSnpCount; i++)
//            {
//                int currSnp = currStartIndex + i;
//                BitSet currSdp = sdpStream.getNextSdp();
//                for(int strainIndex = 0; strainIndex < strainCount; strainIndex++)
//                {
//                    double currCall = currSdp.get(strainIndex) ? 1.0 : 0.0;
//                    genos[currSnp * strainCount + strainIndex] = currCall;
//                }
//            }
//        }
//        
//        return calculateKinship(strainCount, genos);
//    }
//    
//    /**
//     * Perform a scan on the given chromosome using EMMA
//     * @param chrDataSource
//     *          the chromosome data source
//     * @param phenotypeDataSource
//     *          the phenotype data source
//     * @param kinship
//     *          the kinship matrix (if null it's calculated based on data)
//     * @return
//     *          the p-values
//     * @throws IOException
//     */
//    public double[] emmaScan(
//            ChromosomeDataSource chrDataSource,
//            PhenotypeDataSource phenotypeDataSource,
//            double[] kinship)
//            throws IOException
//    {
//        int snpCount = (int)chrDataSource.getSnpPositionInputStream().getSnpCount();
//        Map<String, List<Double>> phenotypeDataMap = phenotypeDataSource.getPhenotypeData();
//        phenotypeDataMap.keySet().retainAll(chrDataSource.getAvailableStrains());
//        String[] commonStrains = phenotypeDataMap.keySet().toArray(new String[0]);
//        Arrays.sort(commonStrains);
//        int strainCount = commonStrains.length;
//        
//        String[] sortedCommonStrains = phenotypeDataMap.keySet().toArray(
//                new String[phenotypeDataMap.size()]);
//        Arrays.sort(sortedCommonStrains);
//        
//        double[] genos = new double[snpCount * strainCount];
//        SdpInputStream sdpStream = chrDataSource.getSdpInputStream(commonStrains); // TODO FIXME
//        for(int snpIndex = 0; snpIndex < snpCount; snpIndex++)
//        {
//            BitSet currSDP = sdpStream.getNextSdp();
//            for(int strainIndex = 0; strainIndex < strainCount; strainIndex++)
//            {
//                double currCall = currSDP.get(strainIndex) ? 1.0 : 0.0;
//                genos[snpIndex * strainCount + strainIndex] = currCall;
//            }
//        }
//        
//        if(kinship == null)
//        {
//            kinship = calculateKinship(strainCount, genos);
//        }
//        
//        double[] phenotypeMeans = new double[strainCount];
//        for(int strainIndex = 0; strainIndex < strainCount; strainIndex++)
//        {
//            phenotypeMeans[strainIndex] = StatisticUtilities.calculateMean(
//                    phenotypeDataMap.get(commonStrains[strainIndex]));
//        }
//        
//        return emmaScan(strainCount, phenotypeMeans, genos, kinship);
//    }
    
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
