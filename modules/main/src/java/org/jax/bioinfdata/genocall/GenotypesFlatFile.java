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

package org.jax.bioinfdata.genocall;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.jax.util.datastructure.SequenceUtilities;
import org.jax.util.io.FlatFileReader;
import org.jax.util.io.FlatFileWriter;
import org.jax.util.io.IllegalFormatException;

/**
 * Class of static functions for reading/writing genotype matrices stored
 * in flast files
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class GenotypesFlatFile
{
    /**
     * Private constructor. This class should not be instantiated.
     */
    private GenotypesFlatFile() {}
    
    /**
     * Convert from a genotype HDF5 file to a flat file
     * @param genoCalls
     *          the genotype calls
     * @param flatFileWriter
     *          the flat file to write to
     * @throws IOException
     *          if anything goes wrong while writing
     */
    public static void writeGenoCallMatrix(
            AbstractGenotypeCallMatrix genoCalls,
            FlatFileWriter flatFileWriter) throws IOException
    {
        int colCount = 0;
        
        int snpIdIndex = -1;
        String[] snpIds = genoCalls.getSnpIds();
        if(snpIds != null)
        {
            snpIdIndex = colCount++;
        }
        
        int aAlleleIndex = -1;
        int bAlleleIndex = -1;
        char[] aAlleles = genoCalls.getAAlleles();
        char[] bAlleles = genoCalls.getBAlleles();
        if(aAlleles != null && bAlleles != null)
        {
            aAlleleIndex = colCount++;
            bAlleleIndex = colCount++;
        }
        
        int chrIndex = -1;
        String[] chrIds = genoCalls.getChrIDs();
        if(chrIds != null)
        {
            chrIndex = colCount++;
        }
        
        int posIndex = -1;
        long[] bpPos = genoCalls.getBpPositions();
        if(bpPos != null)
        {
            posIndex = colCount++;
        }
        
        int callStartIndex = colCount;
        byte[][] callMatrix = genoCalls.getCallMatrix();
        String[] sampleIds = genoCalls.getSampleIds();
        colCount += sampleIds.length;
        
        // start with the header
        String[] currRow = new String[colCount];
        if(snpIds != null)
        {
            currRow[snpIdIndex] = AbstractGenotypeCallMatrix.SNP_IDS_NAME;
        }
        
        if(aAlleles != null && bAlleles != null)
        {
            currRow[aAlleleIndex] = AbstractGenotypeCallMatrix.A_ALLELES_NAME;
            currRow[bAlleleIndex] = AbstractGenotypeCallMatrix.B_ALLELES_NAME;
        }
        
        if(chrIds != null)
        {
            currRow[chrIndex] = AbstractGenotypeCallMatrix.CHR_IDS_NAME;
        }
        
        if(bpPos != null)
        {
            currRow[posIndex] = AbstractGenotypeCallMatrix.BP_POSITIONS_NAME;
        }
        
        for(int i = 0; i < sampleIds.length; i++)
        {
            currRow[i + callStartIndex] = sampleIds[i];
        }
        
        // now fill in everything below the header
        int rowIndex = 0;
        while(currRow != null)
        {
            flatFileWriter.writeRow(currRow);
            if(rowIndex < callMatrix.length)
            {
                if(snpIds != null)
                {
                    currRow[snpIdIndex] = snpIds[rowIndex];
                }
                
                if(aAlleles != null && bAlleles != null)
                {
                    currRow[aAlleleIndex] = Character.toString(aAlleles[rowIndex]);
                    currRow[bAlleleIndex] = Character.toString(bAlleles[rowIndex]);
                }
                
                if(chrIds != null)
                {
                    currRow[chrIndex] = chrIds[rowIndex];
                }
                
                if(bpPos != null)
                {
                    currRow[posIndex] = Long.toString(bpPos[rowIndex]);
                }
                
                for(int i = 0; i < sampleIds.length; i++)
                {
                    currRow[i + callStartIndex] = Byte.toString(callMatrix[rowIndex][i]);
                }
            }
            else
            {
                currRow = null;
            }
            rowIndex++;
        }
    }
    
    private static char[] toChars(List<String> callChars)
    {
        char[] val = new char[callChars.size()];
        int i = 0;
        for(String str : callChars)
        {
            val[i] = str.charAt(0);
            i++;
        }
        return val;
    }

    /**
     * Convert a genotype CSV file to an HDF5 file
     * @param flatFileReaders
     *          the genotype files
     * @param aAlleleColumn
     *          the column index for the A allele
     * @param bAlleleColumn
     *          the column index for the B allele
     * @param snpIdColumn
     *          the SNP ID column (-1 indicates no such column)
     * @param chrColumn
     *          the chromosome ID column (-1 indicates no such column)
     * @param bpPositionColumn
     *          the BP position column (-1 indicates no such column)
     * @param bpBuildId
     *          the base pair build identifier
     * @param firstGenotypeColumn
     *          the 1st genotype column
     * @param lastGenotypeColumnExclusive
     *          the index after the last genotype column. You can use -1 to indicate that
     *          all of the remaining columns after firstGenotypeColumn are
     *          genotype columns
     * @return
     *          the geno matrix read from the flat file
     * @throws IllegalFormatException
     *          if there is a problem with how with how the file is formatted
     * @throws IOException
     *          if there is a problem with file IO while reading the flat file
     */
    public static GenotypeCallMatrix readGenoCallMatrix(
            FlatFileReader[] flatFileReaders,
            int aAlleleColumn,
            int bAlleleColumn,
            int snpIdColumn,
            int chrColumn,
            int bpPositionColumn,
            String bpBuildId,
            int firstGenotypeColumn,
            int lastGenotypeColumnExclusive) throws IllegalFormatException, IOException
    {
        // start with the geno headers
        String[] header = flatFileReaders[0].readRow();
        if(header == null)
        {
            throw new IllegalFormatException("Failed to read the header");
        }
        
        if(lastGenotypeColumnExclusive == -1)
        {
            lastGenotypeColumnExclusive = header.length;
        }
        
        int strainCount = lastGenotypeColumnExclusive - firstGenotypeColumn;
        String[] headerStrains = new String[strainCount];
        for(int i = 0; i < headerStrains.length; i++)
        {
            headerStrains[i] = header[i + firstGenotypeColumn];
        }
        
        // read the genotype data
        List<byte[]> callValues = new ArrayList<byte[]>();
        List<String> aAlleles = new ArrayList<String>();
        List<String> bAlleles = new ArrayList<String>();
        List<String> snpIds = snpIdColumn >= 0 ? new ArrayList<String>() : null;
        List<String> chrs = chrColumn >= 0 ? new ArrayList<String>() : null;
        List<Long> bpPos = bpPositionColumn >= 0 ? new ArrayList<Long>() : null;
        
        // keep track of whether all of the rows are given in sort order or not
        boolean isSorted = chrColumn >= 0 && bpPositionColumn >= 0;
        String prevChr = null;
        long prevPos = -1L;
        
        for(int i = 0; i < flatFileReaders.length; i++)
        {
            String[] currRow = null;
            if(i >= 1)
            {
                currRow = flatFileReaders[i].readRow();
                if(!Arrays.equals(currRow, header))
                {
                    throw new IllegalFormatException(
                            "All file headers must match but\n\n\"" +
                            SequenceUtilities.toString(Arrays.asList(currRow)) +
                            "\"\n\ndoes not match\n\n\"" +
                            SequenceUtilities.toString(Arrays.asList(header)) + "\"");
                }
            }
            
            while((currRow = flatFileReaders[i].readRow()) != null)
            {
                String aAllele = null;
                String bAllele = null;
                if(aAlleleColumn >= 0 && bAlleleColumn >= 0)
                {
                    aAllele = currRow[aAlleleColumn];
                    aAlleles.add(aAllele);
                    bAllele = currRow[bAlleleColumn];
                    bAlleles.add(bAllele);
                }
                
                if(snpIds != null)
                {
                    snpIds.add(currRow[snpIdColumn]);
                }
                
                String currChr = null;
                if(chrs != null)
                {
                    currChr = currRow[chrColumn];
                    chrs.add(currChr);
                }
                
                long currPos = -1L;
                if(bpPos != null)
                {
                    currPos = Long.parseLong(currRow[bpPositionColumn]);
                    bpPos.add(currPos);
                }
                
                String[] callStrings = Arrays.copyOfRange(
                        currRow,
                        firstGenotypeColumn,
                        firstGenotypeColumn + strainCount);
                callValues.add(GenotypeCallMatrix.toCallValues(
                        aAllele,
                        bAllele,
                        callStrings));
                
                if(isSorted)
                {
                    if(currChr == null || currPos < 0)
                    {
                        isSorted = false;
                    }
                    else if(prevChr != null)
                    {
                        if(prevChr.equals(currChr))
                        {
                            isSorted = prevPos <= currPos;
                        }
                        else
                        {
                            int chrComp = CallMatrixSorter.compareChrs(
                                    prevChr,
                                    currChr);
                            if(chrComp > 0)
                            {
                                isSorted = false;
                            }
                            else if(chrComp == 0)
                            {
                                isSorted = prevPos <= currPos;
                            }
                        }
                    }
                }
                
                prevChr = currChr;
                prevPos = currPos;
            }
        }
        
        // write the data
        GenotypeCallMatrix genoCalls = new GenotypeCallMatrix();
        if(aAlleleColumn >= 0 && bAlleleColumn >= 0)
        {
            genoCalls.setAAlleles(toChars(aAlleles));
            genoCalls.setBAlleles(toChars(bAlleles));
        }
        genoCalls.setSampleIds(headerStrains);
        genoCalls.setCallMatrix(callValues.toArray(new byte[callValues.size()][]));
        if(snpIds != null) genoCalls.setSnpIds(snpIds.toArray(new String[snpIds.size()]));
        if(chrs != null) genoCalls.setChrIDs(chrs.toArray(new String[chrs.size()]));
        if(bpPos != null) genoCalls.setBpPositions(SequenceUtilities.toLongArray(bpPos));
        genoCalls.setSortedByPosition(isSorted);
        
        return genoCalls;
    }
    
    /**
     * @param flatFileReader
     *          the alchemy flat file
     * @return
     *          the geno matrix read from the flat file
     * @throws IllegalFormatException
     *          if there is a problem with how with how the file is formatted
     * @throws IOException
     *          if there is a problem with file IO while reading the flat file
     */
    public static GenotypeCallMatrix readAlchemyGenoCalls(FlatFileReader flatFileReader)
    throws IllegalFormatException, IOException
    {
        final int snpIDCol = 0;
        final int sampleIDCol = 1;
        final int abCallCol = 2;
        final int expectedColCount = 14;
        
        // read in the calls and IDs
        List<String> sampleIDs = new ArrayList<String>();
        List<String> snpIDs = new ArrayList<String>();
        List<byte[]> callRows = new ArrayList<byte[]>();
        List<Byte> currCallRow = new ArrayList<Byte>();
        String prevSNPId = null;
        
        String[] currRow = null;
        int snpIndex = 0;
        while((currRow = flatFileReader.readRow()) != null)
        {
            if(currRow.length != expectedColCount)
            {
                throw new IllegalFormatException(
                        "Bad row count. Expected " + expectedColCount +
                        " columns but there were " + currRow.length);
            }
            
            String currSnpID = currRow[snpIDCol];
            String currSampleID = currRow[sampleIDCol];
            String currABCall = currRow[abCallCol];
            
            if(!currSnpID.equals(prevSNPId))
            {
                if(prevSNPId != null)
                {
                    callRows.add(SequenceUtilities.toByteArray(currCallRow));
                    currCallRow.clear();
                    
                    snpIndex++;
                    //if(snpIndex % 1000 == 0)
                    //{
                    //    System.out.println("Processing SNP " + snpIndex);
                    //}
                }
                
                snpIDs.add(currSnpID);
                prevSNPId = currSnpID;
            }
            
            currCallRow.add(alchemyCallToByteCall(currABCall));
            if(snpIndex == 0)
            {
                sampleIDs.add(currSampleID);
            }
        }
        
        if(prevSNPId == null)
        {
            throw new IllegalFormatException(
                    "Failed to convert: the alchemy file appears to be empty");
        }
        
        // pick up the last row that didn't get added in our while loop
        callRows.add(SequenceUtilities.toByteArray(currCallRow));
        
        GenotypeCallMatrix callMat = new GenotypeCallMatrix();
        callMat.setCallMatrix(callRows.toArray(new byte[callRows.size()][]));
        callMat.setSampleIds(sampleIDs.toArray(new String[sampleIDs.size()]));
        callMat.setSnpIds(snpIDs.toArray(new String[snpIDs.size()]));
        return callMat;
    }

    private static byte alchemyCallToByteCall(String abCall)
    {
        if(abCall.equals("AA"))
        {
            return 1;
        }
        else if(abCall.equals("BB"))
        {
            return 2;
        }
        else if(abCall.equals("AB"))
        {
            return 3;
        }
        else
        {
            throw new IllegalArgumentException("Unexpected AB call value: " + abCall);
        }
    }
}
