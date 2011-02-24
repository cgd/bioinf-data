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
import java.util.List;

import org.jax.util.datastructure.SequenceUtilities;
import org.jax.util.io.FlatFileReader;
import org.jax.util.io.FlatFileWriter;
import org.jax.util.io.IllegalFormatException;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class GenotypesFlatFile
{
    /**
     * Convert from a genotype HDF5 file to a flat file
     * @param genoCalls
     *          the genotype calls
     * @param flatFileWriter
     *          the flat file to write to
     * @throws IOException
     *          if anything goes wrong while writing
     */
    public void writeGenoCallMatrix(
            GenotypeCallMatrix genoCalls,
            FlatFileWriter flatFileWriter) throws IOException
    {
        int colCount = 0;
        
        int snpIdIndex = -1;
        String[] snpIds = genoCalls.getSnpIds();
        if(snpIds != null)
        {
            snpIdIndex = colCount++;
        }
        
        int aAlleleIndex = colCount++;
        char[] aAlleles = genoCalls.getAAlleles();
        int bAlleleIndex = colCount++;
        char[] bAlleles = genoCalls.getBAlleles();
        
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
            currRow[snpIdIndex] = GenotypesHDF5.SNP_IDS_PATH;
        }
        
        currRow[aAlleleIndex] = GenotypesHDF5.A_ALLELES_PATH;
        currRow[bAlleleIndex] = GenotypesHDF5.B_ALLELES_PATH;
        
        if(chrIds != null)
        {
            currRow[chrIndex] = GenotypesHDF5.CHR_IDS_PATH;
        }
        
        if(bpPos != null)
        {
            currRow[posIndex] = GenotypesHDF5.BP_POSITIONS_PATH;
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
                
                currRow[aAlleleIndex] = Character.toString(aAlleles[rowIndex]);
                currRow[bAlleleIndex] = Character.toString(bAlleles[rowIndex]);
                
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
    
    private char[] toChars(List<String> callChars)
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
     * @param flatFileReader
     *          the genotype file
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
    public GenotypeCallMatrix readGenoCallMatrix(
            FlatFileReader flatFileReader,
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
        String[] currRow = flatFileReader.readRow();
        if(currRow == null)
        {
            throw new IllegalFormatException("Failed to read the header");
        }
        
        if(lastGenotypeColumnExclusive == -1)
        {
            lastGenotypeColumnExclusive = currRow.length;
        }
        
        int strainCount = lastGenotypeColumnExclusive - firstGenotypeColumn;
        String[] headerStrains = new String[strainCount];
        for(int i = 0; i < headerStrains.length; i++)
        {
            headerStrains[i] = currRow[i + firstGenotypeColumn];
        }
        
        // read the genotype data
        List<byte[]> callValues = new ArrayList<byte[]>();
        List<String> aAlleles = new ArrayList<String>();
        List<String> bAlleles = new ArrayList<String>();
        boolean snpIdValid = snpIdColumn >= 0;
        List<String> snpIds = snpIdValid ? new ArrayList<String>() : null;
        boolean chrValid = chrColumn >= 0;
        List<String> chrs = chrValid ? new ArrayList<String>() : null;
        boolean bpPosValid = bpPositionColumn >= 0;
        List<Long> bpPos = bpPosValid ? new ArrayList<Long>() : null;
        
        while((currRow = flatFileReader.readRow()) != null)
        {
            String aAllele = currRow[aAlleleColumn];
            aAlleles.add(aAllele);
            String bAllele = currRow[bAlleleColumn];
            bAlleles.add(bAllele);
            
            if(snpIdValid)
            {
                snpIds.add(currRow[snpIdColumn]);
            }
            
            if(chrValid)
            {
                chrs.add(currRow[chrColumn]);
            }
            
            if(bpPosValid)
            {
                bpPos.add(Long.parseLong(currRow[bpPositionColumn]));
            }
            
            byte[] currSnpGenos = new byte[strainCount];
            for(int strainIndex = 0; strainIndex < strainCount; strainIndex++)
            {
                currSnpGenos[strainIndex] = GenotypeCallMatrix.toCallValue(
                        aAllele,
                        bAllele,
                        currRow[strainIndex + firstGenotypeColumn]);
            }
            callValues.add(currSnpGenos);
        }
        
        // write the data
        GenotypeCallMatrix genoCalls = new GenotypeCallMatrix();
        genoCalls.setAAlleles(toChars(aAlleles));
        genoCalls.setBAlleles(toChars(bAlleles));
        genoCalls.setSampleIds(headerStrains);
        genoCalls.setCallMatrix(callValues.toArray(new byte[callValues.size()][]));
        genoCalls.setSnpIds(snpIds.toArray(new String[snpIds.size()]));
        genoCalls.setChrIDs(chrs.toArray(new String[chrs.size()]));
        genoCalls.setBpPositions(SequenceUtilities.toLongArray(bpPos));
        
        return genoCalls;
    }
}
