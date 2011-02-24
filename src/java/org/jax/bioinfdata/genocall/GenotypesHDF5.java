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

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;

/**
 * For going to and from HDF5 format for genotype data
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class GenotypesHDF5
{
    public static final String A_ALLELES_PATH = "aAllele";
    public static final String B_ALLELES_PATH = "bAllele";
    public static final String SAMPLE_IDS_PATH = "sampleID";
    public static final String CALL_MATRIX_PATH = "callMatrix";
    public static final String SNP_IDS_PATH = "snpID";
    public static final String CHR_IDS_PATH = "chrID";
    public static final String BP_POSITIONS_PATH = "bpPosition";
    public static final String BUILD_ID_ATTR = "buildID";
    
    private char[] toChars(String[] callChars)
    {
        char[] val = new char[callChars.length];
        for(int i = 0; i < val.length; i++)
        {
            val[i] = callChars[i].charAt(0);
        }
        return val;
    }
    
    private String[] toStrings(char[] callChars)
    {
        String[] val = new String[callChars.length];
        for(int i = 0; i < val.length; i++)
        {
            val[i] = Character.toString(callChars[i]);
        }
        return val;
    }
    
    /**
     * Convert from a genotype HDF5 file to a flat file
     * @param hdf5Reader the HDF5 file
     * @return the call matrix
     */
    public GenotypeCallMatrix readGenoCallMatrix(IHDF5Reader hdf5Reader)
    {
        GenotypeCallMatrix genoCalls = new GenotypeCallMatrix();
        
        if(hdf5Reader.exists(SNP_IDS_PATH))
        {
            genoCalls.setSnpIds(hdf5Reader.readStringArray(SNP_IDS_PATH));
        }
        
        genoCalls.setAAlleles(toChars(hdf5Reader.readStringArray(A_ALLELES_PATH)));
        genoCalls.setBAlleles(toChars(hdf5Reader.readStringArray(B_ALLELES_PATH)));
        
        if(hdf5Reader.exists(CHR_IDS_PATH))
        {
            genoCalls.setChrIDs(hdf5Reader.readStringArray(CHR_IDS_PATH));
        }
        
        if(hdf5Reader.exists(BP_POSITIONS_PATH))
        {
            genoCalls.setBpPositions(hdf5Reader.readLongArray(BP_POSITIONS_PATH));
            if(hdf5Reader.hasAttribute(BP_POSITIONS_PATH, BUILD_ID_ATTR))
            {
                genoCalls.setBuildId(hdf5Reader.getStringAttribute(BP_POSITIONS_PATH, BUILD_ID_ATTR));
            }
        }
        
        genoCalls.setCallMatrix(hdf5Reader.readByteMatrix(CALL_MATRIX_PATH));
        genoCalls.setSampleIds(hdf5Reader.readStringArray(SAMPLE_IDS_PATH));
        
        return genoCalls;
    }

    /**
     * Write the genotype call matrix to the given HDF5 writer
     * @param genoCalls
     *          the calls to write
     * @param hdf5Writer
     *          the writer
     */
    public void writeGenoCallMatrix(GenotypeCallMatrix genoCalls, IHDF5Writer hdf5Writer)
    {
        // write the data
        hdf5Writer.writeStringArray(A_ALLELES_PATH, toStrings(genoCalls.getAAlleles()));
        hdf5Writer.writeStringArray(B_ALLELES_PATH, toStrings(genoCalls.getBAlleles()));
        hdf5Writer.writeStringArray(SAMPLE_IDS_PATH, genoCalls.getSampleIds());
        hdf5Writer.writeByteMatrix(CALL_MATRIX_PATH, genoCalls.getCallMatrix());
        if(genoCalls.getSnpIds() != null)
        {
            hdf5Writer.writeStringArray(SNP_IDS_PATH, genoCalls.getSnpIds());
        }
        if(genoCalls.getChrIDs() != null)
        {
            hdf5Writer.writeStringArray(CHR_IDS_PATH, genoCalls.getChrIDs());
        }
        if(genoCalls.getBpPositions() != null)
        {
            hdf5Writer.writeLongArray(BP_POSITIONS_PATH, genoCalls.getBpPositions());
            if(genoCalls.getBuildId() != null)
            {
                hdf5Writer.setStringAttribute(BP_POSITIONS_PATH, BUILD_ID_ATTR, genoCalls.getBuildId());
            }
        }
    }

    /**
     * Convert from a genotype HDF5 file to a flat file
     * @param hdf5Reader
     *          the HDF5 file
     * @param flatFileWriter
     *          the flat file to write to
     * @throws IOException
     *          if anything goes wrong while reading or writing
     */
    public void genoHDF5ToCSV(
            IHDF5Reader hdf5Reader,
            FlatFileWriter flatFileWriter) throws IOException
    {
        int colCount = 0;
        
        int snpIdIndex = -1;
        String[] snpIds = null;
        if(hdf5Reader.exists(SNP_IDS_PATH))
        {
            snpIdIndex = colCount++;
            snpIds = hdf5Reader.readStringArray(SNP_IDS_PATH);
        }
        
        int aAlleleIndex = colCount++;
        String[] aAlleles = hdf5Reader.readStringArray(A_ALLELES_PATH);
        int bAlleleIndex = colCount++;
        String[] bAlleles = hdf5Reader.readStringArray(B_ALLELES_PATH);
        
        int chrIndex = -1;
        String[] chrIds = null;
        if(hdf5Reader.exists(CHR_IDS_PATH))
        {
            chrIndex = colCount++;
            chrIds = hdf5Reader.readStringArray(CHR_IDS_PATH);
        }
        
        int posIndex = -1;
        long[] bpPos = null;
        if(hdf5Reader.exists(BP_POSITIONS_PATH))
        {
            posIndex = colCount++;
            bpPos = hdf5Reader.readLongArray(BP_POSITIONS_PATH);
        }
        
        int callStartIndex = colCount;
        byte[][] callMatrix = hdf5Reader.readByteMatrix(CALL_MATRIX_PATH);
        String[] sampleIds = hdf5Reader.readStringArray(SAMPLE_IDS_PATH);
        colCount += sampleIds.length;
        
        // start with the header
        String[] currRow = new String[colCount];
        if(snpIds != null)
        {
            currRow[snpIdIndex] = SNP_IDS_PATH;
        }
        
        currRow[aAlleleIndex] = A_ALLELES_PATH;
        currRow[bAlleleIndex] = B_ALLELES_PATH;
        
        if(chrIds != null)
        {
            currRow[chrIndex] = CHR_IDS_PATH;
        }
        
        if(bpPos != null)
        {
            currRow[posIndex] = BP_POSITIONS_PATH;
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
                
                currRow[aAlleleIndex] = aAlleles[rowIndex];
                currRow[bAlleleIndex] = bAlleles[rowIndex];
                
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
     * @param hdf5Writer
     *          the file that we should write HDF5 data to
     * @throws IllegalFormatException
     *          if there is a problem with how with how the file is formatted
     * @throws IOException
     *          if there is a problem with file IO while reading the flat file
     */
    public void genoFlatFileToHDF5(
            FlatFileReader flatFileReader,
            int aAlleleColumn,
            int bAlleleColumn,
            int snpIdColumn,
            int chrColumn,
            int bpPositionColumn,
            String bpBuildId,
            int firstGenotypeColumn,
            int lastGenotypeColumnExclusive,
            IHDF5Writer hdf5Writer) throws IllegalFormatException, IOException
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
        hdf5Writer.writeStringArray(
                A_ALLELES_PATH,
                aAlleles.toArray(new String[aAlleles.size()]));
        hdf5Writer.writeStringArray(
                B_ALLELES_PATH,
                bAlleles.toArray(new String[bAlleles.size()]));
        hdf5Writer.writeStringArray(
                SAMPLE_IDS_PATH,
                headerStrains);
        hdf5Writer.writeByteMatrix(
                CALL_MATRIX_PATH,
                callValues.toArray(new byte[callValues.size()][]));
        if(snpIdValid)
        {
            hdf5Writer.writeStringArray(
                    SNP_IDS_PATH,
                    snpIds.toArray(new String[snpIds.size()]));
        }
        if(chrValid)
        {
            hdf5Writer.writeStringArray(
                    CHR_IDS_PATH,
                    chrs.toArray(new String[chrs.size()]));
        }
        if(bpPosValid)
        {
            hdf5Writer.writeLongArray(
                    BP_POSITIONS_PATH,
                    SequenceUtilities.toLongArray(bpPos));
            if(bpBuildId != null)
            {
                hdf5Writer.setStringAttribute(BP_POSITIONS_PATH, BUILD_ID_ATTR, bpBuildId);
            }
        }
    }
}
