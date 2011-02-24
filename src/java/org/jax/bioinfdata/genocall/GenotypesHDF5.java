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
}
