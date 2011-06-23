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
 * A genotype call matrix backed by an HDF5 file
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class HDF5GenotypeCallMatrix extends AbstractGenotypeCallMatrix
{
    private final IHDF5Reader hdf5Reader;
    private final IHDF5Writer hdf5Writer;
    private HDF5MatrixRowColumnReader hdf5RowReader;
    
    /**
     * Constructor
     * @param hdf5Reader the reader
     * @param hdf5Writer the writer
     */
    public HDF5GenotypeCallMatrix(IHDF5Reader hdf5Reader, IHDF5Writer hdf5Writer)
    {
        this.hdf5Reader = hdf5Reader;
        this.hdf5Writer = hdf5Writer;
    }
    
//    private static String pathCombine(String... pathElements)
//    {
//        StringBuilder sb = new StringBuilder();
//        for(int i = 0; i < pathElements.length; i++)
//        {
//            if(i >= 1)
//            {
//                if(!(pathElements[i - 1].endsWith("/") || pathElements[i].startsWith("/")))
//                {
//                    sb.append('/');
//                }
//            }
//            
//            sb.append(pathElements[i]);
//        }
//        
//        return sb.toString();
//    }
    
    /**
     * Constructor
     * @param hdf5Reader the reader
     */
    public HDF5GenotypeCallMatrix(IHDF5Reader hdf5Reader)
    {
        this(hdf5Reader, null);
    }
    
    /**
     * Constructor
     * @param hdf5Writer the writer
     */
    public HDF5GenotypeCallMatrix(IHDF5Writer hdf5Writer)
    {
        this(null, hdf5Writer);
    }
    
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
     * {@inheritDoc}
     */
    @Override
    public char[] getAAlleles()
    {
        if(this.hdf5Reader.exists(A_ALLELES_NAME))
        {
            return toChars(this.hdf5Reader.readStringArray(A_ALLELES_NAME));
        }
        else
        {
            return null;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setAAlleles(char[] aAlleles)
    {
        if(aAlleles != null)
        {
            this.hdf5Writer.writeStringArray(A_ALLELES_NAME, toStrings(aAlleles));
        }
        else if(this.hdf5Writer.exists(A_ALLELES_NAME))
        {
            this.hdf5Writer.delete(A_ALLELES_NAME);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public char[] getBAlleles()
    {
        if(this.hdf5Reader.exists(B_ALLELES_NAME))
        {
            return toChars(this.hdf5Reader.readStringArray(B_ALLELES_NAME));
        }
        else
        {
            return null;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setBAlleles(char[] bAlleles)
    {
        if(bAlleles != null)
        {
            this.hdf5Writer.writeStringArray(B_ALLELES_NAME, toStrings(bAlleles));
        }
        else if(this.hdf5Writer.exists(B_ALLELES_NAME))
        {
            this.hdf5Writer.delete(B_ALLELES_NAME);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String[] getSampleIds()
    {
        if(this.hdf5Reader.exists(SAMPLE_IDS_NAME))
        {
            return this.hdf5Reader.readStringArray(SAMPLE_IDS_NAME);
        }
        else
        {
            return null;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setSampleIds(String[] sampleIds)
    {
        if(sampleIds != null)
        {
            this.hdf5Writer.writeStringArray(SAMPLE_IDS_NAME, sampleIds);
        }
        else if(this.hdf5Writer.exists(SAMPLE_IDS_NAME))
        {
            this.hdf5Writer.delete(SAMPLE_IDS_NAME);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public byte[][] getCallMatrix()
    {
        if(this.hdf5Reader.exists(CALL_MATRIX_NAME))
        {
            return this.hdf5Reader.readByteMatrix(CALL_MATRIX_NAME);
        }
        else
        {
            return null;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setCallMatrix(byte[][] callMatrix)
    {
        if(callMatrix != null)
        {
            this.hdf5Writer.writeByteMatrix(CALL_MATRIX_NAME, callMatrix);
        }
        else if(this.hdf5Writer.exists(CALL_MATRIX_NAME))
        {
            this.hdf5Writer.delete(CALL_MATRIX_NAME);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public synchronized byte[] getSNPCalls(long snpIndex)
    {
        if(this.hdf5RowReader == null)
        {
            this.hdf5RowReader = new HDF5MatrixRowColumnReader(
                    this.hdf5Reader,
                    CALL_MATRIX_NAME);
        }
        
        return this.hdf5RowReader.getRowBytes(snpIndex);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public synchronized long getSNPCount()
    {
        if(this.hdf5RowReader == null)
        {
            this.hdf5RowReader = new HDF5MatrixRowColumnReader(
                    this.hdf5Reader,
                    CALL_MATRIX_NAME);
        }
        
        return this.hdf5RowReader.getRowCount();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String[] getSnpIds()
    {
        if(this.hdf5Reader.exists(SNP_IDS_NAME))
        {
            return this.hdf5Reader.readStringArray(SNP_IDS_NAME);
        }
        else
        {
            return null;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setSnpIds(String[] snpIds)
    {
        if(snpIds != null)
        {
            this.hdf5Writer.writeStringArray(SNP_IDS_NAME, snpIds);
        }
        else if(this.hdf5Writer.exists(SNP_IDS_NAME))
        {
            this.hdf5Writer.delete(SNP_IDS_NAME);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String[] getChrIDs()
    {
        if(this.hdf5Reader.exists(CHR_IDS_NAME))
        {
            return this.hdf5Reader.readStringArray(CHR_IDS_NAME);
        }
        else
        {
            return null;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setChrIDs(String[] chrIDs)
    {
        if(chrIDs != null)
        {
            this.hdf5Writer.writeStringArray(CHR_IDS_NAME, chrIDs);
        }
        else if(this.hdf5Writer.exists(CHR_IDS_NAME))
        {
            this.hdf5Writer.delete(CHR_IDS_NAME);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public long[] getBpPositions()
    {
        if(this.hdf5Reader.exists(BP_POSITIONS_NAME))
        {
            return this.hdf5Reader.readLongArray(BP_POSITIONS_NAME);
        }
        else
        {
            return null;
        }
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void setBpPositions(long[] bpPositions)
    {
        if(bpPositions != null)
        {
            this.hdf5Writer.writeLongArray(BP_POSITIONS_NAME, bpPositions);
        }
        else if(this.hdf5Writer.exists(BP_POSITIONS_NAME))
        {
            this.hdf5Writer.delete(BP_POSITIONS_NAME);
        }
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public String getBuildId()
    {
        if(this.hdf5Reader.exists(BUILD_ID_NAME))
        {
            return this.hdf5Reader.readString(BUILD_ID_NAME);
        }
        else
        {
            return null;
        }
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void setBuildId(String buildId)
    {
        if(buildId != null)
        {
            this.hdf5Writer.writeString(BUILD_ID_NAME, buildId);
        }
        else if(this.hdf5Writer.exists(BUILD_ID_NAME))
        {
            this.hdf5Writer.delete(BUILD_ID_NAME);
        }
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public boolean isSortedByPosition()
    {
        if(this.hdf5Reader.exists(SORTED_BY_POS_NAME))
        {
            return this.hdf5Reader.readBoolean(SORTED_BY_POS_NAME);
        }
        else
        {
            return false;
        }
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void setSortedByPosition(boolean sortedByPosition)
    {
        this.hdf5Writer.writeBoolean(SORTED_BY_POS_NAME, sortedByPosition);
    }
}
