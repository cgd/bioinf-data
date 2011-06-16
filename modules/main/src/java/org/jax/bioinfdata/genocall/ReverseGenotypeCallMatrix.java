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

import java.util.Arrays;
import java.util.Collections;

/**
 * Provides a reversed view of a given {@link AbstractGenotypeCallMatrix}.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class ReverseGenotypeCallMatrix extends AbstractGenotypeCallMatrix
{
    private final AbstractGenotypeCallMatrix originalCallMatrix;
    private final long snpCount;
    
    /**
     * Constructor
     * @param originalCallMatrix the matrix to reverse
     */
    public ReverseGenotypeCallMatrix(AbstractGenotypeCallMatrix originalCallMatrix)
    {
        this.originalCallMatrix = originalCallMatrix;
        this.snpCount = originalCallMatrix.getSNPCount();
    }
    
    private static char[] reverseCharArray(char[] arr)
    {
        if(arr == null)
        {
            return arr;
        }
        else
        {
            arr = arr.clone();
            int halfLen = arr.length / 2;
            for(int i = 0; i < halfLen; i++)
            {
                int j = arr.length - (i + 1);
                char tmp = arr[i];
                arr[i] = arr[j];
                arr[j] = tmp;
            }
            return arr;
        }
    }

    private static long[] reverseLongArray(long[] arr)
    {
        if(arr == null)
        {
            return arr;
        }
        else
        {
            arr = arr.clone();
            int halfLen = arr.length / 2;
            for(int i = 0; i < halfLen; i++)
            {
                int j = arr.length - (i + 1);
                long tmp = arr[i];
                arr[i] = arr[j];
                arr[j] = tmp;
            }
            return arr;
        }
    }
    
    private static <T> T[] reverseObjectArray(T[] arr)
    {
        if(arr == null)
        {
            return arr;
        }
        else
        {
            arr = arr.clone();
            Collections.reverse(Arrays.asList(arr));
            return arr;
        }
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public char[] getAAlleles()
    {
        return reverseCharArray(this.originalCallMatrix.getAAlleles());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setAAlleles(char[] aAlleles)
    {
        this.originalCallMatrix.setAAlleles(reverseCharArray(aAlleles));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public char[] getBAlleles()
    {
        return reverseCharArray(this.originalCallMatrix.getBAlleles());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setBAlleles(char[] bAlleles)
    {
        this.originalCallMatrix.setBAlleles(reverseCharArray(bAlleles));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String[] getSampleIds()
    {
        return this.originalCallMatrix.getSampleIds();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setSampleIds(String[] sampleIds)
    {
        this.originalCallMatrix.setSampleIds(sampleIds);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public byte[][] getCallMatrix()
    {
        return reverseObjectArray(this.getCallMatrix());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setCallMatrix(byte[][] callMatrix)
    {
        this.originalCallMatrix.setCallMatrix(reverseObjectArray(callMatrix));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public byte[] getSNPCalls(long snpIndex)
    {
        return this.originalCallMatrix.getSNPCalls(this.snpCount - (snpIndex + 1));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public long getSNPCount()
    {
        return this.snpCount;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String[] getSnpIds()
    {
        return reverseObjectArray(this.originalCallMatrix.getSnpIds());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setSnpIds(String[] snpIds)
    {
        this.originalCallMatrix.setSnpIds(reverseObjectArray(snpIds));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String[] getChrIDs()
    {
        return reverseObjectArray(this.originalCallMatrix.getChrIDs());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setChrIDs(String[] chrIDs)
    {
        this.originalCallMatrix.setChrIDs(reverseObjectArray(chrIDs));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public long[] getBpPositions()
    {
        return reverseLongArray(this.originalCallMatrix.getBpPositions());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String getBuildId()
    {
        return this.originalCallMatrix.getBuildId();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setBpPositions(long[] bpPositions)
    {
        this.originalCallMatrix.setBpPositions(reverseLongArray(bpPositions));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setBpPositions(long[] bpPositions, String buildId)
    {
        this.originalCallMatrix.setBpPositions(
                reverseLongArray(bpPositions),
                buildId);
    }
}
