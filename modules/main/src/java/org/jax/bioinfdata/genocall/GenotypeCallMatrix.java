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

/**
 * An object for holding data describing a genotype call matrix. Most uses will
 * require that both sampleIds and callMatrix are set to a non-null value.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class GenotypeCallMatrix extends AbstractGenotypeCallMatrix
{
    private char[] aAlleles;
    private char[] bAlleles;
    private String[] sampleIds;
    private byte[][] callMatrix;
    private String[] snpIds;
    private String[] chrIDs;
    private long[] bpPositions;
    private String buildId;
    private boolean sortedByPosition;
    
    /**
     * {@inheritDoc}
     */
    @Override
    public char[] getAAlleles()
    {
        return this.aAlleles;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void setAAlleles(char[] aAlleles)
    {
        this.aAlleles = aAlleles;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public char[] getBAlleles()
    {
        return this.bAlleles;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void setBAlleles(char[] bAlleles)
    {
        this.bAlleles = bAlleles;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public String[] getSampleIds()
    {
        return this.sampleIds;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void setSampleIds(String[] sampleIds)
    {
        this.sampleIds = sampleIds;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public byte[][] getCallMatrix()
    {
        return this.callMatrix;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void setCallMatrix(byte[][] callMatrix)
    {
        this.callMatrix = callMatrix;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public long getSNPCount()
    {
        return this.callMatrix.length;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public String[] getSnpIds()
    {
        return this.snpIds;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void setSnpIds(String[] snpIds)
    {
        this.snpIds = snpIds;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public String[] getChrIDs()
    {
        return this.chrIDs;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void setChrIDs(String[] chrIDs)
    {
        this.chrIDs = chrIDs;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public long[] getBpPositions()
    {
        return this.bpPositions;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void setBpPositions(long[] bpPositions)
    {
        this.bpPositions = bpPositions;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public String getBuildId()
    {
        return this.buildId;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void setBuildId(String buildId)
    {
        this.buildId = buildId;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void setSortedByPosition(boolean sortedByPosition)
    {
        this.sortedByPosition = sortedByPosition;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public boolean getSortedByPosition()
    {
        return this.sortedByPosition;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public byte[] getSNPCalls(long snpIndex)
    {
        return this.callMatrix[(int)snpIndex];
    }
}
