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

package org.jax.bioinfdata.phylogeny;

import java.io.Serializable;

/**
 * Interval based on indices (rather than physical position). For these blocks
 * to make any sense, you need to know which group of SNPs you're indexing.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class IndexedSnpInterval implements Comparable<IndexedSnpInterval>, Serializable
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = 9076158079661611167L;

    private final int startIndex;
    
    private final int extentInIndices;

    /**
     * Constructor
     * @param startIndex
     *          see {@link #getStartIndex()}
     * @param extentInIndices
     *          see {@link #getExtentInIndices()}
     */
    public IndexedSnpInterval(
            int startIndex,
            int extentInIndices)
    {
        if(startIndex < 0)
        {
            throw new IndexOutOfBoundsException(
                    "Start index should be >= 0: " + startIndex);
        }
        
        if(extentInIndices < 0)
        {
            throw new IndexOutOfBoundsException(
                    "Extent should be >= 0: " + extentInIndices);
        }
        
        this.startIndex = startIndex;
        this.extentInIndices = extentInIndices;
    }
    
    /**
     * Getter for the starting index inclusive
     * @return
     *          the startIndex
     */
    public int getStartIndex()
    {
        return this.startIndex;
    }
    
    /**
     * Getter for the extent in indices (not base pairs!)
     * @return
     *          the extent in indices
     */
    public int getExtentInIndices()
    {
        return this.extentInIndices;
    }
    
    /**
     * Getter for the end index inclusive
     * @return
     *          the end index
     */
    public int getEndIndex()
    {
        return this.startIndex + this.extentInIndices - 1;
    }
    
    /**
     * Determine if this interval intersects the given interval
     * @param otherInterval
     *          the other interval we're checking for intersection
     * @return
     *          true if they intersect
     */
    public boolean intersects(IndexedSnpInterval otherInterval)
    {
        return this.startIndex <= otherInterval.getEndIndex() &&
               this.getEndIndex() >= otherInterval.startIndex;
    }
    
    /**
     * Determine if the this interval fully contains the given interval
     * @param otherInterval
     *          the interval we're comparing against
     * @return
     *          true iff this interval fully contains the other
     */
    public boolean contains(IndexedSnpInterval otherInterval)
    {
        return this.startIndex <= otherInterval.startIndex &&
               this.getEndIndex() >= otherInterval.getEndIndex();
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public String toString()
    {
        return "Start Index: " + this.startIndex +
               ", End Index: " + this.getEndIndex();
    }
    
    /**
     * {@inheritDoc}
     */
    public int compareTo(IndexedSnpInterval otherIndexedInterval)
    {
        int startDifference =
            this.startIndex - otherIndexedInterval.startIndex;
        if(startDifference != 0)
        {
            return startDifference;
        }
        else
        {
            return this.extentInIndices - otherIndexedInterval.extentInIndices;
        }
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(Object otherObject)
    {
        if(otherObject instanceof IndexedSnpInterval)
        {
            IndexedSnpInterval otherIndexedSnpInterval =
                (IndexedSnpInterval)otherObject;
            return this.startIndex == otherIndexedSnpInterval.startIndex &&
                   this.extentInIndices == otherIndexedSnpInterval.extentInIndices;
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
    public int hashCode()
    {
        return this.startIndex | (this.extentInIndices << 16);
    }
}
