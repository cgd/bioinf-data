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

import org.jax.bioinfdata.genocall.CallMatrixSorter;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class GenomeInterval implements Comparable<GenomeInterval>
{
    private final String chrId;
    private final long startPositionBp;
    private final long endPositionBp;
    
    /**
     * Constructor
     * @param chrId             chromosome ID
     * @param startPositionBp   the starting position in base pairs
     * @param endPositionBp     the end position in base pairs
     */
    public GenomeInterval(String chrId, long startPositionBp, long endPositionBp)
    {
        this.chrId = chrId;
        this.startPositionBp = startPositionBp;
        this.endPositionBp = endPositionBp;
    }
    
    /**
     * Getter for the chromosome ID
     * @return the chrId
     */
    public String getChrId()
    {
        return this.chrId;
    }
    
    /**
     * Getter for the start position
     * @return the startPositionBp
     */
    public long getStartPositionBp()
    {
        return this.startPositionBp;
    }
    
    /**
     * Getter for the end position
     * @return the endPositionBp
     */
    public long getEndPositionBp()
    {
        return this.endPositionBp;
    }

    /**
     * {@inheritDoc}
     */
    public int compareTo(GenomeInterval other)
    {
        int comp = CallMatrixSorter.compareChrs(this.chrId, other.chrId);
        if(comp != 0)
        {
            return comp;
        }
        else
        {
            comp = (int)(this.startPositionBp - other.endPositionBp);
            if(comp != 0)
            {
                return comp;
            }
            else
            {
                return (int)(this.endPositionBp - other.endPositionBp);
            }
        }
    }
}
