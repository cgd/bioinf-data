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

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;

import org.jax.util.datastructure.SetUtilities;

/**
 * A recursive structure representing the SDP inclusion hierarchy for a
 * perfect phylogeny.
 * @author Keith Sheppard
 */
/*package-protected*/ class SdpInclusionHierarchy
{
    /**
     * A comparator that only looks at the SDP bits. It's capable of taking
     * an {@link SdpInclusionHierarchy} and/or a {@link BitSet} as input
     * @see SetUtilities#BIT_SET_COMPARATOR
     */
    public static final Comparator<Object> SDP_ONLY_COMPARATOR =
        new Comparator<Object>()
        {
            /**
             * {@inheritDoc}
             */
            public int compare(Object object1, Object object2)
            {
                final BitSet sdpBits1;
                if(object1 instanceof BitSet)
                {
                    sdpBits1 = (BitSet)object1;
                }
                else if(object1 instanceof SdpInclusionHierarchy)
                {
                    SdpInclusionHierarchy inclusionHierarchy1 =
                        (SdpInclusionHierarchy)object1;
                    sdpBits1 = inclusionHierarchy1.getSdpBits();
                }
                else
                {
                    throw new IllegalArgumentException(
                            "argument 1 must be a BitSet or an " +
                            "SdpInclusionHierarchy");
                }
                
                final BitSet sdpBits2;
                if(object2 instanceof BitSet)
                {
                    sdpBits2 = (BitSet)object2;
                }
                else if(object2 instanceof SdpInclusionHierarchy)
                {
                    SdpInclusionHierarchy inclusionHierarchy2 =
                        (SdpInclusionHierarchy)object2;
                    sdpBits2 = inclusionHierarchy2.getSdpBits();
                }
                else
                {
                    throw new IllegalArgumentException(
                            "argument 2 must be a BitSet or an " +
                            "SdpInclusionHierarchy");
                }
                
                return SetUtilities.BIT_SET_COMPARATOR.compare(
                        sdpBits1,
                        sdpBits2);
            }
        };
    
    private final BitSet sdpBits;
    
    private final List<SdpInclusionHierarchy> children;

    /**
     * Constructor
     * @param sdpBits
     *          {@link #getSdpBits()}
     */
    public SdpInclusionHierarchy(BitSet sdpBits)
    {
        this(sdpBits, new ArrayList<SdpInclusionHierarchy>());
    }
    
    /**
     * Constructor
     * @param sdpBits
     *          see {@link #getSdpBits()}
     * @param children
     *          see {@link #getChildren()}
     */
    public SdpInclusionHierarchy(
            BitSet sdpBits,
            List<SdpInclusionHierarchy> children)
    {
        this.sdpBits = sdpBits;
        this.children = children;
    }
    
    /**
     * Getter for the SDP bits
     * @return the sdpBits
     */
    public BitSet getSdpBits()
    {
        return this.sdpBits;
    }
    
    /**
     * Getter for the children
     * @return the children
     */
    public List<SdpInclusionHierarchy> getChildren()
    {
        return this.children;
    }
}
