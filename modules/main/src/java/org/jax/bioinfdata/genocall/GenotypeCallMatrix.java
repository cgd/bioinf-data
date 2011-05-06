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

import org.jax.util.io.IllegalFormatException;

/**
 * An object for holding data describing a genotype call matrix. Most uses will
 * require that both sampleIds and callMatrix are set to a non-null value.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class GenotypeCallMatrix
{
    public enum AlleleCallCode
    {
        ACall
        {
            @Override
            public byte getNumCode()
            {
                return 1;
            }
        },
        
        BCall
        {
            @Override
            public byte getNumCode()
            {
                return 2;
            }
        },
        
        HCall
        {
            @Override
            public byte getNumCode()
            {
                return 3;
            }
        },
        
        NCall
        {
            @Override
            public byte getNumCode()
            {
                return -1;
            }
        };
        
        public static AlleleCallCode numCodeToEnum(byte numCode)
        {
            switch(numCode)
            {
                case 1: return ACall;
                case 2: return BCall;
                case 3: return HCall;
                case -1: return NCall;
                default: throw new IllegalArgumentException(
                        "invalid numCode: " + numCode);
            }
        }
        
        public abstract byte getNumCode();
    }
    
    public static byte toCallValue(String aAllele, String bAllele, String genoCall) throws IllegalFormatException
    {
        if(genoCall.equals("NA"))
        {
            return AlleleCallCode.NCall.getNumCode();
        }
        else if(aAllele == null || bAllele == null)
        {
            return Byte.parseByte(genoCall.trim());
        }
        else
        {
            return toCallCode(aAllele, bAllele, genoCall).getNumCode();
        }
    }
    
    public static AlleleCallCode toCallCode(String aAllele, String bAllele, String genoCall) throws IllegalFormatException
    {
        aAllele = aAllele.toUpperCase();
        bAllele = bAllele.toUpperCase();
        genoCall = genoCall.toUpperCase();
        if(genoCall.equals(aAllele) ||
           genoCall.equals(Byte.toString(AlleleCallCode.ACall.getNumCode())))
        {
            return AlleleCallCode.ACall;
        }
        else if(genoCall.equals(bAllele) ||
                genoCall.equals(Byte.toString(AlleleCallCode.BCall.getNumCode())))
        {
            return AlleleCallCode.BCall;
        }
        else if(genoCall.equals("H") || genoCall.equals("HH") ||
                genoCall.equals(Byte.toString(AlleleCallCode.HCall.getNumCode())))
        {
            return AlleleCallCode.HCall;
        }
        else if(genoCall.length() == 0 || genoCall.equals("N") || genoCall.equals("-") || genoCall.equals("NN") ||
                genoCall.equals(Byte.toString(AlleleCallCode.NCall.getNumCode())))
        {
            return AlleleCallCode.NCall;
        }
        else
        {
            return AlleleCallCode.NCall;
        }
    }

    private char[] aAlleles;
    private char[] bAlleles;
    private String[] sampleIds;
    private byte[][] callMatrix;
    private String[] snpIds;
    private String[] chrIDs;
    private long[] bpPositions;
    private String buildId;
    
    public char[] getAAlleles()
    {
        return this.aAlleles;
    }
    
    public void setAAlleles(char[] aAlleles)
    {
        this.aAlleles = aAlleles;
    }
    
    public char[] getBAlleles()
    {
        return this.bAlleles;
    }
    
    public void setBAlleles(char[] bAlleles)
    {
        this.bAlleles = bAlleles;
    }
    
    public String[] getSampleIds()
    {
        return this.sampleIds;
    }
    
    public void setSampleIds(String[] sampleIds)
    {
        this.sampleIds = sampleIds;
    }
    
    public byte[][] getCallMatrix()
    {
        return this.callMatrix;
    }
    
    public void setCallMatrix(byte[][] callMatrix)
    {
        this.callMatrix = callMatrix;
    }
    
    public String[] getSnpIds()
    {
        return this.snpIds;
    }
    
    public void setSnpIds(String[] snpIds)
    {
        this.snpIds = snpIds;
    }
    
    public String[] getChrIDs()
    {
        return this.chrIDs;
    }
    
    public void setChrIDs(String[] chrIDs)
    {
        this.chrIDs = chrIDs;
    }
    
    public long[] getBpPositions()
    {
        return this.bpPositions;
    }
    
    public void setBpPositions(long[] bpPositions)
    {
        this.bpPositions = bpPositions;
    }
    
    public String getBuildId()
    {
        return this.buildId;
    }
    
    public void setBuildId(String buildId)
    {
        this.buildId = buildId;
    }
}
