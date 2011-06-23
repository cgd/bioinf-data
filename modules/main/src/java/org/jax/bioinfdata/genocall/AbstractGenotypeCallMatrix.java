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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Abstract representation of a genotype matrix.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public abstract class AbstractGenotypeCallMatrix
{
    public static final String A_ALLELES_NAME = "aAllele";
    public static final String B_ALLELES_NAME = "bAllele";
    public static final String SAMPLE_IDS_NAME = "sampleID";
    public static final String CALL_MATRIX_NAME = "callMatrix";
    public static final String SNP_IDS_NAME = "snpID";
    public static final String CHR_IDS_NAME = "chrID";
    public static final String BP_POSITIONS_NAME = "bpPosition";
    public static final String BUILD_ID_NAME = "buildID";
    public static final String SORTED_BY_POS_NAME = "sortedByPosition";
    
    public static final byte A_CALL_CODE = 1;
    public static final byte B_CALL_CODE = 2;
    public static final byte H_CALL_CODE = 3;
    public static final byte N_CALL_CODE = -1;
    
    private static final List<String> NUC_CODES = Arrays.asList("G", "A", "T", "C");
    
    /**
     * Convert the given genoCalls strings into a byte array.
     * @param aAllele
     *          call representing an A allele for the SNPs in question (can be null)
     * @param bAllele
     *          call representing a B allele for the SNPs in question (can be null)
     * @param genoCalls
     *          the call string
     * @return
     *          the byte representation of the call
     */
    public static byte[] toCallValues(String aAllele, String bAllele, String[] genoCalls)
    {
        if(aAllele == null || bAllele == null)
        {
            // we don't have A or B alleles to match against. First see if
            // everything is encoded according to it's byte code
            try
            {
                return toCallValuesWithoutAB(genoCalls);
            }
            catch(NumberFormatException ex)
            {
                // we failed to encode using a byte code so we will fall
                // back on trying to encode with nucleotide values
                Set<String> nucCodesInSNP = new HashSet<String>(4);
                for(int i = 0; i < genoCalls.length; i++)
                {
                    String genoCall = genoCalls[i].toUpperCase();
                    if(NUC_CODES.contains(genoCall))
                    {
                        nucCodesInSNP.add(genoCall);
                    }
                }
                
                if(nucCodesInSNP.size() == 2)
                {
                    // we can't do anything with this SNP unless we have
                    // exactly 2 calls
                    String[] codes = nucCodesInSNP.toArray(new String[2]);
                    return toCallValuesWithAB(codes[0], codes[1], genoCalls);
                }
                else
                {
                    byte[] callValues = new byte[genoCalls.length];
                    Arrays.fill(callValues, N_CALL_CODE);
                    return callValues;
                }
            }
        }
        else
        {
            return toCallValuesWithAB(aAllele, bAllele, genoCalls);
        }
    }
    
    private static byte[] toCallValuesWithAB(String aAllele, String bAllele, String[] genoCalls)
    {
        byte[] callValues = new byte[genoCalls.length];
        Arrays.fill(callValues, N_CALL_CODE);
        
        aAllele = aAllele.toUpperCase();
        bAllele = bAllele.toUpperCase();
        for(int i = 0; i < genoCalls.length; i++)
        {
            String genoCall = genoCalls[i].toUpperCase();
            if(genoCall.equals(aAllele) ||
               genoCall.equals(Byte.toString(A_CALL_CODE)))
            {
                callValues[i] = A_CALL_CODE;
            }
            else if(genoCall.equals(bAllele) ||
                    genoCall.equals(Byte.toString(B_CALL_CODE)))
            {
                callValues[i] = B_CALL_CODE;
            }
            else if(genoCall.equals("H") || genoCall.equals("HH") ||
                    genoCall.equals(Byte.toString(H_CALL_CODE)))
            {
                callValues[i] = H_CALL_CODE;
            }
        }
        
        return callValues;
    }
    
    private static byte[] toCallValuesWithoutAB(String[] genoCalls)
    {
        byte[] callValues = new byte[genoCalls.length];
        Arrays.fill(callValues, N_CALL_CODE);
        
        // everything is encoded according to it's byte code
        for(int i = 0; i < genoCalls.length; i++)
        {
            callValues[i] = Byte.parseByte(genoCalls[i].trim());
        }
        
        return callValues;
    }
    
    /**
     * Copy values from one genotype matrix to the other
     * @param fromMat   the source matrix
     * @param toMat     the destination matrix
     */
    public static void copyGenoMatrix(
            AbstractGenotypeCallMatrix fromMat,
            AbstractGenotypeCallMatrix toMat)
    {
        toMat.setAAlleles(fromMat.getAAlleles());
        toMat.setBAlleles(fromMat.getBAlleles());
        toMat.setSampleIds(fromMat.getSampleIds());
        toMat.setCallMatrix(fromMat.getCallMatrix());
        toMat.setSnpIds(fromMat.getSnpIds());
        toMat.setChrIDs(fromMat.getChrIDs());
        toMat.setBpPositions(fromMat.getBpPositions());
        toMat.setBuildId(fromMat.getBuildId());
        toMat.setSortedByPosition(fromMat.isSortedByPosition());
    }
    
    /**
     * Convert the following calls into a BitSet. This function expects all
     * snpCalls to be A or B and will throw an exception otherwise
     * @param snpCalls
     *          the SNP calls to convert
     * @param minNorm
     *          should the calls be "minority normalized". If true then the
     *          minority allele will be 1 or in the case where both alleles
     *          occur with the same frequency the allele that occurs in the
     *          first position will be 1. If false the A allele will be 1 and
     *          the B allele will be 0
     * @return
     *          the bit set
     * @throws IllegalArgumentException
     *          in the case where the snpCalls are not all A or B
     */
    public static BitSet toBitSet(byte[] snpCalls, boolean minNorm) throws IllegalArgumentException
    {
        BitSet snpBits = new BitSet(snpCalls.length);
        for(int i = 0; i < snpCalls.length; i++)
        {
            if(snpCalls[i] == AbstractGenotypeCallMatrix.A_CALL_CODE)
            {
                snpBits.set(i);
            }
            else if(snpCalls[i] != B_CALL_CODE)
            {
                throw new IllegalArgumentException(
                        "the toBitSet function is only valid for SNPs where " +
                        "all calls are A or B");
            }
        }
        
        if(minNorm)
        {
            // if there are more 1's than 0's OR if the 1's count is equal to
            // the 0's count AND the 1st bit is on, then we need to flip all
            // of the bits to normalize
            int doubleOnesCount = snpBits.cardinality() * 2;
            if(doubleOnesCount > snpCalls.length ||
               (doubleOnesCount == snpCalls.length && snpBits.get(0)))
            {
                snpBits.flip(0, snpCalls.length);
            }
        }
        
        return snpBits;
    }

    /**
     * getter for the A alleles
     * @return  the A alleles
     */
    public abstract char[] getAAlleles();
    
    /**
     * setter for the A alleles
     * @param aAlleles  the A alleles
     */
    public abstract void setAAlleles(char[] aAlleles);
    
    /**
     * getter for the B alleles
     * @return  the B alleles
     */
    public abstract char[] getBAlleles();
    
    /**
     * setter for the B alleles
     * @param bAlleles  the B alleles
     */
    public abstract void setBAlleles(char[] bAlleles);
    
    /**
     * getter for the sample IDs
     * @return  the sample IDs
     */
    public abstract String[] getSampleIds();
    
    /**
     * setter for the sample IDs
     * @param sampleIds the sample IDs
     */
    public abstract void setSampleIds(String[] sampleIds);
    
    /**
     * getter for the call matrix
     * @return  the call matrix
     */
    public abstract byte[][] getCallMatrix();
    
    /**
     * setter for the call matrix
     * @param callMatrix    the call matrix
     */
    public abstract void setCallMatrix(byte[][] callMatrix);
    
    /**
     * Extracts a row from the {@link #getCallMatrix()}. Note that if the
     * implementation is {@link HDF5GenotypeCallMatrix} using this function will be more
     * memory efficient than {@link #getCallMatrix()} but may be less CPU
     * efficient (depending on how many SNP rows you need to access)
     * @param snpIndex
     *      this is the same as the row in the call matrix
     * @return
     *      the call values for the given SNP
     */
    public abstract byte[] getSNPCalls(long snpIndex);
    
    /**
     * Getter for the SNP count
     * @return  the SNP count
     */
    public abstract long getSNPCount();
    
    /**
     * getter for SNP IDs
     * @return  the SNP IDs
     */
    public abstract String[] getSnpIds();
    
    /**
     * setter for the SNP IDs
     * @param snpIds    the SNP IDs
     */
    public abstract void setSnpIds(String[] snpIds);
    
    /**
     * getter for the chromosome IDs
     * @return  the chromosome IDs
     */
    public abstract String[] getChrIDs();
    
    /**
     * setter for chromosome IDs
     * @param chrIDs    the chromosome IDs
     */
    public abstract void setChrIDs(String[] chrIDs);
    
    /**
     * getter for the BP positions
     * @return  the BP positions
     */
    public abstract long[] getBpPositions();
    
    /**
     * setter for the BP positions
     * @param bpPositions   the bp positions
     */
    public abstract void setBpPositions(long[] bpPositions);
    
    /**
     * Getter for the build ID (eg "build 37")
     * @return  the build ID
     */
    public abstract String getBuildId();
    
    /**
     * Setter for the build ID
     * @param buildId   the build ID
     */
    public abstract void setBuildId(String buildId);
    
    /**
     * Determines if this is sorted by position
     * @return  true iff this is sorted by position
     */
    public abstract boolean isSortedByPosition();
    
    /**
     * Setter indicating if the entries are sorted by position or not
     * @param sortedByPosition  indicates if this is sorted or not
     */
    public abstract void setSortedByPosition(boolean sortedByPosition);
    
    /**
     * Uses {@link #getSnpSubsetView(int, int)} to extract views for all of
     * the regions where there is a contiguous chromosome value
     * @return  the contiguous chromosome views
     * @throws NullPointerException if {@link #getChrIDs()} returns null
     */
    public List<AbstractGenotypeCallMatrix> getChromosomeViews()
    {
        List<AbstractGenotypeCallMatrix> chrViews =
            new ArrayList<AbstractGenotypeCallMatrix>();
        
        String[] chrIds = this.getChrIDs();
        if(chrIds.length >= 1)
        {
            int currStart = 0;
            String currContigChr = chrIds[0];
            for(int i = 1; i < chrIds.length; i++)
            {
                if(!currContigChr.equals(chrIds[i]))
                {
                    chrViews.add(this.getSnpSubsetView(currStart, i - currStart));
                    currStart = i;
                    currContigChr = chrIds[i];
                }
            }
            
            if(chrViews.size() == 0)
            {
                chrViews.add(this);
            }
            else
            {
                chrViews.add(this.getSnpSubsetView(currStart, chrIds.length - currStart));
            }
        }
        
        return chrViews;
    }
    
    /**
     * Get a read-only subset view of this call matrix
     * @param startSnpIndex the starting index
     * @param extent        the extent
     * @return              the subset view
     */
    public AbstractGenotypeCallMatrix getSnpSubsetView(
            final int startSnpIndex,
            final int extent)
    {
        return new AbstractGenotypeCallMatrix()
        {
            @Override
            public void setSortedByPosition(boolean sortedByPosition)
            {
                throw new UnsupportedOperationException(
                        "This function is not supported in a subset view");
            }
            
            @Override
            public void setSnpIds(String[] snpIds)
            {
                throw new UnsupportedOperationException(
                        "This function is not supported in a subset view");
            }
            
            @Override
            public void setSampleIds(String[] sampleIds)
            {
                throw new UnsupportedOperationException(
                        "This function is not supported in a subset view");
            }
            
            @Override
            public void setChrIDs(String[] chrIDs)
            {
                throw new UnsupportedOperationException(
                        "This function is not supported in a subset view");
            }
            
            @Override
            public void setCallMatrix(byte[][] callMatrix)
            {
                throw new UnsupportedOperationException(
                        "This function is not supported in a subset view");
            }
            
            @Override
            public void setBuildId(String buildId)
            {
                throw new UnsupportedOperationException(
                        "This function is not supported in a subset view");
            }
            
            @Override
            public void setBpPositions(long[] bpPositions)
            {
                throw new UnsupportedOperationException(
                        "This function is not supported in a subset view");
            }
            
            @Override
            public void setBAlleles(char[] bAlleles)
            {
                throw new UnsupportedOperationException(
                        "This function is not supported in a subset view");
            }
            
            @Override
            public void setAAlleles(char[] aAlleles)
            {
                throw new UnsupportedOperationException(
                        "This function is not supported in a subset view");
            }
            
            @Override
            public boolean isSortedByPosition()
            {
                return AbstractGenotypeCallMatrix.this.isSortedByPosition();
            }
            
            @Override
            public String[] getSnpIds()
            {
                String[] snpIds = AbstractGenotypeCallMatrix.this.getSnpIds();
                return Arrays.copyOfRange(snpIds, startSnpIndex, startSnpIndex + extent);
            }
            
            @Override
            public String[] getSampleIds()
            {
                return AbstractGenotypeCallMatrix.this.getSampleIds();
            }
            
            @Override
            public long getSNPCount()
            {
                return extent;
            }
            
            @Override
            public byte[] getSNPCalls(long snpIndex)
            {
                return AbstractGenotypeCallMatrix.this.getSNPCalls(snpIndex + startSnpIndex);
            }
            
            @Override
            public String[] getChrIDs()
            {
                String[] chrIds = AbstractGenotypeCallMatrix.this.getChrIDs();
                return Arrays.copyOfRange(chrIds, startSnpIndex, startSnpIndex + extent);
            }
            
            @Override
            public byte[][] getCallMatrix()
            {
                byte[][] callMat = AbstractGenotypeCallMatrix.this.getCallMatrix();
                return Arrays.copyOfRange(callMat, startSnpIndex, startSnpIndex + extent);
            }
            
            @Override
            public String getBuildId()
            {
                return AbstractGenotypeCallMatrix.this.getBuildId();
            }
            
            @Override
            public long[] getBpPositions()
            {
                long[] bpPos = AbstractGenotypeCallMatrix.this.getBpPositions();
                return Arrays.copyOfRange(bpPos, startSnpIndex, startSnpIndex + extent);
            }
            
            @Override
            public char[] getBAlleles()
            {
                char[] bAlleles = AbstractGenotypeCallMatrix.this.getBAlleles();
                return Arrays.copyOfRange(bAlleles, startSnpIndex, startSnpIndex + extent);
            }
            
            @Override
            public char[] getAAlleles()
            {
                char[] aAlleles = AbstractGenotypeCallMatrix.this.getAAlleles();
                return Arrays.copyOfRange(aAlleles, startSnpIndex, startSnpIndex + extent);
            }
        };
    }
}
