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
import java.util.Comparator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A class for sorting call matrix by genomic position.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class CallMatrixSorter
{
    // after matching this pattern against a valid chromosome ID group 2
    // should contain the chromosome number (or in the case or X, Y and M
    // the chromosome letter)
    private static final Pattern CHR_PATTERN = Pattern.compile(
            "^(chromosome|chr)?\\s*(\\S+)$", Pattern.CASE_INSENSITIVE);

    private static class IndexComparator implements Comparator<Integer>
    {
        private final String[] chromosomeIds;
        private final long[] bpPositions;
        
        public IndexComparator(String[] chromosomeIds, long[] bpPositions)
        {
            this.chromosomeIds = chromosomeIds;
            this.bpPositions = bpPositions;
        }
        
        /**
         * {@inheritDoc}
         */
        public int compare(Integer index1, Integer index2)
        {
            int i1 = index1.intValue();
            int i2 = index2.intValue();
            int chrComp = compareChrs(
                    this.chromosomeIds[i1],
                    this.chromosomeIds[i2]);
            if(chrComp != 0)
            {
                return chrComp;
            }
            else
            {
                return (int)(this.bpPositions[i1] - this.bpPositions[i2]);
            }
        }
    }
    
    private static <T> T[] reorder(T[] arr, Integer[] order)
    {
        T[] newArr = arr.clone();
        for(int i = 0; i < newArr.length; i++)
        {
            newArr[i] = arr[order[i]];
        }
        return newArr;
    }
    
    private static char[] reorderChar(char[] arr, Integer[] order)
    {
        char[] newArr = arr.clone();
        for(int i = 0; i < newArr.length; i++)
        {
            newArr[i] = arr[order[i]];
        }
        return newArr;
    }
    
    private static long[] reorderLong(long[] arr, Integer[] order)
    {
        long[] newArr = arr.clone();
        for(int i = 0; i < newArr.length; i++)
        {
            newArr[i] = arr[order[i]];
        }
        return newArr;
    }
    
    /**
     * Compare two chromosome strings following the rules set by
     * {@link Comparator}. Numbered chromosomes will be compared
     * according to their natural number order. The named chromosomes are
     * considered greater than any numbered chromosome and will be ordered as:
     * X, Y and finally M. This comparison will tolerate a proceeding "chr" or
     * "chromosome" before the chromosome ID as in "chrX" or "chromosome y" but
     * other unrecognized formats will case a {@link IllegalArgumentException}
     * to be thrown.
     * @param chrName1  the 1st chromosome string
     * @param chrName2  the 2nd chromosome string
     * @return  the comparison value
     */
    public static int compareChrs(String chrName1, String chrName2)
    {
        return chrToInt(chrName1) - chrToInt(chrName2);
    }
    
    private static int chrToInt(String chr)
    {
        Matcher matcher = CHR_PATTERN.matcher(chr);
        if(matcher.matches())
        {
            String matchedChr = matcher.group(2);
            // numbers come first followed by X, Y and M
            try
            {
                return Integer.parseInt(matchedChr);
            }
            catch(NumberFormatException ex)
            {
                matchedChr = matchedChr.toUpperCase();
                if(matchedChr.equals("X"))
                {
                    return Integer.MAX_VALUE - 2;
                }
                else if(matchedChr.equals("Y"))
                {
                    return Integer.MAX_VALUE - 1;
                }
                else if(matchedChr.equals("M"))
                {
                    return Integer.MAX_VALUE;
                }
                else
                {
                    throw new IllegalArgumentException(
                            "\"" + chr + "\" is not a valid chromosome name");
                }
            }
        }
        else
        {
            throw new IllegalArgumentException(
                    "\"" + chr + "\" is not a valid chromosome name");
        }
    }
    
    /**
     * Sort the given call matrix by genomic position.
     * @param callMat the call matrix to sort
     */
    public static void sortCallMatrix(AbstractGenotypeCallMatrix callMat)
    {
        IndexComparator c = new IndexComparator(
                callMat.getChrIDs(),
                callMat.getBpPositions());
        Integer[] indices = new Integer[(int)callMat.getSNPCount()];
        for(int i = 0; i < indices.length; i++)
        {
            indices[i] = i;
        }
        Arrays.sort(indices, c);
        
        callMat.setAAlleles(reorderChar(callMat.getAAlleles(), indices));
        callMat.setBAlleles(reorderChar(callMat.getBAlleles(), indices));
        callMat.setBpPositions(reorderLong(callMat.getBpPositions(), indices));
        callMat.setChrIDs(reorder(callMat.getChrIDs(), indices));
        callMat.setSnpIds(reorder(callMat.getSnpIds(), indices));
        callMat.setCallMatrix(reorder(callMat.getCallMatrix(), indices));
        
        callMat.setSortedByPosition(true);
    }
}
