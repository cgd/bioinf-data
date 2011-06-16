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
import java.util.List;

import org.jax.bioinfdata.genocall.AbstractGenotypeCallMatrix;
import org.jax.bioinfdata.genocall.ReverseGenotypeCallMatrix;
import org.jax.util.datastructure.ExposedArrayList;
import org.jax.util.datastructure.SequenceUtilities;

/**
 * A class for building compatible intervals where a compatible interval
 * is a contiguous region of the genome where all SNP value differences can
 * be explained by mutations using the "infinite sites" assumption. One
 * nicety of compatible intervals is that you can always build perfect
 * phylogenies from them.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class IntervalScanner
{
    /**
     * Class for grouping an SDP with its index
     */
    private static class SdpIndexPair
    {
        private final byte[] sdpBits;
        
        private final int index;

        /**
         * Constructor
         * @param sdpBits
         *          see {@link #getSdpBits()}
         * @param index
         *          see {@link #getIndex()}
         */
        public SdpIndexPair(byte[] sdpBits, int index)
        {
            this.sdpBits = sdpBits;
            this.index = index;
        }
        
        /**
         * Getter for the index
         * @return the index
         */
        public int getIndex()
        {
            return this.index;
        }
        
        /**
         * getter for the SDP
         * @return the SDP
         */
        public byte[] getSdpBits()
        {
            return this.sdpBits;
        }
    }
    
    /**
     * Private constructor. This class should not be instantiated since all
     * functions and members are static.
     */
    private IntervalScanner() {}
    
    /**
     * Do a max-k scan which involves doing a lot of other scans.
     * @param callMatrix
     *          the calls to scan
     * @return
     *          the max-k interval
     */
    public static List<IndexedSnpInterval> maxKScan(AbstractGenotypeCallMatrix callMatrix)
    {
        List<IndexedSnpInterval> forwardIntervals = greedyScan(callMatrix);
        List<IndexedSnpInterval> reverseIntervals = greedyScan(new ReverseGenotypeCallMatrix(callMatrix));
        reverseIndexedIntervals(reverseIntervals, (int)callMatrix.getSNPCount());
        List<IndexedSnpInterval> uberIntervals = uberScan(callMatrix);
        List<IndexedSnpInterval> coreIntervals = createCoreIntervals(forwardIntervals, reverseIntervals);
        List<List<IndexedSnpInterval>> uberCores = createUberCores(uberIntervals, coreIntervals);
        List<IndexedSnpInterval> maxKIntervals = createMaxKIntervals(uberCores);
        
        return maxKIntervals;
    }
    
    /**
     * Do an uber-scan looking for every possible maximal compatible interval
     * @param callMatrix
     *          the calls to scan
     * @return
     *          the uber list of compatible intervals
     */
    public static List<IndexedSnpInterval> uberScan(AbstractGenotypeCallMatrix callMatrix)
    {
        ArrayList<IndexedSnpInterval> intervals = new ArrayList<IndexedSnpInterval>();
        
        long sdpCount = callMatrix.getSNPCount();
        int startIndex = 0;
        ExposedArrayList<SdpIndexPair> intervalSdps = new ExposedArrayList<SdpIndexPair>();
        int nearestIncompatibleIndex = -1;
        while(startIndex < sdpCount)
        {
            nearestIncompatibleIndex = testCompatibleAndUberAdd(
                    intervalSdps,
                    callMatrix.getSNPCalls(startIndex),
                    startIndex);
            assert nearestIncompatibleIndex == -1;
            int nextIndex = startIndex + 1;
            while(nextIndex < sdpCount)
            {
                nearestIncompatibleIndex = testCompatibleAndUberAdd(
                        intervalSdps,
                        callMatrix.getSNPCalls(nextIndex),
                        nextIndex);
                if(nearestIncompatibleIndex == -1)
                {
                    nextIndex++;
                }
                else
                {
                    break;
                }
            }
            intervals.add(new IndexedSnpInterval(
                    startIndex,
                    nextIndex - startIndex));
            if(nearestIncompatibleIndex != -1)
            {
                startIndex = nearestIncompatibleIndex + 1;
            }
            else
            {
                assert nextIndex == sdpCount;
                break;
            }
        }
        
        intervals.trimToSize();
        
        return intervals;
    }
    
    /**
     * Create a set of core intervals from the given forward and reverse
     * greedy scan results
     * @param forwardGreedyIntervals
     *          the forward scan intervals
     * @param reverseGreedyIntervals
     *          the reverse scan intervals
     * @return
     *          the cores
     */
    public static List<IndexedSnpInterval> createCoreIntervals(
            List<IndexedSnpInterval> forwardGreedyIntervals,
            List<IndexedSnpInterval> reverseGreedyIntervals)
    {
        if(forwardGreedyIntervals.size() != reverseGreedyIntervals.size())
        {
            throw new IllegalArgumentException(
                    "the reverse and forward interval lists should be the " +
                    "same size");
        }
        
        int intervalCount = forwardGreedyIntervals.size();
        List<IndexedSnpInterval> coreIntervals =
            new ArrayList<IndexedSnpInterval>(intervalCount);
        for(int intervalIndex = 0; intervalIndex < intervalCount; intervalIndex++)
        {
            IndexedSnpInterval currForwardInterval =
                forwardGreedyIntervals.get(intervalIndex);
            IndexedSnpInterval currReverseInterval =
                reverseGreedyIntervals.get(intervalIndex);
            
            int coreIntervalStart = currForwardInterval.getStartIndex();
            int coreIntervalEnd = currReverseInterval.getEndIndex();
            assert coreIntervalStart <= coreIntervalEnd;
            
            coreIntervals.add(new IndexedSnpInterval(
                    coreIntervalStart,
                    1 + coreIntervalEnd - coreIntervalStart));
        }
        
        return coreIntervals;
    }
    
    /**
     * Prune the uber intervals in preparation for calculating the max-k
     * intervals. We can do this by throwing out any uber intervals that do
     * not intersect exactly one core interval since we know that max-k
     * intervals should intersect one and only one core interval
     * @param uberIntervals
     *          the uber set of intervals (unmodified by this function)
     * @param coreIntervals
     *          the core set of intervals (unmodified by this function)
     * @return
     *          a new list of intervals which is a subset of the uber set
     *          following some rules that we know apply to the max-k set
     *          of intervals
     */
    public static List<List<IndexedSnpInterval>> createUberCores(
            List<IndexedSnpInterval> uberIntervals,
            List<IndexedSnpInterval> coreIntervals)
    {
        if(uberIntervals.size() < coreIntervals.size())
        {
            throw new IllegalArgumentException(
                    "the list of uber intervals should be at least as big as " +
                    "the list of core intervals");
        }
        
        int uberSize = uberIntervals.size();
        int coreSize = coreIntervals.size();
        List<List<IndexedSnpInterval>> uberIntervalsSubset =
            new ArrayList<List<IndexedSnpInterval>>(coreSize);
        
        if(coreSize >= 1)
        {
            // initialize the data that we'll use in our subsetting loop
            int coreIndex = 0;
            IndexedSnpInterval prevCore = null;
            IndexedSnpInterval currCore = coreIntervals.get(coreIndex);
            IndexedSnpInterval nextCore = null;
            if(coreIndex + 1 < coreSize)
            {
                nextCore = coreIntervals.get(coreIndex + 1);
            }
            
            ArrayList<IndexedSnpInterval> currCoreUberIntervals =
                new ArrayList<IndexedSnpInterval>();
            
            // OK, take care of the subsetting
            for(int uberIndex = 0; uberIndex < uberSize && currCore != null; uberIndex++)
            {
                IndexedSnpInterval currUberInterval =
                    uberIntervals.get(uberIndex);
                
                int uberStart = currUberInterval.getStartIndex();
                
                // do we need to move on to the next core?
                if(uberStart > currCore.getEndIndex())
                {
                    assert !currCoreUberIntervals.isEmpty();
                    assert SequenceUtilities.isSorted(currCoreUberIntervals);
                    
                    // move on to the next core
                    coreIndex++;
                    prevCore = currCore;
                    currCore = nextCore;
                    nextCore = null;
                    if(coreIndex + 1 < coreSize)
                    {
                        nextCore = coreIntervals.get(coreIndex + 1);
                    }
                    uberIntervalsSubset.add(currCoreUberIntervals);
                    currCoreUberIntervals.trimToSize();
                    currCoreUberIntervals =
                        new ArrayList<IndexedSnpInterval>(1);
                }
                
                if(currCore != null)
                {
                    assert uberStart <= currCore.getEndIndex();
                    
                    // we have the correct core, now see if the uber interval
                    // has any chance to become a max-k interval (for this it
                    // must cover only the current core, without intersecting
                    // with the previous or next cores)
                    if(currUberInterval.contains(currCore) &&
                       (prevCore == null || !currUberInterval.intersects(prevCore)) &&
                       (nextCore == null || !currUberInterval.intersects(nextCore)))
                    {
                        // pass this one through
                        currCoreUberIntervals.add(currUberInterval);
                    }
                }
            }
            
            // clean-up
            if(!currCoreUberIntervals.isEmpty())
            {
                assert SequenceUtilities.isSorted(currCoreUberIntervals);
                uberIntervalsSubset.add(currCoreUberIntervals);
            }
            
            assert uberIntervalsSubset.size() == coreSize;
        }
        
        return uberIntervalsSubset;
    }

    /**
     * Add to the interval SDPs assuming that we're working toward the
     * Uber set of compatible intervals. This function modifies the list
     * in different ways depending on whether it finds a match or a conflict,
     * but you shouldn't really have to care about that.
     * <br/><br/>
     * If this function finds a conflict, it returns the the index of the
     * nearest conflicting SNP.
     * @param intervalSdps
     *          the list of SDPs for the current interval. this list is "owned"
     *          by this function and it's probably best not to touch it or try
     *          to make sense of it outside of this function
     * @param sdpToAdd
     *          the SDP that we're trying to introduce to the interval
     * @param sdpIndex
     *          the index of the SDP
     * @return
     *          -1 if the given SDP is compatible with current interval SDPs
     *          or the index of the nearest incompatibility if we find one
     */
    private static int testCompatibleAndUberAdd(
            ExposedArrayList<SdpIndexPair> intervalSdps,
            byte[] sdpToAdd,
            int sdpIndex)
    {
        int intervalSdpCount = intervalSdps.size();
        for(int i = intervalSdpCount - 1; i >= 0 ; i--)
        {
            SdpIndexPair currPair = intervalSdps.get(i);
            byte[] currSdp = currPair.getSdpBits();
            if(sdpToAdd.equals(currSdp))
            {
                // this SDP was already added, so return compatible after
                // moving the snp to the end of the list
                intervalSdps.remove(i);
                intervalSdps.add(new SdpIndexPair(sdpToAdd, sdpIndex));
                return -1;
            }
            else if(!areSdpsCompatible(sdpToAdd, currSdp))
            {
                // found an incompatibility. remove everything before the
                // incompatibility and return the index of the sdp that
                // caused the trouble
                intervalSdps.removeRange(
                        0,
                        i + 1);
                intervalSdps.add(new SdpIndexPair(sdpToAdd, sdpIndex));
                
                return currPair.getIndex();
            }
        }
        
        // we're compatible and didn't find any exact matches for the SDP
        // which means that we need to add it to the list
        intervalSdps.add(new SdpIndexPair(sdpToAdd, sdpIndex));
        return -1;
    }
    
    /**
     * Perform a greedy scan on the given streams. If the stream we're given
     * reads in the reverse direction then we'll do a
     * {@link #reverseIndexedIntervals(List, int)} on the intervals before
     * returning. This means that the indices should be easily comparable with
     * intervals read in the forward direction.
     * @param callMatrix
     *          the genotypes to scan
     * @return
     *          the list of SNP intervals. these intervals will cover the
     *          entire genome with no overlap
     */
    public static List<IndexedSnpInterval> greedyScan(AbstractGenotypeCallMatrix callMatrix)
    {
        ArrayList<IndexedSnpInterval> intervals = new ArrayList<IndexedSnpInterval>();
        
        List<byte[]> intervalSdps = new ArrayList<byte[]>();
        long sdpCount = callMatrix.getSNPCount();
        int startIndex = 0;
        while(startIndex < sdpCount)
        {
            intervalSdps.add(callMatrix.getSNPCalls(startIndex));
            int incompatIndex = startIndex + 1;
            while(incompatIndex < sdpCount &&
                  checkCompatibilityAndAddSdp(intervalSdps, callMatrix.getSNPCalls(incompatIndex)))
            {
                incompatIndex++;
            }
            intervals.add(new IndexedSnpInterval(
                    startIndex,
                    incompatIndex - startIndex));
            
            intervalSdps.clear();
            startIndex = incompatIndex;
        }
        
        intervals.trimToSize();
        return intervals;
    }

    /**
     * Check the compatibility of the SDP and if it's compatible, add it
     * @param intervalSdps  the interval SDPs to compare against
     * @param sdpToAdd      the SDP we'll try to add to the interval
     * @return  true if the SDP is fully compatible, false otherwise
     */
    private static boolean checkCompatibilityAndAddSdp(
            List<byte[]> intervalSdps,
            byte[] sdpToAdd)
    {
        int intervalSdpCount = intervalSdps.size();
        for(int i = 0; i < intervalSdpCount; i++)
        {
            byte[] currSdp = intervalSdps.get(i);
            if(sdpToAdd.equals(currSdp))
            {
                // this SDP was already added, so return compatible
                // without adding SDP
                return true;
            }
            else if(!areSdpsCompatible(sdpToAdd, currSdp))
            {
                // found an incompatibility, so return incompatible
                // without adding SDP
                return false;
            }
        }
        
        // we're compatible and didn't find any exact matches for the SDP
        // which means that we need to add it to the list
        intervalSdps.add(sdpToAdd);
        return true;
    }
    
    /**
     * Test if the two SDPs are compatible. Both SDPs must be "minority
     * normalized" meaning that the minority allele bits will be set and
     * the majority allele bits will be unset
     * @param minorityNormalizedSdp1    the 1st SDP
     * @param minorityNormalizedSdp2    the 2nd SDP
     * @return  true if the given SDPs are compatible
     */
    public static boolean areMinorityNormalizedSdpsCompatible(
            BitSet minorityNormalizedSdp1,
            BitSet minorityNormalizedSdp2)
    {
        if(!minorityNormalizedSdp1.intersects(minorityNormalizedSdp2))
        {
            // disjoint SDPs are compatible
            return true;
        }
        else
        {
            BitSet intersection = (BitSet)minorityNormalizedSdp1.clone();
            intersection.and(minorityNormalizedSdp2);
            
            if(intersection.equals(minorityNormalizedSdp1))
            {
                // SDP1 is a subset of SDP2 indicating compatible
                return true;
            }
            else if(intersection.equals(minorityNormalizedSdp2))
            {
                // SDP2 is a subset of SDP1 indicating compatible
                return true;
            }
            else
            {
                // they're incompatible since they're not disjoint and they
                // don't have a subset/superset relationship
                return false;
            }
        }
    }
    
    /**
     * Test if the two SDPs are compatible. If we observe all four gametes they
     * are not compatible.
     * @param sdp1  the 1st SDP
     * @param sdp2  the 2nd SDP
     * @return      true if the given SDPs are compatible
     */
    public static boolean areSdpsCompatible(byte[] sdp1, byte[] sdp2)
    {
        int len = sdp1.length;
        assert sdp2.length == len;
        
        boolean observedAA = false;
        boolean observedAB = false;
        boolean observedBA = false;
        boolean observedBB = false;
        
        for(int i = 0; i < len; i++)
        {
            byte call1 = sdp1[i];
            byte call2 = sdp2[i];
            boolean call1IsA = call1 == AbstractGenotypeCallMatrix.A_CALL_CODE;
            boolean call1IsB = call1 == AbstractGenotypeCallMatrix.B_CALL_CODE;
            boolean call2IsA = call2 == AbstractGenotypeCallMatrix.A_CALL_CODE;
            boolean call2IsB = call2 == AbstractGenotypeCallMatrix.B_CALL_CODE;
            
            if(call1IsA)
            {
                if(call2IsA)
                {
                    observedAA = true;
                }
                else if(call2IsB)
                {
                    observedAB = true;
                }
            }
            else if(call1IsB)
            {
                if(call2IsA)
                {
                    observedBA = true;
                }
                else if(call2IsB)
                {
                    observedBB = true;
                }
            }
        }
        
        boolean observedFourGametes = observedAA && observedAB && observedBA && observedBB;
        assert !observedFourGametes == areMinorityNormalizedSdpsCompatible(
                        AbstractGenotypeCallMatrix.toBitSet(sdp1, true),
                        AbstractGenotypeCallMatrix.toBitSet(sdp2, true));
        return !observedFourGametes;
    }
    
    /**
     * Performs {@link #reverseIndexedInterval(IndexedSnpInterval, int)} on
     * every interval in the list and also reverses the ordering of intervals
     * in the list
     * @param intervalsToReverse
     *          the list to reverse
     * @param totalSnpCount
     *          the total snp count (range) that we're flipping the intervals
     *          over
     */
    private static void reverseIndexedIntervals(
            List<IndexedSnpInterval> intervalsToReverse,
            int totalSnpCount)
    {
        int intervalCount = intervalsToReverse.size();
        int halfIntervalCount = intervalCount / 2;
        
        for(int index1 = 0; index1 < halfIntervalCount; index1++)
        {
            int index2 = (intervalCount - index1) - 1;
            
            IndexedSnpInterval interval1 = intervalsToReverse.get(index1);
            IndexedSnpInterval interval2 = intervalsToReverse.get(index2);
            
            intervalsToReverse.set(
                    index1,
                    reverseIndexedInterval(interval2, totalSnpCount));
            intervalsToReverse.set(
                    index2,
                    reverseIndexedInterval(interval1, totalSnpCount));
        }
        
        // if we have an odd count we still need to reverse the middle interval
        if(intervalCount % 2 == 1)
        {
            int middleIndex = halfIntervalCount;
            IndexedSnpInterval middleInterval =
                intervalsToReverse.get(middleIndex);
            intervalsToReverse.set(
                    middleIndex,
                    reverseIndexedInterval(middleInterval, totalSnpCount));
        }
    }
    
    /**
     * Reverse the indecies assuming the given total index count
     * @param intervalToReverse
     *          the interval 
     * @param totalSnpCount
     *          the total # of SNPs. we need to know this in order to be able
     *          to flip the interval
     * @return
     *          the interval which will be reveresed (or "flipped") along the
     *          range indicated by the total count
     */
    private static IndexedSnpInterval reverseIndexedInterval(
            IndexedSnpInterval intervalToReverse,
            int totalSnpCount)
    {
        int extent = intervalToReverse.getExtentInIndices();
        int newStart = (totalSnpCount - intervalToReverse.getStartIndex()) - extent;
        return new IndexedSnpInterval(
                newStart,
                extent);
    }
    
    /**
     * Scan the given uber cores to come up with the max-k data set
     * @param uberCores
     *          the uber cores to search. see
     *          {@link #createUberCores(List, List)}
     * @return
     *          the max-k interval list
     */
    public static List<IndexedSnpInterval> createMaxKIntervals(List<List<IndexedSnpInterval>> uberCores)
    {
        int coreCount = uberCores.size();
        List<IndexedSnpInterval> maxKIntervals = new ArrayList<IndexedSnpInterval>(
                coreCount);
        
        if(coreCount >= 1)
        {
            // initialize data pre-scan
            int[][] forwardPointers = new int[uberCores.size() - 1][];
            List<IndexedSnpInterval> prevCoreGroup = uberCores.get(
                    uberCores.size() - 1);
            long[] cumulativeExtents = new long[prevCoreGroup.size()];
            for(int i = 0; i < cumulativeExtents.length; i++)
            {
                cumulativeExtents[i] = prevCoreGroup.get(i).getExtentInIndices();
            }
            
            // perform the scan by working your way back from the tail and by
            // constructing forward pointers as you go
            for(int i = coreCount - 2; i >= 0; i--)
            {
                List<IndexedSnpInterval> currCoreGroup = uberCores.get(i);
                int currCoreGroupSize = currCoreGroup.size();
                int[] currForwardPointers = new int[currCoreGroupSize];
                long[] currCumulativeExtents = new long[currCoreGroupSize];
                
                for(int j = 0; j < currCoreGroupSize; j++)
                {
                    IndexedSnpInterval currGroupInterval = currCoreGroup.get(j);
                    int prevCoreGroupSize = prevCoreGroup.size();
                    long maxCumulativeExtent = 0L;
                    for(int k = 0; k < prevCoreGroupSize; k++)
                    {
                        IndexedSnpInterval prevGroupInterval = prevCoreGroup.get(k);
                        long currCumulativeExtent =
                            cumulativeExtents[k] + currGroupInterval.getExtentInIndices();
                        if(currCumulativeExtent > maxCumulativeExtent &&
                           currGroupInterval.getEndIndex() >= prevGroupInterval.getStartIndex() - 1)
                        {
                            maxCumulativeExtent = currCumulativeExtent;
                            currCumulativeExtents[j] = currCumulativeExtent;
                            currForwardPointers[j] = k;
                        }
                    }
                    assert maxCumulativeExtent > 0L;
                }
                
                forwardPointers[i] = currForwardPointers;
                cumulativeExtents = currCumulativeExtents;
                prevCoreGroup = currCoreGroup;
            }
            
            // now we need to work our way forward through the pointers.
            // initialize the 1st pointer using the max cumulative extent
            // start with the max pointer
            int currPointer = 0;
            for(int i = 0; i < cumulativeExtents.length; i++)
            {
                if(cumulativeExtents[i] > cumulativeExtents[currPointer])
                {
                    currPointer = i;
                }
            }
            maxKIntervals.add(uberCores.get(0).get(currPointer));
            
            // the rest is easy. just hop through the forward pointers
            for(int i = 0; i < forwardPointers.length; i++)
            {
                currPointer = forwardPointers[i][currPointer];
                IndexedSnpInterval currMaxKInterval =
                    uberCores.get(i + 1).get(currPointer);
                maxKIntervals.add(currMaxKInterval);
            }
        }
        
        return maxKIntervals;
    }
}
