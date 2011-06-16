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
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.jax.util.ObjectUtil;
import org.jax.util.datastructure.ListComparator;
import org.jax.util.io.IllegalFormatException;
/**
 * A node in a perfect phylogeny tree
 * @see PhylogenyTreeEdge
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class PhylogenyTreeNode implements Cloneable, Serializable
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = 601593203889453930L;

    private static final Logger LOG = Logger.getLogger(
            PhylogenyTreeNode.class.getName());
    
    // Group 1: the recursive part of the newick string
    private static final Pattern NEWICK_OUTER_PATTERN = Pattern.compile(
            "^(.*);$");
    
    // Group 2: the edge list (if null then leaf node)
    // Group 3: the node name (if null, anonymous)
    private static final Pattern NEWICK_NODE_PATTERN = Pattern.compile(
            "^(\\((.+)\\))?([^\\(\\)]*)$");
    private static final int EDGE_LIST_GROUP_NUMBER = 2;
    private static final int NODE_NAME_GROUP_NUMBER = 3;
    
    // Group 1: the node part
    // Group 2: the edge distance
    private static final Pattern NEWICK_EDGE_PATTERN = Pattern.compile(
            "^(.*):(.*)$");
    private static final int NODE_PART_GROUP_NUMBER = 1;
    private static final int EDGE_DISTANCE_GROUP_NUMBER = 2;
    
    private List<PhylogenyTreeEdge> childEdges;
    
    private List<String> strains;
    
    private static final String STRAIN_SEPARATOR = "|";
    private static final char NEWICK_CHILD_SEPARATOR = ',';
    
    /**
     * Constructor
     */
    public PhylogenyTreeNode()
    {
        this(new ArrayList<PhylogenyTreeEdge>(), new ArrayList<String>());
    }
    
    /**
     * Constructor
     * @param childEdges
     *          edges pointing to child nodes
     * @param strains
     *          the strains that exist in this node of the phylogeny
     */
    public PhylogenyTreeNode(
            List<PhylogenyTreeEdge> childEdges,
            List<String> strains)
    {
        this.childEdges = childEdges;
        this.strains = strains;
    }
    
    /**
     * Getter for the edges pointing to children
     * @return the childEdges
     */
    public List<PhylogenyTreeEdge> getChildEdges()
    {
        return this.childEdges;
    }
    
    /**
     * Setter for edges pointing to children
     * @param childEdges the childEdges to set
     */
    public void setChildEdges(List<PhylogenyTreeEdge> childEdges)
    {
        this.childEdges = childEdges;
    }
    
    /**
     * Getter for the strains at this node
     * @return the strains
     */
    public List<String> getStrains()
    {
        return this.strains;
    }
    
    /**
     * Recusively get all strains at or below this node
     * @return
     *          the strains
     */
    public List<String> getAllStrains()
    {
        List<String> allStrains = new ArrayList<String>(this.strains);
        
        for(PhylogenyTreeEdge childEdge: this.childEdges)
        {
            allStrains.addAll(childEdge.getNode().getAllStrains());
        }
        
        return allStrains;
    }
    
    /**
     * Find all of the leaf nodes under this node (including this node if
     * its a leaf)
     * @return
     *          the leaf nodes
     */
    public List<PhylogenyTreeNode> getAllLeafNodes()
    {
        List<PhylogenyTreeNode> leaves = new ArrayList<PhylogenyTreeNode>();
        
        this.getAllLeafNodesRecursive(leaves);
        
        return leaves;
    }
    
    /**
     * Find all of the leaf nodes under this node and add them to the
     * leaves list
     * @param leaves
     *          the leaves
     */
    private void getAllLeafNodesRecursive(List<PhylogenyTreeNode> leaves)
    {
        if(this.isLeafNode())
        {
            leaves.add(this);
        }
        else
        {
            for(PhylogenyTreeEdge childEdge: this.getChildEdges())
            {
                PhylogenyTreeNode childNode = childEdge.getNode();
                childNode.getAllLeafNodesRecursive(leaves);
            }
        }
    }

    /**
     * Setter for the strains at this node
     * @param strains the strains to set
     */
    public void setStrains(List<String> strains)
    {
        this.strains = strains;
    }
    
    /**
     * Convenience function for determining if this is a leaf (no children)
     * @return
     *          true iff there are no children
     */
    public boolean isLeafNode()
    {
        return this.childEdges.isEmpty();
    }
    
    /**
     * Convert a newick formatted string into a phylogeny tree
     * @param newickFormattedTree
     *          the newick string
     * @return
     *          the node
     * @throws IllegalFormatException
     *          if bad string formating is detected
     */
    public static PhylogenyTreeNode fromNewickFormat(String newickFormattedTree)
    throws IllegalFormatException
    {
        Matcher outerTreeMatcher = NEWICK_OUTER_PATTERN.matcher(newickFormattedTree);
        if(outerTreeMatcher.matches())
        {
            String nodePart = outerTreeMatcher.group(NODE_PART_GROUP_NUMBER);
            return PhylogenyTreeNode.fromNewickNode(nodePart);
        }
        else
        {
            throw new IllegalFormatException(
                    "can't convert " + newickFormattedTree +
                    " to phylogeny tree");
        }
    }
    
    /**
     * Convert a newick node fragment into a phylogeny node
     * @param newickNodeString
     *          the newick formatted string
     * @return
     *          the node
     * @throws IllegalFormatException
     *          if we detect bad formatting
     */
    private static PhylogenyTreeNode fromNewickNode(String newickNodeString)
            throws IllegalFormatException
    {
        Matcher nodeMatcher = NEWICK_NODE_PATTERN.matcher(newickNodeString);
        if(nodeMatcher.matches())
        {
            String nodeNamePart = nodeMatcher.group(NODE_NAME_GROUP_NUMBER);
            String edgeListPart = nodeMatcher.group(EDGE_LIST_GROUP_NUMBER);
            
            PhylogenyTreeNode node = new PhylogenyTreeNode();
            if(nodeNamePart != null && nodeNamePart.length() >= 1)
            {
                node.setStrains(PhylogenyTreeNode.newickNodeNameToStrainList(
                        nodeNamePart));
            }
            
            if(edgeListPart != null && edgeListPart.length() >= 1)
            {
                node.setChildEdges(PhylogenyTreeNode.newickEdgesToEdgeList(
                        edgeListPart));
            }
            return node;
        }
        else
        {
            throw new IllegalFormatException(
                    "failed to parse " + newickNodeString + " as a node");
        }
    }

    /**
     * Convert the comma separated newick edge list into
     * @param edgeListPart
     *          the newick formatted edge list fragment
     * @return
     *          the edges
     * @throws IllegalFormatException
     *          if we detect bad formatting
     */
    private static List<PhylogenyTreeEdge> newickEdgesToEdgeList(
            String edgeListPart)
            throws IllegalFormatException
    {
        List<PhylogenyTreeEdge> edges = new ArrayList<PhylogenyTreeEdge>();
        int currEdgeStartIndex = 0;
        int unmatchedParenCount = 0;
        for(int i = 0; i < edgeListPart.length(); i++)
        {
            char currChar = edgeListPart.charAt(i);
            if(currChar == '(')
            {
                unmatchedParenCount++;
            }
            else if(currChar == ')')
            {
                unmatchedParenCount--;
            }
            else if(unmatchedParenCount == 0 && currChar == NEWICK_CHILD_SEPARATOR)
            {
                String edgeString = edgeListPart.substring(
                        currEdgeStartIndex,
                        i);
                edges.add(PhylogenyTreeNode.newickEdgeToPhlylogenyEdge(
                        edgeString));
                currEdgeStartIndex = i + 1;
            }
        }
        
        if(unmatchedParenCount != 0)
        {
            throw new IllegalFormatException(
                    "unbalanced parentheses in edge list: " + edgeListPart);
        }
        else
        {
            String edgeString = edgeListPart.substring(
                    currEdgeStartIndex);
            edges.add(PhylogenyTreeNode.newickEdgeToPhlylogenyEdge(
                    edgeString));
            return edges;
        }
    }

    /**
     * Convert a single edge string
     * @param edgeString
     *          the newick edge string fragment
     * @return
     *          the edge
     * @throws IllegalFormatException
     *          if we detect bad formatting
     */
    private static PhylogenyTreeEdge newickEdgeToPhlylogenyEdge(
            String edgeString)
            throws IllegalFormatException
    {
        Matcher edgeMatcher = NEWICK_EDGE_PATTERN.matcher(edgeString);
        if(edgeMatcher.matches())
        {
            String edgeLengthPart = edgeMatcher.group(EDGE_DISTANCE_GROUP_NUMBER).trim();
            String edgeNodePart = edgeMatcher.group(NODE_PART_GROUP_NUMBER).trim();
            
            PhylogenyTreeEdge edge = new PhylogenyTreeEdge();
            edge.setEdgeLength(Double.parseDouble(edgeLengthPart));
            edge.setNode(PhylogenyTreeNode.fromNewickNode(
                    edgeNodePart));
            return edge;
        }
        else
        {
            throw new IllegalFormatException(
                    "Can't parse edge: " + edgeString);
        }
    }

    /**
     * Convert the newick node name into a strain list using
     * {@link #STRAIN_SEPARATOR} as the separator value
     * @param newickNodeName
     *          the node name to turn into a strain list
     * @return
     *          the list
     */
    private static List<String> newickNodeNameToStrainList(String newickNodeName)
    {
        if(newickNodeName == null)
        {
            // use an empty list
            return new ArrayList<String>();
        }
        else
        {
            StringTokenizer strainTokenizer = new StringTokenizer(
                    newickNodeName,
                    STRAIN_SEPARATOR);
            List<String> strainList = new ArrayList<String>(
                    strainTokenizer.countTokens());
            while(strainTokenizer.hasMoreTokens())
            {
                strainList.add(strainTokenizer.nextToken());
            }
            return strainList;
        }
    }
    
    /**
     * Create a newick formatted string representing this tree
     * @return
     *          the newick string representation
     */
    public String toNewickFormat()
    {
        return this.toNewickFormatRecursive() + ';';
    }
    
    /**
     * The recursive part of the newick formatting logic
     * @return
     *          the recursive part of the format (everything but the ';')
     */
    private String toNewickFormatRecursive()
    {
        StringBuffer sb = new StringBuffer();
        
        // generate children recursively. separate them with commas
        if(!this.isLeafNode())
        {
            sb.append('(');
            boolean firstChildEdge = true;
            for(PhylogenyTreeEdge childEdge: this.childEdges)
            {
                if(firstChildEdge)
                {
                    firstChildEdge = false;
                }
                else
                {
                    sb.append(NEWICK_CHILD_SEPARATOR);
                }
                sb.append(childEdge.getNode().toNewickFormatRecursive());
                sb.append(':');
                sb.append(childEdge.getEdgeLength());
            }
            sb.append(')');
        }
        
        // add the node name which is all of the strains separated with a bar
        boolean firstStrain = true;
        for(String strainName: this.strains)
        {
            if(firstStrain)
            {
                firstStrain = false;
            }
            else
            {
                sb.append(STRAIN_SEPARATOR);
            }
            sb.append(strainName);
        }
        
        return sb.toString();
    }
    
    /**
     * For a interior nodes with any strains or leaf nodes with more than one
     * strain, push all of those strains out onto new leaf nodes whose edges
     * are length epsilon
     * @param epsilon
     *          a (usually) small value used to push the nodes out onto
     *          individual leaves
     * @return
     *          the tree (this tree is not modified)
     */
    public PhylogenyTreeNode resolveToSingleStrainLeafNodes(
            double epsilon)
    {
        List<PhylogenyTreeEdge> newEdges = new ArrayList<PhylogenyTreeEdge>();
        
        // recursively resolve children
        for(PhylogenyTreeEdge edge: this.childEdges)
        {
            newEdges.add(new PhylogenyTreeEdge(
                    edge.getSdpBits(),
                    edge.getNode().resolveToSingleStrainLeafNodes(epsilon),
                    edge.getEdgeLength()));
        }
        
        // now resolve this node
        int strainCount = this.strains.size();
        if(strainCount >= 2 || !this.isLeafNode() && strainCount >= 1)
        {
            // we need to resolve the strains to new leaf nodes
            for(String strainName: this.strains)
            {
                List<PhylogenyTreeEdge> emptyEdges = Collections.emptyList();
                PhylogenyTreeNode newLeafNode = new PhylogenyTreeNode(
                        emptyEdges,
                        Collections.singletonList(strainName));
                newEdges.add(new PhylogenyTreeEdge(
                        new BitSet(0),
                        newLeafNode,
                        epsilon));
            }
        }
        
        // build the new node
        if(newEdges.isEmpty())
        {
            return new PhylogenyTreeNode(
                    newEdges,
                    new ArrayList<String>(this.strains));
        }
        else
        {
            List<String> emptyStrains = Collections.emptyList();
            return new PhylogenyTreeNode(
                    newEdges,
                    emptyStrains);
        }
    }
    
    /**
     * Creates a new tree with all of the non-branching interior nodes
     * removed. The when an interior node is removed, the newly created edge
     * length will equal (removed parent edge + removed child edge) so the
     * overall size of the tree will not change. This function does not
     * modify this tree.
     * @return
     *          the new tree
     */
    public PhylogenyTreeNode removeNonBranchingInteriorNodes()
    {
        List<PhylogenyTreeEdge> newEdges = new ArrayList<PhylogenyTreeEdge>();
        
        for(PhylogenyTreeEdge edge: this.childEdges)
        {
            // recursively resolve child
            PhylogenyTreeEdge newEdge = new PhylogenyTreeEdge(
                    edge.getSdpBits(),
                    edge.getNode().removeNonBranchingInteriorNodes(),
                    edge.getEdgeLength());
            
            // now collapse child edge if we need to
            if(newEdge.getNode().getChildEdges().size() == 1)
            {
                PhylogenyTreeEdge grandChildEdge =
                    newEdge.getNode().getChildEdges().get(0);
                newEdge = new PhylogenyTreeEdge(
                        newEdge.getSdpBits(),
                        grandChildEdge.getNode(),
                        newEdge.getEdgeLength() + grandChildEdge.getEdgeLength());
            }
            
            newEdges.add(newEdge);
        }
        
        return new PhylogenyTreeNode(
                newEdges,
                new ArrayList<String>(this.strains));
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(Object obj)
    {
        if(obj instanceof PhylogenyTreeNode)
        {
            PhylogenyTreeNode otherNode = (PhylogenyTreeNode)obj;
            
            return ObjectUtil.areEqual(this.strains, otherNode.strains) &&
                   ObjectUtil.areEqual(this.childEdges, otherNode.childEdges);
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
        return ObjectUtil.hashObject(this.strains);
    }
    
    /**
     * Sets all SDPs to null recursively
     * @return
     *          the resulting tree
     */
    public PhylogenyTreeNode createSdpFreeTree()
    {
        PhylogenyTreeNode zeroedTree = this.clone();
        for(PhylogenyTreeEdge childEdge: zeroedTree.childEdges)
        {
            childEdge.setSdpBits(null);
            childEdge.getNode().createSdpFreeTree();
        }
        return zeroedTree;
    }
    
    /**
     * Normalize the tree. Any trees w/c are "unrooted equivalent" will
     * normalize to the same rooted tree for this function
     * @return
     *          the newly normalized tree
     */
    public PhylogenyTreeNode createNormalizedTree()
    {
        // make a copy so that none of our changes affect the original
        PhylogenyTreeNode thisCopy = this.clone();
        
        // recursively sort strain lists so that the ordering is consistent
        thisCopy.sortStrainsRecursive();
        
        // normalize around a deterministic root
        List<PhylogenyTreeNode> allNodes = new ArrayList<PhylogenyTreeNode>();
        thisCopy.getAllNodesWithStrainsRecursive(allNodes);
        if(allNodes.isEmpty())
        {
            LOG.severe("failed to find any nodes containing strains!");
            return null;
        }
        else
        {
            ListComparator<String> strListComp = new ListComparator<String>();
            PhylogenyTreeNode deterministicRoot = allNodes.get(0);
            int allNodesCount = allNodes.size();
            for(int i = 0; i < allNodesCount; i++)
            {
                PhylogenyTreeNode node = allNodes.get(i);
                if(deterministicRoot == null ||
                   strListComp.compare(node.strains, deterministicRoot.strains) < 0)
                {
                    deterministicRoot = node;
                }
            }
            
            if(thisCopy == deterministicRoot)
            {
                return thisCopy;
            }
            else
            {
                if(thisCopy.reRootRecursive(deterministicRoot))
                {
                    deterministicRoot.sortChildEdgesRecursive();
                    return deterministicRoot;
                }
                else
                {
                    LOG.severe("failed to reroot tree during normalization");
                    return null;
                }
            }
        }
    }
    
    /**
     * Create a copy of this tree that has the given strains removed. In
     * addition, any branches that are empty after strain removal are pruned
     * @param strainsToRetain
     *          the strains that we should keep
     * @return
     *          the pruned tree
     */
    public PhylogenyTreeNode createStrainPrunedTree(Set<String> strainsToRetain)
    {
        PhylogenyTreeNode thisCopy = this.clone();
        thisCopy.createStrainPrunedTreeRecursive(strainsToRetain);
        
        return thisCopy;
    }

    /**
     * Recursively prune
     * @param strainsToRetain
     *          the strains we should keep
     * @return
     *          true iff there are any strains under this branch
     */
    private boolean createStrainPrunedTreeRecursive(Set<String> strainsToRetain)
    {
        boolean anyChildren = false;
        Iterator<PhylogenyTreeEdge> childEdgeIter = this.childEdges.iterator();
        while(childEdgeIter.hasNext())
        {
            PhylogenyTreeEdge childEdge = childEdgeIter.next();
            if(childEdge.getNode().createStrainPrunedTreeRecursive(strainsToRetain))
            {
                anyChildren = true;
            }
            else
            {
                childEdgeIter.remove();
            }
        }
        
        this.strains.retainAll(strainsToRetain);
        
        return anyChildren || !this.strains.isEmpty();
    }

    /**
     * Sort the child edges for tree normalization
     */
    private void sortChildEdgesRecursive()
    {
        // strain lists should all be mutually exclusive so this should be
        // a consistent sort
        List<String> allStrains = new ArrayList<String>(); 
        this.getAllStrainsRecursive(allStrains);
        Collections.sort(allStrains);
    }
    
    /**
     * Get all strains in this node and children
     * @param allStrains
     *          the list to fill
     */
    private void getAllStrainsRecursive(List<String> allStrains)
    {
        allStrains.addAll(this.strains);
        for(PhylogenyTreeEdge childEdge: this.childEdges)
        {
            childEdge.getNode().getAllStrainsRecursive(allStrains);
        }
    }
    
    /**
     * Reroot this tree using the given root node. This function creates
     * a new tree without modifying this tree
     * @param newRoot
     *          the new root node
     * @return
     *          the re-rooted tree or null if this tree doesn't contain
     *          the given new root
     */
    public PhylogenyTreeNode reRoot(PhylogenyTreeNode newRoot)
    {
        PhylogenyTreeNode thisCopy = this.clone();
        PhylogenyTreeNode matchingNewRoot = this.findMatchingCopyNode(
                thisCopy,
                this,
                newRoot);
        if(thisCopy.reRootRecursive(matchingNewRoot))
        {
            return matchingNewRoot;
        }
        else
        {
            return null;
        }
    }

    /**
     * This function finds a match for "nodeToMatch" in the copy tree.
     * This is done assuming that the structure of the copy exactly matches
     * the structure of the original
     * @param copyTree
     *          the copy
     * @param originalTree
     *          the original
     * @param nodeToMatch
     *          the node that we're tring to find in the original so that we
     *          can return the copy's match
     * @return
     *          the "structural" match from the copy tree (no equality
     *          value checks are performed)
     */
    private PhylogenyTreeNode findMatchingCopyNode(
            PhylogenyTreeNode copyTree,
            PhylogenyTreeNode originalTree,
            PhylogenyTreeNode nodeToMatch)
    {
        if(originalTree == nodeToMatch)
        {
            return copyTree;
        }
        else
        {
            int edgeCount = originalTree.getChildEdges().size();
            for(int i = 0; i < edgeCount; i++)
            {
                PhylogenyTreeNode match = this.findMatchingCopyNode(
                        copyTree.getChildEdges().get(i).getNode(),
                        originalTree.getChildEdges().get(i).getNode(),
                        nodeToMatch);
                
                if(match != null)
                {
                    return match;
                }
            }
            
            return null;
        }
    }

    /**
     * Reroot this [sub]tree recursively using {@code newRoot} as the
     * new tree root
     * @param newRoot
     *          the new root node
     * @return
     *          true if re-rooting succeeded
     */
    private boolean reRootRecursive(PhylogenyTreeNode newRoot)
    {
        if(this == newRoot)
        {
            return true;
        }
        else
        {
            Iterator<PhylogenyTreeEdge> childIter = this.childEdges.iterator();
            while(childIter.hasNext())
            {
                PhylogenyTreeEdge nextChild = childIter.next();
                if(nextChild.getNode().reRootRecursive(newRoot))
                {
                    // reverse the edge
                    childIter.remove();
                    nextChild.getNode().childEdges.add(new PhylogenyTreeEdge(
                            nextChild.getSdpBits(),
                            this,
                            nextChild.getEdgeLength()));
                    return true;
                }
            }
            
            return false;
        }
    }

    /**
     * Recursively grab any nodes that have associated strains and throw them
     * into {@code nodes}
     * @param nodes
     *          the nodes list to fill up
     */
    private void getAllNodesWithStrainsRecursive(List<PhylogenyTreeNode> nodes)
    {
        if(!this.strains.isEmpty())
        {
            nodes.add(this);
        }
        
        for(PhylogenyTreeEdge childEdge: this.childEdges)
        {
            childEdge.getNode().getAllNodesWithStrainsRecursive(nodes);
        }
    }

    /**
     * Recursively modifies the strains lists so that they're all sorted.
     */
    private void sortStrainsRecursive()
    {
        Collections.sort(this.strains);
        
        for(PhylogenyTreeEdge childEdge: this.childEdges)
        {
            childEdge.getNode().sortStrainsRecursive();
        }
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public PhylogenyTreeNode clone()
    {
        try
        {
            PhylogenyTreeNode newNode = (PhylogenyTreeNode)super.clone();
            
            newNode.strains = new ArrayList<String>(this.strains);
            
            newNode.childEdges = new ArrayList<PhylogenyTreeEdge>(this.childEdges.size());
            for(PhylogenyTreeEdge child: this.childEdges)
            {
                    newNode.childEdges.add(child.clone());
            }
            
            return newNode;
        }
        catch(CloneNotSupportedException ex)
        {
            LOG.log(Level.SEVERE,
                    "Clone failed. This should never happen",
                    ex);
            return null;
        }
    }
}
