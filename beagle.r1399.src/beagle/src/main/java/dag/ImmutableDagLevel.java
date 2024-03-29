/*
 * Copyright (C) 2014 Brian L. Browning
 *
 * This file is part of Beagle
 *
 * Beagle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Beagle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package dag;

import blbutil.Const;
import blbutil.IndexSet;
import java.util.Arrays;
import vcf.Marker;

/**
 * <p>Class {@code ImmutableDagLevel} represents a level of a leveled
 * directed acyclic graph (DAG).
 * </p>
 * Instances of {@code ImmutableDagLevel} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class ImmutableDagLevel implements DagLevel {

/*
 * The k-th edge parent node index is stored in {@code this.parentNodes[k]}.
 * The k-th edge child node index is stored in {@code this.childNodes[k]}.
 * The k-th edge symbol is stored in {@code this.symbols[k]}.
 * The k-th edge count is stored in {@code this.edgeCounts[k]}.
 * The k-th edge conditional edge probability is stored in
 * {@code this.condEdgeProbs[k]}, and is defined to be the
 * k-th edge count divided by the k-th edge's parent node count.
 * The k-th node count is stored in {@code this.nodeCounts[k]}.
 *
 * The outgoing edges indices of the k-th parent node are stored in consecutive
 * entries of {@code this.parents} beginning with
 * {@code this.parentIndices[k]} (inclusive) and ending with
 * {@code this.parentIndices[k+1]} (exclusive).
 *
 * The ingoing edges indices of the k-th child node are stored in consecutive
 * entries of {@code this.children} beginning with
 * {@code this.childIndices[k]} (inclusive) and ending with
 * {@code this.childIndices[k+1]} (exclusive).
 */
    private final Marker marker;
    private final float count;
    private final char[] parentNodes;
    private final char[] childNodes;
    private final char[] parentIndices;
    private final char[] parents;
    private final char[] childIndices;
    private final char[] children;
    private final byte[] symbols;
    private final float[] edgeCounts;
    private final float[] condEdgeProbs;
    private final float[] nodeCounts;

    /**
     * Constructs a new {@code ImmutableDagLevel} instance.
     *
     * @param marker the marker corresponding to this level of the DAG.
     * @param parentNodes an array mapping edge index to parent node index.
     * @param childNodes an array mapping edge index to child node index.
     * @param symbols an array mapping edge index to the symbol labeling the
     * edge.
     * @param counts an array mapping edge index to edge count.
     *
     * @throws IllegalArgumentException if the specified arrays do not all
     * have the same length
     * @throws IllegalArgumentException if any array has length greater than
     * {@code Character.MAX_VALUE}
     * @throws IllegalArgumentException if any element of the symbols array
     * is negative or equal to {@code Byte.MAX_VALUE}
     * @throws IllegalArgumentException if any two edges have the same
     * parent node and are both labeled with the same symbol
     * @throws IllegalArgumentException if the set of values of the
     * {@code parentNodes} array is not equal to {@code {0, 1, 2, ..., k}} for
     * some {@code k}.
     * @throws IllegalArgumentException if the set of values of the
     * {@code childNodes} array is not equal to {@code {0, 1, 2, ..., k}}
     * for some {@code k}.
     *
     * @throws NullPointerException if any parameter is {@code null}.
     */
    public ImmutableDagLevel(Marker marker, char[] parentNodes,
            char[] childNodes, byte[] symbols, float[] counts) {
        int nEdges = checkLengths(parentNodes, childNodes, symbols, counts);
        checkSymbols(marker, symbols);
        this.marker = marker;
        this.parentIndices = getIndicesArray(parentNodes);
        this.childIndices = getIndicesArray(childNodes);
        this.parentNodes = parentNodes.clone();
        this.childNodes = childNodes.clone();
        this.symbols = symbols.clone();
        this.edgeCounts = counts.clone();
        this.condEdgeProbs = new float[nEdges];
        this.parents = new char[nEdges];
        this.children = new char[nEdges];

        char[] pIndices = Arrays.copyOfRange(parentIndices, 0,
                parentIndices.length-1);
        char[] cIndices = Arrays.copyOfRange(childIndices, 0,
                childIndices.length-1);
        this.nodeCounts = nodeCnts(parentNodes, counts, pIndices.length);
        this.count = sum(this.nodeCounts);

        for (char j=0; j<nEdges; ++j) {
            char p = parentNodes[j];
            char c = childNodes[j];
            this.parents[pIndices[p]++] = j;
            this.children[cIndices[c]++] = j;
            this.condEdgeProbs[j] = counts[j] / nodeCounts[p];
        }
        checkForDuplicateOutEdges(parentIndices, parents, symbols,
                marker.nAlleles());
    }

    private static int checkLengths(char[] parentNodes, char[] childNodes,
            byte[] symbols, float[] counts) {
        if (parentNodes.length > Character.MAX_VALUE) {
            String s = "parentNodes.length>Character.MAX_VALUE";
            throw new IllegalArgumentException(s);
        }
        if ( ((parentNodes.length != childNodes.length)
                || (parentNodes.length != symbols.length))
                || (parentNodes.length != counts.length) ) {
            throw new IllegalArgumentException("inconsistent arrays");
        }
        return parentNodes.length;
    }

    private static void checkSymbols(Marker marker, byte[] symbols) {
        int n = marker.nAlleles();
        for (byte b : symbols) {
            if (b<0 || b>=n) {
                throw new IllegalArgumentException("symbol: " + b);
            }
        }
    }

    private static void checkForDuplicateOutEdges(char[] parentIndices,
            char[] parents, byte[] symbols, int nSymbols) {
        IndexSet indexSet = new IndexSet(nSymbols);
        for (int j=1; j<parentIndices.length; ++j) {
            indexSet.clear();
            for (int k=parentIndices[j-1], n=parentIndices[j]; k<n; ++k) {
                int edge = parents[k];
                if (indexSet.add(symbols[edge])==false) {
                    throw new IllegalArgumentException("duplicate edge");
                }
            }
        }
    }

    private static char[] getIndicesArray(char[] nodes) {
        int[] countArray = elementCounts(nodes);
        char[] indicesArray = new char[countArray.length + 1];
        for (int j=1; j<indicesArray.length; ++j) {
            assert countArray[j-1]>0;
            int x = indicesArray[j-1] + countArray[j-1];
            assert x <= Character.MAX_VALUE;
            indicesArray[j] = (char) x;
        }
        return indicesArray;
    }

    /*
     * Returns an array of length {@code max(nodes) + 1}
     * whose {@code j}-th element is the number of
     * elements of the specified array that have value {@code j}.
     *
     * @param nodes an array of non-negative values.
     * @return an array whose {@code j}-th element is the number of
     * elements of the specified array that have value {@code j}.
     * @throws IllegalArgumenException if set of elements of the
     * specified array is not equal to {@code {0, 1, 2, ..., k}} for some
     * {@code k}.
     */
    private static int[] elementCounts(char[] array) {
        int maxNode = max(array);
        int[] nodeCounts = new int[maxNode + 1];
        for (char c : array) {
            ++nodeCounts[c];
        }
        for (int j=0; j<nodeCounts.length; ++j) {
            if (nodeCounts[j]==0) {
                throw new IllegalArgumentException("no element with value " + j);
            }
        }
        return nodeCounts;
    }

    private static int max(char[] ca) {
        char max = 0;
        for (char c : ca) {
            if (c>max) {
                max=c;
            }
        }
        return max;
    }

    private float sum(float[] fa) {
        float sum = 0.0f;
        for (float f : fa) {
            sum += f;
        }
        return sum;
    }

    private float[] nodeCnts(char[] parentNodes, float[] counts, int nNodes) {
        float[] nodeCnts = new float[nNodes];
        for (int j=0; j<condEdgeProbs.length; ++j) {
            char p = parentNodes[j];
            nodeCnts[p] += counts[j];
        }
        return nodeCnts;
    }

    @Override
    public Marker marker() {
        return marker;
    }

    @Override
    public int nEdges() {
        return condEdgeProbs.length;
    }

    @Override
    public char nParentNodes() {
        return (char) (parentIndices.length - 1);
    }

    @Override
    public int nChildNodes() {
        return (childIndices.length - 1);
    }

    @Override
    public char parentNode(int edge) {
        return parentNodes[edge];
    }

    @Override
    public char childNode(int edge) {
        return childNodes[edge];
    }

    @Override
    public byte symbol(int edge) {
        return symbols[edge];
    }

    @Override
    public float edgeCnt(int edge) {
        return edgeCounts[edge];
    }

    @Override
    public float nodeCnt(int parentNode) {
        return nodeCounts[parentNode];
    }

    @Override
    public float condEdgeProb(int edge) {
        return condEdgeProbs[edge];
    }

    @Override
    public float edgeProb(int edge) {
        return (edgeCounts[edge] / count);
    }

    @Override
    public float nodeProb(int node) {
        return (nodeCounts[node] / count);
    }

    @Override
    public char nOutEdges(int parentNode) {
        return (char) (parentIndices[parentNode+1] - parentIndices[parentNode]);
    }

    @Override
    public char outEdge(int parentNode, int outEdgeIndex) {
        if (outEdgeIndex<0 || outEdgeIndex>=nOutEdges(parentNode)) {
            throw new IndexOutOfBoundsException(String.valueOf(outEdgeIndex));
        }
        return this.parents[parentIndices[parentNode] + outEdgeIndex];
    }

    @Override
    public char outEdgeBySymbol(int parentNode, byte symbol) {
        if (symbol < 0 || symbol >= marker.nAlleles()) {
            throw new IndexOutOfBoundsException(String.valueOf(symbol));
        }
        char start = parentIndices[parentNode];
        char end = parentIndices[parentNode+1];
        for (char j=start; j<end; ++j) {
            char edgeIndex = parents[j];
            if (symbols[edgeIndex]==symbol) {
                return edgeIndex;
            }
        }
        return Character.MAX_VALUE;
    }

    @Override
    public char nInEdges(int childNode) {
        return (char) (childIndices[childNode+1] - childIndices[childNode]);
    }

    @Override
    public char inEdge(int childNode, int inEdgeIndex) {
        if (inEdgeIndex<0 || inEdgeIndex>=nInEdges(childNode)) {
            throw new IndexOutOfBoundsException(String.valueOf(inEdgeIndex));
        }
        return this.children[childIndices[childNode] + inEdgeIndex];
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(500);
        sb.append(Const.nl);
        sb.append("parentNodes=");
        sb.append(charArrayToString(parentNodes));
        sb.append(Const.nl);
        sb.append("childNodes=");
        sb.append(charArrayToString(childNodes));
        sb.append(Const.nl);
        sb.append("symbols=");
        sb.append(Arrays.toString(symbols));
        sb.append(Const.nl);
        sb.append("condEdgeProbs=");
        sb.append(Arrays.toString(condEdgeProbs));
        sb.append(Const.nl);
        sb.append("edgeCounts=");
        sb.append(Arrays.toString(condEdgeProbs));
        sb.append(Const.nl);
        sb.append("nodeCounts=");
        sb.append(Arrays.toString(edgeCounts));
        sb.append(Const.nl);
        sb.append("parentIndices=");
        sb.append(charArrayToString(parentIndices));
        sb.append(Const.nl);
        sb.append("parents=");
        sb.append(charArrayToString(parents));
        sb.append(Const.nl);
        sb.append("childIndices=");
        sb.append(charArrayToString(childIndices));
        sb.append(Const.nl);
        sb.append("children=");
        sb.append(charArrayToString(children));
        return sb.toString();
    }

    /*
     * Returns a string representation of the specified character array
     * with the character elements converted to integers for printing.
     * The representation has the format
     * "[a[1], a[2], a[3], ..., a[length-1]]".
     *
     * @param a a character array.
     *
     * @return a string representation of the specified character array
     * with the character elements converted to integers for printing.
     */
    private static String charArrayToString(char[] a) {
        StringBuilder sb = new StringBuilder(a.length * 4 + 10);
        sb.append("[");
        sb.append((int) a[0]);
        for (int j = 1; j < a.length; ++j) {
            sb.append(", ");
            sb.append((int) a[j]);
        }
        sb.append("]");
        return sb.toString();
    }
}
