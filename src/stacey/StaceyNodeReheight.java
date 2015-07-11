/*

        This file is part of STACEY.

        STACEY is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        STACEY is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with STACEY.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 * Derived from BEAST's NodeReheight.
 */

package stacey;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import stacey.debugtune.Checks;
import stacey.util.Bindings;
import stacey.util.BitUnion;
import stacey.util.UnionArrays;

import java.util.*;

@Description("Faster implementation of BEAST's NodeReheight move. " +
        "It implements the method of Newton, Mau, and Larget (1999) on the SMC-tree while maintaining compatibuility with gene trees.")
public class StaceyNodeReheight extends Operator {

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<Tree> smcTreeInput =
            new Input<Tree>("smcTree",
                    "The species tree or minimal clusters tree", Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<List<Tree>> geneTreesInput =
            new Input<List<Tree>>("geneTree",
                    "All gene trees",
                    new ArrayList<Tree>());

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<Long> delayInput =
            new Input<>("delay",
                    "Number of times the operator is disabled.");

    private Tree sTree;
    private List<Tree> gTrees;
    private long delay = 0;
    private int callCount = 0;
    private boolean sTreeTooSmall;
    private UnionArrays unionArrays;

    private final boolean debugFlag = Boolean.valueOf(System.getProperty("stacey.debug"));
    private int numberofdebugchecks = 0;
    private final static int maxnumberofdebugchecks = 100000;



    @Override
    public void initAndValidate() throws Exception {
        sTree = smcTreeInput.get();
        gTrees = geneTreesInput.get();
        sTreeTooSmall = (sTree.getLeafNodeCount() < 3);
        if (delayInput.get() != null) {
            delay = delayInput.get().longValue();
        }
        Bindings bindings = Bindings.initialise(sTree, gTrees);
        unionArrays = UnionArrays.initialise(sTree, gTrees, bindings);
    }



    @Override
    public double proposal() {

        // TODO Input says use get(this). See TODO in CoordinatedPruneRegraft

        if (sTreeTooSmall) {
            return Double.NEGATIVE_INFINITY;
        }
        callCount++;
        if (callCount < delay) {
            return Double.NEGATIVE_INFINITY;
        }

        if (debugFlag  &&  numberofdebugchecks < maxnumberofdebugchecks) {
            Checks.allTreesAndCompatibility(sTree, gTrees, "StaceyNodeReheight", "before move");
            numberofdebugchecks++;
        }

        unionArrays.update();
        double logHR = doNodeReheightMove();        //  The business
        unionArrays.reset();

        if (debugFlag  &&  numberofdebugchecks < maxnumberofdebugchecks) {
            Checks.allTreesAndCompatibility(sTree, gTrees, "StaceyNodeReheight", "after move");
            numberofdebugchecks++;
        }
        return logHR;
    }




/***********************************************************************************/


    private double doNodeReheightMove() {
        Node [] sNodes = sTree.getNodesAsArray();
        // randomly change left/right order
        sTree.startEditing(this);  // we change the tree
        reorder(sTree.getRoot());
        // collect heights
        final double[] fHeights = new double[sNodes.length];
        final int[] iReverseOrder = new int[sNodes.length];
        collectHeights(sTree.getRoot(), fHeights, iReverseOrder, 0);
        // change height of an internal node
        int iNode = Randomizer.nextInt(fHeights.length);
        while (sTree.getNodesAsArray()[iReverseOrder[iNode]].isLeaf()) {
            iNode = Randomizer.nextInt(fHeights.length);
        }
        final double fMaxHeight = calcMaxHeight(iReverseOrder, iNode);
        fHeights[iNode] = Randomizer.nextDouble() * fMaxHeight;
        sNodes[iReverseOrder[iNode]].setHeight(fHeights[iNode]);
        // reconstruct tree from heights
        final Node root = reconstructTree(fHeights, iReverseOrder, 0,
                      fHeights.length, new boolean[fHeights.length]);
        assert checkConsistency(root, new boolean[fHeights.length]) ;
        root.setParent(null);
        sTree.setRoot(root);
        return 0.0;
    }





    private boolean checkConsistency(final Node node, final boolean[] bUsed) {
        if (bUsed[node.getNr()]) {
            // used twice? tha's bad
            return false;
        }
        bUsed[node.getNr()] = true;
        if ( node.isLeaf() ) {
            return true;
        }
        return checkConsistency(node.getLeft(), bUsed) && checkConsistency(node.getRight(), bUsed);
    }

    /**
     * calculate maximum height that node iNode can become restricted
     * by nodes on the left and right
     */
    // This is speed-critical. I've seen the BEAST version take 40% of overall runtime. GRJ.
    private double calcMaxHeight(final int[] iReverseOrder, final int iNode) {
        // find species (or minimal clusters) on the left of selected node
        final BitUnion sppL = new BitUnion(sTree.getLeafNodeCount());
        final Node[] nodes = sTree.getNodesAsArray();
        for (int i = 0; i < iNode; i++) {
            final Node node = nodes[iReverseOrder[i]];
            if (node.isLeaf()) {
                sppL.insert(node.getNr());
            }
        }
        // find species (or minimal clusters) on the right of selected node
        BitUnion sppR = new BitUnion(sTree.getLeafNodeCount());
        for (int i = iNode + 1; i < nodes.length; i++) {
            final Node node = nodes[iReverseOrder[i]];
            if (node.isLeaf()) {
                sppR.insert(node.getNr());
            }
        }
        // for each gene tree, get the straddlers, ie those nodes whose
        // set of species (or minimal clusters) overlaps both sppL and sppR
        double maxHeight = Double.POSITIVE_INFINITY;
        for (int j = 0; j < gTrees.size(); j++) {
            ArrayList<Node> straddlers = unionArrays.getStraddlers(j, sppL, sppR);
            for (Node straddler : straddlers) {
                maxHeight = Math.min(maxHeight, straddler.getHeight());
            }
        }
        return maxHeight;
    }



    /**
     * construct tree top down by joining heighest left and right nodes *
     */
    private Node reconstructTree(final double[] fHeights, final int[] iReverseOrder,
                                 final int iFrom, final int iTo, final boolean[] bHasParent) {
        Node [] sNodes = sTree.getNodesAsArray();
        //iNode = maxIndex(fHeights, 0, fHeights.length);
        int iNode = -1;
        double fMax = Double.NEGATIVE_INFINITY;
        for (int j = iFrom; j < iTo; j++) {
            if (fMax < fHeights[j] && !sNodes[iReverseOrder[j]].isLeaf()) {
                fMax = fHeights[j];
                iNode = j;
            }
        }
        if (iNode < 0) {
            return null;
        }
        final Node node = sNodes[iReverseOrder[iNode]];

        //int iLeft = maxIndex(fHeights, 0, iNode);
        int iLeft = -1;
        fMax = Double.NEGATIVE_INFINITY;
        for (int j = iFrom; j < iNode; j++) {
            if (fMax < fHeights[j] && !bHasParent[j]) {
                fMax = fHeights[j];
                iLeft = j;
            }
        }

        //int iRight = maxIndex(fHeights, iNode+1, fHeights.length);
        int iRight = -1;
        fMax = Double.NEGATIVE_INFINITY;
        for (int j = iNode + 1; j < iTo; j++) {
            if (fMax < fHeights[j] && !bHasParent[j]) {
                fMax = fHeights[j];
                iRight = j;
            }
        }

        node.setLeft(sNodes[iReverseOrder[iLeft]]);
        node.getLeft().setParent(node);
        node.setRight(sNodes[iReverseOrder[iRight]]);
        node.getRight().setParent(node);
        if (node.getLeft().isLeaf()) {
            fHeights[iLeft] = Double.NEGATIVE_INFINITY;
        }
        if (node.getRight().isLeaf()) {
            fHeights[iRight] = Double.NEGATIVE_INFINITY;
        }
        bHasParent[iLeft] = true;
        bHasParent[iRight] = true;
        fHeights[iNode] = Double.NEGATIVE_INFINITY;


        reconstructTree(fHeights, iReverseOrder, iFrom, iNode, bHasParent);
        reconstructTree(fHeights, iReverseOrder, iNode, iTo, bHasParent);
        return node;
    }


    /**
     ** gather height of each node, and the node index associated with the height.*
     **/
    private int collectHeights(final Node node, final double[] fHeights, final int[] iReverseOrder, int iCurrent) {
        if (node.isLeaf()) {
            fHeights[iCurrent] = node.getHeight();
            iReverseOrder[iCurrent] = node.getNr();
            iCurrent++;
        } else {
            iCurrent = collectHeights(node.getLeft(), fHeights, iReverseOrder, iCurrent);
            fHeights[iCurrent] = node.getHeight();
            iReverseOrder[iCurrent] = node.getNr();
            iCurrent++;
            iCurrent = collectHeights(node.getRight(), fHeights, iReverseOrder, iCurrent);
        }
        return iCurrent;
    }

    /**
     * randomly changes left and right children in every internal node *
     */
    private void reorder(final Node node) {
        if (!node.isLeaf()) {
            if (Randomizer.nextBoolean()) {
                final Node tmp = node.getLeft();
                node.setLeft(node.getRight());
                node.setRight(tmp);
            }
            reorder(node.getLeft());
            reorder(node.getRight());
        }
    }



}
