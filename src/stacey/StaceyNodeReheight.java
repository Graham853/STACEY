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
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import stacey.debugtune.Checks;
import stacey.util.Bindings;
import stacey.util.BitUnion;
import stacey.util.UnionArrays;

import java.util.*;

@Description("Faster implementation of BEAST's NodeReheight move. " +
        "It implements the method of Newton, Mau, and Larget (1999) on the SMC-tree while maintaining compatibility with gene trees.")
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

    public Input<RealParameter> popSFInput =
            new Input<>("popSF",
                    "The population scaling factor for the STACEY coalescent", Input.Validate.REQUIRED);

    public Input<Double> propUniformInput =
            new Input<>("proportionUniform",
                    "The fraction of times which the operator uses a uniform density for sampling new heights. " +
                            "It must be between 0.0 and 1.0. " +
                            "The rest of the time the operator uses a density skewed towards the maximum compatible height.");

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

    private double propUniform;

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
        propUniform = 0.1;
        if (propUniformInput.get() != null) {
            propUniform = propUniformInput.get().doubleValue();
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
        final double maxHeight = calcMaxHeight(iReverseOrder, iNode);

        double newHeight;
        double logHR;
        if (Randomizer.nextDouble() > propUniform) {
            /*
            Non-uniform sampling of new heights. NodeReheight is rarely accepted for a node when there
            are large number of loci (unless the node has a very small height.) After burnin, only small
            changes in heights will work. So seems here new heights are concentrated near the max.
            This requires HRs to counter the bias of course, but should still help.

            The cdf we're going to sample from is defined in [0,h] as
            F(x) = (log(h+a) - log(h+a-x)) / (log(h+a)-log(a))
            where h is maxHeight and a is hgtS = popSF / gTrees.size();
            The inverse is
            x = G(y) = h + a - exp(  log(h+a) (1-y)  +  log(a) y  )
            The pdf is
            f(x) = 1 / ( (log(h+a)-log(a)) (h+a-x) )
            The median is G(.5) = h+a - sqrt(a(h+a))
            Eg h = 0.001, a = 0.00001, G(.5) =  0.00101-sqrt(.00001*.00101) = 0.0009095
            Eg h = 0.001, a = 0.0000001, G(.5) =  0.0010001-sqrt(.0000001*.0010001) = 0.0009900995
            For large h/a, median ~= h(1 - 1/sqrt(h/a))
            */

            double popSF = popSFInput.get().getValue();
            double hgtS = 0.1 * popSF / gTrees.size();
            double oldHeight = fHeights[iNode];
            newHeight = newHeightSample(maxHeight, hgtS);
            double oldDensity = newHeightPDF(oldHeight, maxHeight, hgtS);
            double newDensity = newHeightPDF(newHeight, maxHeight, hgtS);
            logHR = Math.log(oldDensity/newDensity);
        } else {
            newHeight = Randomizer.nextDouble() * maxHeight;
            logHR = 0.0;
        }
        fHeights[iNode] = newHeight;
        sNodes[iReverseOrder[iNode]].setHeight(fHeights[iNode]);
        // reconstruct tree from heights
        final Node root = reconstructTree(fHeights, iReverseOrder, 0,
                fHeights.length, new boolean[fHeights.length]);
        assert checkConsistency(root, new boolean[fHeights.length]) ;
        root.setParent(null);
        sTree.setRoot(root);
        return logHR;



    }





    private double newHeightPDF(double x, double h, double a) {
        double d = 1.0 / (Math.log(h + a) - Math.log(a));
        d /= (h + a - x);
        return d;
    }

    // not used but maybe for debugging
    /*private double newHeightCDF(double x, double h, double a) {
        double p = 1.0 / (Math.log(h + a) - Math.log(a));
        p *= (Math.log(h + a) - Math.log(h + a - x));
        return p;
    }*/

    private double newHeightSample(double h, double a) {
        double y = Randomizer.nextDouble();
        double q = h + a - Math.exp(Math.log(h + a) * (1 - y) + Math.log(a) * y);
        q = Math.max(q, 0.0);
        q = Math.min(q, h);
        return q;
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
